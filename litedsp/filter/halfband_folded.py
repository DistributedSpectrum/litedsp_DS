#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

"""Folded symmetric half-band decimator (decimate by 2).

The non-zero taps of a half-band are the centre plus symmetric pairs (c-o, c+o), o odd: P unique
products per output and rail, each a DSP48E1 pre-add (a + d) * coef. The fold factor F spreads
them over ceil(P/F) multiply-accumulate lanes per rail, each lane taking F clocks per output, so:

    DSP48s       = 2 * ceil(P / F)                (I and Q)
    max input    = 2 / F samples per clock        (one output per F clocks)

F = 1 is fully parallel. Input faster than 2/F per clock is back-pressured. Pick F = floor(f_clk / f_out) for the cheapest design that keeps up,
e.g. 245.76 MHz clock and 122.88 MSps in (61.44 MSps out) -> F = 4.

The input window is snapshotted at each launch so the delay line can keep shifting while the
lanes work. Arithmetic: full-precision sum of integer coefficients (signed, cw bits, Q1.(cw-1)),
then round half-up by cw-1 bits and saturate to data_width, as litedsp's FIR decimator does.
Coefficients are writable at runtime (coeff_we/coeff_adr/coeff_dat over the P unique taps).

Measured on xc7a200tsbg484-3 at 245.76 MHz, 39 taps (70% kept, 82.6 dB, 18-bit coefficients):
F=1 22 DSP48, F=2 12, F=4 6, F=8 4, F=11 2; all close timing (268-329 MHz).
"""

from math import ceil

from migen import *
from litex.gen import *
from litex.soc.interconnect import stream

from litedsp.common import iq_layout, rounded, saturated


class LiteDSPHalfbandDecimatorFolded(LiteXModule):
    def __init__(self, coefficients, fold=1, data_width=16, coeff_width=18, phase=0):
        n_taps = len(coefficients)
        assert (n_taps - 3) % 4 == 0, "half-band length must be 4K+3"
        W, CW = data_width, coeff_width
        c     = (n_taps - 1) // 2
        offs  = list(range(1, c + 1, 2))
        for o in offs:
            assert coefficients[c - o] == coefficients[c + o], "half-band must be symmetric"
        uniq  = [coefficients[c]] + [coefficients[c - o] for o in offs]
        pairs = [(c, None)] + [(c - o, c + o) for o in offs]
        P     = len(uniq)
        F     = max(1, min(fold, P))
        L     = ceil(P / F)                                   # MAC lanes per rail.
        self.n_products, self.fold, self.lanes = P, F, L
        self.dsp48 = 2 * L

        self.sink   = stream.Endpoint(iq_layout(W))
        self.source = stream.Endpoint(iq_layout(W))

        # Runtime-writable unique coefficients; slot (lane, t) holds product lane*F + t.
        self.coeff_we  = Signal()
        self.coeff_adr = Signal(max=max(2, P))
        self.coeff_dat = Signal((CW, True))
        coef = [Signal((CW, True), reset=v, name=f"coef{k}") for k, v in enumerate(uniq)]
        self.sync += If(self.coeff_we, Case(self.coeff_adr, {k: coef[k].eq(self.coeff_dat) for k in range(P)}))

        tree_levels = (L - 1).bit_length()
        # snap, operand select, pre-add, multiply, F accumulate beats, tree, round, saturate.
        latency      = 1 + 1 + 1 + 1 + F + tree_levels + 2
        self.latency = latency

        self.fifo = fifo = stream.SyncFIFO(iq_layout(W), depth=8 + latency)   # Room for every output in flight.
        self.comb += fifo.source.connect(self.source)
        accept = Signal()
        # Launch one output every second accepted sample. A launching sample is held off while the
        # lanes are still busy (until their last slot), so input is limited to 2/F per clock.
        ph, launch = Signal(reset=phase), Signal()
        busy = Signal(); t = Signal(max=max(2, F))
        self.comb += [self.sink.ready.eq((fifo.level < 4) & ~(ph & busy & (t != F - 1))),
                      accept.eq(self.sink.valid & self.sink.ready)]
        self.sync += If(accept, ph.eq(~ph))
        self.comb += launch.eq(accept & ph)

        # Slot sequencer: t counts 0..F-1 after each launch.
        if F > 1:
            self.sync += If(launch, busy.eq(1), t.eq(0)).Elif(busy, If(t == F - 1, busy.eq(0)).Else(t.eq(t + 1)))
        else:
            self.sync += busy.eq(launch)
        # Control pipeline matching the datapath: sel -> preadd -> mult -> accumulate.
        first = Signal(); last = Signal()
        self.comb += [first.eq(busy & (t == 0)), last.eq(busy & (t == F - 1))]
        first_p = Signal(3); last_p = Signal(3); act_p = Signal(3)
        self.sync += [first_p.eq(Cat(first, first_p[:-1])), last_p.eq(Cat(last, last_p[:-1])), act_p.eq(Cat(busy, act_p[:-1]))]
        done = Signal()                       # lane accumulators hold a complete sum.
        self.sync += done.eq(last_p[2])
        done_pipe = Signal(tree_levels + 3)                # hold, tree levels, round, saturate.
        self.sync += done_pipe.eq(Cat(done, done_pipe[:-1]))

        for rail in ("i", "q"):
            x = [Signal((W, True), name=f"x_{rail}_{n}") for n in range(n_taps)]
            self.sync += If(accept, x[0].eq(getattr(self.sink, rail)), *[x[n].eq(x[n - 1]) for n in range(1, n_taps)])
            win  = [getattr(self.sink, rail)] + x[:-1]          # Window incl. the sample arriving now.
            snap = [Signal((W, True), name=f"snap_{rail}_{n}") for n in range(n_taps)]
            self.sync += If(launch, *[snap[n].eq(win[n]) for n in range(n_taps)])

            accs = []
            for lane in range(L):
                a, d   = Signal((W, True), name=f"a_{rail}_{lane}"), Signal((W, True), name=f"d_{rail}_{lane}")
                b      = Signal((CW, True), name=f"b_{rail}_{lane}")
                ad     = Signal((W + 1, True), name=f"ad_{rail}_{lane}")
                b2     = Signal((CW, True), name=f"b2_{rail}_{lane}")
                m      = Signal((W + 1 + CW, True), name=f"m_{rail}_{lane}")
                acc    = Signal((W + 1 + CW + (F - 1).bit_length() + 1, True), name=f"acc_{rail}_{lane}")
                cases_a, cases_d, cases_b = {}, {}, {}
                for s in range(F):
                    k = lane * F + s
                    lo, hi = pairs[k] if k < P else (c, None)
                    cases_a[s] = a.eq(snap[lo])
                    cases_d[s] = d.eq(snap[hi] if hi is not None else 0)
                    cases_b[s] = b.eq(coef[k] if k < P else 0)
                self.sync += [
                    Case(t, cases_a), Case(t, cases_d), Case(t, cases_b),  # operand select (A/D/B regs)
                    ad.eq(a + d), b2.eq(b),                                  # pre-add (AD reg)
                    m.eq(ad * b2),                                           # multiply (M reg)
                    If(first_p[2], acc.eq(m)).Elif(act_p[2], acc.eq(acc + m)),  # accumulate (P reg)
                ]
                accs.append(acc)

            # Registered adder tree over the lanes, sampled once the accumulators are complete.
            hold = []
            for j, s in enumerate(accs):
                h = Signal.like(s, name=f"h_{rail}_{j}")
                self.sync += If(done, h.eq(s))
                hold.append(h)
            level, lv = hold, 0
            while len(level) > 1:
                nxt = []
                for j in range(0, len(level), 2):
                    s = Signal((len(level[0]) + 1, True), name=f"s_{rail}_{lv}_{j}")
                    self.sync += s.eq(level[j] + level[j + 1] if j + 1 < len(level) else level[j])
                    nxt.append(s)
                level, lv = nxt, lv + 1
            total = level[0]
            rnd   = Signal((len(total) - (CW - 1) + 1, True), name=f"r_{rail}")
            out   = Signal((W, True), name=f"y_{rail}")
            self.sync += [rnd.eq(rounded(total, CW - 1)), out.eq(saturated(rnd, W))]
            setattr(self, f"_y_{rail}", out)

        self.comb += [
            fifo.sink.valid.eq(done_pipe[-1]),
            fifo.sink.i.eq(self._y_i),
            fifo.sink.q.eq(self._y_q),
        ]

# Design -------------------------------------------------------------------------------------------

def halfband_kaiser_coefficients(usable=0.70, atten_db=80.0, coeff_width=18, max_taps=799):
    """Shortest 4K+3 Kaiser half-band keeping ``usable`` of the OUTPUT rate (two-sided) with at
    least ``atten_db`` stopband after quantization to ``coeff_width``-bit Q1.(cw-1) integers.

    Returns ``(coefficients, achieved_stopband_db)``. The stopband edge is 0.5 - usable/4 of the
    input rate, i.e. exactly the band that folds onto the kept passband when decimating by 2.
    (70%, 80 dB, 18-bit -> 39 taps, 82.6 dB.)
    """
    import numpy as np
    from scipy.signal import freqz
    from scipy.signal.windows import kaiser
    fp = usable / 4
    for n_taps in range(7, max_taps + 1, 4):
        n = np.arange(n_taps) - (n_taps - 1)//2
        for a in np.arange(atten_db, atten_db + 40, 1.0):
            h = 0.5*np.sinc(n/2.0)*kaiser(n_taps, 0.1102*(a - 8.7))
            h[(n % 2 == 0) & (n != 0)] = 0.0; h[n == 0] = 0.5
            q = [int(v) for v in np.round(h * 2**(coeff_width - 1))]
            w, H = freqz(np.array(q) / 2**(coeff_width - 1), worN=65536, fs=1.0)
            sb = -20*np.log10(np.abs(H[w >= 0.5 - fp]).max() + 1e-15)
            if sb >= atten_db:
                return q, sb
    raise ValueError("no half-band up to max_taps meets the specification")
