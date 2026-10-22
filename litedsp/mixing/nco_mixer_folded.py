#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

"""Folded NCO + complex down-mixer: output = input * conj(NCO), with the multiplies time-shared.

Down-conversion needs four real products per sample:
    I = a.i*c + a.q*s        Q = a.q*c - a.i*s        (c, s = NCO cos/sin)
The fold factor F (1, 2 or 4) runs them on 4/F DSP48 multipliers, each taking F consecutive
products and accumulating them in its P register:

    F = 1 : 4 DSP48, 1 sample per clock
    F = 2 : 2 DSP48 (one builds I, one builds Q), 1 sample every 2 clocks
    F = 4 : 1 DSP48 (I then Q),                    1 sample every 4 clocks

Pick F = floor(f_clk / f_in); e.g. 245.76 MHz with 122.88 MSps in -> F = 2.
The NCO is litedsp's LiteDSPNCO, advanced once per accepted input sample. Arithmetic matches
litedsp's LiteDSPMixer (down mode): full-precision sum, round half-up by data_width-1, saturate.

Measured on xc7a200tsbg484-3 at 245.76 MHz: F=1 4 DSP48 (263 MHz), F=2 2 DSP48 (283 MHz),
F=4 1 DSP48 (285 MHz); litedsp's LiteDSPNCO + LiteDSPMixer reaches 147 MHz.
"""

from migen import *
from litex.gen import *
from litex.soc.interconnect import stream

from litedsp.common import iq_layout, rounded, saturated
from litedsp.generation.nco import LiteDSPNCO

# Product schedule: (rail, input field, NCO field, sign). Rail 0 = I, 1 = Q.
PRODUCTS = [(0, "i", "i", +1), (0, "q", "q", +1), (1, "q", "i", +1), (1, "i", "q", -1)]


class LiteDSPFoldedNCOMixer(LiteXModule):
    def __init__(self, fold=2, data_width=16, phase_bits=32, lut_depth=1024, quarter_wave=True):
        assert fold in (1, 2, 4)
        W, F = data_width, fold
        L = 4 // F                                            # DSP48 lanes.
        self.fold, self.dsp48 = F, L
        self.sink      = stream.Endpoint(iq_layout(W))
        self.source    = stream.Endpoint(iq_layout(W))
        self.phase_inc = Signal(phase_bits)

        self.nco = nco = LiteDSPNCO(phase_bits=phase_bits, data_width=W, lut_depth=lut_depth,
                                    quarter_wave=quarter_wave, with_csr=False)
        self.comb += nco.phase_inc.eq(self.phase_inc)

        self.fifo = fifo = stream.SyncFIFO(iq_layout(W), depth=16)
        self.comb += fifo.source.connect(self.source)

        # Accept one sample every F clocks (and only while the output FIFO has room).
        gap    = Signal(max=max(2, F), reset=0)
        accept = Signal()
        room   = Signal()
        self.comb += [
            room.eq((gap == 0) & (fifo.level < 8)),
            self.sink.ready.eq(room & nco.source.valid),
            nco.source.ready.eq(room & self.sink.valid),
            accept.eq(self.sink.valid & self.sink.ready),
        ]
        if F > 1:
            self.sync += If(accept, gap.eq(F - 1)).Elif(gap != 0, gap.eq(gap - 1))

        # Snapshot the sample and its NCO value; slot counter t walks 0..F-1 after each launch.
        a = {f: Signal((W, True), name=f"a_{f}") for f in ("i", "q")}
        b = {f: Signal((W, True), name=f"b_{f}") for f in ("i", "q")}
        self.sync += If(accept, a["i"].eq(self.sink.i), a["q"].eq(self.sink.q),
                                b["i"].eq(nco.source.i), b["q"].eq(nco.source.q))
        busy = Signal(); t = Signal(max=max(2, F))
        if F > 1:
            self.sync += If(accept, busy.eq(1), t.eq(0)).Elif(busy, If(t == F - 1, busy.eq(0)).Else(t.eq(t + 1)))
        else:
            self.sync += busy.eq(accept)

        PW = 2 * W + 2                                        # Product + accumulation headroom.
        part = [[Signal((PW, True), name=f"part_{j}_{r}") for r in range(2)] for j in range(L)]
        for j in range(L):
            x, y = Signal((W, True), name=f"x_{j}"), Signal((W, True), name=f"y_{j}")
            m    = Signal((2 * W, True), name=f"m_{j}")
            acc  = Signal((PW, True), name=f"acc_{j}")
            # Per-slot control, computed at select time and delayed with the data.
            first, last, neg, rail = Signal(), Signal(), Signal(), Signal()
            sel = {}
            for s in range(F):
                k = j * F + s
                r, fa, fb, sg = PRODUCTS[k]
                is_first = (s == 0) or PRODUCTS[k - 1][0] != r
                is_last  = (s == F - 1) or k == 3 or PRODUCTS[k + 1][0] != r
                sel[s] = [x.eq(a[fa]), y.eq(b[fb]), first.eq(int(is_first)), last.eq(int(is_last) & 1),
                          neg.eq(int(sg < 0)), rail.eq(r)]
            v1 = Signal()
            v2, f2, l2, n2, r2 = Signal(), Signal(), Signal(), Signal(), Signal()
            acc_new = Signal((PW, True), name=f"accn_{j}")
            self.sync += [
                If(busy, Case(t, sel)), v1.eq(busy),                      # select (A/B regs)
                m.eq(x * y), v2.eq(v1), f2.eq(first), l2.eq(last), n2.eq(neg), r2.eq(rail),  # multiply (M reg); flags ride along
            ]
            self.comb += acc_new.eq(Mux(f2, 0, acc) + Mux(n2, -m, m))
            self.sync += If(v2, acc.eq(acc_new),                          # accumulate (P reg)
                            If(l2, If(r2, part[j][1].eq(acc_new)).Else(part[j][0].eq(acc_new))))

        # All partials are final F+3 clocks after a launch; then sum lanes per rail, round, saturate.
        done = Signal(F + 4)
        self.sync += done.eq(Cat(accept, done[:-1]))
        res = {}
        for r, name in ((0, "i"), (1, "q")):
            lanes = [j for j in range(L) if any(PRODUCTS[j * F + s][0] == r for s in range(F))]
            tot = Signal((PW + 1, True), name=f"tot_{name}")
            rnd = Signal((PW - W + 3, True), name=f"rnd_{name}")
            out = Signal((W, True), name=f"out_{name}")
            self.sync += [
                If(done[F + 2], tot.eq(sum(part[j][r] for j in lanes))),
                rnd.eq(rounded(tot, W - 1)),
                out.eq(saturated(rnd, W)),
            ]
            res[name] = out
        valid = Signal(2)
        self.sync += valid.eq(Cat(done[F + 2], valid[:-1]))
        out_valid = Signal()
        self.sync += out_valid.eq(valid[-1])
        self.comb += [fifo.sink.valid.eq(out_valid), fifo.sink.i.eq(res["i"]), fifo.sink.q.eq(res["q"])]
