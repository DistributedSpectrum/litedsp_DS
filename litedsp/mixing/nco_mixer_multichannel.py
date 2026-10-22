#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

"""Multi-channel low-spur NCO + complex down-mixer: C channels of one stream share one table.

For an input with F clocks per sample (F >= 4; 61.44 MSps at 245.76 MHz -> F = 4), up to C = F
channels share one NCO datapath:

- One block-RAM quarter-wave table (1024 x 18, one RAMB18). The NCO pipeline is time-multiplexed:
  channel k is looked up in clock slot k after each accepted sample, the RAM's two read ports
  giving its sin and cos together.
- The first-order Taylor correction, built in LUTs. The correction offset d fits in 9 bits, and
  cos/sin are cut to ``corr_bits`` (10) bits before multiplying, which keeps the worst NCO spur
  near -110 dBc (full-width correction: -109 dBc). Each product is split in two halves so it
  closes timing.
- Each channel then has its own mixer: one DSP48 doing the four products over four clocks.

A group of 4 channels at 61.44 MSps therefore costs 4 DSP48 and one RAMB18. Each channel has its
own ``phase_inc`` and output stream. All channels advance together on each accepted sample, and
input pauses while any channel's output FIFO is filling (they share back-pressure).

Arithmetic per channel matches LiteDSPLowSpurNCOMixer except for the correction products, which
are round(d * (c >> s) / 2^(KW - s)) with s = B + 1 - corr_bits.

Measured on xc7a200tsbg484-3 at 245.76 MHz, 4 channels, F = 4: 4 DSP48, 1 RAMB18, 849 LUT
(109 as FIFO memory), 1556 FF, 282 MHz. Latency 19 + k clocks for channel k.
"""

import math
from functools import reduce
from operator import and_

from migen import *
from litex.gen import *
from litex.soc.interconnect import stream

from litedsp.common import iq_layout, rounded, saturated
from litedsp.mixing.nco_mixer_lowspur import PRODUCTS, nco_table


class LiteDSPMultiChannelNCOMixer(LiteXModule):
    def __init__(self, channels=4, fold=4, data_width=16, phase_bits=32, A=12, D=10, B=18, KW=18, corr_bits=10):
        C, F, W = channels, fold, data_width
        assert F >= 4, "each channel's mixer needs four clocks per sample"
        assert 1 <= C <= F, "the table ports give one lookup per clock: at most F channels"
        Q  = 1 << (A - 2)
        K  = int(round(2 * math.pi / (1 << A) * (1 << KW)))   # radians per table step, x 2^KW
        DW = ((K + 1) // 2).bit_length() + 1                  # |d| <= K/2
        SH = B + 1 - corr_bits                                # cos/sin bits dropped before the correction
        RS = KW - SH                                          # rounding shift of the narrowed product
        LB = corr_bits // 2                                   # split of the narrowed operand
        assert RS >= 1
        self.channels, self.fold, self.dsp48 = C, F, C
        self.sink       = stream.Endpoint(iq_layout(W))
        self.sources    = [stream.Endpoint(iq_layout(W)) for _ in range(C)]
        self.phase_incs = [Signal(phase_bits, name=f"phase_inc{k}") for k in range(C)]

        # Per-channel output FIFOs. A sample takes 19 + k clocks to reach channel k's output, so at
        # most ceil((19 + F - 1) / F) are in flight per channel (6 at F = 4): gating input at
        # level < 8 leaves room in 16.
        fifos = [stream.SyncFIFO(iq_layout(W), depth=16) for _ in range(C)]
        self.submodules += fifos
        for fifo, source in zip(fifos, self.sources):
            self.comb += fifo.source.connect(source)
        gap, accept = Signal(max=F), Signal()
        self.comb += [
            self.sink.ready.eq((gap == 0) & reduce(and_, [fifo.level < 8 for fifo in fifos])),
            accept.eq(self.sink.valid & self.sink.ready),
        ]
        self.sync += If(accept, gap.eq(F - 1)).Elif(gap != 0, gap.eq(gap - 1))

        # ---- Phase accumulators: every channel steps on each accepted sample ----
        phs = []
        for k in range(C):
            phase, ph = Signal(phase_bits, name=f"phase{k}"), Signal(phase_bits, name=f"ph{k}")
            self.sync += If(accept, phase.eq(phase + self.phase_incs[k]), ph.eq(phase + self.phase_incs[k]))
            phs.append(ph)
        xi0, xq0 = Signal((W, True)), Signal((W, True))
        self.sync += If(accept, xi0.eq(self.sink.i), xq0.eq(self.sink.q))

        # ---- Time-multiplexed NCO pipeline: channel `slot` enters one clock after another ----
        # The sample and the lookup phases hold until the next accept (at least F >= C clocks later).
        issuing, slot = Signal(), Signal(max=max(2, C))
        if C > 1:
            self.sync += If(accept, issuing.eq(1), slot.eq(0)).Elif(issuing,
                If(slot == C - 1, issuing.eq(0)).Else(slot.eq(slot + 1)))
        else:
            self.sync += issuing.eq(accept)
        ph_sel = Signal(phase_bits)
        self.comb += Case(slot, {k: ph_sel.eq(phs[k]) for k in range(C)})

        def stage(n, **sigs):
            # Registers carrying the slot's valid, channel and sample to stage n.
            return {name: Signal.like(s, name=f"{name}{n}") for name, s in sigs.items()}

        # Stage 1: table addresses, quadrant signs, correction offset.
        a   = ph_sel[phase_bits - A:]
        qs  = a[A - 2:]; idx = a[:A - 2]
        qc  = Signal(2); self.comb += qc.eq(qs + 1)
        fr  = ph_sel[phase_bits - A - D: phase_bits - A]
        adr_s, adr_c = Signal(A - 2), Signal(A - 2)
        neg_s1, neg_c1, frs1 = Signal(), Signal(), Signal((D + 1, True))
        p1 = stage(1, v=issuing, k=slot, xi=xi0, xq=xq0)
        self.sync += [p1["v"].eq(issuing), p1["k"].eq(slot), p1["xi"].eq(xi0), p1["xq"].eq(xq0),
            adr_s.eq(Mux(qs[0], Q - 1 - idx, idx)), neg_s1.eq(qs[1]),
            adr_c.eq(Mux(qc[0], Q - 1 - idx, idx)), neg_c1.eq(qc[1]),
            frs1.eq(fr - (1 << (D - 1)))]

        # Stage 2: dual-port table read (RAMB18), offset scaled to radians in LUTs.
        mem = Memory(B, Q, init=[x & ((1 << B) - 1) for x in nco_table(A, B)], name="nco_tab")
        ps, pc = mem.get_port(), mem.get_port()
        self.specials += mem, ps, pc
        self.comb += [ps.adr.eq(adr_s), pc.adr.eq(adr_c)]
        dq2 = Signal((DW, True), attr={("use_dsp", "no")})
        neg_s2, neg_c2 = Signal(), Signal()
        p2 = stage(2, **p1)
        self.sync += [*[p2[n].eq(p1[n]) for n in p1], neg_s2.eq(neg_s1), neg_c2.eq(neg_c1), dq2.eq((frs1 * K) >> D)]

        # Stage 3: RAM output register (maps onto the RAMB18 DO register).
        sraw, craw = Signal((B, True)), Signal((B, True))
        neg_s3, neg_c3, dq3 = Signal(), Signal(), Signal.like(dq2)
        p3 = stage(3, **p2)
        self.sync += [*[p3[n].eq(p2[n]) for n in p2], sraw.eq(ps.dat_r), craw.eq(pc.dat_r),
            neg_s3.eq(neg_s2), neg_c3.eq(neg_c2), dq3.eq(dq2)]

        # Stage 4: signed coarse values.
        s4, c4, dq4 = Signal((B + 1, True)), Signal((B + 1, True)), Signal.like(dq2)
        p4 = stage(4, **p3)
        self.sync += [*[p4[n].eq(p3[n]) for n in p3],
            s4.eq(Mux(neg_s3, -sraw, sraw)), c4.eq(Mux(neg_c3, -craw, craw)), dq4.eq(dq3)]

        # Stages 5-6: correction products d * (c >> SH) and d * (s >> SH) in LUTs, each as two
        # registered half products (low half unsigned, high half signed), summed the clock after.
        p5, p6 = stage(5, **p4), stage(6, **p4)
        s5, c5, s6, c6 = Signal.like(s4), Signal.like(c4), Signal.like(s4), Signal.like(c4)
        self.sync += [*[p5[n].eq(p4[n]) for n in p4], s5.eq(s4), c5.eq(c4),
                      *[p6[n].eq(p5[n]) for n in p5], s6.eq(s5), c6.eq(c5)]
        corr = {}
        for name, val in (("c", c4), ("s", s4)):
            t  = Signal((corr_bits, True), name=f"corr_{name}_t")
            lo = Signal((LB + 1, True), name=f"corr_{name}_lo")
            hi = Signal((corr_bits - LB, True), name=f"corr_{name}_hi")
            self.comb += [t.eq(val >> SH), lo.eq(t[:LB]), hi.eq(t[LB:])]
            plo = Signal((DW + LB + 1, True), name=f"corr_{name}_plo", attr={("use_dsp", "no")})
            phi = Signal((DW + corr_bits - LB, True), name=f"corr_{name}_phi", attr={("use_dsp", "no")})
            prod = Signal((DW + corr_bits + 1, True), name=f"corr_{name}_prod")
            self.sync += [plo.eq(dq4 * lo), phi.eq(dq4 * hi), prod.eq((phi << LB) + plo)]
            corr[name] = prod

        # Stage 7: corrected sin/cos, ready for the channel's mixer.
        p7 = stage(7, **p6)
        sn, cn = Signal((B + 1, True), name="nco_s"), Signal((B + 1, True), name="nco_c")
        self.sync += [*[p7[n].eq(p6[n]) for n in p6],
            sn.eq(s6 + rounded(corr["c"], RS)), cn.eq(c6 - rounded(corr["s"], RS))]

        # ---- Per-channel mixers: one DSP48 each, four products over four clocks ----
        FM, PW = 4, W + B + 3
        for k, fifo in enumerate(fifos):
            launch = Signal(name=f"launch{k}")
            self.comb += launch.eq(p7["v"] & (p7["k"] == k))
            av = {"i": Signal((W, True), name=f"a{k}_i"), "q": Signal((W, True), name=f"a{k}_q")}
            bv = {"c": Signal((B + 1, True), name=f"b{k}_c"), "s": Signal((B + 1, True), name=f"b{k}_s")}
            self.sync += If(launch, av["i"].eq(p7["xi"]), av["q"].eq(p7["xq"]), bv["c"].eq(cn), bv["s"].eq(sn))
            busy, t = Signal(name=f"busy{k}"), Signal(2, name=f"t{k}")
            self.sync += If(launch, busy.eq(1), t.eq(0)).Elif(busy, If(t == FM - 1, busy.eq(0)).Else(t.eq(t + 1)))
            part = [Signal((PW, True), name=f"part{k}_{r}") for r in range(2)]
            x, y = Signal((W, True), name=f"x{k}"), Signal((B + 1, True), name=f"y{k}")
            m, acc = Signal((W + B + 1, True), name=f"m{k}"), Signal((PW, True), name=f"acc{k}")
            first, last, neg, rail = Signal(), Signal(), Signal(), Signal()
            sel = {}
            for s in range(FM):
                r, fa, fb, sg = PRODUCTS[s]
                is_first = (s == 0) or PRODUCTS[s - 1][0] != r
                is_last  = (s == FM - 1) or PRODUCTS[s + 1][0] != r
                sel[s] = [x.eq(av[fa]), y.eq(bv[fb]), first.eq(int(is_first)), last.eq(int(is_last)),
                          neg.eq(int(sg < 0)), rail.eq(r)]
            v1, v2, f2, l2, n2, r2 = Signal(), Signal(), Signal(), Signal(), Signal(), Signal()
            acc_new = Signal((PW, True), name=f"accn{k}")
            self.sync += [If(busy, Case(t, sel)), v1.eq(busy),
                          m.eq(x * y), v2.eq(v1), f2.eq(first), l2.eq(last), n2.eq(neg), r2.eq(rail)]
            self.comb += acc_new.eq(Mux(f2, 0, acc) + Mux(n2, -m, m))
            self.sync += If(v2, acc.eq(acc_new), If(l2, If(r2, part[1].eq(acc_new)).Else(part[0].eq(acc_new))))
            done = Signal(FM + 4, name=f"done{k}")
            self.sync += done.eq(Cat(launch, done[:-1]))
            res = {}
            for r, name in ((0, "i"), (1, "q")):
                tot = Signal((PW + 1, True), name=f"tot{k}_{name}")
                rnd = Signal((PW - B + 3, True), name=f"rnd{k}_{name}")
                out = Signal((W, True), name=f"out{k}_{name}")
                self.sync += [If(done[FM + 2], tot.eq(part[r])), rnd.eq(rounded(tot, B - 1)), out.eq(saturated(rnd, W))]
                res[name] = out
            valid = Signal(3, name=f"valid{k}")
            self.sync += valid.eq(Cat(done[FM + 2], valid[:-1]))
            self.comb += [fifo.sink.valid.eq(valid[2]), fifo.sink.i.eq(res["i"]), fifo.sink.q.eq(res["q"])]
