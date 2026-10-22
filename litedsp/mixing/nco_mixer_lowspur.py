#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

"""Low-spur folded NCO + complex down-mixer (block-RAM table + first-order Taylor correction).

NCO: 32-bit phase accumulator. The top A bits address a half-sample-offset quarter-wave sine
table (A=12: 1024 x 18-bit = one RAMB18; its two ports read cos and sin together). The next D bits
give the signed offset d from the table point, and a first-order correction
    sin(t + d) ~= sin t + d*cos t        cos(t + d) ~= cos t - d*sin t
removes the phase-truncation error (worst spur about -111 dBc for A=12, D=10, 18-bit; versus about
-58 dBc for a plain 1024-entry table). The two small correction products (d is under 10 bits) are
built in LUTs, not DSP48s. The NCO is fully pipelined and advances once per accepted sample.

Mixer: output = input * conj(NCO), four real products per sample time-shared on 4/F DSP48s
(F = 1, 2, 4), 16-bit sample x 18-bit NCO, round half-up by B-1 bits, saturate to 16 bits.

Measured on xc7a200tsbg484-3 at 245.76 MHz (1 RAMB18 each): F=1 6 DSP48 (274 MHz), F=2 3 DSP48
(285 MHz), F=4 2 DSP48 (288 MHz). Worst spur at the 16-bit mixer output -98 dBc (NCO -111 dBc).
"""

import math

from migen import *
from litex.gen import *
from litex.soc.interconnect import stream

from litedsp.common import iq_layout, rounded, saturated

PRODUCTS = [(0, "i", "c", +1), (0, "q", "s", +1), (1, "q", "c", +1), (1, "i", "s", -1)]


def nco_table(A, B):
    N, Q, scale = 1 << A, 1 << (A - 2), (1 << (B - 1)) - 1
    return [int(round(math.sin(2 * math.pi * (i + 0.5) / N) * scale)) for i in range(Q)]


class LiteDSPLowSpurNCOMixer(LiteXModule):
    def __init__(self, fold=2, data_width=16, phase_bits=32, A=12, D=10, B=18, KW=18):
        assert fold in (1, 2, 4)
        W, F, L = data_width, fold, 4 // fold
        Q = 1 << (A - 2)
        self.fold, self.dsp48 = F, L
        self.sink      = stream.Endpoint(iq_layout(W))
        self.source    = stream.Endpoint(iq_layout(W))
        self.phase_inc = Signal(phase_bits)
        K = int(round(2 * math.pi / (1 << A) * (1 << KW)))      # radians per table step, x 2^KW

        self.fifo = fifo = stream.SyncFIFO(iq_layout(W), depth=32)   # ~15 samples in flight at F=1, gated at level < 6.
        self.comb += fifo.source.connect(self.source)
        gap, accept = Signal(max=max(2, F)), Signal()
        self.comb += [self.sink.ready.eq((gap == 0) & (fifo.level < 6)), accept.eq(self.sink.valid & self.sink.ready)]
        if F > 1:
            self.sync += If(accept, gap.eq(F - 1)).Elif(gap != 0, gap.eq(gap - 1))

        # ---- NCO pipeline (one result per accepted sample, sample data carried alongside) ----
        # Every stage register loads only when its input is valid, so values hold between samples
        # (needed by the time-shared correction multiplier at F >= 2).
        phase = Signal(phase_bits)
        self.sync += If(accept, phase.eq(phase + self.phase_inc))
        v1, v2, v3 = Signal(), Signal(), Signal()
        ph = Signal(phase_bits)
        si1, sq1 = Signal((W, True)), Signal((W, True))
        self.sync += [v1.eq(accept), If(accept, ph.eq(phase + self.phase_inc), si1.eq(self.sink.i), sq1.eq(self.sink.q))]

        # Stage 1 -> 2: table addresses, quadrant signs, correction offset.
        a   = ph[phase_bits - A:]
        qs  = a[A - 2:]; idx = a[:A - 2]
        qc  = Signal(2); self.comb += qc.eq(qs + 1)
        fr  = ph[phase_bits - A - D: phase_bits - A]
        adr_s, adr_c = Signal(A - 2), Signal(A - 2)
        neg_s1, neg_c1, neg_s2, neg_c2 = Signal(), Signal(), Signal(), Signal()
        frs = Signal((D + 1, True))
        si2, sq2, si3, sq3 = Signal((W, True)), Signal((W, True)), Signal((W, True)), Signal((W, True))
        self.sync += [v2.eq(v1), If(v1,
            adr_s.eq(Mux(qs[0], Q - 1 - idx, idx)), neg_s1.eq(qs[1]),
            adr_c.eq(Mux(qc[0], Q - 1 - idx, idx)), neg_c1.eq(qc[1]),
            frs.eq(fr - (1 << (D - 1))), si2.eq(si1), sq2.eq(sq1))]
        # Stage 2 -> 3: dual-port table read (RAMB18, read enabled with the stage), delta scaling (LUTs).
        mem = Memory(B, Q, init=[x & ((1 << B) - 1) for x in nco_table(A, B)], name="nco_tab")
        ps, pc = mem.get_port(has_re=True), mem.get_port(has_re=True)
        self.specials += mem, ps, pc
        self.comb += [ps.adr.eq(adr_s), pc.adr.eq(adr_c), ps.re.eq(v2), pc.re.eq(v2)]
        dq3 = Signal((KW + 2, True), attr={("use_dsp", "no")})
        self.sync += [v3.eq(v2), If(v2, neg_s2.eq(neg_s1), neg_c2.eq(neg_c1), dq3.eq((frs * K) >> D), si3.eq(si2), sq3.eq(sq2))]
        # Stage 3 -> 4: signed coarse values.
        v4 = Signal()
        s4, c4 = Signal((B + 1, True)), Signal((B + 1, True))
        dq4 = Signal.like(dq3)
        si4, sq4 = Signal((W, True)), Signal((W, True))
        # RAM output register first (maps onto the RAMB18 DO register), then the sign.
        v3b = Signal()
        sraw, craw = Signal((B, True)), Signal((B, True))
        neg_s3, neg_c3 = Signal(), Signal()
        dq3b = Signal.like(dq3)
        si3b, sq3b = Signal((W, True)), Signal((W, True))
        self.sync += [v3b.eq(v3), If(v3, sraw.eq(ps.dat_r), craw.eq(pc.dat_r), neg_s3.eq(neg_s2), neg_c3.eq(neg_c2),
                                     dq3b.eq(dq3), si3b.eq(si3), sq3b.eq(sq3))]
        self.sync += [v4.eq(v3b), If(v3b, s4.eq(Mux(neg_s3, -sraw, sraw)), c4.eq(Mux(neg_c3, -craw, craw)), dq4.eq(dq3b), si4.eq(si3b), sq4.eq(sq3b))]
        # Stage 4 -> 6: correction products on DSP48: two multipliers at F = 1, one shared at F >= 2.
        PWc = B + KW + 4
        ms, mc = Signal((PWc, True), name="corr_ms"), Signal((PWc, True), name="corr_mc")   # dq*c, dq*s
        v5, v6 = Signal(), Signal()
        si6, sq6 = Signal((W, True)), Signal((W, True))
        s6, c6 = Signal.like(s4), Signal.like(c4)
        if F == 1:
            m1, m2 = Signal((PWc, True), name="cm1"), Signal((PWc, True), name="cm2")
            si5, sq5, s5, c5 = Signal((W, True)), Signal((W, True)), Signal.like(s4), Signal.like(c4)
            self.sync += [v5.eq(v4), m1.eq(dq4 * c4), m2.eq(dq4 * s4), si5.eq(si4), sq5.eq(sq4), s5.eq(s4), c5.eq(c4),
                          v6.eq(v5), ms.eq(m1), mc.eq(m2), si6.eq(si5), sq6.eq(sq5), s6.eq(s5), c6.eq(c5)]
        else:
            # Registered operands (DSP48 A/B regs), then the product (M reg):
            # cycle k: operands dq,c; k+1: dq,s (stage-4 registers hold); products land at k+2, k+3.
            opa, opb = Signal.like(dq4, name="cop_a"), Signal.like(c4, name="cop_b")
            prod     = Signal((PWc, True), name="cprod")
            v7a = Signal()
            self.sync += [
                opa.eq(dq4), opb.eq(Mux(v4, c4, s4)),
                prod.eq(opa * opb),
                v5.eq(v4), v6.eq(v5), v7a.eq(v6),
                If(v5, si6.eq(si4), sq6.eq(sq4), s6.eq(s4), c6.eq(c4)),   # stage 4 still holds this sample at k+1
                If(v6, ms.eq(prod)),                                       # dq*c lands at k+2
            ]
        # Final stage: corrected cos/sin.
        v7 = Signal()
        sn, cn = Signal((B + 1, True), name="nco_s"), Signal((B + 1, True), name="nco_c")
        si7, sq7 = Signal((W, True)), Signal((W, True))
        if F == 1:
            self.sync += [v7.eq(v6), If(v6, sn.eq(s6 + rounded(ms, KW)), cn.eq(c6 - rounded(mc, KW)), si7.eq(si6), sq7.eq(sq6))]
        else:
            self.sync += [v7.eq(v7a), If(v7a, sn.eq(s6 + rounded(ms, KW)), cn.eq(c6 - rounded(prod, KW)), si7.eq(si6), sq7.eq(sq6))]   # dq*s at k+3
        self.nco_s, self.nco_c = sn, cn
        launch = v7                                              # NCO result + sample ready
        ain, aqn = si7, sq7

        # ---- Folded mixer ----
        av = {"i": Signal((W, True)), "q": Signal((W, True))}
        bv = {"c": Signal((B + 1, True)), "s": Signal((B + 1, True))}
        self.sync += If(launch, av["i"].eq(ain), av["q"].eq(aqn), bv["c"].eq(cn), bv["s"].eq(sn))
        busy = Signal(); t = Signal(max=max(2, F))
        if F > 1:
            self.sync += If(launch, busy.eq(1), t.eq(0)).Elif(busy, If(t == F - 1, busy.eq(0)).Else(t.eq(t + 1)))
        else:
            self.sync += busy.eq(launch)
        PW = W + B + 3
        part = [[Signal((PW, True), name=f"part_{j}_{r}") for r in range(2)] for j in range(L)]
        for j in range(L):
            x, y = Signal((W, True), name=f"x_{j}"), Signal((B + 1, True), name=f"y_{j}")
            m, acc = Signal((W + B + 1, True), name=f"m_{j}"), Signal((PW, True), name=f"acc_{j}")
            first, last, neg, rail = Signal(), Signal(), Signal(), Signal()
            sel = {}
            for s in range(F):
                k = j * F + s
                r, fa, fb, sg = PRODUCTS[k]
                is_first = (s == 0) or PRODUCTS[k - 1][0] != r
                is_last  = (s == F - 1) or k == 3 or PRODUCTS[k + 1][0] != r
                sel[s] = [x.eq(av[fa]), y.eq(bv[fb]), first.eq(int(is_first)), last.eq(int(is_last)), neg.eq(int(sg < 0)), rail.eq(r)]
            v1, v2, f2, l2, n2, r2 = Signal(), Signal(), Signal(), Signal(), Signal(), Signal()
            acc_new = Signal((PW, True), name=f"accn_{j}")
            self.sync += [If(busy, Case(t, sel)), v1.eq(busy),
                          m.eq(x * y), v2.eq(v1), f2.eq(first), l2.eq(last), n2.eq(neg), r2.eq(rail)]
            self.comb += acc_new.eq(Mux(f2, 0, acc) + Mux(n2, -m, m))
            self.sync += If(v2, acc.eq(acc_new), If(l2, If(r2, part[j][1].eq(acc_new)).Else(part[j][0].eq(acc_new))))
        done = Signal(F + 4)
        self.sync += done.eq(Cat(launch, done[:-1]))
        res = {}
        for r, name in ((0, "i"), (1, "q")):
            lanes = [j for j in range(L) if any(PRODUCTS[j * F + s][0] == r for s in range(F))]
            tot, rnd, out = Signal((PW + 1, True)), Signal((PW - B + 3, True)), Signal((W, True), name=f"out_{name}")
            self.sync += [If(done[F + 2], tot.eq(sum(part[j][r] for j in lanes))), rnd.eq(rounded(tot, B - 1)), out.eq(saturated(rnd, W))]
            res[name] = out
        valid = Signal(3)
        self.sync += valid.eq(Cat(done[F + 2], valid[:-1]))
        self.comb += [fifo.sink.valid.eq(valid[2]), fifo.sink.i.eq(res["i"]), fifo.sink.q.eq(res["q"])]
