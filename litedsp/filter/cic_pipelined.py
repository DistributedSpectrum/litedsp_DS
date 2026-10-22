#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

"""Pipelined CIC decimator with a runtime rate, for full-rate (1 sample/clock) input.

Hogenauer structure with every stage registered: N integrators at the input rate (each adds the
previous stage's registered value, one adder per clock), a rate counter set at runtime, N comb
stages as a registered pipeline at the output rate, then rounding and the runtime rescale shift
in separate registered steps. Datapath sized for r_max: growth = ceil(N * log2(r_max * M)) bits.
The host sets ``rate`` and ``shift = ceil(N*log2(rate*M))`` together (same contract as litedsp's
CICDecimatorRuntime). Output = saturate(round_half_up(comb_out >> shift)).

Measured on xc7a200tsbg484-3 at 245.76 MHz, r_max=512: N=4 1434 LUT / 1576 FF / 0 DSP (288 MHz),
N=5 1921 LUT / 2205 FF (272 MHz). (The unpipelined CICDecimatorRuntime reaches about 60 MHz.)
"""

from math import ceil, log2

from migen import *
from litex.gen import *
from litex.soc.interconnect import stream

from litedsp.common import iq_layout, saturated


class LiteDSPCICDecimatorPipelined(LiteXModule):
    def __init__(self, data_width=16, r_max=512, N=4, M=1):
        W      = data_width
        growth = ceil(N * log2(r_max * M))
        Wi     = W + growth
        self.growth = growth
        self.sink   = stream.Endpoint(iq_layout(W))
        self.source = stream.Endpoint(iq_layout(W))
        self.rate   = Signal(max=r_max + 1, reset=r_max)
        self.shift  = Signal(max=growth + 1, reset=growth)

        self.fifo = fifo = stream.SyncFIFO(iq_layout(W), depth=16)
        self.comb += fifo.source.connect(self.source)
        accept = Signal()
        self.comb += [self.sink.ready.eq(fifo.level < 8), accept.eq(self.sink.valid & self.sink.ready)]

        # Rate counter: an output strobe every `rate` accepted samples.
        cnt  = Signal(max=r_max + 1)
        take = Signal()
        self.sync += If(accept, If(cnt >= self.rate - 1, cnt.eq(0)).Else(cnt.eq(cnt + 1)))
        self.comb += take.eq(accept & (cnt >= self.rate - 1))

        # Integrator pipeline: stage k adds stage k-1's registered value, so the strobe that
        # samples the last stage is delayed by N-1 accepted samples relative to `take`.
        take_d = Signal(N)
        self.sync += If(accept, take_d.eq(Cat(take, take_d[:-1])))
        strobe = Signal()
        self.comb += strobe.eq(accept & (take_d[N - 2] if N > 1 else take))

        # Rescale constants, registered so the shift is off the data path.
        bias = Signal(Wi)
        self.sync += bias.eq(Mux(self.shift == 0, 0, (1 << self.shift) >> 1))

        n_out = N + 3                       # comb stages + bias add + shift + saturate
        vpipe = Signal(n_out)
        self.sync += vpipe.eq(Cat(strobe, vpipe[:-1]))

        for rail in ("i", "q"):
            x = getattr(self.sink, rail)
            integ = [Signal((Wi, True), name=f"int_{rail}_{k}") for k in range(N)]
            upd = [integ[0].eq(integ[0] + x)] + [integ[k].eq(integ[k] + integ[k - 1]) for k in range(1, N)]
            self.sync += If(accept, *upd)

            # Comb pipeline at the output rate (each stage one subtract per clock).
            stage = integ[-1]
            for k in range(N):
                prevs = [Signal((Wi, True), name=f"cprev_{rail}_{k}_{m}") for m in range(M)]
                out   = Signal((Wi, True), name=f"comb_{rail}_{k}")
                vk    = vpipe[k - 1] if k > 0 else strobe
                self.sync += If(vk, out.eq(stage - prevs[-1]), prevs[0].eq(stage),
                                *[prevs[m].eq(prevs[m - 1]) for m in range(1, M)])
                stage = out
            biased  = Signal((Wi + 1, True), name=f"biased_{rail}")
            shifted = Signal((Wi + 1, True), name=f"shifted_{rail}")
            res     = Signal((W, True), name=f"res_{rail}")
            self.sync += [biased.eq(stage + bias), shifted.eq(biased >> self.shift), res.eq(saturated(shifted, W))]
            setattr(self, f"_res_{rail}", res)

        self.comb += [fifo.sink.valid.eq(vpipe[-1]), fifo.sink.i.eq(self._res_i), fifo.sink.q.eq(self._res_q)]
