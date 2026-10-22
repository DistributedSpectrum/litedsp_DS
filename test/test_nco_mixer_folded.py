#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

import unittest

import numpy as np
from migen import *
from litex.gen import LiteXModule

from litedsp.generation.nco import LiteDSPNCO
from litedsp.mixing.mixer import LiteDSPMixer
from litedsp.mixing.nco_mixer_folded import LiteDSPFoldedNCOMixer

from test.common import run_stream, column

class _Reference(LiteXModule):
    """litedsp's own NCO + down-mixer, the arithmetic the folded block must reproduce."""
    def __init__(self, phase_inc):
        self.nco = LiteDSPNCO(quarter_wave=True, with_csr=False)
        self.mix = LiteDSPMixer(with_csr=False)
        self.comb += [self.nco.source.connect(self.mix.sink_b), self.mix.mode.eq(0), self.nco.phase_inc.eq(phase_inc)]
        self.sink, self.source = self.mix.sink_a, self.mix.source

class TestFoldedNCOMixer(unittest.TestCase):
    def test_bit_exact_vs_litedsp(self):
        rng = np.random.RandomState(3)
        x_i = rng.randint(-32768, 32768, 200); x_i[::17] = -32768
        x_q = rng.randint(-32768, 32768, 200)
        samples = [{"i": int(i), "q": int(q)} for i, q in zip(x_i, x_q)]
        for pinc in (0x12345678, 0x80000001):
            ref = run_stream(_Reference(pinc), samples, len(samples), ["i", "q"], ["i", "q"], sink_throttle=0.0, source_ready_rate=1.0)
            for fold, dsp48 in ((1, 4), (2, 2), (4, 1)):
                with self.subTest(phase_inc=hex(pinc), fold=fold):
                    dut = LiteDSPFoldedNCOMixer(fold=fold)
                    self.assertEqual(dut.dsp48, dsp48)
                    dut.comb += dut.phase_inc.eq(pinc)
                    cap = run_stream(dut, samples, len(samples), ["i", "q"], ["i", "q"], sink_throttle=0.2, source_ready_rate=0.6)
                    np.testing.assert_array_equal(column(cap, "i", 16), column(ref, "i", 16))
                    np.testing.assert_array_equal(column(cap, "q", 16), column(ref, "q", 16))

if __name__ == "__main__":
    unittest.main()
