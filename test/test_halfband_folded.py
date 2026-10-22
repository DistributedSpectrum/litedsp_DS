#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

import unittest

import numpy as np

from litedsp.filter.halfband_folded import LiteDSPHalfbandDecimatorFolded, halfband_kaiser_coefficients

from test.common import run_stream, column
from test.models import fir_decimator_model

COEFFS, STOPBAND = halfband_kaiser_coefficients(usable=0.70, atten_db=80.0, coeff_width=18)

class TestHalfbandFolded(unittest.TestCase):
    def test_design(self):
        self.assertEqual(len(COEFFS), 39)
        self.assertGreaterEqual(STOPBAND, 80.0)

    def test_bit_exact_every_fold(self):
        rng = np.random.RandomState(7)
        x_i = rng.randint(-32768, 32768, 160)
        x_q = rng.randint(-32768, 32768, 160)
        samples = [{"i": int(i), "q": int(q)} for i, q in zip(x_i, x_q)]
        for fold, dsp48 in ((1, 22), (2, 12), (4, 6), (8, 4), (11, 2)):
            with self.subTest(fold=fold):
                dut = LiteDSPHalfbandDecimatorFolded(COEFFS, fold=fold, coeff_width=18)
                self.assertEqual(dut.dsp48, dsp48)
                cap = run_stream(dut, samples, len(samples)//2, ["i", "q"], ["i", "q"],
                    sink_throttle=0.2, source_ready_rate=0.6)
                np.testing.assert_array_equal(column(cap, "i", 16), fir_decimator_model(x_i, COEFFS, 2, shift=17))
                np.testing.assert_array_equal(column(cap, "q", 16), fir_decimator_model(x_q, COEFFS, 2, shift=17))

    def test_m2sdr_front_halfband(self):
        # 122.88 -> 61.44 MSps keeping 50 MHz: 59 taps, 8 DSP48 at fold 4 (245.76 MHz clock).
        coeffs, stopband = halfband_kaiser_coefficients(usable=50/61.44, atten_db=80.0, coeff_width=18)
        self.assertEqual(len(coeffs), 59)
        self.assertGreaterEqual(stopband, 80.0)
        rng = np.random.RandomState(3)
        x_i = rng.randint(-32768, 32768, 160)
        x_q = rng.randint(-32768, 32768, 160)
        dut = LiteDSPHalfbandDecimatorFolded(coeffs, fold=4, coeff_width=18)
        self.assertEqual(dut.dsp48, 8)
        cap = run_stream(dut, [{"i": int(i), "q": int(q)} for i, q in zip(x_i, x_q)], 80, ["i", "q"], ["i", "q"],
            sink_throttle=0.2, source_ready_rate=0.6)
        np.testing.assert_array_equal(column(cap, "i", 16), fir_decimator_model(x_i, coeffs, 2, shift=17))
        np.testing.assert_array_equal(column(cap, "q", 16), fir_decimator_model(x_q, coeffs, 2, shift=17))

if __name__ == "__main__":
    unittest.main()
