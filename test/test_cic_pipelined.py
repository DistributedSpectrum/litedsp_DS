#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

import unittest
from math import ceil, log2

import numpy as np

from litedsp.filter.cic_pipelined import LiteDSPCICDecimatorPipelined

from test.common import run_stream, column

def cic_model(x, R, N, shift):
    """Ideal CIC: N-fold length-R boxcar, decimated, round half-up by `shift`, saturated to 16 bits."""
    k = np.ones(R, dtype=object)
    for _ in range(N - 1):
        k = np.convolve(k, np.ones(R, dtype=object))
    y = np.convolve(np.array(x, dtype=object), k)
    bias = (1 << shift) >> 1 if shift else 0
    return np.array([max(-32768, min(32767, (int(v) + bias) >> shift)) for v in y], dtype=np.int64)

class TestCICPipelined(unittest.TestCase):
    def test_bit_exact_runtime_rates(self):
        for N, R in ((4, 8), (4, 64), (5, 32)):
            with self.subTest(N=N, R=R):
                rng = np.random.RandomState(R + N)
                n_out = 12
                x_i = rng.randint(-32768, 32768, R * (n_out + 2))
                x_q = rng.randint(-32768, 32768, R * (n_out + 2))
                samples = [{"i": int(i), "q": int(q)} for i, q in zip(x_i, x_q)]
                dut = LiteDSPCICDecimatorPipelined(r_max=512, N=N)
                shift = ceil(N * log2(R))
                dut.comb += [dut.rate.eq(R), dut.shift.eq(shift)]
                cap = run_stream(dut, samples, n_out, ["i", "q"], ["i", "q"],
                    sink_throttle=0.3, source_ready_rate=0.5)
                # The pipelined integrators put output k at input index (k+1)*R - 2.
                np.testing.assert_array_equal(column(cap, "i", 16), cic_model(x_i, R, N, shift)[R - 2::R][:n_out])
                np.testing.assert_array_equal(column(cap, "q", 16), cic_model(x_q, R, N, shift)[R - 2::R][:n_out])

if __name__ == "__main__":
    unittest.main()
