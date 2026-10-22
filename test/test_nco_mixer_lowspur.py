#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

import math
import unittest

import numpy as np
from scipy.signal.windows import kaiser

from litedsp.mixing.nco_mixer_lowspur import LiteDSPLowSpurNCOMixer, nco_table

from test.common import run_stream, column

A, D, B, KW = 12, 10, 18, 18
TAB = nco_table(A, B); Q = 1 << (A - 2); K = int(round(2*math.pi/(1 << A)*(1 << KW)))

def nco_model(phase):
    """Integer model of the low-spur NCO: (cos, sin) for a 32-bit phase."""
    a  = phase >> (32 - A)
    fr = ((phase >> (32 - A - D)) & ((1 << D) - 1)) - (1 << (D - 1))
    def sin_of(a):
        q, i = a >> (A - 2), a & (Q - 1)
        v = TAB[Q - 1 - i] if q & 1 else TAB[i]
        return -v if q & 2 else v
    s, c = sin_of(a), sin_of((a + Q) & ((1 << A) - 1))
    dq = (fr * K) >> D
    return c - ((dq * s + (1 << (KW - 1))) >> KW), s + ((dq * c + (1 << (KW - 1))) >> KW)

def mixer_model(x_i, x_q, phase_inc):
    sat = lambda v: max(-32768, min(32767, (v + (1 << (B - 2))) >> (B - 1)))
    out_i, out_q = [], []
    for k in range(len(x_i)):
        c, s = nco_model(((k + 1) * phase_inc) & 0xFFFFFFFF)
        out_i.append(sat(int(x_i[k])*c + int(x_q[k])*s)); out_q.append(sat(int(x_q[k])*c - int(x_i[k])*s))
    return np.array(out_i), np.array(out_q)

def worst_spur_dbc(x):
    X = np.abs(np.fft.fft(x * kaiser(len(x), 38)))**2; n = len(X); k = np.argmax(X)
    d = np.minimum((np.arange(n) - k) % n, (k - np.arange(n)) % n)
    return 10*np.log10(X[d > 40].max() / X[k])

PHASE_INCS = (0x12345678, 0x0F0F0F0F, 0x2AAAAAAB, 0x19999999, 0x0000F001)

class TestLowSpurNCOMixer(unittest.TestCase):
    def test_bit_exact_vs_model(self):
        rng = np.random.RandomState(5)
        x_i = rng.randint(-32768, 32768, 200); x_i[::13] = -32768
        x_q = rng.randint(-32768, 32768, 200)
        samples = [{"i": int(i), "q": int(q)} for i, q in zip(x_i, x_q)]
        for pinc in (0x12345678, 0x80000001, 0x0000F001):
            ref_i, ref_q = mixer_model(x_i, x_q, pinc)
            for fold, dsp48 in ((1, 4), (2, 2), (4, 1)):
                with self.subTest(phase_inc=hex(pinc), fold=fold):
                    dut = LiteDSPLowSpurNCOMixer(fold=fold)
                    self.assertEqual(dut.dsp48, dsp48)          # mixer lanes (corrections add 2 / 1 / 1)
                    dut.comb += dut.phase_inc.eq(pinc)
                    cap = run_stream(dut, samples, len(samples), ["i", "q"], ["i", "q"], sink_throttle=0.2, source_ready_rate=0.6)
                    np.testing.assert_array_equal(column(cap, "i", 16), ref_i)
                    np.testing.assert_array_equal(column(cap, "q", 16), ref_q)

    def test_spurs_below_16_bit(self):
        n = 1 << 14
        nco = [worst_spur_dbc(np.array([complex(*nco_model(((k + 1) * p) & 0xFFFFFFFF)) for k in range(n)])) for p in PHASE_INCS]
        self.assertLess(max(nco), -105.0)                       # NCO alone (measured about -111 dBc)
        out = []
        for p in PHASE_INCS:
            i, q = mixer_model([23170] * n, [0] * n, p)         # -3 dBFS constant input -> output is the NCO tone
            out.append(worst_spur_dbc(i + 1j*q))
        self.assertLess(max(out), -95.0)                        # 16-bit output (measured about -98 dBc)

if __name__ == "__main__":
    unittest.main()
