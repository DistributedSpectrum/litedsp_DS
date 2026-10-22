#
# This file is part of LiteDSP.
#
# Copyright (c) 2026 DistributedSpectrum
# SPDX-License-Identifier: BSD-2-Clause

import math
import unittest

import numpy as np
from scipy.signal.windows import kaiser

from migen import run_simulation

from litedsp.mixing.nco_mixer_lowspur import nco_table
from litedsp.mixing.nco_mixer_multichannel import LiteDSPMultiChannelNCOMixer

from test.common import stream_driver, stream_capture, column

A, D, B, KW, CB = 12, 10, 18, 18, 10
TAB = nco_table(A, B); Q = 1 << (A - 2); K = int(round(2*math.pi/(1 << A)*(1 << KW)))
SH = B + 1 - CB; RS = KW - SH

def nco_model(phase):
    """Integer model of one channel's NCO: (cos, sin), correction on cos/sin cut to CB bits."""
    a  = phase >> (32 - A)
    fr = ((phase >> (32 - A - D)) & ((1 << D) - 1)) - (1 << (D - 1))
    def sin_of(a):
        q, i = a >> (A - 2), a & (Q - 1)
        v = TAB[Q - 1 - i] if q & 1 else TAB[i]
        return -v if q & 2 else v
    s, c = sin_of(a), sin_of((a + Q) & ((1 << A) - 1))
    dq = (fr * K) >> D
    return c - ((dq * (s >> SH) + (1 << (RS - 1))) >> RS), s + ((dq * (c >> SH) + (1 << (RS - 1))) >> RS)

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

PHASE_INCS = (0x12345678, 0x80000001, 0x0000F001, 0xD6555555)

class TestMultiChannelNCOMixer(unittest.TestCase):
    def test_bit_exact_per_channel(self):
        rng = np.random.RandomState(7)
        x_i = rng.randint(-32768, 32768, 160); x_i[::11] = -32768
        x_q = rng.randint(-32768, 32768, 160)
        samples = [{"i": int(i), "q": int(q)} for i, q in zip(x_i, x_q)]
        for channels, fold in ((4, 4), (3, 4), (1, 4), (4, 8)):
            with self.subTest(channels=channels, fold=fold):
                dut = LiteDSPMultiChannelNCOMixer(channels=channels, fold=fold)
                self.assertEqual(dut.dsp48, channels)
                for k in range(channels):
                    dut.comb += dut.phase_incs[k].eq(PHASE_INCS[k])
                caps = [[] for _ in range(channels)]
                gens = [stream_driver(dut.sink, samples, ["i", "q"], seed=1, throttle=0.2)]
                gens += [stream_capture(dut.sources[k], caps[k], len(samples), ["i", "q"], seed=20 + k, ready_rate=0.6)
                         for k in range(channels)]
                run_simulation(dut, gens)
                for k in range(channels):
                    ref_i, ref_q = mixer_model(x_i, x_q, PHASE_INCS[k])
                    np.testing.assert_array_equal(column(caps[k], "i", 16), ref_i, err_msg=f"channel {k} I")
                    np.testing.assert_array_equal(column(caps[k], "q", 16), ref_q, err_msg=f"channel {k} Q")

    def test_spurs_below_16_bit(self):
        n = 1 << 14
        incs = PHASE_INCS + (0x0F0F0F0F, 0x2AAAAAAB, 0x19999999)
        nco = [worst_spur_dbc(np.array([complex(*nco_model(((k + 1) * p) & 0xFFFFFFFF)) for k in range(n)])) for p in incs]
        self.assertLess(max(nco), -105.0)                       # NCO alone (measured about -110 dBc)
        out = []
        for p in incs:
            i, q = mixer_model([23170] * n, [0] * n, p)         # -3 dBFS constant input -> output is the NCO tone
            out.append(worst_spur_dbc(i + 1j*q))
        self.assertLess(max(out), -95.0)                        # 16-bit output (about -98 dBc)

if __name__ == "__main__":
    unittest.main()
