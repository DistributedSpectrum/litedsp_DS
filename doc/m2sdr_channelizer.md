# M2SDR Channelizer: Layout, Costs and Next Steps

Status as of 2026-10-22, branch `m2sdr-channelizer-blocks` of litedsp_DS. Written as the
hand-off for the next stage: packetizing the channel streams as VITA 49 and sending them out
over PCIe.

- **Built and tested:** every block below exists in this branch, bit-exact against integer
  models under random input gaps and output stalls. Each block was placed and routed **on its
  own** on xc7a200tsbg484-3 at 245.76 MHz.
- **Not built yet:** the channelizer itself, meaning the top level that wires these blocks into
  the eight channels and connects them to the AD9361 stream.
- **Not integrated:** nothing is in the litex_m2sdr gateware. `litex_m2sdr_stack_DS` still pins
  litedsp at `8babdd4`, before these blocks.
- All costs below are sums of standalone place-and-route results, not a full build.

## Signal flow

```
AD9361 RX, 122.88 MSps (1T1R wide mode, rfic_clk 245.76 MHz)
  |
  v
Front halfband /2: 59 taps, keeps 50 MHz, 8 DSP48 (fold 4), shared by every channel
  |  61.44 MSps, 1 sample per 4 clocks
  v
DDC groups: LiteDSPMultiChannelNCOMixer, 4 channels per group, 1 RAMB18 per group
  |  every channel: its own phase_inc, still at 61.44 MSps
  +--> ch0  no further filtering                          61.44 MSps
  +--> ch1  HB39                                          30.72 MSps
  +--> ch2  HB19 > HB39                                   15.36 MSps
  +--> ch3  HB15 > HB19 > HB39                             7.68 MSps
  +--> ch4  HB15 > HB15 > HB19 > HB39                      3.84 MSps
  +--> ch5  HB15 > HB15 > HB15 > HB19 > HB39               1.92 MSps
  +--> ch6  CIC /16 > LPF53 /4                              960 ksps
  +--> ch7  CIC /32 > LPF53 /4                              480 ksps
```

HBn is an n-tap halfband (decimate by 2). LPF53 is the 53-tap low-pass that decimates by 4.

Rules behind the layout:

- **Halfbands only down to 1.92 MSps; below that, a CIC then a low-pass.**
- **Each channel has a fixed output rate.** So only its last halfband needs the sharp 39-tap
  filter (70 % of its output rate kept). Earlier halfbands only have to keep aliases off the
  final band, so they shrink to 19 and 15 taps.
- **The CIC goes first, never last.** Followed by the ÷4 low-pass, it rejects aliases by
  ≥ 80 dB. A CIC at the end of a chain gives only about 22 dB at 70 % usable bandwidth.
- **The two DDC groups have exactly 8 slots.** A ninth channel needs a third group (4 more
  DSP48 slots for mixers, 1 more RAMB18).

## Blocks

| Block | File | Parameters used here |
|---|---|---|
| Front halfband | `litedsp/filter/halfband_folded.py` `LiteDSPHalfbandDecimatorFolded` | `halfband_kaiser_coefficients(usable=50/61.44)` (59 taps, 80.2 dB), `fold=4`, `coeff_width=18` |
| DDC group | `litedsp/mixing/nco_mixer_multichannel.py` `LiteDSPMultiChannelNCOMixer` | `channels=4, fold=4` |
| HB39 / HB19 / HB15 | `LiteDSPHalfbandDecimatorFolded` | `halfband_kaiser_coefficients(usable=u)` with u = 0.70 (39 taps), 0.35 (19), ≤ 0.175 (15); fold = 245.76 MHz / output rate (capped at the number of distinct products) |
| CIC | `litedsp/filter/cic_pipelined.py` `LiteDSPCICDecimatorPipelined` | `N=4, r_max=16` or `32`. `rate` and `shift` are runtime registers. |
| LPF53 /4 | `litedsp/filter/fir_poly.py` `LiteDSPFIRDecimator` | `decimating_lowpass_coefficients(4)` from `litedsp/filter/design.py` (53 taps, 82.2 dB); `data_width=18` (samples sign-extended), `architecture="pipelined"`, `with_csr=False` |

Other blocks in this branch, not used in the layout: `LiteDSPLowSpurNCOMixer` (one channel,
folds 1/2/4) and `LiteDSPFoldedNCOMixer`. Tests are in `test/test_*` for each file above.

## Channels and usable bandwidth

"Usable" means the band kept with at least 80 dB of rejection for everything that folds onto
it. Group delay is that of the channel's own filters after the DDC. Add 236.0 ns for the front
halfband (29 samples at 122.88 MSps).

| Ch | Output rate | Usable bandwidth | Alias rejection | Passband | Filter group delay |
|---|---|---|---|---|---|
| 0 | 61.44 MSps | 50.0 MHz (front halfband) | 80.2 dB | ±0.001 dB | 0 |
| 1 | 30.72 MSps | 21.50 MHz | 82.6 dB | flat | 309.2 ns |
| 2 | 15.36 MSps | 10.75 MHz | ≥ 80 dB | flat | 765.0 ns |
| 3 | 7.68 MSps | 5.376 MHz | ≥ 80 dB | flat | 1,643.9 ns |
| 4 | 3.84 MSps | 2.688 MHz | ≥ 80 dB | flat | 3,401.7 ns |
| 5 | 1.92 MSps | 1.344 MHz | ≥ 80 dB | flat | 6,917.3 ns |
| 6 | 960 ksps | 672 kHz | 81.3 dB | −0.45 dB at edge (CIC droop) | 7,259.1 ns |
| 7 | 480 ksps | 336 kHz | 81.4 dB | −0.45 dB at edge (CIC droop) | 14,550.8 ns |

- Group delays assume linear-phase FIRs: (N − 1)/2 input samples each. The CIC adds
  N(R − 1)/2 input samples.
- On top of that, each block has a fixed pipeline latency of a few to a few tens of 245.76 MHz
  clocks. For example, the DDC takes 19 + k clocks for channel slot k.
- The low-pass could correct the CIC droop at no extra cost; `cic_comp_coefficients` in
  `design.py` is the starting point.

## Tuning

Each channel shifts the signal at f_off (relative to the RF centre) down to 0 Hz:

```
phase_inc = round(f_off * 2^32 / 61.44e6) mod 2^32
```

- Resolution is 0.0143 Hz.
- f_off must stay within ±25 MHz, the front halfband's 50 MHz.
- The first sample after reset uses phase (n + 1)·phase_inc.
- A retune is phase-continuous and takes effect on the next accepted sample.
- The oscillator's worst spur is about −110 dBc; at the 16-bit output it is about −98 dBc.

## Costs

Every block was placed and routed alone on xc7a200tsbg484-3 at 245.76 MHz.

| Block | DSP48 | LUT | FF | BRAM | Fmax |
|---|---|---|---|---|---|
| Front halfband (59 taps, fold 4) | 8 | 584 | 3,478 | 0 | 265.8 MHz |
| DDC group (4 channels) | 4 | 849 | 1,556 | 1 RAMB18 | 282.2 MHz |
| HB39, fold 8 | 4 | 390 | 2,210 | 0 | 271.4 MHz |
| HB39, fold ≥ 11 | 2 | 317 | 2,055 | 0 | 277.0 MHz |
| HB19 (fold 8) | 2 | 223 | 1,090 | 0 | 294.8 MHz |
| HB15 (fold 8; assumed the same at larger folds) | 2 | 205 | 898 | 0 | 303.8 MHz |
| CIC /16 | 0 | 846 | 991 | 0 | 286.8 MHz |
| CIC /32 | 0 | 941 | 1,104 | 0 | 292.3 MHz |
| LPF53 /4, pipelined, 18-bit | 2 | 463 | 190 | 0 | 232.0 MHz |

Per channel, with the DDC counted as a quarter of its group:

| Ch | Output rate | DSP48 | LUT | FF |
|---|---|---|---|---|
| 0 | 61.44 MSps | 1 | 212 | 389 |
| 1 | 30.72 MSps | 5 | 602 | 2,599 |
| 2 | 15.36 MSps | 5 | 752 | 3,534 |
| 3 | 7.68 MSps | 7 | 957 | 4,432 |
| 4 | 3.84 MSps | 9 | 1,162 | 5,330 |
| 5 | 1.92 MSps | 11 | 1,367 | 6,228 |
| 6 | 960 ksps | 3 | 1,521 | 1,570 |
| 7 | 480 ksps | 3 | 1,616 | 1,683 |
| | **8 channels** | **44** | **8,189** | **25,765** |
| | **+ front halfband** | **52** | **8,773** | **29,243** |

Block RAM is 2 RAMB18, one per DDC group.

How it fits beside the PCIe x4 base build. The base uses 0 DSP48, 12,111 LUT, 11,383 FF and
76.5 RAMB36 (2026-10-20 build).

| Resource | Base + channelizer | Device | Share |
|---|---|---|---|
| DSP48 | 52 | 740 | 7.0 % |
| LUT | 20,884 | 133,800 | 15.6 % |
| FF | 40,626 | 267,600 | 15.2 % |
| RAMB36 (equivalent) | 77.5 | 365 | 21.2 % |

Flip-flops limit how many copies of this 8-channel set fit: 6 sets (48 channels) at a 70 % fill (the
usual routability ceiling) and 9 sets (72 channels) using the whole chip. Most of each channel's
flip-flops are its halfband delay lines.

## Timing status

- **Every block meets 245.76 MHz except LPF53.** The low-pass reaches 232.0 MHz. Its critical
  path runs from the tap counter, through the delay-line address and the LUT-RAM read, into
  the DSP48 A input. Two fixes:
  1. Run both LPF53s on the 125 MHz system clock behind a small clock-domain-crossing FIFO.
     The channel output goes to the DMA in that domain anyway, and at 125 MHz there are still
     at least 130 clocks per output; the filter needs 59.
  2. Add one register on the read address in `LiteDSPFIRDecimator`. This costs one more cycle
     per output (60 of the 256 or 512 available).
- **Not yet checked in a full design.** All of these numbers are for blocks on their own;
  timing has to be confirmed once the channelizer shares the chip with the PCIe design.

## Integration notes for the next stage (VITA 49 over PCIe)

### Clock domain: decide this first

- **All costs above assume the processing runs at 245.76 MHz.** That is `rfic_clk` in the
  122.88 MSps 1T1R wide-mode build.
- **litex_m2sdr's processing hook is on the system clock.** `add_rfic_datapath()` in
  `litex_m2sdr.py` hands over the RX samples on the 125 MHz system clock, as 64-bit
  `dma_layout` words, after the AD9361 RX clock-domain crossing in
  `litex_m2sdr/gateware/ad9361/core.py`.
- **Option A: tap the samples in the `rfic` domain before that crossing.** The layout and costs
  above then hold.
- **Option B: rebuild the layout for 125 MHz.** There are about half as many clocks per
  sample, so every fold halves and DSP48 counts roughly double. For example, the front
  halfband becomes about fold 2, 16 DSP48.

### Channel outputs

- **Format:** each output is a `stream.Endpoint(iq_layout(16))`, with signed 16-bit `i` and
  `q`, valid/ready handshaking, and no `first`/`last` framing. The packetizer adds framing.
- **The ADC never pauses, so every output must be drained at its average rate.** Within one
  DDC group, all four channels accept a sample together. If any of the four output FIFOs (16
  deep) holds 8 or more samples, input to the whole group stops. Halfbands, CIC and low-pass
  also apply back-pressure. A packetizer that holds a channel for a whole packet needs a FIFO
  of at least one packet per channel ahead of it.
- **Total payload:** all eight channels together carry 122.40 MSps × 32 bits = 3.92 Gb/s. That
  is the same as the raw 122.88 MSps sc16 stream the Gen2 x4 build already carries without
  overflow. Gen2 x1 tops out at about 3.32 Gb/s, so on x1 drop channel 0 (the remaining seven
  carry 60.96 MSps = 1.95 Gb/s) or use 8-bit samples. VITA headers and context packets add to
  this.

### Existing pieces to reuse

- **`litex_m2sdr/gateware/vrt.py`.** `VRTSignalPacketInserter` and
  `VRTSignalPacketStreamer` already build VRT (VITA 49) signal-data packets, with packet types
  and TSI/TSF timestamp enums. Today they are used only for the Ethernet UDP streamer
  (`with_eth_vrt`); the time inputs come from `time_s_vrt`/`time_ps_vrt` in `litex_m2sdr.py`.
- **`litex_m2sdr/gateware/header.py`.** The current PCIe DMA framing: a 16-byte header
  (sync word plus a timestamp in ns) with the RX timestamp taken from `time_gen`.
- **More DMA channels.** `pcie_dmas > 1` (`--pcie-dmas`) adds `self.pcie_dma1`, … after the
  hook, if the channels go out on separate DMAs.
- **litedsp `doc/timestamps.md` and `litedsp/frontend/packet.py`.** A time core, a
  timestamper and a packetizer header modelled on VITA 49 but explicitly not wire-compliant.
  Useful for the timestamp arithmetic.

### Timestamps and context packets

- **Timestamps.** Latch the time where samples enter the channelizer, at the front halfband
  input. An output sample's time is then that ingress time, plus its sample index scaled by
  the channel's decimation, minus the filter group delay in the channels table, plus the fixed
  pipeline latency (constant per channel; measure once in simulation).
- **Decimation phase.** Each halfband's `phase` argument sets which input sample an output
  lines up with. Keep it at the default so the alignment stays deterministic.
- **Context packets** need, per channel:
  - RF reference frequency: AD9361 LO + f_off, with f_off from `phase_inc`.
  - Bandwidth: the usable bandwidth in the channels table.
  - Sample rate: the output rate.
  - Stream ID: the channel number.

## Possible future change: CIC + low-pass for more rates

A CIC followed by LPF53 is far cheaper in flip-flops than a halfband chain: about 1,441 FF
against 6,228 FF at 1.92 MSps. The low-pass needs 59 clocks per output, so it fits wherever
the output rate leaves at least 59 clocks. Measured alias rejection with the same 53-tap
low-pass and a CIC N = 4 in front:

| Output | Chain | Alias rejection | Clocks per output (needs 59) | Verdict |
|---|---|---|---|---|
| 7.68 MSps | CIC /2 > LPF53 /4 | 68.7 dB | 32 | No: too few clocks, too little rejection |
| 3.84 MSps | CIC /4 > LPF53 /4 | 78.5 dB | 64 | Possible: 1.5 dB short of 80 dB; a longer low-pass needs a 2-MAC FIR |
| 1.92 MSps | CIC /8 > LPF53 /4 | 80.7 dB | 128 | Yes |
| 960 ksps | CIC /16 > LPF53 /4 | 81.3 dB | 256 | In use |
| 480 ksps | CIC /32 > LPF53 /4 | 81.4 dB | 512 | In use |

Moving 1.92 MSps over would cost about 3 DSP48, 1,362 LUT and 1,441 FF,
instead of 11 DSP48, 1,367 LUT and 6,228 FF. That saves 8 DSP48 and about 4,800 FF per set.
Further options:

- **Lower rates.** The same pair reaches 240 ksps (CIC /64) or below, with alias rejection
  above 81 dB and growing room per output.
- **CIC droop.** A droop-compensated low-pass (`cic_comp_coefficients`) would flatten the
  −0.45 dB edge at no extra cost.
- **One CIC per DDC group.** A single CIC with several runtime rates, fed by a time-shared
  low-pass, could serve all four channels of a group at different low rates. That would need
  a multi-channel CIC; none exists yet.
