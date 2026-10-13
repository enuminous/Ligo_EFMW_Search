# GW170817 EFMW precursor search

<!-- ENUMINOUS-NETWORK:START -->
**eNuminous network:** [All repositories](https://enuminous.github.io/EFMW/repositories.html) · [EFMW](https://enuminous.github.io/EFMW/) · [102 equations](https://github.com/enuminous/Monolithic_102_EFMW) · [165 triplets](https://enuminous.github.io/FieldSpace/) · [Zoo](https://enuminous.github.io/Tortoise/) · [Lean](https://enuminous.github.io/Aristotle-102-Monolithic-Lean/) · [Engine](https://enuminous.github.io/Archimedes-Engine/) · [Papers](https://enuminous.github.io/medium/papers-essays-index.html) · [Audit](https://github.com/enuminous/EFMW_Post156_Zoo_Audit/blob/main/portfolio/INTERLOCK_AUDIT.md)

[Repository](https://github.com/enuminous/Ligo_EFMW_Search) · This repository is linked through its source; GitHub Pages is not enabled.

<details>
<summary>Recent repositories (9)</summary>

- [Oroboros](https://github.com/enuminous/Oroboros) · [README](https://github.com/enuminous/Oroboros/blob/main/README.md) · [Page](https://enuminous.github.io/Oroboros/)
- [P-NP](https://github.com/enuminous/P-NP) · [README](https://github.com/enuminous/P-NP/blob/main/README.md) · [Page](https://enuminous.github.io/P-NP/)
- [Research-Program](https://github.com/enuminous/Research-Program) · [README](https://github.com/enuminous/Research-Program/blob/main/README.md) · [Page](https://enuminous.github.io/Research-Program/)
- [I-Have-Guts](https://github.com/enuminous/I-Have-Guts) · [README](https://github.com/enuminous/I-Have-Guts/blob/main/README.md) · [Page](https://enuminous.github.io/I-Have-Guts/)
- [Wrights-Law](https://github.com/enuminous/Wrights-Law) · [README](https://github.com/enuminous/Wrights-Law/blob/main/README.md) · [Page](https://enuminous.github.io/Wrights-Law/)
- [Aristotle-eNuminous](https://github.com/enuminous/Aristotle-eNuminous) · [README](https://github.com/enuminous/Aristotle-eNuminous/blob/main/README.md) · [Page](https://enuminous.github.io/Aristotle-eNuminous/)
- [Multiverse-Simulator](https://github.com/enuminous/Multiverse-Simulator) · [README](https://github.com/enuminous/Multiverse-Simulator/blob/main/README.md) · [Page](https://enuminous.github.io/Multiverse-Simulator/)
- [accelerando-8192](https://github.com/enuminous/accelerando-8192) · [README](https://github.com/enuminous/accelerando-8192/blob/main/README.md) · [Page source](https://github.com/enuminous/accelerando-8192/blob/main/docs/index.html)
- [Codex-888](https://github.com/enuminous/Codex-888) · [README](https://github.com/enuminous/Codex-888/blob/main/README.md) · [Page](https://enuminous.github.io/Codex-888/)

</details>
<!-- ENUMINOUS-NETWORK:END -->

This folder contains a real, limited search for an early signal before the
binary-neutron-star merger GW170817.

## The question

EFMW ME-002 predicts that a scalar signal could travel faster than light. If
GW170817 produced such a signal, it should look like an early, weaker copy of
the published inspiral chirp.

This run tested that idea as **Channel A**: a time-shifted, linear breathing-mode
copy of the tensor waveform. It did not test Channel B, and it did not treat a
day-long signal as a chirp.

## What was done

- Public 4096 Hz GWOSC strain was downloaded for H1, L1, and V1.
- The search used the published GW170817 low-spin parameters and a Phenom
  NRTidal waveform.
- Real nearby O2 noise was used for 125 injections before the search windows
  were opened.
- The windows were opened in the requested order: W1, W2, W3, W0, W4.
- A candidate needed network SNR at least 8, at least two detectors with SNR at
  least 4, ordinary inter-site timing, acceptable chi-square/autocorrelation
  tests, and clean category-1 data quality.

## Main result

No window produced a candidate that passed all of those tests.

| Window | Early time tested | Maximum network SNR | Result |
|---|---:|---:|---|
| W0 | 0.1 s | 12.205 | Rejected: detector timing disagreed by 91.8 ms and L1 failed both vetoes |
| W1 | 1.7 s | 7.209 | Rejected: below network threshold and timing disagreed by 542.7 ms |
| W2 | 10 s | 6.638 | Rejected: no detector reached single-detector SNR 4 |
| W3 | 100 s | 7.011 | Rejected: below network threshold and timing disagreed by 288.8 ms |
| W4 | 3600 s | 7.371 | Rejected: below network threshold and timing disagreed by 412.8 ms |

The W0 SNR of 12.2 is not a detection. It is a noise/template feature that
failed the required coincidence and waveform-consistency tests.

## Injection check

The injection strength (r) is the scalar amplitude divided by the tensor
GW170817 amplitude. Each point below used five real-noise trials. The number in
parentheses is the median network SNR.

| Window | r = 1.0 | r = 0.3 | r = 0.1 | r = 0.03 | r = 0.01 |
|---|---:|---:|---:|---:|---:|
| W0 | 80% (32.04) | 20% (10.68) | 0% (7.08) | 0% (7.06) | 0% (7.06) |
| W1 | 40% (32.39) | 0% (11.16) | 0% (7.38) | 0% (7.36) | 0% (7.36) |
| W2 | 60% (31.67) | 0% (10.80) | 0% (7.05) | 0% (6.82) | 0% (6.82) |
| W3 | 40% (32.03) | 20% (10.67) | 0% (6.88) | 0% (6.88) | 0% (6.88) |
| W4 | 60% (32.25) | 40% (11.01) | 0% (6.80) | 0% (6.80) | 0% (6.80) |

These percentages use the SNR, two-detector, and timing cuts. The injection
branch did not rerun the chi-square/autocorrelation vetoes for every injection,
so this is not a full statistical detection-efficiency study.

## What the amplitude number means

For this template and noise estimate, an amplitude of 1 produced network SNR
7.3368. Using the brief's GW170817 tensor SNR of about 32 gives

```text
A = 4.3616 r
r expected to give SNR 8 = 0.25
A expected to give SNR 8 = 1.0904
```

Those are threshold scales, not 90% or 95% confidence limits. Five injections
per point are not enough to make a calibrated upper-limit statement.

## What this does not show

This run does not confirm ME-002, and a null result does not rule out EFMW or
the parameter alpha. EFMW fixes the timing shift, but it does not specify the
scalar amplitude. Under Channel A, the result can constrain only the amplitude
coupling (for example a product such as xi*kappa/M-star) after that coupling
model is supplied.

No Channel-B explanation was used to rescue a non-detection. No feature was
credited to EFMW unless it passed the frozen tests.

## Reproducibility

The exact scripts and all row-level results are included here:

- `gw170817_efmw_search.py` — data handling, template, injections, and initial search
- `gw170817_efmw_veto_search.py` — final veto-enabled search
- `gw170817_efmw_results.json` — measured injection and search results

The final search windows opened at `2026-09-08T14:42:37.238763+00:00`.

The injection random seed was `230908`. The JSON file records every injection
time, SNR, detector result, data file name, and checksum.

