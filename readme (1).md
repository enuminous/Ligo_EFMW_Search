# GW170817 EFMW precursor search

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
  least 4, an ordinary-light-travel **pair**, acceptable
  chi-square/autocorrelation tests, and clean category-1 data quality.
- A separate common-time coherent diagnostic used the signed scalar antenna
  factors and one shared phase.

## Main result

No window produced a candidate that passed all of those tests.

| Window | Early time tested | Maximum network SNR | Result |
|---|---:|---:|---|
| W0 | 0.1 s | 12.205 | Rejected: detector timing disagreed by 95.1 ms and L1 failed both vetoes |
| W1 | 1.7 s | 7.209 | Rejected: below network threshold and best pair disagreed by 90.9 ms |
| W2 | 10 s | 6.638 | Rejected: no detector reached single-detector SNR 4 |
| W3 | 100 s | 7.011 | Rejected: below network threshold and best pair disagreed by 285.6 ms |
| W4 | 3600 s | 7.371 | Rejected: below network threshold and best pair disagreed by 436.5 ms |

The W0 SNR of 12.2 is not a detection. It is a mostly single-detector
L1/template feature that failed the required pair timing and waveform
consistency tests. The common-time coherent value was 10.39, but no second
detector reached SNR 4.

The first implementation incorrectly required every detector above SNR 4 to
agree. The corrected run requires one valid pair, as the frozen experiment
specifies.

## Injection check

The injection strength (r) is the scalar amplitude divided by the tensor
GW170817 amplitude. Each point below used five real-noise trials. The number in
parentheses is the median network SNR.

| Window | r = 1.0 | r = 0.3 | r = 0.1 | r = 0.03 | r = 0.01 |
|---|---:|---:|---:|---:|---:|
| W0 | 100% (32.04) | 40% (10.68) | 20% (7.08) | 0% (7.06) | 0% (7.06) |
| W1 | 100% (32.39) | 40% (11.16) | 0% (7.38) | 0% (7.36) | 0% (7.36) |
| W2 | 100% (31.67) | 0% (10.80) | 0% (7.05) | 0% (6.82) | 0% (6.82) |
| W3 | 100% (32.03) | 20% (10.67) | 0% (6.88) | 0% (6.88) | 0% (6.88) |
| W4 | 100% (32.25) | 60% (11.01) | 0% (6.80) | 0% (6.80) | 0% (6.80) |

These percentages use the SNR, two-detector, and corrected pairwise timing cuts.
The injection branch did not rerun the chi-square/autocorrelation vetoes for
every injection, so this is not a full statistical detection-efficiency study.

## What the amplitude number means

For this template and noise estimate, an amplitude of 1 produced network SNR
7.3368. Using the brief's GW170817 tensor SNR of about 32 gives

```text
A = 4.3616 r
r expected to give SNR 8 = 0.25
A expected to give SNR 8 = 1.0904
```

Those are threshold scales, not 90% or 95% confidence limits. In this sky
geometry, a signal needs roughly A=1.59 (r=0.36) to reach SNR 4 in both H1 and
L1 from signal power alone. Five injections per point are not enough to make a
calibrated upper-limit statement.

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
- `gw170817_efmw_corrected_run.py` — corrected pairwise injection-first run
- `gw170817_efmw_coherent_check.py` — common-time coherent diagnostic
- `gw170817_efmw_corrected_results.json` — corrected measured injection/search results
- `gw170817_efmw_coherent_results.json` — coherent diagnostic results
- `gw170817_efmw_survey_requirements.md` — requirements for a publication-grade survey

The corrected injection campaign ran before the corrected search windows. The
corrected windows opened at `2026-09-08T15:38:32.358789+00:00`.

The injection random seed was `230908`. The JSON file records every injection
time, SNR, detector result, data file name, and checksum.
