# EFMW/GW170817 precursor-survey requirements

This document separates what EFMW actually predicts from the additional
assumptions and analysis choices needed to search LIGO/Virgo data.

## 1. The precise EFMW claim

ME-002 fixes a propagation speed and therefore a time advance:

```text
v = c / sqrt(1 - alpha^2)
Delta-t(alpha) = T * (1 - sqrt(1 - alpha^2))
```

For GW170817, use the frozen 40 Mpc distance and the supplied W0-W4 time
shifts. The shift is a source-clock prediction. It is not an allowance to
slide the signal by an arbitrary amount between detectors.

EFMW does **not** fix the scalar amplitude. A result must therefore be reported
as an amplitude limit or sensitivity curve, `A(Delta-t)`, or as the relative
amplitude `r = A / h_tensor`. A null cannot by itself be reported as a bound on
alpha.

## 2. Channel-A assumptions

The search tested only this additional model:

```text
h_i(t) = A * F_b,i * C(t + Delta-t)
          * cos(Psi(t + Delta-t) + phi0)
```

The same inspiral clock is reused in every detector. The scalar breathing
antenna factor and the published sky position determine the relative detector
response. Channel B, involving quadratic scalar stress in vacuum, is a
separate null and must not be used to reinterpret a Channel-A non-detection.

## 3. Data requirements

Use public calibrated strain for H1, L1, and V1 at 4096 Hz. The data span must
cover every template sample, including the approximately 66-second inspiral
before each proposed scalar coalescence time. Screen the whole analyzed segment
with the DATA and CBC category-1 quality bits.

Do not begin with a 30-day pull. Start with the frozen 4096-second event data
and add only the adjacent data needed for W4 and off-source injections.

## 4. Waveform requirements

The production search should use a posterior bank, not only one median
waveform. Draw or tile over the published GW170817 masses, mass ratio, spins,
tidal parameters, sky position, and the waveform approximant choice
(Phenom-NRTidal or SEOBNRT). Retain the inspiral clock and truncate at a
predefined coalescence point.

The nuisance phase `phi0` must be handled by two explicitly defined quadrature
templates or an equivalent validated complex statistic. The template and the
amplitude normalization must be frozen before the windows are opened.

## 5. Network-statistic requirements

The primary statistic should be a coherent scalar network likelihood. It should

1. evaluate one shared geocentric coalescence time;
2. maximize over one common scalar phase;
3. use the signed breathing antenna factors and detector PSDs;
4. retain the frozen network-SNR threshold of 8; and
5. require at least one detector pair with SNR at least 4 and ordinary
   light-travel timing.

An unrelated third-detector noise excursion must not veto a valid pair. The
pairwise timing test is applied after correcting detector times back to the
geocenter. Delays of order the EFMW time advance are never permitted between
sites.

## 6. Injection requirements

Injections must be completed before opening any search window. Use real nearby
noise, not the exact search stretch. For each W0-W4, use the frozen amplitudes

```text
r = 1.0, 0.3, 0.1, 0.03, 0.01
```

and enough independent trials to measure an efficiency curve. Randomize the
phase and sample waveform/sky nuisance parameters from the published posterior
unless the experiment is explicitly defined as a single-template test.

The recovery definition must include every cut used for the search: network
SNR, two-detector threshold, pairwise timing, data quality, chi-square,
autocorrelation, and any ranking statistic. Five trials per point are useful as
a smoke test but cannot establish a 90% or 95% upper limit.

## 7. Veto, background, and calibration requirements

Declare the exact power-chi-square and autocorrelation implementation and
thresholds before unblinding. Estimate the off-source background with many
time slides or equivalent trials that preserve ordinary inter-site delays.
Account for the five tested windows, waveform-bank trials, calibration
uncertainty, and PSD nonstationarity.

The W4 search needs a local or demonstrably stationary PSD. A single PSD
segment near the merger should not automatically be assumed valid several
thousand seconds earlier. Glitch gating and calibration validity must be
recorded per detector and per window.

## 8. Corrected work completed here

The first implementation had a stricter-than-requested timing rule: it required
all detectors above SNR 4 to agree, rather than requiring one valid pair. That
was corrected to a geocentrically corrected pairwise test. The corrected runner
also records the pair, residual, detector times, and injection components.

The corrected injection run recovered every `r=1.0` injection in all five
windows. The earlier 40%-80% result was therefore partly a pipeline artifact,
not evidence that strong signals were intrinsically missed.

A separate common-time coherent diagnostic was validated on an `r=1.0`
injection and then applied to the search windows. It found no coherent SNR-8
precursor. The W0 coherent excess remained dominated by L1 and did not have a
second detector at SNR 4.

## 9. What is still required before a final null

The corrected run is a stronger and more honest smoke test, but not yet a
publication-grade upper-limit analysis. Before calling the result a final null:

- run the full posterior waveform bank;
- use the coherent statistic as the prespecified decision statistic;
- rerun injections through the complete veto pipeline with randomized phase;
- build an off-source background and trials correction;
- use local PSDs and explicit glitch/calibration checks; and
- report `A_max(Delta-t)` with a stated confidence level.

If that complete analysis finds no surviving candidate, the valid conclusion is
a Channel-A amplitude bound at the tested delays. It is not confirmation of
ME-002, and it is not a standalone exclusion of alpha.

