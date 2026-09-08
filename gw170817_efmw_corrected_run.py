#!/usr/bin/env python3
"""Corrected Channel-A GW170817 run.

This reruns the frozen injection campaign and the frozen search order after
fixing the coincidence rule: the event needs one ordinary-light-travel pair,
not every detector that happens to cross the single-detector threshold.

The injection branch exploits linearity of the matched filter.  For each
noise segment and detector it computes the real-noise filter once and the
unit-amplitude injection response once, then evaluates all five amplitudes.
This is algebraically equivalent to filtering each injected segment while
avoiding unnecessary repeated FFTs.
"""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pycbc.filter import matched_filter
from pycbc.types import TimeSeries

import gw170817_efmw_search as g
import gw170817_efmw_veto_search as veto


OUT = Path("gw170817_efmw_corrected_results.json")
SEARCH_ORDER = ["W1", "W2", "W3", "W0", "W4"]


def filter_template(data: TimeSeries, tmpl: g.Template, psd) -> TimeSeries:
    h = tmpl.hp.copy()
    h.resize(len(data))
    return matched_filter(
        h,
        data,
        psd=g.psd_for_data(psd, len(data)),
        low_frequency_cutoff=30.0,
    )


def raw_layout(tmpl: g.Template, det: str, t_geo: float) -> tuple[int, int]:
    t_arr = t_geo + g.delay(det, t_geo)
    start = t_arr - tmpl.peak_offset - g.PRE_PAD
    end = t_arr + g.POST_PAD
    raw = g.read_detector(det, start, end)
    i0 = int(round((t_arr - tmpl.peak_offset - float(raw.start_time)) * g.FS))
    return len(raw), i0


def unit_response(
    tmpl: g.Template,
    det: str,
    psd,
    n: int,
    i0: int,
    epoch: float,
) -> TimeSeries:
    """Matched-filter response for A=1, with the detector's F_b included."""
    x = np.zeros(n, dtype=np.float64)
    sig = g.fb(det) * np.asarray(tmpl.hp)
    x[i0:i0 + len(sig)] += sig
    raw = TimeSeries(x, delta_t=1.0 / g.FS, epoch=epoch)
    return filter_template(g.clean(raw), tmpl, psd)


def noise_and_unit(
    tmpl: g.Template,
    det: str,
    t_geo: float,
    psd,
    response_cache: dict,
) -> tuple[TimeSeries, TimeSeries]:
    t_arr = t_geo + g.delay(det, t_geo)
    start = t_arr - tmpl.peak_offset - g.PRE_PAD
    end = t_arr + g.POST_PAD
    raw = g.read_detector(det, start, end)
    noise = g.clean(raw)
    i0 = int(round((t_arr - tmpl.peak_offset - float(raw.start_time)) * g.FS))
    key = (det, len(raw), i0)
    if key not in response_cache:
        # The unit response is translation invariant for these equal-length
        # segments.  Reset its epoch below to the current cleaned segment.
        response_cache[key] = np.asarray(
            unit_response(tmpl, det, psd, len(raw), i0, 0.0), dtype=np.complex128
        )
    z_noise = filter_template(noise, tmpl, psd)
    z_unit = TimeSeries(
        response_cache[key], delta_t=z_noise.delta_t, epoch=z_noise.start_time
    )
    return z_noise, z_unit


def local_complex_peak(z: TimeSeries, center: float, half_width: float = 0.5):
    lo = max(0, int((center - half_width - float(z.start_time)) * g.FS))
    hi = min(len(z), int((center + half_width - float(z.start_time)) * g.FS) + 1)
    if hi <= lo:
        raise ValueError("empty local search interval")
    arr = np.asarray(z)
    j = lo + int(np.argmax(np.abs(arr[lo:hi])))
    return float(abs(arr[j])), complex(arr[j]), float(z.start_time + j / g.FS)


def injection_campaign(tmpl, psds, norms):
    rng = np.random.default_rng(g.SEED)
    inj_times = {
        name: sorted(
            rng.uniform(g.T0 - 3400.0, g.T0 - 3000.0, g.N_INJECTION_REPS).tolist()
        )
        for name, _, _, _ in g.WINDOWS
    }
    response_cache = {}
    rows = []
    started = datetime.now(timezone.utc).isoformat()
    for name, dt, _, _ in g.WINDOWS:
        for rep, t_geo in enumerate(inj_times[name]):
            zbase = {}
            for det in ["H1", "L1", "V1"]:
                zbase[det] = noise_and_unit(
                    tmpl, det, t_geo, psds[det], response_cache
                )
            for r in g.R_LEVELS:
                A = r * 32.0 / norms["network"]
                dscores = {}
                for det in ["H1", "L1", "V1"]:
                    zn, zu = zbase[det]
                    z = TimeSeries(
                        np.asarray(zn) + A * np.asarray(zu),
                        delta_t=zn.delta_t,
                        epoch=zn.start_time,
                    )
                    center = t_geo + g.delay(det, t_geo) - tmpl.peak_offset
                    peak, _, tp = local_complex_peak(z, center, 0.5)
                    dscores[det] = {"snr": peak, "time": tp}
                active = [d for d in dscores if dscores[d]["snr"] >= 4.0]
                net = math.sqrt(sum(dscores[d]["snr"] ** 2 for d in dscores))
                coincidence = g.pairwise_coincidence(
                    {d: dscores[d]["time"] for d in dscores}, active, t_geo
                )
                quality = all(
                    g.quality_pass(
                        d,
                        t_geo + g.delay(d, t_geo)
                        - tmpl.peak_offset
                        - g.PRE_PAD,
                        t_geo + g.delay(d, t_geo) + g.POST_PAD,
                    )
                    for d in ["H1", "L1", "V1"]
                )
                rows.append(
                    {
                        "window": name,
                        "delta_t_s": dt,
                        "rep": rep,
                        "t_geo": t_geo,
                        "r": r,
                        "A": A,
                        "single_snr": {
                            d: dscores[d]["snr"] for d in dscores
                        },
                        "times": {d: dscores[d]["time"] for d in dscores},
                        "network_snr": net,
                        "active_ge4": active,
                        "two_detector_pass": len(active) >= 2,
                        "coincidence_pass": coincidence["pass"],
                        "coincidence_pair": coincidence["pair"],
                        "coincidence_residual_ms": coincidence["residual_ms"],
                        "pass": bool(net >= 8.0 and coincidence["pass"]),
                        "quality_pass": quality,
                    }
                )
    finished = datetime.now(timezone.utc).isoformat()
    return rows, inj_times, response_cache, started, finished


def search_campaign(tmpl, psds, windows_opened_utc):
    rows = []
    for target in SEARCH_ORDER:
        name, dt, t_geo, alpha = next(w for w in g.WINDOWS if w[0] == target)
        data = {
            d: g.analysis_data(tmpl, d, t_geo) for d in ["H1", "L1", "V1"]
        }
        single = {
            d: veto.detector_stat(data[d], tmpl, psds[d], d, t_geo)
            for d in ["H1", "L1", "V1"]
        }
        times = {d: single[d]["time"] for d in single}
        active = [d for d in single if single[d]["snr"] >= 4.0]
        net = math.sqrt(sum(single[d]["snr"] ** 2 for d in single))
        coincidence = g.pairwise_coincidence(times, active, t_geo)
        pair = coincidence["pair"] or []
        pair_veto = bool(pair) and all(single[d]["veto_pass"] for d in pair)
        all_active_veto = all(single[d]["veto_pass"] for d in active)
        quality = all(
            g.quality_pass(
                d,
                t_geo + g.delay(d, t_geo) - tmpl.peak_offset - g.PRE_PAD,
                t_geo + g.delay(d, t_geo) + g.POST_PAD,
            )
            for d in ["H1", "L1", "V1"]
        )
        rows.append(
            {
                "window": name,
                "delta_t_s": dt,
                "gps_nominal": t_geo,
                "alpha_abs": alpha,
                "single": single,
                "times": times,
                "network": net,
                "active_ge4": active,
                "network_pass": net >= 8.0,
                "two_detector_pass": len(active) >= 2,
                "coincidence_pass": coincidence["pass"],
                "coincidence_pair": coincidence["pair"],
                "coincidence_residual_ms": coincidence["residual_ms"],
                "coincidence_pairs": coincidence["pairs"],
                "chi2_autochisq_veto_pass": pair_veto,
                "all_active_chi2_autochisq_veto_pass": all_active_veto,
                "candidate_pass": bool(net >= 8.0 and coincidence["pass"] and pair_veto),
                "quality_pass": quality,
                "best_delta_t_s": g.T0
                - np.mean(
                    [
                        times[d] + tmpl.peak_offset - g.delay(d, t_geo)
                        for d in times
                    ]
                ),
            }
        )
    return rows


def main():
    tmpl = g.make_template()
    psds = {d: g.prepare_psd(d)[0] for d in ["H1", "L1", "V1"]}
    probe = g.analysis_data(tmpl, "H1", g.T0 - 3300.0)
    norms = g.scalar_norms(tmpl, psds, len(probe))
    injections, inj_times, response_cache, inj_started, inj_finished = injection_campaign(
        tmpl, psds, norms
    )
    windows_opened = datetime.now(timezone.utc).isoformat()
    search = search_campaign(tmpl, psds, windows_opened)
    meta = {
        "template": {
            "approximant": "IMRPhenomPv2_NRTidalv2",
            "mass1_detector_Msun": 1.4790595881585529,
            "mass2_detector_Msun": 1.2807580711111188,
            "spin1z": 0.004953281544706123,
            "spin2z": 0.0028389018809797492,
            "lambda1": 269.0083686566036,
            "lambda2": 445.770904970051,
            "distance_Mpc": 39.95566849878731,
            "inclination_rad": 2.5596915857494036,
            "f_lower_Hz": 30.0,
            "samples": len(tmpl.hp),
            "duration_s": len(tmpl.hp) / g.FS,
            "peak_offset_s": tmpl.peak_offset,
        },
        "antenna_fb": {d: g.fb(d) for d in ["H1", "L1", "V1"]},
        "files": {
            d: {
                "event": str(g.EVENT_FILES[d]),
                "event_sha256": g.sha256(g.EVENT_FILES[d]),
                "pre": str(g.PRE_FILES[d]),
                "pre_sha256": g.sha256(g.PRE_FILES[d]),
            }
            for d in ["H1", "L1", "V1"]
        },
    }
    result = {
        "metadata": meta,
        "scalar_norms": norms,
        "injection_times": inj_times,
        "injections": injections,
        "injection_started_utc": inj_started,
        "injection_finished_utc": inj_finished,
        "search": search,
        "search_order": SEARCH_ORDER,
        "windows_opened_utc": windows_opened,
        "seed": g.SEED,
        "corrections": {
            "coincidence": "At least one active detector pair must agree after geocentric delay correction within 30 ms.",
            "injection_filtering": "Cached exact linear unit-injection matched-filter response; same signal and PSD as original branch.",
            "veto_pair_scope": "Chi-square/autocorrelation veto applied to the pair establishing coincidence.",
            "coherent_network_statistic": "Not used as a decision statistic in this run; remains a required next-stage improvement.",
            "injection_efficiency_scope": "Network SNR, two-detector, pairwise timing, and data quality; vetoes not rerun for every injection.",
        },
        "veto_settings": {
            "power_chisq_bins": veto.CHI2_BINS,
            "chi2_reduced_max": veto.CHI2_REDUCED_MAX,
            "autochisq_stride_samples": veto.AUTO_STRIDE,
            "autochisq_points": veto.AUTO_POINTS,
            "autochisq_reduced_max": veto.AUTO_REDUCED_MAX,
        },
        "response_cache_keys": [list(k) for k in sorted(response_cache)],
    }
    OUT.write_text(json.dumps(result, indent=2))
    print(json.dumps({
        "output": str(OUT),
        "injection_started_utc": inj_started,
        "injection_finished_utc": inj_finished,
        "windows_opened_utc": windows_opened,
        "search": [
            {
                "window": r["window"],
                "network": r["network"],
                "active": r["active_ge4"],
                "pair": r["coincidence_pair"],
                "pair_residual_ms": r["coincidence_residual_ms"],
                "candidate": r["candidate_pass"],
            }
            for r in search
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
