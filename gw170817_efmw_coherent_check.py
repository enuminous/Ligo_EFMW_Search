#!/usr/bin/env python3
"""Common-time coherent diagnostic for the corrected frozen windows.

This is deliberately reported as a diagnostic, not silently substituted for
the frozen network-SNR decision rule.  It evaluates one shared geocentric
coalescence offset and maximizes over the common scalar phase using the signed
breathing antenna weights.
"""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pycbc.filter import matched_filter

import gw170817_efmw_search as g


OUT = Path("gw170817_efmw_coherent_results.json")
ORDER = ["W1", "W2", "W3", "W0", "W4"]


def complex_filter(data, tmpl, psd):
    h = tmpl.hp.copy()
    h.resize(len(data))
    return matched_filter(
        h,
        data,
        psd=g.psd_for_data(psd, len(data)),
        low_frequency_cutoff=30.0,
    )


def sample(z, target):
    idx = int(round((target - float(z.start_time)) * g.FS))
    if idx < 0 or idx >= len(z):
        return 0j
    return complex(np.asarray(z)[idx])


def coherent_scan(zs, tmpl, psds, norms, t_geo, is_w0=False):
    # A half-second local nuisance offset is searched at native sample
    # resolution.  W0 is capped at t0-0.05 s in geocentric coalescence time.
    lo = -0.5
    hi = 0.5
    if is_w0:
        hi = min(hi, g.T0 - 0.05 - t_geo)
    n = int(round((hi - lo) * g.FS))
    deltas = lo + np.arange(n + 1, dtype=np.float64) / g.FS
    q = {d: math.copysign(norms[d], g.fb(d)) for d in ["H1", "L1", "V1"]}
    denom = norms["network"]
    best = None
    for delta in deltas:
        vals = {
            d: sample(
                zs[d],
                t_geo + float(delta) + g.delay(d, t_geo) - tmpl.peak_offset,
            )
            for d in ["H1", "L1", "V1"]
        }
        coh = abs(sum(q[d] * vals[d] for d in vals)) / denom
        net = math.sqrt(sum(abs(vals[d]) ** 2 for d in vals))
        if best is None or coh > best["coherent_snr"]:
            best = {
                "delta_t_offset_s": float(delta),
                "coherent_snr": float(coh),
                "common_time_network_snr": float(net),
                "complex_snr": {
                    d: [float(vals[d].real), float(vals[d].imag)] for d in vals
                },
                "single_snr_at_common_time": {
                    d: float(abs(vals[d])) for d in vals
                },
            }
    return best


def main():
    tmpl = g.make_template()
    psds = {d: g.prepare_psd(d)[0] for d in ["H1", "L1", "V1"]}
    probe = g.analysis_data(tmpl, "H1", g.T0 - 3300.0)
    norms = g.scalar_norms(tmpl, psds, len(probe))
    opened = datetime.now(timezone.utc).isoformat()
    rows = []
    for target in ORDER:
        name, dt, t_geo, alpha = next(w for w in g.WINDOWS if w[0] == target)
        data = {d: g.analysis_data(tmpl, d, t_geo) for d in ["H1", "L1", "V1"]}
        zs = {d: complex_filter(data[d], tmpl, psds[d]) for d in ["H1", "L1", "V1"]}
        best = coherent_scan(zs, tmpl, psds, norms, t_geo, is_w0=(name == "W0"))
        best["window"] = name
        best["delta_t_s"] = dt
        best["gps_nominal"] = t_geo
        best["alpha_abs"] = alpha
        best["best_delta_t_s"] = g.T0 - best["delta_t_offset_s"] - t_geo
        rows.append(best)
        print(json.dumps(best))
    result = {
        "windows_opened_utc": opened,
        "search_order": ORDER,
        "scalar_norms": norms,
        "rows": rows,
        "statistic": "max over a common geocentric offset of |sum(signed F_b sigma_i z_i)| / sqrt(sum((F_b sigma_i)^2)); phase-maximized.",
        "decision_scope": "Diagnostic only; the corrected frozen candidate table remains the independent network-SNR plus pairwise timing/veto result.",
    }
    OUT.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
