#!/usr/bin/env python3
"""Final frozen-window pass with PyCBC chi-square and autocorrelation vetoes."""
import json
import math
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from pycbc.filter import matched_filter, matched_filter_core
from pycbc.vetoes.chisq import power_chisq
from pycbc.vetoes.autochisq import SingleDetAutoChisq
import gw170817_efmw_search as g

CHI2_BINS = 16
CHI2_REDUCED_MAX = 5.0
AUTO_STRIDE = 16
AUTO_POINTS = 8
AUTO_REDUCED_MAX = 5.0


def detector_stat(data, tmpl, psd, det, t_geo):
    h = tmpl.hp.copy()
    h.resize(len(data))
    pp = g.psd_for_data(psd, len(data))
    z = abs(matched_filter(h, data, psd=pp, low_frequency_cutoff=30.0))
    center = t_geo + g.delay(det, t_geo) - tmpl.peak_offset
    upper = g.T0 - 0.05 - g.delay(det, t_geo) - tmpl.peak_offset
    lo = center - 0.5
    hi = min(center + 0.5, upper)
    peak, tp = g.peak_between(z, lo, hi)
    idx = int(round((tp - float(z.start_time)) * g.FS))
    chi = power_chisq(h, data, CHI2_BINS, pp, low_frequency_cutoff=30.0)
    chi_red = float(chi[idx]) / (2 * CHI2_BINS - 2)
    raw, _, norm = matched_filter_core(h, data, psd=pp, low_frequency_cutoff=30.0)
    auto = SingleDetAutoChisq(
        stride=AUTO_STRIDE,
        num_points=AUTO_POINTS,
        onesided=None,
        twophase=True,
        take_maximum_value=False,
    )
    aval, adof = auto.values(
        raw, np.array([idx]), h, pp, norm, low_frequency_cutoff=30.0
    )
    auto_red = float(aval[0]) / adof
    return {
        "snr": peak,
        "time": tp,
        "chi2": float(chi[idx]),
        "chi2_reduced": chi_red,
        "autochisq": float(aval[0]),
        "autochisq_dof": int(adof),
        "autochisq_reduced": auto_red,
        "chi2_pass": chi_red <= CHI2_REDUCED_MAX,
        "autochisq_pass": auto_red <= AUTO_REDUCED_MAX,
        "veto_pass": chi_red <= CHI2_REDUCED_MAX and auto_red <= AUTO_REDUCED_MAX,
    }


def main():
    tmpl = g.make_template()
    psds = {d: g.prepare_psd(d)[0] for d in ["H1", "L1", "V1"]}
    rows = []
    # Persist the moment the frozen search windows are opened.  This is kept
    # separate from the later result-file completion time for auditability.
    windows_opened_utc = datetime.now(timezone.utc).isoformat()
    for target in ["W1", "W2", "W3", "W0", "W4"]:
        name, dt, t_geo, alpha = next(w for w in g.WINDOWS if w[0] == target)
        data = {d: g.analysis_data(tmpl, d, t_geo) for d in ["H1", "L1", "V1"]}
        single = {d: detector_stat(data[d], tmpl, psds[d], d, t_geo)
                  for d in ["H1", "L1", "V1"]}
        times = {d: single[d]["time"] for d in single}
        active = [d for d in single if single[d]["snr"] >= 4.0]
        net = math.sqrt(sum(single[d]["snr"] ** 2 for d in single))
        ordinary = len(active) >= 2 and max(times[d] for d in active) - min(
            times[d] for d in active
        ) <= 0.030
        veto = all(single[d]["veto_pass"] for d in active)
        row = {
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
            "coincidence_pass": ordinary,
            "chi2_autochisq_veto_pass": veto,
            "candidate_pass": net >= 8.0 and len(active) >= 2 and ordinary and veto,
            "quality_pass": all(
                g.quality_pass(
                    d,
                    t_geo + g.delay(d, t_geo) - tmpl.peak_offset - g.PRE_PAD,
                    t_geo + g.delay(d, t_geo) + g.POST_PAD,
                )
                for d in ["H1", "L1", "V1"]
            ),
            "best_delta_t_s": g.T0 - np.mean([
                times[d] + tmpl.peak_offset - g.delay(d, t_geo)
                for d in times
            ]),
        }
        rows.append(row)
        print(json.dumps(row))
    p = Path("gw170817_efmw_results.json")
    old = json.loads(p.read_text())
    old["search"] = rows
    old["search_order"] = ["W1", "W2", "W3", "W0", "W4"]
    old["windows_opened_utc"] = windows_opened_utc
    old["veto_settings"] = {
        "power_chisq_bins": CHI2_BINS,
        "chi2_reduced_max": CHI2_REDUCED_MAX,
        "autochisq_stride_samples": AUTO_STRIDE,
        "autochisq_points": AUTO_POINTS,
        "autochisq_reduced_max": AUTO_REDUCED_MAX,
    }
    p.write_text(json.dumps(old, indent=2))


if __name__ == "__main__":
    main()
