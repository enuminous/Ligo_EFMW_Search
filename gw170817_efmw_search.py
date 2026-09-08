#!/usr/bin/env python3
"""Frozen Channel-A precursor search scaffold for GW170817.

This script is intentionally conservative: every result is marked measured only
after the corresponding injection/search branch completes. It uses the public
GWOSC HDF5 strain files and PyCBC/LALSuite. The source template is a published
GW170817-like IMRPhenomPv2_NRTidal waveform (GWTC-1 low-spin best values).
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from pycbc.detector import Detector
from pycbc.filter import highpass_fir, matched_filter, sigmasq
from pycbc.psd import interpolate, inverse_spectrum_truncation, welch
from pycbc.types import TimeSeries
from pycbc.waveform import get_td_waveform

FS = 4096
T0 = 1187008882.4457
# Fixed to the GW170817 low-spin posterior medians in the public
# GWTC-1 PE file (RA=3.44616 rad, Dec=-0.408084 rad).
RA = 3.44616
DEC = -0.408084
GPS_DATA_START = 1187004416.0
GPS_EVENT_START = 1187006835.0
DATA_DIR = Path(__file__).parent / "gw170817_data"

EVENT_FILES = {
    "H1": DATA_DIR / "H-H1_GWOSC_4KHZ_R1-1187006835-4096.hdf5",
    "L1": DATA_DIR / "L-L1_GWOSC_4KHZ_R1-1187006835-4096.full.hdf5",
    "V1": DATA_DIR / "V-V1_GWOSC_4KHZ_R1-1187006835-4096.full.hdf5",
}
PRE_FILES = {
    "H1": DATA_DIR / "pre/H-H1_GWOSC_O2_4KHZ_R1-1187004416-4096.hdf5",
    "L1": DATA_DIR / "pre/L-L1_GWOSC_O2_4KHZ_R1-1187004416-4096.hdf5",
    "V1": DATA_DIR / "pre/V-V1_GWOSC_O2_4KHZ_R1-1187004416-4096.hdf5",
}

WINDOWS = [
    ("W0", 0.1, 1187008882.3457, 7.0e-9),
    ("W1", 1.7, 1187008880.7457, 2.9e-8),
    ("W2", 10.0, 1187008872.4457, 7.0e-8),
    ("W3", 100.0, 1187008782.4457, 2.2e-7),
    ("W4", 3600.0, 1187005282.4457, 1.3e-6),
]
R_LEVELS = [1.0, 0.3, 0.1, 0.03, 0.01]
SEED = 230908
N_INJECTION_REPS = 5
PRE_PAD = 20.0
POST_PAD = 60.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_hdf(path: Path, start: float, end: float) -> TimeSeries:
    with h5py.File(path, "r") as f:
        ds = f["strain/Strain"]
        x0 = float(ds.attrs["Xstart"])
        dt = float(ds.attrs["Xspacing"])
        i0 = max(0, int(math.floor((start - x0) / dt)))
        i1 = min(len(ds), int(math.ceil((end - x0) / dt)))
        arr = np.asarray(ds[i0:i1], dtype=np.float64)
        return TimeSeries(arr, delta_t=dt, epoch=x0 + i0 * dt)


def read_detector(det: str, start: float, end: float) -> TimeSeries:
    # Prefer the event-release file whenever it contains the requested interval.
    p = EVENT_FILES[det]
    if start >= GPS_EVENT_START and end <= GPS_EVENT_START + 4096:
        return read_hdf(p, start, end)
    p = PRE_FILES[det]
    if start >= GPS_DATA_START and end <= GPS_DATA_START + 4096:
        return read_hdf(p, start, end)
    raise ValueError(f"no downloaded file covers {det} [{start}, {end}]")


def clean(x: TimeSeries) -> TimeSeries:
    # 8-second FIR transient padding on each side.
    y = highpass_fir(x, 15.0, 8 * FS)
    return y.crop(8.0, 8.0)


@dataclass
class Template:
    hp: TimeSeries
    hc: TimeSeries
    peak_offset: float


def make_template() -> Template:
    hp, hc = get_td_waveform(
        approximant="IMRPhenomPv2_NRTidalv2",
        mass1=1.4790595881585529,
        mass2=1.2807580711111188,
        spin1z=0.004953281544706123,
        spin2z=0.0028389018809797492,
        lambda1=269.0083686566036,
        lambda2=445.770904970051,
        distance=39.95566849878731,
        inclination=2.5596915857494036,
        delta_t=1.0 / FS,
        f_lower=30.0,
    )
    # Keep the inspiral through the strain-amplitude maximum. This defines the
    # template's coalescence timestamp and excludes post-merger/ringdown samples.
    i = int(np.argmax(np.abs(np.asarray(hp))))
    hp = TimeSeries(np.asarray(hp[: i + 1]), delta_t=1.0 / FS, epoch=0.0)
    hc = TimeSeries(np.asarray(hc[: i + 1]), delta_t=1.0 / FS, epoch=0.0)
    return Template(hp=hp, hc=hc, peak_offset=(len(hp) - 1) / FS)


def fb(det: str) -> float:
    return float(Detector(det).antenna_pattern(RA, DEC, 0.0, T0, polarization_type="scalar")[0])


def delay(det: str, t_geo: float) -> float:
    return float(Detector(det).time_delay_from_earth_center(RA, DEC, t_geo))


def prepare_psd(det: str) -> tuple:
    # Stationary reference segment, outside the event and outside search windows.
    x = clean(read_detector(det, T0 - 900.0, T0 - 380.0))
    p = welch(x, seg_len=8 * FS, seg_stride=4 * FS, avg_method="median")
    return p, x


def psd_for_data(psd, n: int) -> object:
    p = interpolate(psd, 1.0 / (n / FS))
    return inverse_spectrum_truncation(p, 8 * FS, low_frequency_cutoff=30.0)


def quadrature_snr(data: TimeSeries, tmpl: Template, psd) -> tuple:
    # PyCBC's complex matched-filter output already maximizes over the unknown
    # phase of a real inspiral template. The scalar phase phi0 is therefore
    # represented by abs(z), avoiding double-counting two correlated filters.
    n = len(data)
    p = psd_for_data(psd, n)
    hp = tmpl.hp.copy()
    hp.resize(n)
    z0 = matched_filter(hp, data, psd=p, low_frequency_cutoff=30.0)
    return z0, None, np.abs(z0)


def local_peak(z: TimeSeries, center: float, half_width: float = 0.5) -> tuple:
    lo = max(0, int((center - half_width - float(z.start_time)) * FS))
    hi = min(len(z), int((center + half_width - float(z.start_time)) * FS) + 1)
    if hi <= lo:
        raise ValueError("empty local search interval")
    j = lo + int(np.argmax(np.asarray(z)[lo:hi]))
    return float(z[j]), float(z.start_time + j / FS)


def peak_between(z: TimeSeries, lo_time: float, hi_time: float) -> tuple:
    lo = max(0, int(math.ceil((lo_time - float(z.start_time)) * FS)))
    hi = min(len(z), int(math.floor((hi_time - float(z.start_time)) * FS)) + 1)
    if hi <= lo:
        raise ValueError("empty bounded search interval")
    j = lo + int(np.argmax(np.asarray(z)[lo:hi]))
    return float(z[j]), float(z.start_time + j / FS)


def window_score(det_data: dict, tmpl: Template, psds: dict, t_geo: float) -> dict:
    snrs = {}
    times = {}
    for d in ["H1", "L1", "V1"]:
        _, _, z = quadrature_snr(det_data[d], tmpl, psds[d])
        # The output timestamp is the padded template start. Search at the
        # detector's ordinary-light-travel arrival time, not at a Delta-t-sized
        # inter-site offset.
        score_center = t_geo + delay(d, t_geo) - tmpl.peak_offset
        # W0 credit is capped in coalescence time, not template-start time.
        upper_coal = T0 - 0.05
        lo = score_center - 0.5
        hi = min(score_center + 0.5, upper_coal - delay(d, t_geo) - tmpl.peak_offset)
        peak, t = peak_between(z, lo, hi)
        snrs[d] = peak
        times[d] = t
    active = [d for d in snrs if snrs[d] >= 4.0]
    net = float(math.sqrt(sum(v * v for v in snrs.values())))
    ordinary = (
        len(active) >= 2
        and max(times[d] for d in active) - min(times[d] for d in active) <= 0.030
    )
    return {
        "single": snrs,
        "times": times,
        "network": net,
        "active_ge4": active,
        "network_pass": net >= 8.0,
        "two_detector_pass": len(active) >= 2,
        "coincidence_pass": ordinary,
    }


def analysis_data(tmpl: Template, det: str, t_geo: float, injected_scale: float = 0.0) -> TimeSeries:
    """Read a short real-noise segment and optionally add a frozen-channel-A injection."""
    t_arr = t_geo + delay(det, t_geo)
    start = t_arr - tmpl.peak_offset - PRE_PAD
    end = t_arr + POST_PAD
    raw = read_detector(det, start, end)
    if injected_scale:
        x = np.asarray(raw, dtype=np.float64).copy()
        i0 = int(round((t_arr - tmpl.peak_offset - float(raw.start_time)) * FS))
        sig = injected_scale * fb(det) * np.asarray(tmpl.hp)
        i1 = i0 + len(sig)
        if i0 < 0 or i1 > len(x):
            raise ValueError(f"injection does not fit {det}: {i0}:{i1} / {len(x)}")
        x[i0:i1] += sig
        raw = TimeSeries(x, delta_t=1.0 / FS, epoch=raw.start_time)
    return clean(raw)


def scalar_norms(tmpl: Template, psds: dict, n: int) -> dict:
    h = tmpl.hp.copy()
    h.resize(n)
    hf = h.to_frequencyseries()
    per = {}
    for d in ["H1", "L1", "V1"]:
        norm = math.sqrt(float(sigmasq(hf, psd=psd_for_data(psds[d], n),
                                       low_frequency_cutoff=30.0)))
        per[d] = abs(fb(d)) * norm
    per["network"] = math.sqrt(sum(per[d] ** 2 for d in ["H1", "L1", "V1"]))
    return per


def quality_pass(det: str, start: float, end: float) -> bool:
    """Require DATA and CBC_CAT1 pass bits throughout a segment."""
    p = EVENT_FILES[det]
    if not (start >= GPS_EVENT_START and end <= GPS_EVENT_START + 4096):
        p = PRE_FILES[det]
    with h5py.File(p, "r") as f:
        ds = f["quality/simple/DQmask"]
        x0 = float(ds.attrs["Xstart"])
        i0 = max(0, int(math.floor(start - x0)))
        i1 = min(len(ds), int(math.ceil(end - x0)))
        return bool(np.all((np.asarray(ds[i0:i1]) & 3) == 3))


def single_detector_filter(data: TimeSeries, tmpl: Template, psd) -> TimeSeries:
    h = tmpl.hp.copy()
    h.resize(len(data))
    return matched_filter(h, data, psd=psd_for_data(psd, len(data)),
                          low_frequency_cutoff=30.0)


def main():
    tmpl = make_template()
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
            "duration_s": len(tmpl.hp) / FS,
            "peak_offset_s": tmpl.peak_offset,
        },
        "antenna_fb": {d: fb(d) for d in ["H1", "L1", "V1"]},
        "files": {
            d: {
                "event": str(EVENT_FILES[d]),
                "event_sha256": sha256(EVENT_FILES[d]),
                "pre": str(PRE_FILES[d]),
                "pre_sha256": sha256(PRE_FILES[d]),
            }
            for d in ["H1", "L1", "V1"]
        },
    }
    print("META " + json.dumps(meta))

    psds = {d: prepare_psd(d)[0] for d in ["H1", "L1", "V1"]}
    # Calibrate the dimensionless intrinsic amplitude A to the requested
    # network-equivalent relative amplitude r using rho_tensor ~= 32.
    probe = analysis_data(tmpl, "H1", T0 - 3300.0)
    norms = scalar_norms(tmpl, psds, len(probe))
    rho_tensor = 32.0
    print("NORMS " + json.dumps({"scalar_unit_snr": norms, "rho_tensor_reference": rho_tensor,
                                  "A_per_r": rho_tensor / norms["network"],
                                  "A_at_network_8": 8.0 / norms["network"]}))

    # Frozen injection campaign: real nearby O2 noise, five deterministic
    # coalescence times per window, before opening any search window.
    rng = np.random.default_rng(SEED)
    inj_times = {
        name: sorted(rng.uniform(T0 - 3400.0, T0 - 3000.0, N_INJECTION_REPS).tolist())
        for name, _, _, _ in WINDOWS
    }
    inj_rows = []
    for name, dt, _, _ in WINDOWS:
        for rep, t_geo in enumerate(inj_times[name]):
            for r in R_LEVELS:
                A = r * rho_tensor / norms["network"]
                dscores = {}
                for d in ["H1", "L1", "V1"]:
                    injected = analysis_data(tmpl, d, t_geo, injected_scale=A)
                    z = abs(single_detector_filter(injected, tmpl, psds[d]))
                    center = t_geo + delay(d, t_geo) - tmpl.peak_offset
                    peak, tp = local_peak(z, center, 0.5)
                    dscores[d] = {"snr": peak, "time": tp}
                active = [d for d in dscores if dscores[d]["snr"] >= 4.0]
                net = math.sqrt(sum(dscores[d]["snr"] ** 2 for d in dscores))
                coinc = len(active) >= 2 and max(dscores[d]["time"] for d in active) - min(
                    dscores[d]["time"] for d in active
                ) <= 0.030
                inj_rows.append({
                    "window": name, "delta_t_s": dt, "rep": rep, "t_geo": t_geo,
                    "r": r, "A": A, "single_snr": {d: dscores[d]["snr"] for d in dscores},
                    "network_snr": net, "pass": net >= 8 and len(active) >= 2 and coinc,
                    "quality_pass": all(quality_pass(d, t_geo - tmpl.peak_offset - PRE_PAD,
                                                      t_geo + POST_PAD) for d in ["H1", "L1", "V1"]),
                })
    print("INJECTIONS " + json.dumps(inj_rows))

    # Search order is explicitly W1, W2, W3, W0, W4.
    search_rows = []
    for target in ["W1", "W2", "W3", "W0", "W4"]:
        name, dt, t_geo, alpha = next(w for w in WINDOWS if w[0] == target)
        data = {d: analysis_data(tmpl, d, t_geo) for d in ["H1", "L1", "V1"]}
        row = window_score(data, tmpl, psds, t_geo)
        row.update({
            "window": name, "delta_t_s": dt, "gps_nominal": t_geo, "alpha_abs": alpha,
            "quality_pass": all(quality_pass(d, t_geo - tmpl.peak_offset - PRE_PAD,
                                             t_geo + POST_PAD) for d in ["H1", "L1", "V1"]),
            "best_delta_t_s": T0 - np.mean([
                row["times"][d] + tmpl.peak_offset - delay(d, t_geo)
                for d in row["times"]
            ]),
        })
        search_rows.append(row)
        print("SEARCH " + json.dumps(row))

    result = {
        "metadata": meta,
        "scalar_norms": norms,
        "injection_times": inj_times,
        "injections": inj_rows,
        "search": search_rows,
        "seed": SEED,
    }
    Path("gw170817_efmw_results.json").write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
