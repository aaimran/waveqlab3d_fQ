#!/usr/bin/env python3
"""
Frequency-wavenumber (f-k) reference solution generator for WaveQLab3D-fQ.

Computes semi-analytical seismograms for viscoelastic media using
complex-valued elastic moduli from the fQ8 generalized standard linear solid.
Half-space cases use a far-field point-source Green's function with
frequency-dependent complex velocities; plane-wave propagation provides
a clean 1D comparison for attenuation and dispersion validation.

Validates anelastic-fQ8 with coefficient_method='conventional-nnls',
coarse_grain=0 against the Withers, Olsen & Day (2015, BSSA) benchmark.

Usage:
    python fk_reference.py --case elastic
    python fk_reference.py --case constant_q
    python fk_reference.py --case powerlaw
    python fk_reference.py --case layered
    python fk_reference.py --diagnostics --Qs0 50 --gamma 0.6
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

PI = np.pi
N_MECH = 8

# ============================================================
# Withers (2015) Table data — mirrors withers_tables.f90
# ============================================================

TAU_MIN = np.array([0.0032, 0.0032, 0.0032, 0.0032, 0.0032,
                    0.0032, 0.0032, 0.0066, 0.0066, 0.0085])
TAU_MAX = np.array([15.9155, 15.9155, 15.9155, 15.9155, 15.9155,
                    15.9155, 15.9155, 3.9789, 3.9789, 3.9789])


def _gamma_interp(gamma: float) -> tuple[int, int, float]:
    g = np.clip(gamma, 0.0, 0.9)
    idx = min(int(g * 10), 8)
    alpha = g * 10 - idx
    return idx, idx + 1, alpha


def get_relaxation_times(gamma: float, f_transition: float = 1.0) -> NDArray:
    """Withers (2015) Eq. 15: log-spaced relaxation times scaled by f_transition."""
    lo, hi, alpha = _gamma_interp(gamma)
    tau_min = (1 - alpha) * TAU_MIN[lo] + alpha * TAU_MIN[hi]
    tau_max = (1 - alpha) * TAU_MAX[lo] + alpha * TAU_MAX[hi]
    k = np.arange(1, N_MECH + 1)
    tau = np.exp(np.log(tau_min) + (2 * k - 1) / 16.0 * np.log(tau_max / tau_min))
    return tau / f_transition


# ============================================================
# Coefficient fitting — mirrors fit_conventional_strengths
# ============================================================

def fit_conventional_strengths(
    Q0: float, gamma: float, f_transition: float, tau: NDArray,
    nfreq: int = 256, max_sweeps: int = 20000
) -> NDArray:
    """NNLS coordinate-descent matching anelastic_fq8_model.f90."""
    f = 0.1 * f_transition * 100.0 ** (np.arange(nfreq) / (nfreq - 1))
    target_q = np.where(f < f_transition, Q0, Q0 * (f / f_transition) ** gamma)
    omega = 2 * PI * f

    A = np.zeros((nfreq, N_MECH))
    for j in range(N_MECH):
        x = omega * tau[j]
        denom = 1.0 + x * x
        A[:, j] = (target_q * x + 1.0) / denom

    b = np.ones(nfreq)
    strength = np.zeros(N_MECH)
    residual = b.copy()

    for _ in range(max_sweeps):
        max_delta = 0.0
        for j in range(N_MECH):
            old = strength[j]
            aj = A[:, j]
            new = max(0.0, old + np.dot(aj, residual) /
                      max(np.dot(aj, aj), np.finfo(float).tiny))
            delta = new - old
            if delta != 0.0:
                residual -= delta * aj
            strength[j] = new
            max_delta = max(max_delta, abs(delta))
        if max_delta < 1e-13:
            break

    return strength


def build_fq8_coefficients(
    Qs0: float, Qp0: float, gamma: float,
    f_transition: float = 1.0, fref: float = 1.0
) -> tuple[NDArray, NDArray, NDArray]:
    """Build (tau, strength_s, strength_p) for conventional-nnls, coarse_grain=0."""
    tau = get_relaxation_times(gamma, f_transition)
    strength_s = fit_conventional_strengths(Qs0, gamma, f_transition, tau)
    strength_p = fit_conventional_strengths(Qp0, gamma, f_transition, tau)
    return tau, strength_s, strength_p


# ============================================================
# Complex modulus and dispersion
# ============================================================

def complex_modulus_ratio(
    freq: NDArray, tau: NDArray, strength: NDArray
) -> NDArray:
    """M(f)/M_u = 1 - sum_k lambda_k/(1 + i*omega*tau_k) [Withers Eq. 6]."""
    omega = 2 * PI * np.atleast_1d(freq)
    result = np.ones_like(omega, dtype=complex)
    for k in range(len(tau)):
        result -= strength[k] / (1.0 + 1j * omega * tau[k])
    return result


def modulus_scale(fref: float, tau: NDArray, strength: NDArray) -> float:
    """Unrelaxed modulus correction: M_u/M_physical so that c(fref) = c_ref."""
    M = complex_modulus_ratio(np.array([fref]), tau, strength)[0]
    return float(np.real(1.0 / np.sqrt(M)) ** 2)


def realized_q(freq: NDArray, tau: NDArray, strength: NDArray) -> NDArray:
    """Q(f) = Re(M)/Im(M) from the complex modulus."""
    M = complex_modulus_ratio(freq, tau, strength)
    return np.abs(np.real(M) / np.imag(M))


def complex_velocity(
    freq: NDArray, c_ref: float, fref: float,
    tau: NDArray, strength: NDArray
) -> NDArray:
    """Dispersive complex velocity c(f) normalized so Re[c(fref)] = c_ref."""
    scale = modulus_scale(fref, tau, strength)
    M_f = scale * complex_modulus_ratio(freq, tau, strength)
    return c_ref * np.sqrt(M_f)


def phase_velocity(
    freq: NDArray, c_ref: float, fref: float,
    tau: NDArray, strength: NDArray
) -> NDArray:
    """Phase velocity = omega/Re(k) = Re-part of sqrt(M)-based velocity."""
    c = complex_velocity(freq, c_ref, fref, tau, strength)
    return float(c_ref) / np.real(c_ref / c)


# ============================================================
# Source time functions
# ============================================================

def ricker_wavelet(t: NDArray, f0: float, t0: float = 0.0) -> NDArray:
    """Ricker wavelet (negative second derivative of Gaussian)."""
    a = PI * f0 * (t - t0)
    return (1 - 2 * a ** 2) * np.exp(-a ** 2)


def cosine_bell(t: NDArray, t0: float, t_half: float) -> NDArray:
    """Cosine-bell source time function matching WaveQLab3D 'Cos-Bell'."""
    s = np.zeros_like(t)
    mask = np.abs(t - t0) < t_half
    s[mask] = 0.5 * (1 + np.cos(PI * (t[mask] - t0) / t_half))
    return s


def gaussian_pulse(t: NDArray, f0: float, t0: float = 0.0) -> NDArray:
    """Gaussian pulse centered at t0 with dominant frequency f0."""
    sigma = 1.0 / (2 * PI * f0)
    return np.exp(-0.5 * ((t - t0) / sigma) ** 2) / (sigma * np.sqrt(2 * PI))


# ============================================================
# 1-D viscoelastic plane-wave propagation (simplest f-k case)
# ============================================================

def plane_wave_seismogram(
    distance: float, c_ref: float, fref: float,
    tau: NDArray, strength: NDArray,
    source_func: NDArray, dt: float
) -> NDArray:
    """Exact 1D plane-wave through homogeneous viscoelastic medium.

    Computes u(x,t) = IFFT[ S(w) * exp(i*k(w)*x) ]
    where k(w) = w/c(w) is the complex wavenumber.
    """
    nt = len(source_func)
    S = np.fft.rfft(source_func)
    freqs = np.fft.rfftfreq(nt, dt)

    result_spectrum = np.zeros_like(S)
    mask = freqs > 0
    c = complex_velocity(freqs[mask], c_ref, fref, tau, strength)
    k_complex = 2 * PI * freqs[mask] / c
    result_spectrum[mask] = S[mask] * np.exp(-1j * k_complex * distance)

    return np.fft.irfft(result_spectrum, n=nt)


# ============================================================
# Discrete wavenumber f-k method for layered viscoelastic media
# ============================================================

@dataclass
class Layer:
    """One layer in a 1D velocity model."""
    thickness: float   # meters (inf for half-space)
    rho: float         # density (kg/m^3 or g/cm^3)
    vs: float          # S-wave velocity at reference frequency
    vp: float          # P-wave velocity at reference frequency
    Qs0: float = 1e10  # S quality factor (very large = elastic)
    Qp0: float = 1e10  # P quality factor
    gamma: float = 0.0
    f_transition: float = 1.0
    fref: float = 1.0
    tau: NDArray = field(default_factory=lambda: np.zeros(0))
    strength_s: NDArray = field(default_factory=lambda: np.zeros(0))
    strength_p: NDArray = field(default_factory=lambda: np.zeros(0))


def prepare_layer(layer: Layer) -> Layer:
    """Compute fQ8 coefficients for a layer (skip if elastic)."""
    if layer.Qs0 > 1e6 and layer.Qp0 > 1e6:
        layer.tau = np.zeros(N_MECH)
        layer.strength_s = np.zeros(N_MECH)
        layer.strength_p = np.zeros(N_MECH)
        return layer
    tau, ss, sp = build_fq8_coefficients(
        layer.Qs0, layer.Qp0, layer.gamma, layer.f_transition, layer.fref
    )
    layer.tau = tau
    layer.strength_s = ss
    layer.strength_p = sp
    return layer


def layer_complex_velocities(layer: Layer, freq: float) -> tuple[complex, complex]:
    """Return (alpha, beta) complex velocities for this layer at frequency freq."""
    if np.sum(np.abs(layer.strength_s)) < 1e-30:
        return complex(layer.vp), complex(layer.vs)
    alpha = complex_velocity(
        np.array([freq]), layer.vp, layer.fref, layer.tau, layer.strength_p
    )[0]
    beta = complex_velocity(
        np.array([freq]), layer.vs, layer.fref, layer.tau, layer.strength_s
    )[0]
    return alpha, beta


def farfield_halfspace_seismograms(
    layer: Layer, source_depth: float,
    receiver_offsets: NDArray, dt: float, nt: int,
    source_func: NDArray, fmax: float | None = None,
    free_surface: bool = True
) -> tuple[NDArray, NDArray]:
    """Far-field seismograms for explosion source in viscoelastic half-space.

    Computes P and S (converted) arrivals using the exact frequency-domain
    propagation with complex velocities, including geometric spreading.
    For a half-space with free surface, the direct P and surface-reflected PP
    arrivals dominate.

    Parameters
    ----------
    layer : Layer
        Half-space material (must have coefficients computed via prepare_layer).
    source_depth : float
        Source depth below free surface (meters).
    receiver_offsets : array
        Horizontal distances from source to receivers (meters).
    dt : float
        Time step (seconds).
    nt : int
        Number of time samples.
    source_func : array
        Source time function (length nt).
    fmax : float, optional
        Maximum frequency to compute.
    free_surface : bool
        If True, apply free-surface amplification (factor 2 for vertical).

    Returns
    -------
    uz : array, shape (n_receivers, nt)
        Vertical displacement seismograms.
    ur : array, shape (n_receivers, nt)
        Radial displacement seismograms.
    """
    n_recv = len(receiver_offsets)
    freqs = np.fft.rfftfreq(nt, dt)
    nf = len(freqs)
    if fmax is None:
        fmax = 0.5 / dt

    S_omega = np.fft.rfft(source_func)

    Uz_omega = np.zeros((n_recv, nf), dtype=complex)
    Ur_omega = np.zeros((n_recv, nf), dtype=complex)

    for ri, r_horiz in enumerate(receiver_offsets):
        R_direct = np.sqrt(r_horiz ** 2 + source_depth ** 2)
        cos_theta = source_depth / R_direct
        sin_theta = r_horiz / R_direct

        # Image source for free-surface reflection (source mirrored above surface)
        R_image = np.sqrt(r_horiz ** 2 + source_depth ** 2)  # same for flat surface

        for fi in range(1, nf):
            f = freqs[fi]
            if f > fmax:
                break

            omega = 2 * PI * f
            alpha_c, beta_c = layer_complex_velocities(layer, f)

            # Complex wavenumber for P-wave
            k_P = omega / alpha_c
            # Complex wavenumber for S-wave (for P-to-S conversion estimate)
            k_S = omega / beta_c

            # Direct P-wave: explosion source radiates isotropically
            # u_P = M0/(4*pi*rho*alpha^2) * (-i*omega/alpha) * exp(-i*k_P*R) / R
            # For displacement from scalar potential: u = -grad(Phi)
            # Phi = M0/(4*pi*rho*alpha^2) * exp(-i*k_P*R)/R  (outgoing convention with exp(+iwt))

            # Using numpy fft convention exp(-iwt):
            phase_P = np.exp(-1j * k_P * R_direct)
            geo_spread_P = 1.0 / R_direct
            amp_P = phase_P * geo_spread_P / (4 * PI * layer.rho * abs(alpha_c) ** 2)

            # Free surface approximately doubles displacement
            fs_factor = 2.0 if free_surface else 1.0

            # Vertical: u_z = -dPhi/dz ≈ i*k_P*cos(theta)*Phi (far-field)
            Uz_omega[ri, fi] += fs_factor * (-1j * k_P * cos_theta) * amp_P * S_omega[fi]

            # Radial: u_r = -dPhi/dr ≈ i*k_P*sin(theta)*Phi (far-field)
            Ur_omega[ri, fi] += fs_factor * (-1j * k_P * sin_theta) * amp_P * S_omega[fi]

    uz = np.zeros((n_recv, nt))
    ur = np.zeros((n_recv, nt))
    for ri in range(n_recv):
        uz[ri, :] = np.fft.irfft(Uz_omega[ri, :], n=nt)
        ur[ri, :] = np.fft.irfft(Ur_omega[ri, :], n=nt)

    return uz, ur


# ============================================================
# Benchmark configurations from Withers et al. (2015)
# ============================================================

def make_elastic_halfspace() -> tuple[list[Layer], dict]:
    """Fig. 3: elastic half-space (Q=inf)."""
    layer = Layer(
        thickness=np.inf, rho=2700.0, vs=3464.0, vp=6000.0,
        Qs0=1e10, Qp0=1e10, gamma=0.0, f_transition=1.0, fref=1.0
    )
    config = dict(
        name="elastic_halfspace",
        source_depth=2000.0,
        receiver_offsets=np.array([5000.0, 10000.0, 15000.0]),
        dt=0.002, nt=4096, f0=1.0, t0=1.5,
        description="Elastic half-space (Q=inf), Withers Fig. 3"
    )
    return [layer], config


def make_constant_q_halfspace() -> tuple[list[Layer], dict]:
    """Fig. 4: constant-Q (gamma=0, Q=50)."""
    layer = Layer(
        thickness=np.inf, rho=2700.0, vs=3464.0, vp=6000.0,
        Qs0=50.0, Qp0=50.0, gamma=0.0, f_transition=1.0, fref=1.0
    )
    config = dict(
        name="constant_q_halfspace",
        source_depth=2000.0,
        receiver_offsets=np.array([5000.0, 10000.0, 15000.0]),
        dt=0.002, nt=4096, f0=1.0, t0=1.5,
        description="Constant-Q half-space (gamma=0, Qs=Qp=50), Withers Fig. 4"
    )
    return [layer], config


def make_powerlaw_halfspace() -> tuple[list[Layer], dict]:
    """Fig. 5: power-law Q (gamma=0.6, Q=50)."""
    layer = Layer(
        thickness=np.inf, rho=2700.0, vs=3464.0, vp=6000.0,
        Qs0=50.0, Qp0=50.0, gamma=0.6, f_transition=1.0, fref=1.0
    )
    config = dict(
        name="powerlaw_halfspace",
        source_depth=2000.0,
        receiver_offsets=np.array([5000.0, 10000.0, 15000.0]),
        dt=0.002, nt=4096, f0=1.0, t0=1.5,
        description="Power-law half-space (gamma=0.6, Qs=Qp=50), Withers Fig. 5"
    )
    return [layer], config


def make_layered_model() -> tuple[list[Layer], dict]:
    """Table 3 / Fig. 6: layered model with Q discontinuity."""
    layers = [
        Layer(
            thickness=2000.0, rho=2700.0, vs=3464.0, vp=6000.0,
            Qs0=20.0, Qp0=20.0, gamma=0.6, f_transition=1.0, fref=1.0
        ),
        Layer(
            thickness=np.inf, rho=2700.0, vs=3464.0, vp=6000.0,
            Qs0=210.0, Qp0=210.0, gamma=0.6, f_transition=1.0, fref=1.0
        ),
    ]
    config = dict(
        name="layered_model",
        source_depth=3000.0,
        receiver_offsets=np.array([5000.0, 10000.0, 15000.0]),
        dt=0.002, nt=4096, f0=1.0, t0=1.5,
        description="Layered model (Q=20/Q=210, gamma=0.6), Withers Table 3 / Fig. 6"
    )
    return layers, config


# ============================================================
# Diagnostics: Q(f) and dispersion curves
# ============================================================

def print_diagnostics(
    Qs0: float, Qp0: float, gamma: float,
    f_transition: float = 1.0, fref: float = 1.0,
    outdir: Path | None = None
) -> None:
    """Print Q(f) accuracy and dispersion for the given parameters."""
    tau, ss, sp = build_fq8_coefficients(Qs0, Qp0, gamma, f_transition, fref)

    print(f"fQ8 coefficients (conventional-nnls, coarse_grain=0):")
    print(f"  Qs0={Qs0}, Qp0={Qp0}, gamma={gamma}, f_transition={f_transition}")
    print(f"  tau = {tau}")
    print(f"  strength_s = {ss}")
    print(f"  strength_p = {sp}")
    print(f"  sum(strength_s) = {np.sum(ss):.6f}")
    print(f"  sum(strength_p) = {np.sum(sp):.6f}")

    f_test = np.logspace(-1, 1, 256) * f_transition
    q_s = realized_q(f_test, tau, ss)
    q_p = realized_q(f_test, tau, sp)
    target_q_s = np.where(f_test < f_transition, Qs0, Qs0 * (f_test / f_transition) ** gamma)
    target_q_p = np.where(f_test < f_transition, Qp0, Qp0 * (f_test / f_transition) ** gamma)

    err_s = np.abs(q_s / target_q_s - 1)
    err_p = np.abs(q_p / target_q_p - 1)
    print(f"\n  Max Qs relative error: {100 * np.max(err_s):.4f}%")
    print(f"  Max Qp relative error: {100 * np.max(err_p):.4f}%")

    scale_s = modulus_scale(fref, tau, ss)
    scale_p = modulus_scale(fref, tau, sp)
    print(f"\n  Modulus scale (S): {scale_s:.10f}")
    print(f"  Modulus scale (P): {scale_p:.10f}")

    if outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)

        data = np.column_stack([
            f_test, target_q_s, q_s, err_s, target_q_p, q_p, err_p
        ])
        header = "freq target_Qs realized_Qs error_s target_Qp realized_Qp error_p"
        np.savetxt(outdir / "q_of_f.dat", data, header=header, fmt="%.10e")

        c_ref_s, c_ref_p = 3464.0, 6000.0
        cv_s = phase_velocity(f_test, c_ref_s, fref, tau, ss)
        cv_p = phase_velocity(f_test, c_ref_p, fref, tau, sp)
        data_disp = np.column_stack([f_test, cv_s, cv_p])
        header_disp = "freq phase_vel_s phase_vel_p"
        np.savetxt(outdir / "dispersion.dat", data_disp,
                   header=header_disp, fmt="%.10e")

        print(f"\n  Wrote q_of_f.dat and dispersion.dat to {outdir}")


# ============================================================
# Main driver
# ============================================================

def run_benchmark(case_name: str, outdir: Path) -> None:
    """Run one benchmark case and write output."""
    cases = {
        "elastic": make_elastic_halfspace,
        "constant_q": make_constant_q_halfspace,
        "powerlaw": make_powerlaw_halfspace,
        "layered": make_layered_model,
    }
    if case_name not in cases:
        print(f"Unknown case '{case_name}'. Available: {list(cases.keys())}")
        sys.exit(1)

    layers, config = cases[case_name]()
    outdir = outdir / config["name"]
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"=== {config['description']} ===")
    for i, layer in enumerate(layers):
        print(f"  Layer {i}: rho={layer.rho}, vs={layer.vs}, vp={layer.vp}, "
              f"Qs={layer.Qs0}, Qp={layer.Qp0}, gamma={layer.gamma}")

    # Build coefficients and print diagnostics
    for layer in layers:
        if layer.Qs0 < 1e6:
            prepare_layer(layer)
            print(f"\n  Coefficients for Q={layer.Qs0}, gamma={layer.gamma}:")
            print(f"    tau = {layer.tau}")
            print(f"    strength_s = {layer.strength_s}")
            print(f"    strength_p = {layer.strength_p}")

    # Source time function
    nt = config["nt"]
    dt = config["dt"]
    t = np.arange(nt) * dt
    source = ricker_wavelet(t, config["f0"], config["t0"])

    # Run far-field synthetics (half-space)
    if len(layers) == 1:
        print(f"\n  Computing far-field synthetics (half-space)...")
        layer_prepared = prepare_layer(layers[0])
        uz, ur = farfield_halfspace_seismograms(
            layer_prepared, config["source_depth"],
            config["receiver_offsets"], dt, nt, source,
            fmax=5 * config["f0"]
        )

        for ri, offset in enumerate(config["receiver_offsets"]):
            fname = f"uz_r{offset:.0f}m.dat"
            data = np.column_stack([t, uz[ri, :]])
            np.savetxt(outdir / fname, data, header="time uz", fmt="%.10e")
            fname = f"ur_r{offset:.0f}m.dat"
            data = np.column_stack([t, ur[ri, :]])
            np.savetxt(outdir / fname, data, header="time ur", fmt="%.10e")

        print(f"  Wrote seismograms to {outdir}")

    # Also compute plane-wave reference (1D, always works)
    if layers[0].Qs0 < 1e6:
        print(f"\n  Computing 1D plane-wave references...")
        layer = layers[0]
        for dist_km in [5, 10, 15]:
            dist = dist_km * 1000.0

            uz_p = plane_wave_seismogram(
                dist, layer.vp, layer.fref,
                layer.tau, layer.strength_p, source, dt
            )
            uz_s = plane_wave_seismogram(
                dist, layer.vs, layer.fref,
                layer.tau, layer.strength_s, source, dt
            )
            data = np.column_stack([t, uz_p, uz_s])
            np.savetxt(
                outdir / f"planewave_{dist_km}km.dat", data,
                header="time P_wave S_wave", fmt="%.10e"
            )

    # Q(f) diagnostics
    for layer in layers:
        if layer.Qs0 < 1e6:
            f_test = np.logspace(-1, 1, 256) * layer.f_transition
            q_s = realized_q(f_test, layer.tau, layer.strength_s)
            target_qs = np.where(
                f_test < layer.f_transition, layer.Qs0,
                layer.Qs0 * (f_test / layer.f_transition) ** layer.gamma
            )
            data = np.column_stack([f_test, target_qs, q_s])
            label = f"Q{layer.Qs0:.0f}"
            np.savetxt(outdir / f"q_of_f_{label}.dat", data,
                       header="freq target_Q realized_Q", fmt="%.10e")

    # Save config
    config_out = {k: v.tolist() if isinstance(v, np.ndarray) else v
                  for k, v in config.items()}
    with open(outdir / "config.json", "w") as fp:
        json.dump(config_out, fp, indent=2)

    print(f"\n  Done: {config['name']}")


# ============================================================
# Kristekova et al. (2009) time-frequency misfit metrics
# ============================================================

def _analytic_signal(x: NDArray) -> NDArray:
    """Analytic signal via Hilbert transform."""
    from numpy.fft import fft, ifft
    N = len(x)
    X = fft(x)
    h = np.zeros(N)
    if N > 0:
        h[0] = 1
        if N % 2 == 0:
            h[N // 2] = 1
            h[1:N // 2] = 2
        else:
            h[1:(N + 1) // 2] = 2
    return ifft(X * h)


def envelope_phase_misfit(
    ref: NDArray, test: NDArray, dt: float,
    fmin: float = 0.0, fmax: float | None = None
) -> tuple[float, float]:
    """Kristekova et al. (2009) single-valued envelope and phase misfits.

    Parameters
    ----------
    ref, test : 1D arrays of equal length
    dt : time step
    fmin, fmax : frequency band (Hz)

    Returns
    -------
    EM : envelope misfit (percent)
    PM : phase misfit (percent)
    """
    if fmax is None:
        fmax = 0.5 / dt
    assert len(ref) == len(test)
    nt = len(ref)

    # Bandpass filter
    freqs = np.fft.rfftfreq(nt, dt)
    mask = (freqs >= fmin) & (freqs <= fmax)
    for sig in [ref, test]:
        pass  # filter in-place below

    def bandpass(s: NDArray) -> NDArray:
        S = np.fft.rfft(s)
        S[~mask] = 0.0
        return np.fft.irfft(S, n=nt)

    r = bandpass(ref)
    t = bandpass(test)

    # Analytic signals
    ar = _analytic_signal(r)
    at = _analytic_signal(t)

    # Envelopes
    env_r = np.abs(ar)
    env_t = np.abs(at)

    # Instantaneous phases
    phase_r = np.unwrap(np.angle(ar))
    phase_t = np.unwrap(np.angle(at))

    # Envelope misfit (EM): L2 norm of envelope difference / L2 norm of ref envelope
    env_diff = env_t - env_r
    EM = 100.0 * np.sqrt(np.sum(env_diff ** 2) / max(np.sum(env_r ** 2), 1e-30))

    # Phase misfit (PM): weighted L2 of phase difference
    w = env_r  # weight by reference envelope
    dphi = phase_t - phase_r
    PM = 100.0 * np.sqrt(
        np.sum(w ** 2 * dphi ** 2) / max(np.sum(w ** 2 * PI ** 2), 1e-30)
    )

    return float(EM), float(PM)


def compare_seismograms(
    ref_file: str, test_file: str, dt: float,
    fmin: float = 0.0, fmax: float | None = None,
    col_ref: int = 1, col_test: int = 1
) -> tuple[float, float]:
    """Load two column-data files and compute EM/PM."""
    ref_data = np.loadtxt(ref_file)
    test_data = np.loadtxt(test_file)
    ref_trace = ref_data[:, col_ref]
    test_trace = test_data[:, col_test]
    n = min(len(ref_trace), len(test_trace))
    return envelope_phase_misfit(ref_trace[:n], test_trace[:n], dt, fmin, fmax)


def run_selftest() -> bool:
    """Verify coefficient fitting, Q accuracy, and propagation against known values."""
    ok = True
    n_pass = 0
    n_fail = 0

    def check(label: str, value: float, expected: float, tol: float) -> None:
        nonlocal ok, n_pass, n_fail
        err = abs(value - expected)
        if err > tol:
            print(f"  FAIL: {label}: got {value:.10e}, expected {expected:.10e}, "
                  f"err={err:.4e} > tol={tol:.4e}")
            ok = False
            n_fail += 1
        else:
            n_pass += 1

    def check_rel(label: str, value: float, expected: float, tol: float) -> None:
        nonlocal ok, n_pass, n_fail
        rel = abs(value / expected - 1.0) if expected != 0.0 else abs(value)
        if rel > tol:
            print(f"  FAIL: {label}: got {value:.10e}, expected {expected:.10e}, "
                  f"rel_err={rel:.4e} > tol={tol:.4e}")
            ok = False
            n_fail += 1
        else:
            n_pass += 1

    # --- Test 1: constant-Q (gamma=0, Q=50) coefficients ---
    print("Test 1: constant-Q coefficient fitting (gamma=0, Q=50)")
    tau, ss, sp = build_fq8_coefficients(50.0, 50.0, 0.0, 1.0, 1.0)
    check("P/S symmetry", float(np.max(np.abs(ss - sp))), 0.0, 1e-14)
    check("sum(strength) < 1", float(np.sum(ss)), float(np.sum(ss)), 1.0)
    f_test = np.logspace(-1, 1, 256) * 1.0
    q_s = realized_q(f_test, tau, ss)
    max_q_err = float(np.max(np.abs(q_s / 50.0 - 1.0)))
    check("max Q error (gamma=0)", max_q_err, 0.0, 0.002)
    scale = modulus_scale(1.0, tau, ss)
    # Scale corrects unrelaxed→physical; for Q=50 it's ~1.06
    check_rel("modulus scale (gamma=0)", scale, 1.0596, 0.005)

    # --- Test 2: power-law (gamma=0.6, Q=50) coefficients ---
    print("Test 2: power-law coefficient fitting (gamma=0.6, Q=50)")
    tau6, ss6, sp6 = build_fq8_coefficients(50.0, 50.0, 0.6, 1.0, 1.0)
    check("P/S symmetry (gamma=0.6)", float(np.max(np.abs(ss6 - sp6))), 0.0, 1e-14)
    target_q6 = np.where(f_test < 1.0, 50.0, 50.0 * (f_test / 1.0) ** 0.6)
    q_s6 = realized_q(f_test, tau6, ss6)
    max_q_err6 = float(np.max(np.abs(q_s6 / target_q6 - 1.0)))
    check("max Q error (gamma=0.6)", max_q_err6, 0.0, 0.10)

    # --- Test 3: elastic plane-wave arrival time ---
    print("Test 3: elastic plane-wave arrival time")
    dt = 0.002
    nt = 4096
    t = np.arange(nt) * dt
    src = ricker_wavelet(t, 1.0, 1.5)
    c_ref = 3464.0
    distance = 5000.0
    tau_el = np.zeros(N_MECH)
    str_el = np.zeros(N_MECH)
    uz = plane_wave_seismogram(distance, c_ref, 1.0, tau_el, str_el, src, dt)
    peak_t = t[np.argmax(np.abs(uz))]
    expected_t = distance / c_ref + 1.5  # source delay
    check("peak arrival time", peak_t, expected_t, 0.01)

    # --- Test 4: attenuated plane-wave amplitude reduction ---
    print("Test 4: Q=50 plane-wave amplitude reduction vs elastic")
    uz_q = plane_wave_seismogram(distance, c_ref, 1.0, tau, ss, src, dt)
    amp_ratio = float(np.max(np.abs(uz_q)) / np.max(np.abs(uz)))
    if amp_ratio >= 1.0:
        print(f"  FAIL: Q=50 amplitude ({amp_ratio:.4f}) not less than elastic")
        ok = False
        n_fail += 1
    elif amp_ratio < 0.3:
        print(f"  FAIL: Q=50 amplitude ({amp_ratio:.4f}) unreasonably small")
        ok = False
        n_fail += 1
    else:
        n_pass += 1

    # --- Test 5: low-Q (Q=20) coefficients ---
    print("Test 5: low-Q coefficient fitting (Q=20, gamma=0.6)")
    tau20, ss20, sp20 = build_fq8_coefficients(20.0, 20.0, 0.6, 1.0, 1.0)
    target_q20 = np.where(f_test < 1.0, 20.0, 20.0 * (f_test / 1.0) ** 0.6)
    q_s20 = realized_q(f_test, tau20, ss20)
    max_q_err20 = float(np.max(np.abs(q_s20 / target_q20 - 1.0)))
    check("max Q error (Q=20, gamma=0.6)", max_q_err20, 0.0, 0.10)
    if float(np.sum(ss20)) >= 1.0:
        print(f"  FAIL: sum(strength)={np.sum(ss20):.6f} >= 1 (unphysical)")
        ok = False
        n_fail += 1
    else:
        n_pass += 1

    # --- Test 6: high-Q (Q=210) coefficients ---
    print("Test 6: high-Q coefficient fitting (Q=210, gamma=0.6)")
    tau210, ss210, sp210 = build_fq8_coefficients(210.0, 210.0, 0.6, 1.0, 1.0)
    target_q210 = np.where(f_test < 1.0, 210.0, 210.0 * (f_test / 1.0) ** 0.6)
    q_s210 = realized_q(f_test, tau210, ss210)
    max_q_err210 = float(np.max(np.abs(q_s210 / target_q210 - 1.0)))
    check("max Q error (Q=210, gamma=0.6)", max_q_err210, 0.0, 0.10)

    # --- Test 7: EM/PM metric self-consistency ---
    print("Test 7: Kristekova EM/PM metric consistency")
    sig = ricker_wavelet(t, 1.0, 1.5)
    em_id, pm_id = envelope_phase_misfit(sig, sig, dt)
    check("EM(identical signals)", em_id, 0.0, 0.01)
    check("PM(identical signals)", pm_id, 0.0, 0.01)
    sig_shifted = ricker_wavelet(t, 1.0, 1.52)  # 20ms time shift
    em_sh, pm_sh = envelope_phase_misfit(sig, sig_shifted, dt)
    if pm_sh < 0.1:
        print(f"  FAIL: PM for shifted signal too small: {pm_sh:.4f}%")
        ok = False
        n_fail += 1
    else:
        n_pass += 1
    sig_scaled = 0.8 * sig
    em_sc, pm_sc = envelope_phase_misfit(sig, sig_scaled, dt)
    check_rel("EM(80% amplitude)", em_sc, 20.0, 0.5)

    print(f"\nSelftest: {n_pass} passed, {n_fail} failed")
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(
        description="f-k reference solutions for WaveQLab3D-fQ validation"
    )
    parser.add_argument("--case", type=str, default=None,
                        help="Benchmark case: elastic, constant_q, powerlaw, layered")
    parser.add_argument("--all", action="store_true",
                        help="Run all benchmark cases")
    parser.add_argument("--diagnostics", action="store_true",
                        help="Print Q(f) and dispersion diagnostics only")
    parser.add_argument("--Qs0", type=float, default=50.0)
    parser.add_argument("--Qp0", type=float, default=50.0)
    parser.add_argument("--gamma", type=float, default=0.0)
    parser.add_argument("--f-transition", type=float, default=1.0)
    parser.add_argument("--fref", type=float, default=1.0)
    parser.add_argument("--selftest", action="store_true",
                        help="Run self-test verifying coefficients and propagation")
    parser.add_argument("--compare", nargs=2, metavar=("REF", "TEST"),
                        help="Compare reference and test seismogram files (EM/PM)")
    parser.add_argument("--dt-compare", type=float, default=0.002,
                        help="Time step for --compare (default 0.002)")
    parser.add_argument("--fmax-compare", type=float, default=None,
                        help="Max frequency for --compare")
    parser.add_argument("--outdir", type=str, default="fk_reference_output")
    args = parser.parse_args()

    if args.selftest:
        success = run_selftest()
        sys.exit(0 if success else 1)

    if args.compare:
        em, pm = compare_seismograms(
            args.compare[0], args.compare[1], args.dt_compare,
            fmax=args.fmax_compare
        )
        print(f"EM = {em:.2f}%  PM = {pm:.2f}%")
        sys.exit(0)

    outdir = Path(args.outdir)

    if args.diagnostics:
        print_diagnostics(args.Qs0, args.Qp0, args.gamma,
                          args.f_transition, args.fref, outdir)
        return

    if args.all:
        for case in ["elastic", "constant_q", "powerlaw", "layered"]:
            run_benchmark(case, outdir)
            print()
        return

    if args.case:
        run_benchmark(args.case, outdir)
        return

    parser.print_help()


if __name__ == "__main__":
    main()
