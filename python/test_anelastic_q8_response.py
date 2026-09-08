#!/usr/bin/env python3
"""Static constitutive regression for the fixed anelastic-Q8 spectrum."""

from __future__ import annotations

import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "src" / "anelastic_q8_model.f90"
NMECH = 8
TARGET_Q = 50.0
FMIN = 0.08
FMAX = 15.0
NFREQ = 256


def read_weights() -> list[float]:
    text = MODEL.read_text(encoding="utf-8")
    match = re.search(r"weight\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if match is None:
        raise AssertionError("fixed Q8 weight array not found")
    values = re.findall(r"([-+0-9.]+(?:[eEdD][-+0-9]+)?)_wp", match.group(1))
    if len(values) != NMECH:
        raise AssertionError(f"expected {NMECH} Q8 weights, found {len(values)}")
    return [float(value.replace("d", "e").replace("D", "E")) for value in values]


def relaxation_times() -> list[float]:
    tau_min = 1.0 / (2.0 * math.pi * 20.0)
    tau_max = 1.0 / (2.0 * math.pi * 0.05)
    return [
        math.exp(
            math.log(tau_min)
            + (2.0 * k - 1.0)
            / (2.0 * NMECH)
            * math.log(tau_max / tau_min)
        )
        for k in range(1, NMECH + 1)
    ]


def realized_q(frequency: float, weights: list[float], tau: list[float]) -> float:
    omega = 2.0 * math.pi * frequency
    modulus_real = 1.0
    modulus_imag = 0.0
    for weight, relaxation_time in zip(weights, tau):
        x = omega * relaxation_time
        modulus_real -= (weight / TARGET_Q) / (1.0 + x * x)
        modulus_imag += (weight / TARGET_Q) * x / (1.0 + x * x)
    return abs(modulus_real / modulus_imag)


def main() -> None:
    model_text = MODEL.read_text(encoding="utf-8")
    assert "response anelastic-Q8 requires a valid &anelastic_Q8_list" in model_text
    assert "anelastic-Q8 requires finite, positive Qs0 and Qp0" in model_text
    assert "fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz" in model_text

    weights = read_weights()
    tau = relaxation_times()
    assert all(weight >= 0.0 for weight in weights)
    assert all(value > 0.0 for value in tau)

    max_error = 0.0
    for index in range(NFREQ):
        fraction = index / (NFREQ - 1)
        frequency = FMIN * (FMAX / FMIN) ** fraction
        q_value = realized_q(frequency, weights, tau)
        max_error = max(max_error, abs(q_value / TARGET_Q - 1.0))

    print(f"anelastic-Q8 maximum Q error: {100.0 * max_error:.4f}%")
    if max_error > 0.05:
        raise AssertionError("anelastic-Q8 maximum Q error exceeds 5%")


if __name__ == "__main__":
    main()
