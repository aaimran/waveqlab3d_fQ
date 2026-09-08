#!/usr/bin/env python3
"""Check a WaveQLab3D-Q input file the same way src/input_preflight.f90 does.

Usage:
    python check_input.py fname.in np

Reports namelist/preflight errors, the MPI Cartesian decomposition for ``np``
ranks, and the timestep WaveQLab would take from CFL and the block domains.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path


FLOAT_PATTERN = re.compile(
    r"[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eEdD][+-]?\d+)?"
)
IDENT_ASSIGN = re.compile(
    r"(?is)([A-Za-z][A-Za-z0-9_]*)"
    r"(?:\s*\(\s*(\d+)\s*\))?"
    r"(?:\s*%\s*([A-Za-z][A-Za-z0-9_]*)(?:\s*\(\s*(\d+)\s*\))?)?"
    r"\s*="
)

SUPPORTED_RESPONSES = {
    "elastic",
    "plastic",
    "anelastic",
    "low-pass",
    "anelastic-Q",
    "anelastic-Q4",
    "anelastic-Q8",
    "anelastic-cQ8-b2",
    "anelastic-cQ",
    "anelastic-Qf",
    "anelastic-fQ8",
    "constant-Q-4M",
    "constant-Q-8M",
    "frequency-Q-4M",
    "frequency-Q-8M",
}

PROBLEM_FIELDS = {
    "name",
    "problem",
    "response",
    "plastic_model",
    "nblocks",
    "nt",
    "cfl",
    "coupling",
    "fd_type",
    "order",
    "t_final",
    "mesh_source",
    "type_of_mesh",
    "material_source",
    "interpol",
    "w_stride",
    "w_fault",
    "use_topography",
    "topo",
    "mollify_source",
}

BLOCK_FIELDS = {
    "nqrs",
    "aqrs",
    "bqrs",
    "lc",
    "rc",
    "rho_s_p",
    "mu_beta_eta",
    "lqrs",
    "rqrs",
    "profile_type",
    "profile_path",
    "material_path",
    "topography_type",
    "topography_path",
    "pml_lqrs",
    "pml_rqrs",
    "npml",
    "faultsize",
}

OUTPUT_FIELDS = {
    "output_exact_moment",
    "output_seismograms",
    "output_station_info",
    "output_station_mapping",
    "output_fault_topo",
    "output_fields_block1",
    "output_fields_block2",
    "stride_fields",
    "station_xyz_index",
    "station_list",
    "station_list_file",
    "station_file_directory",
    "station_output_order",
    "station_number_in_list",
    "station_number_in_filename",
    "station_use_block_subdirectories",
    "common_stations_blocks",
    "station_add_header",
    "station_add_metadata",
}

Q4_FIELDS = {"weight_method", "qs0", "qp0", "fref", "fmin", "fmax"}
Q8_FIELDS = Q4_FIELDS
CQ8_B2_FIELDS = Q4_FIELDS
CQ_FIELDS = {
    "qs0",
    "qp0",
    "fref",
    "fmin",
    "fmax",
    "n_mechanisms",
    "nnls_samples",
    "coefficient_policy",
    "nnls_objective",
    "nnls_tolerance",
    "max_fit_error",
}
FQ8_FIELDS = {
    "coefficient_method",
    "weight_policy",
    "coarse_grain",
    "qs0",
    "qp0",
    "gamma",
    "f_transition",
    "fref",
}

CONSERVATIVE_MINIMUM_OWNED = 20
MACHINE_EPS = sys.float_info.epsilon
PI = 3.141592653589793

WITHERS_GAMMA = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
WITHERS_TAU_MIN = [
    0.0032, 0.0032, 0.0032, 0.0032, 0.0032,
    0.0032, 0.0032, 0.0066, 0.0066, 0.0085,
]
WITHERS_TAU_MAX = [
    15.9155, 15.9155, 15.9155, 15.9155, 15.9155,
    15.9155, 15.9155, 3.9789, 3.9789, 3.9789,
]


class InputError(ValueError):
    """Raised when the file cannot be parsed far enough to run preflight."""


@dataclass
class Diagnostic:
    severity: str
    code: str
    message: str
    section: str = ""
    field: str = ""
    suggestion: str = ""
    block_id: int = 0


@dataclass
class StencilRequirements:
    halo_width: int = 0
    boundary_width: int = 0
    minimum_global_points: int = 0
    minimum_owned_points: int = 0
    supported: bool = False
    operator_name: str = ""


@dataclass
class Block:
    block_id: int
    nqrs: tuple[int, int, int] = (0, 0, 0)
    aqrs: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bqrs: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rho_s_p: tuple[float, float, float] = (0.0, 0.0, 0.0)
    lqrs: tuple[int, int, int] = (0, 0, 0)
    rqrs: tuple[int, int, int] = (0, 0, 0)
    pml_lqrs: tuple[bool, bool, bool] = (False, False, False)
    pml_rqrs: tuple[bool, bool, bool] = (False, False, False)
    npml: int = 0

    @property
    def spacing(self) -> tuple[float, float, float]:
        return tuple(
            (upper - lower) / (count - 1)
            for lower, upper, count in zip(self.aqrs, self.bqrs, self.nqrs)
        )

    @property
    def h(self) -> float:
        return min(self.spacing)

    @property
    def wave_speed_factor(self) -> float:
        vs, vp = self.rho_s_p[1], self.rho_s_p[2]
        return math.sqrt(vs * vs + vp * vp)

    def elastic_dt(self, cfl: float) -> float:
        return cfl * self.h / self.wave_speed_factor


@dataclass
class QParams:
    kind: str = ""
    weight_method: str = "fixed-q50"
    qs0: list[float] = field(default_factory=lambda: [-1.0])
    qp0: list[float] = field(default_factory=lambda: [-1.0])
    fref: float = 1.0
    fmin: float = 0.05
    fmax: float = 20.0
    n_mechanisms: int = 8
    nnls_samples: int = 256
    coefficient_policy: str = "nnls-block-ps"
    nnls_objective: str = "relative-q"
    nnls_tolerance: float = 1.0e-10
    max_fit_error: float = 0.10
    coefficient_method: str = "conventional-nnls"
    weight_policy: str = "table-exact"
    coarse_grain: int = -1
    gamma: float = -1.0
    f_transition: float = -1.0


@dataclass
class Problem:
    name: str = "default"
    problem: str = "TPV5"
    response: str = "elastic"
    plastic_model: str = "default"
    nblocks: int = 2
    nt: int = 0
    cfl: float = 0.5
    coupling: str = "locked"
    fd_type: str = "traditional"
    order: int = 5
    t_final: float = 0.0
    mesh_source: str = "compute"
    type_of_mesh: str = "cartesian"
    material_source: str = "hardcode"
    interpol: bool = False
    w_stride: int = 1
    w_fault: bool = True
    use_topography: bool = False
    topo: float = 1.0
    mollify_source: bool = False


@dataclass
class OutputConfig:
    output_seismograms: bool = False
    station_list: str = "infile"
    station_list_file: str = ""
    station_number_in_list: bool = False


@dataclass
class Decomposition:
    process_dims: list[tuple[int, int, int]]
    rank_begin: list[int]
    rank_end: list[int]
    serial_shared_blocks: bool = False


def strip_fortran_comments(text: str) -> str:
    output: list[str] = []
    for line in text.splitlines():
        quote: str | None = None
        kept: list[str] = []
        index = 0
        while index < len(line):
            char = line[index]
            if quote:
                kept.append(char)
                if char == quote:
                    if index + 1 < len(line) and line[index + 1] == quote:
                        kept.append(line[index + 1])
                        index += 1
                    else:
                        quote = None
            elif char in ("'", '"'):
                quote = char
                kept.append(char)
            elif char == "!":
                break
            else:
                kept.append(char)
            index += 1
        output.append("".join(kept))
    return "\n".join(output)


def unquoted_slash(text: str, start: int) -> int:
    quote: str | None = None
    index = start
    while index < len(text):
        char = text[index]
        if quote:
            if char == quote:
                if index + 1 < len(text) and text[index + 1] == quote:
                    index += 1
                else:
                    quote = None
        elif char in ("'", '"'):
            quote = char
        elif char == "/":
            return index
        index += 1
    return -1


def namelist_body(text: str, name: str) -> str | None:
    match = re.search(rf"(?is)&\s*{re.escape(name)}\b", text)
    if not match:
        return None
    end = unquoted_slash(text, match.end())
    if end < 0:
        raise InputError(f"&{name} is missing a terminating /")
    return text[match.end() : end]


def quote_aware_assign_spans(body: str) -> list[tuple[int, re.Match[str]]]:
    spans: list[tuple[int, re.Match[str]]] = []
    quote: str | None = None
    index = 0
    while index < len(body):
        char = body[index]
        if quote:
            if char == quote:
                if index + 1 < len(body) and body[index + 1] == quote:
                    index += 1
                else:
                    quote = None
            index += 1
            continue
        if char in ("'", '"'):
            quote = char
            index += 1
            continue
        match = IDENT_ASSIGN.match(body, index)
        if match:
            spans.append((index, match))
            index = match.end()
            continue
        index += 1
    return spans


def parse_assignments(body: str) -> dict[str, str]:
    spans = quote_aware_assign_spans(body)
    values: dict[str, str] = {}
    for item, (start, match) in enumerate(spans):
        name = match.group(1).lower()
        index = match.group(2)
        component = match.group(3).lower() if match.group(3) else None
        if component:
            key = f"{name}({index})%{component}" if index else f"{name}%{component}"
        else:
            key = f"{name}({index})" if index else name
        value_start = match.end()
        value_end = spans[item + 1][0] if item + 1 < len(spans) else len(body)
        raw = body[value_start:value_end].strip()
        if raw.endswith(","):
            raw = raw[:-1].rstrip()
        values[key] = raw
    return values


def split_tokens(raw: str) -> list[str]:
    tokens: list[str] = []
    quote: str | None = None
    current: list[str] = []
    index = 0
    while index < len(raw):
        char = raw[index]
        if quote:
            current.append(char)
            if char == quote:
                if index + 1 < len(raw) and raw[index + 1] == quote:
                    current.append(raw[index + 1])
                    index += 1
                else:
                    quote = None
        elif char in ("'", '"'):
            quote = char
            current.append(char)
        elif char == ",":
            token = "".join(current).strip()
            if token:
                tokens.append(token)
            current = []
        else:
            current.append(char)
        index += 1
    token = "".join(current).strip()
    if token:
        tokens.append(token)
    return tokens


def fortran_float(value: str) -> float:
    token = value.strip().strip("()").replace("D", "e").replace("d", "e")
    match = FLOAT_PATTERN.fullmatch(token)
    if not match:
        raise InputError(f"invalid floating-point value: {value!r}")
    return float(token)


def fortran_int(value: str) -> int:
    number = fortran_float(value)
    integer = int(number)
    if float(integer) != number:
        raise InputError(f"expected an integer value, got {value!r}")
    return integer


def fortran_bool(value: str) -> bool:
    token = value.strip().strip(".").lower()
    if token in {"t", "true"}:
        return True
    if token in {"f", "false"}:
        return False
    raise InputError(f"invalid logical value: {value!r}")


def fortran_string(value: str) -> str:
    token = value.strip()
    if len(token) >= 2 and token[0] in "'\"" and token[-1] == token[0]:
        quote = token[0]
        return token[1:-1].replace(quote * 2, quote)
    return token


def optional_assign(assigns: dict[str, str], name: str) -> str | None:
    return assigns.get(name.lower())


def require_floats(raw: str, count: int, label: str) -> tuple[float, ...]:
    tokens = FLOAT_PATTERN.findall(raw)
    if len(tokens) != count:
        raise InputError(f"{label} has {len(tokens)} values; expected {count}")
    return tuple(fortran_float(token) for token in tokens)


def require_ints(raw: str, count: int, label: str) -> tuple[int, ...]:
    values = require_floats(raw, count, label)
    integers = tuple(int(value) for value in values)
    if any(float(integer) != value for integer, value in zip(integers, values)):
        raise InputError(f"{label} must contain integers")
    return integers


def require_bools(raw: str, count: int, label: str) -> tuple[bool, ...]:
    tokens = split_tokens(raw)
    if len(tokens) != count:
        raise InputError(f"{label} has {len(tokens)} values; expected {count}")
    return tuple(fortran_bool(token) for token in tokens)


def unknown_fields(assigns: dict[str, str], allowed: set[str], section: str) -> list[str]:
    unknown: list[str] = []
    for key in assigns:
        base = key.split("%", 1)[-1]
        base = re.sub(r"\(\d+\)$", "", base)
        if base.lower() not in {item.lower() for item in allowed}:
            unknown.append(key)
    return unknown


def get_stencil_requirements(fd_type: str, order: int) -> StencilRequirements:
    req = StencilRequirements()
    kind = fd_type.strip()
    if kind == "traditional":
        req.halo_width = 3
        req.boundary_width = 6
        req.operator_name = "traditional-6"
        req.supported = True
    elif kind == "upwind":
        if order in (2, 3):
            req.halo_width, req.boundary_width = 3, 2
        elif order in (4, 5):
            req.halo_width, req.boundary_width = 3, 4
        elif order in (6, 7):
            req.halo_width, req.boundary_width = 4, 6
        elif order in (8, 9):
            req.halo_width, req.boundary_width = 5, 8
        else:
            return req
        req.operator_name = f"upwind-{order}"
        req.supported = True
    elif kind == "upwind_drp":
        if order == 3:
            req.halo_width, req.boundary_width = 3, 4
        elif order in (4, 5, 66):
            req.halo_width, req.boundary_width = 4, 6
        elif order in (6, 7, 679):
            req.halo_width, req.boundary_width = 5, 8
        else:
            return req
        req.operator_name = f"upwind-drp-{order}"
        req.supported = True
    else:
        return req
    req.minimum_global_points = 2 * req.boundary_width
    req.minimum_owned_points = max(
        CONSERVATIVE_MINIMUM_OWNED,
        2 * req.halo_width + 1,
        2 * req.boundary_width,
    )
    return req


def topology_fits(
    points: tuple[int, int, int], dims: tuple[int, int, int], required_owned: int
) -> bool:
    for count, dim in zip(points, dims):
        if dim < 1 or dim > count:
            return False
        if dim > 1 and count // dim < required_owned:
            return False
    return True


def select_single_block_topology(
    points: tuple[int, int, int], nranks: int, required_owned: int
) -> tuple[int, int, int] | None:
    best_score = math.inf
    best: tuple[int, int, int] | None = None
    for p in range(1, nranks + 1):
        if nranks % p != 0:
            continue
        for q in range(1, nranks // p + 1):
            if (nranks // p) % q != 0:
                continue
            r = nranks // (p * q)
            dims = (p, q, r)
            if not topology_fits(points, dims, required_owned):
                continue
            local = [count / dim for count, dim in zip(points, dims)]
            score = (
                1.0 / local[0]
                + 1.0 / local[1]
                + 1.0 / local[2]
                + 100.0 * MACHINE_EPS * ((q - 1) + 2 * (r - 1))
            )
            if score < best_score:
                best_score = score
                best = dims
    return best


def resolve_decomposition(
    blocks: list[Block], world_size: int, minimum_owned: int
) -> Decomposition | None:
    nblocks = len(blocks)
    if nblocks == 1:
        dims = select_single_block_topology(blocks[0].nqrs, world_size, minimum_owned)
        if dims is None:
            return None
        return Decomposition([dims], [0], [world_size - 1], False)

    if world_size == 1:
        return Decomposition([(1, 1, 1), (1, 1, 1)], [0, 0], [0, 0], True)

    work = [float(math.prod(block.nqrs)) for block in blocks]
    best_score = math.inf
    best_dims = [(0, 0, 0), (0, 0, 0)]
    best_sizes = [0, 0]
    found = False
    for pr in range(1, world_size + 1):
        for ps in range(1, world_size // pr + 1):
            tangential = pr * ps
            if world_size % tangential != 0:
                continue
            qsum = world_size // tangential
            if qsum < 2:
                continue
            for pq1 in range(1, qsum):
                pq2 = qsum - pq1
                dims1 = (pq1, pr, ps)
                dims2 = (pq2, pr, ps)
                if not topology_fits(blocks[0].nqrs, dims1, minimum_owned):
                    continue
                if not topology_fits(blocks[1].nqrs, dims2, minimum_owned):
                    continue
                load1 = work[0] / float(pq1 * tangential)
                load2 = work[1] / float(pq2 * tangential)
                score = abs(load1 - load2) / max(load1, load2)
                if score < best_score:
                    best_score = score
                    best_dims = [dims1, dims2]
                    best_sizes = [pq1 * tangential, pq2 * tangential]
                    found = True
    if not found:
        return None
    return Decomposition(
        best_dims,
        [0, best_sizes[0]],
        [best_sizes[0] - 1, world_size - 1],
        False,
    )


def interior_tau(fmin: float, fmax: float, n_mechanisms: int) -> list[float]:
    taumin = 1.0 / (2.0 * PI * fmax)
    taumax = 1.0 / (2.0 * PI * fmin)
    ratio = math.log(taumax / taumin)
    return [
        math.exp(math.log(taumin) + (2.0 * k - 1.0) / (2.0 * n_mechanisms) * ratio)
        for k in range(1, n_mechanisms + 1)
    ]


def withers_tau(gamma: float, f_transition: float) -> list[float]:
    clamped = max(0.0, min(0.9, gamma))
    if clamped <= WITHERS_GAMMA[0]:
        idx_low = idx_high = 0
        alpha = 0.0
    elif clamped >= WITHERS_GAMMA[-1]:
        idx_low = idx_high = len(WITHERS_GAMMA) - 1
        alpha = 0.0
    else:
        idx_low = idx_high = 0
        alpha = 0.0
        for i in range(len(WITHERS_GAMMA) - 1):
            if WITHERS_GAMMA[i] <= clamped <= WITHERS_GAMMA[i + 1]:
                idx_low = i
                idx_high = i + 1
                alpha = (clamped - WITHERS_GAMMA[i]) / (
                    WITHERS_GAMMA[i + 1] - WITHERS_GAMMA[i]
                )
                break
    taumin = (1.0 - alpha) * WITHERS_TAU_MIN[idx_low] + alpha * WITHERS_TAU_MIN[idx_high]
    taumax = (1.0 - alpha) * WITHERS_TAU_MAX[idx_low] + alpha * WITHERS_TAU_MAX[idx_high]
    ratio = math.log(taumax) - math.log(taumin)
    tau = [
        math.exp(math.log(taumin) + (2 * k - 1) / 16.0 * ratio) for k in range(1, 9)
    ]
    return [value / f_transition for value in tau]


def q4_relaxation_dt_limit(params: QParams) -> float:
    return 2.0 * min(interior_tau(params.fmin, params.fmax, 4))


def q8_relaxation_dt_limit(params: QParams) -> float:
    return 2.0 * min(interior_tau(params.fmin, params.fmax, 8))


def cq_relaxation_dt_limit(params: QParams) -> float:
    taumin = 1.0 / (2.0 * PI * params.fmax)
    taumax = 1.0 / (2.0 * PI * params.fmin)
    if params.coefficient_policy == "fixed-q50":
        return 2.0 * math.exp(
            math.log(taumin)
            + 1.0 / (2.0 * params.n_mechanisms) * math.log(taumax / taumin)
        )
    return 2.0 * taumin


def fq8_relaxation_dt_limit(params: QParams) -> float:
    return 2.0 * min(withers_tau(params.gamma, params.f_transition))


def finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0.0


class Preflight:
    def __init__(self, path: Path, np: int) -> None:
        self.path = path
        self.np = np
        self.issues: list[Diagnostic] = []
        self.raw = path.read_text(encoding="utf-8")
        self.text = strip_fortran_comments(self.raw)
        self.problem = Problem()
        self.blocks: list[Block] = []
        self.output = OutputConfig()
        self.q: QParams | None = None
        self.stencil = StencilRequirements()
        self.decomp: Decomposition | None = None

    def add(
        self,
        severity: str,
        code: str,
        message: str,
        section: str = "",
        field: str = "",
        suggestion: str = "",
        block_id: int = 0,
    ) -> None:
        self.issues.append(
            Diagnostic(severity, code, message, section, field, suggestion, block_id)
        )

    def has_errors(self) -> bool:
        return any(item.severity == "ERROR" for item in self.issues)

    def error_count(self) -> int:
        return sum(item.severity == "ERROR" for item in self.issues)

    def warning_count(self) -> int:
        return sum(item.severity == "WARNING" for item in self.issues)

    def run(self) -> None:
        self.read_problem()
        if not self.has_errors():
            self.read_blocks()
        blocks_ok = bool(self.blocks) and not any(
            item.severity == "ERROR"
            and item.code.startswith(("CFG-BLOCK", "CFG-STENCIL", "CFG-PML", "CFG-BOUNDARY", "CFG-IFACE"))
            for item in self.issues
        )
        if blocks_ok:
            self.read_anelastic()
        if not self.has_errors():
            self.read_output()
        if blocks_ok and self.stencil.supported:
            self.resolve_decomp()

    def read_problem(self) -> None:
        try:
            body = namelist_body(self.text, "problem_list")
        except InputError as error:
            self.add("ERROR", "CFG-PROBLEM-001", str(error), section="problem_list")
            return
        if body is None:
            self.add(
                "ERROR",
                "CFG-PROBLEM-001",
                "Cannot parse &problem_list: namelist not found",
                section="problem_list",
                suggestion="Fix namelist syntax or unknown fields.",
            )
            return
        try:
            assigns = parse_assignments(body)
        except InputError as error:
            self.add("ERROR", "CFG-PROBLEM-001", f"Cannot parse &problem_list: {error}",
                     section="problem_list")
            return
        for key in unknown_fields(assigns, PROBLEM_FIELDS, "problem_list"):
            self.add(
                "ERROR",
                "CFG-PROBLEM-001",
                f"Cannot parse &problem_list: unknown field {key}",
                section="problem_list",
                field=key,
                suggestion="Fix namelist syntax or unknown fields.",
            )
            return
        try:
            self._fill_problem(assigns)
        except InputError as error:
            self.add("ERROR", "CFG-PROBLEM-001", f"Cannot parse &problem_list: {error}",
                     section="problem_list")
            return
        self.validate_problem()

    def _fill_problem(self, assigns: dict[str, str]) -> None:
        p = self.problem
        if raw := optional_assign(assigns, "name"):
            p.name = fortran_string(raw)
        if raw := optional_assign(assigns, "problem"):
            p.problem = fortran_string(raw)
        if raw := optional_assign(assigns, "response"):
            p.response = fortran_string(raw).strip()
        if raw := optional_assign(assigns, "plastic_model"):
            p.plastic_model = fortran_string(raw)
        if raw := optional_assign(assigns, "nblocks"):
            p.nblocks = fortran_int(raw)
        if raw := optional_assign(assigns, "nt"):
            p.nt = fortran_int(raw)
        if raw := optional_assign(assigns, "cfl"):
            p.cfl = fortran_float(raw)
        if raw := optional_assign(assigns, "coupling"):
            p.coupling = fortran_string(raw)
        if raw := optional_assign(assigns, "fd_type"):
            p.fd_type = fortran_string(raw).strip()
        if raw := optional_assign(assigns, "order"):
            p.order = fortran_int(raw)
        if raw := optional_assign(assigns, "t_final"):
            p.t_final = fortran_float(raw)
        if raw := optional_assign(assigns, "mesh_source"):
            p.mesh_source = fortran_string(raw)
        if raw := optional_assign(assigns, "type_of_mesh"):
            p.type_of_mesh = fortran_string(raw)
        if raw := optional_assign(assigns, "material_source"):
            p.material_source = fortran_string(raw)
        if raw := optional_assign(assigns, "interpol"):
            p.interpol = fortran_bool(raw)
        if raw := optional_assign(assigns, "w_stride"):
            p.w_stride = fortran_int(raw)
        if raw := optional_assign(assigns, "w_fault"):
            p.w_fault = fortran_bool(raw)
        if raw := optional_assign(assigns, "use_topography"):
            p.use_topography = fortran_bool(raw)
        if raw := optional_assign(assigns, "topo"):
            p.topo = fortran_float(raw)
        if raw := optional_assign(assigns, "mollify_source"):
            p.mollify_source = fortran_bool(raw)

    def validate_problem(self) -> None:
        p = self.problem
        if p.response not in SUPPORTED_RESPONSES:
            self.add(
                "ERROR",
                "CFG-PROBLEM-002",
                f"Unsupported response: {p.response}",
                section="problem_list",
                field="response",
            )
        if p.nblocks < 1 or p.nblocks > 2:
            self.add(
                "ERROR",
                "CFG-PROBLEM-003",
                "nblocks must be 1 or 2.",
                section="problem_list",
                field="nblocks",
            )
        if p.cfl <= 0.0:
            self.add(
                "ERROR",
                "CFG-PROBLEM-004",
                "CFL must be positive.",
                section="problem_list",
                field="CFL",
            )
        if p.t_final < 0.0:
            self.add(
                "ERROR",
                "CFG-PROBLEM-005",
                "t_final cannot be negative.",
                section="problem_list",
                field="t_final",
            )
        if p.w_stride < 1:
            self.add(
                "ERROR",
                "CFG-PROBLEM-006",
                "w_stride must be at least 1.",
                section="problem_list",
                field="w_stride",
            )
        if p.t_final == 0.0:
            self.add(
                "WARNING",
                "CFG-PROBLEM-007",
                "t_final is zero; initialization will run but no time steps are expected.",
                section="problem_list",
                field="t_final",
            )
        self.stencil = get_stencil_requirements(p.fd_type, p.order)
        if not self.stencil.supported:
            self.add(
                "ERROR",
                "CFG-FD-001",
                f"Unsupported finite-difference type/order combination: {p.fd_type}",
                section="problem_list",
                field="fd_type/order",
                suggestion="Use traditional, upwind order 2-9, or a supported upwind_drp order.",
            )
        if p.response == "anelastic-Q":
            self.add(
                "WARNING",
                "CFG-Q4-DEP-001",
                "Response anelastic-Q is deprecated and is normalized to anelastic-Q4.",
                section="problem_list",
                field="response",
                suggestion="Set response='anelastic-Q4' and use &anelastic_Q4_list.",
            )
            p.response = "anelastic-Q4"
        if p.response == "frequency-Q-8M":
            self.add(
                "WARNING",
                "CFG-FQ8-DEP-001",
                "Response frequency-Q-8M is deprecated and is normalized to anelastic-fQ8.",
                section="problem_list",
                field="response",
                suggestion="Set response='anelastic-fQ8' and use &anelastic_fQ8_list.",
            )
            p.response = "anelastic-fQ8"

    def read_blocks(self) -> None:
        try:
            body = namelist_body(self.text, "block_list")
        except InputError as error:
            self.add("ERROR", "CFG-BLOCK-001", str(error), section="block_list")
            return
        if body is None:
            self.add(
                "ERROR",
                "CFG-BLOCK-001",
                "Cannot parse &block_list: namelist not found",
                section="block_list",
                suggestion="Define every active block and fix namelist syntax.",
            )
            return
        try:
            assigns = parse_assignments(body)
        except InputError as error:
            self.add("ERROR", "CFG-BLOCK-001", f"Cannot parse &block_list: {error}",
                     section="block_list")
            return
        for key in assigns:
            match = re.fullmatch(r"btp\((\d+)\)%([a-z0-9_]+)", key)
            if not match:
                self.add(
                    "ERROR",
                    "CFG-BLOCK-001",
                    f"Cannot parse &block_list: unexpected field {key}",
                    section="block_list",
                    field=key,
                    suggestion="Define every active block and fix namelist syntax.",
                )
                return
            if int(match.group(1)) not in (1, 2) or match.group(2) not in BLOCK_FIELDS:
                self.add(
                    "ERROR",
                    "CFG-BLOCK-001",
                    f"Cannot parse &block_list: unknown field {key}",
                    section="block_list",
                    field=key,
                    suggestion="Define every active block and fix namelist syntax.",
                )
                return
        try:
            self.blocks = [
                self._fill_block(assigns, block_id)
                for block_id in range(1, self.problem.nblocks + 1)
            ]
        except InputError as error:
            self.add("ERROR", "CFG-BLOCK-001", f"Cannot parse &block_list: {error}",
                     section="block_list")
            return
        self.validate_blocks()
        self.validate_exterior_boundaries()
        if self.problem.nblocks == 2:
            self.validate_two_block_interface()

    def resolve_decomp(self) -> None:
        self.decomp = resolve_decomposition(
            self.blocks, self.np, self.stencil.minimum_owned_points
        )
        if self.decomp is None:
            if self.problem.nblocks == 1:
                self.add(
                    "ERROR",
                    "CFG-DECOMP-002",
                    "No Cartesian factorization satisfies the owned-subdomain bound.",
                    section="decomposition",
                    block_id=1,
                    suggestion="Use fewer ranks or enlarge the grid.",
                )
            else:
                self.add(
                    "ERROR",
                    "CFG-DECOMP-101",
                    "No asymmetric two-block topology satisfies the shared Y-Z topology and 20-point bound.",
                    section="decomposition",
                    suggestion="Use fewer ranks, enlarge a block, or choose a factorable world size.",
                )

    def _fill_block(self, assigns: dict[str, str], block_id: int) -> Block:
        block = Block(block_id)
        prefix = f"btp({block_id})%"
        if raw := assigns.get(prefix + "nqrs"):
            block.nqrs = require_ints(raw, 3, prefix + "nqrs")
        if raw := assigns.get(prefix + "aqrs"):
            block.aqrs = require_floats(raw, 3, prefix + "aqrs")
        if raw := assigns.get(prefix + "bqrs"):
            block.bqrs = require_floats(raw, 3, prefix + "bqrs")
        if raw := assigns.get(prefix + "rho_s_p"):
            block.rho_s_p = require_floats(raw, 3, prefix + "rho_s_p")
        if raw := assigns.get(prefix + "lqrs"):
            block.lqrs = require_ints(raw, 3, prefix + "lqrs")
        if raw := assigns.get(prefix + "rqrs"):
            block.rqrs = require_ints(raw, 3, prefix + "rqrs")
        if raw := assigns.get(prefix + "pml_lqrs"):
            block.pml_lqrs = require_bools(raw, 3, prefix + "pml_lqrs")
        if raw := assigns.get(prefix + "pml_rqrs"):
            block.pml_rqrs = require_bools(raw, 3, prefix + "pml_rqrs")
        if raw := assigns.get(prefix + "npml"):
            block.npml = fortran_int(raw)
        return block

    def validate_blocks(self) -> None:
        req = self.stencil
        for block in self.blocks:
            if any(count < 2 for count in block.nqrs):
                self.add(
                    "ERROR",
                    "CFG-BLOCK-002",
                    "Every active grid dimension must contain at least 2 points.",
                    section="block_list",
                    field="nqrs",
                    block_id=block.block_id,
                )
            if any(upper <= lower for lower, upper in zip(block.aqrs, block.bqrs)):
                self.add(
                    "ERROR",
                    "CFG-BLOCK-003",
                    "Every bqrs coordinate must exceed the corresponding aqrs coordinate.",
                    section="block_list",
                    field="aqrs/bqrs",
                    block_id=block.block_id,
                )
            if any(value <= 0.0 for value in block.rho_s_p):
                self.add(
                    "ERROR",
                    "CFG-BLOCK-004",
                    "Density, Vs, and Vp must be positive.",
                    section="block_list",
                    field="rho_s_p",
                    block_id=block.block_id,
                )
            if block.npml < 0:
                self.add(
                    "ERROR",
                    "CFG-BLOCK-005",
                    "npml cannot be negative.",
                    section="block_list",
                    field="npml",
                    block_id=block.block_id,
                )
            if req.supported and any(
                count < req.minimum_global_points for count in block.nqrs
            ):
                self.add(
                    "ERROR",
                    "CFG-STENCIL-001",
                    "Grid is too small for nonoverlapping physical-boundary closures.",
                    section="block_list",
                    field="nqrs",
                    block_id=block.block_id,
                    suggestion="Increase every grid dimension to the operator global minimum.",
                )
            if block.npml > 0:
                for count, left, right in zip(block.nqrs, block.pml_lqrs, block.pml_rqrs):
                    occupied = block.npml * (int(left) + int(right))
                    if count <= occupied:
                        self.add(
                            "ERROR",
                            "CFG-PML-001",
                            "PML layers leave no non-PML core in at least one direction.",
                            section="block_list",
                            field="npml/pml_lqrs/pml_rqrs",
                            block_id=block.block_id,
                        )
                        break

    def validate_exterior_boundaries(self) -> None:
        if self.problem.nblocks == 1:
            block = self.blocks[0]
            if any(code == 0 for code in block.lqrs) or any(code == 0 for code in block.rqrs):
                self.add(
                    "ERROR",
                    "CFG-BOUNDARY-001",
                    "Boundary code 0 is reserved for an internal block interface; all single-block faces are exterior.",
                    section="block_list",
                    field="btp(1)%lqrs/rqrs",
                    block_id=1,
                    suggestion="Use boundary code 1 (characteristic) or 2 (free surface) on every exterior face.",
                )
        elif self.problem.nblocks == 2:
            b1, b2 = self.blocks
            if any(code == 0 for code in b1.lqrs) or any(code == 0 for code in b1.rqrs[1:]):
                self.add(
                    "ERROR",
                    "CFG-BOUNDARY-001",
                    "Block 1 uses interface boundary code 0 on an exterior face.",
                    section="block_list",
                    field="btp(1)%lqrs or btp(1)%rqrs(2:3)",
                    block_id=1,
                    suggestion="Only btp(1)%rqrs(1) may be 0 in the supported two-block topology.",
                )
            if any(code == 0 for code in b2.rqrs) or any(code == 0 for code in b2.lqrs[1:]):
                self.add(
                    "ERROR",
                    "CFG-BOUNDARY-001",
                    "Block 2 uses interface boundary code 0 on an exterior face.",
                    section="block_list",
                    field="btp(2)%rqrs or btp(2)%lqrs(2:3)",
                    block_id=2,
                    suggestion="Only btp(2)%lqrs(1) may be 0 in the supported two-block topology.",
                )

    def validate_two_block_interface(self) -> None:
        b1, b2 = self.blocks
        coords = list(b1.aqrs) + list(b1.bqrs) + list(b2.aqrs) + list(b2.bqrs)
        scale = max(1.0, max(abs(value) for value in coords))
        tolerance = 1000.0 * MACHINE_EPS * scale
        if b1.nqrs[1:] != b2.nqrs[1:]:
            self.add(
                "ERROR",
                "CFG-IFACE-001",
                "The two blocks must have identical r/s (Y-Z) grid counts.",
                section="block_list",
                field="btp%nqrs(2:3)",
                suggestion="Keep Y and Z counts equal; only the q/X count may differ.",
            )
        if any(abs(a - b) > tolerance for a, b in zip(b1.aqrs[1:], b2.aqrs[1:])) or any(
            abs(a - b) > tolerance for a, b in zip(b1.bqrs[1:], b2.bqrs[1:])
        ):
            self.add(
                "ERROR",
                "CFG-IFACE-002",
                "The two blocks must have identical r/s (Y-Z) physical extents.",
                section="block_list",
                field="btp%aqrs/bqrs(2:3)",
            )
        if abs(b1.bqrs[0] - b2.aqrs[0]) > tolerance:
            self.add(
                "ERROR",
                "CFG-IFACE-003",
                "Block 1 right-q coordinate must equal block 2 left-q coordinate.",
                section="block_list",
                field="q interface coordinate",
            )
        if b1.rqrs[0] != 0 or b2.lqrs[0] != 0:
            self.add(
                "ERROR",
                "CFG-IFACE-004",
                "Internal q faces must both use interface boundary code 0.",
                section="block_list",
                field="btp(1)%rqrs(1), btp(2)%lqrs(1)",
            )
        if b1.pml_rqrs[0] or b2.pml_lqrs[0]:
            self.add(
                "ERROR",
                "CFG-IFACE-005",
                "PML cannot be enabled on either side of the internal interface.",
                section="block_list",
                field="internal q-face PML",
            )

    def read_anelastic(self) -> None:
        response = self.problem.response
        if response == "anelastic-Q4":
            self.q = self.read_q4()
        elif response == "anelastic-Q8":
            self.q = self.read_q8()
        elif response == "anelastic-cQ8-b2":
            if self.problem.nblocks != 2:
                self.add(
                    "ERROR",
                    "CFG-CQ8-B2-002",
                    "Response anelastic-cQ8-b2 requires exactly two blocks.",
                    section="problem_list",
                    field="nblocks",
                    suggestion="Set nblocks=2 or use response='anelastic-Q8'.",
                )
            else:
                self.q = self.read_cq8_b2()
        elif response == "anelastic-cQ":
            if self.problem.nblocks != 2:
                self.add(
                    "ERROR",
                    "CFG-CQ-002",
                    "Response anelastic-cQ requires exactly two blocks.",
                    section="problem_list",
                    field="nblocks",
                    suggestion="Set nblocks=2.",
                )
            else:
                self.q = self.read_cq()
        elif response == "anelastic-fQ8":
            self.q = self.read_fq8()

    def _q_assigns(self, namelist: str, code: str, allowed: set[str], message: str) -> dict[str, str] | None:
        try:
            body = namelist_body(self.text, namelist)
        except InputError as error:
            self.add("ERROR", code, str(error), section=namelist)
            return None
        if body is None:
            self.add("ERROR", code, message, section=namelist)
            return None
        try:
            assigns = parse_assignments(body)
        except InputError as error:
            self.add("ERROR", code, str(error), section=namelist)
            return None
        unknown = unknown_fields(assigns, allowed, namelist)
        if unknown:
            self.add(
                "ERROR",
                code,
                f"unknown field {unknown[0]} in &{namelist}",
                section=namelist,
                field=unknown[0],
            )
            return None
        return assigns

    def read_q4(self) -> QParams | None:
        assigns = self._q_assigns(
            "anelastic_Q4_list",
            "CFG-Q4-001",
            Q4_FIELDS,
            "response anelastic-Q4 requires a valid &anelastic_Q4_list namelist",
        )
        if assigns is None:
            return None
        params = QParams(kind="Q4")
        try:
            self._fill_shared_q(params, assigns)
        except InputError as error:
            self.add("ERROR", "CFG-Q4-001", str(error), section="anelastic_Q4_list")
            return None
        if not all(finite_positive(value) for value in params.qs0 + params.qp0):
            self.add(
                "ERROR",
                "CFG-Q4-001",
                "anelastic-Q4 requires finite, positive Qs0 and Qp0",
                section="anelastic_Q4_list",
                suggestion="Provide positive Qs0/Qp0 in &anelastic_Q4_list; legacy c is unsupported.",
            )
            return None
        if not finite_positive(params.fref):
            self.add("ERROR", "CFG-Q4-001", "anelastic-Q4 fref must be finite and positive",
                     section="anelastic_Q4_list")
            return None
        if not (math.isfinite(params.fmin) and math.isfinite(params.fmax) and 0.0 < params.fmin < params.fmax):
            self.add(
                "ERROR",
                "CFG-Q4-001",
                "anelastic-Q4 requires finite frequencies with 0 < fmin < fmax",
                section="anelastic_Q4_list",
            )
            return None
        if params.weight_method == "fixed-q50":
            if abs(params.fmin - 0.05) > 100.0 * MACHINE_EPS or abs(params.fmax - 20.0) > 100.0 * MACHINE_EPS:
                self.add(
                    "ERROR",
                    "CFG-Q4-001",
                    "fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz",
                    section="anelastic_Q4_list",
                )
                return None
        else:
            self.add(
                "ERROR",
                "CFG-Q4-001",
                f"unsupported anelastic-Q4 weight_method: {params.weight_method}",
                section="anelastic_Q4_list",
            )
            return None
        return params

    def read_q8(self) -> QParams | None:
        assigns = self._q_assigns(
            "anelastic_Q8_list",
            "CFG-Q8-001",
            Q8_FIELDS,
            "response anelastic-Q8 requires a valid &anelastic_Q8_list namelist",
        )
        if assigns is None:
            return None
        params = QParams(kind="Q8")
        try:
            self._fill_shared_q(params, assigns)
        except InputError as error:
            self.add("ERROR", "CFG-Q8-001", str(error), section="anelastic_Q8_list")
            return None
        if not all(finite_positive(value) for value in params.qs0 + params.qp0):
            self.add(
                "ERROR",
                "CFG-Q8-001",
                "anelastic-Q8 requires finite, positive Qs0 and Qp0",
                section="anelastic_Q8_list",
                suggestion="Provide positive Qs0/Qp0 and a supported coefficient setup.",
            )
            return None
        if not finite_positive(params.fref):
            self.add("ERROR", "CFG-Q8-001", "anelastic-Q8 fref must be finite and positive",
                     section="anelastic_Q8_list")
            return None
        if not (math.isfinite(params.fmin) and math.isfinite(params.fmax) and 0.0 < params.fmin < params.fmax):
            self.add(
                "ERROR",
                "CFG-Q8-001",
                "anelastic-Q8 requires finite frequencies with 0 < fmin < fmax",
                section="anelastic_Q8_list",
            )
            return None
        if params.weight_method == "fixed-q50":
            if abs(params.fmin - 0.05) > 100.0 * MACHINE_EPS or abs(params.fmax - 20.0) > 100.0 * MACHINE_EPS:
                self.add(
                    "ERROR",
                    "CFG-Q8-001",
                    "fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz",
                    section="anelastic_Q8_list",
                )
                return None
        else:
            self.add(
                "ERROR",
                "CFG-Q8-001",
                f"unsupported anelastic-Q8 weight_method in this increment: {params.weight_method}",
                section="anelastic_Q8_list",
            )
            return None
        return params

    def read_cq8_b2(self) -> QParams | None:
        assigns = self._q_assigns(
            "anelastic_cQ8_b2_list",
            "CFG-CQ8-B2-001",
            CQ8_B2_FIELDS,
            "response anelastic-cQ8-b2 requires a valid &anelastic_cQ8_b2_list namelist",
        )
        if assigns is None:
            return None
        params = QParams(kind="cQ8-b2", qs0=[-1.0, -1.0], qp0=[-1.0, -1.0])
        try:
            self._fill_shared_q(params, assigns, pair=True)
        except InputError as error:
            self.add("ERROR", "CFG-CQ8-B2-001", str(error), section="anelastic_cQ8_b2_list")
            return None
        if not all(finite_positive(value) for value in params.qs0 + params.qp0):
            self.add(
                "ERROR",
                "CFG-CQ8-B2-001",
                "anelastic-cQ8-b2 requires two finite, positive Qs0 and Qp0 values",
                section="anelastic_cQ8_b2_list",
                suggestion="Provide two positive Qs0/Qp0 values and a supported coefficient setup.",
            )
            return None
        if not finite_positive(params.fref):
            self.add(
                "ERROR",
                "CFG-CQ8-B2-001",
                "anelastic-cQ8-b2 fref must be finite and positive",
                section="anelastic_cQ8_b2_list",
            )
            return None
        if not (math.isfinite(params.fmin) and math.isfinite(params.fmax) and 0.0 < params.fmin < params.fmax):
            self.add(
                "ERROR",
                "CFG-CQ8-B2-001",
                "anelastic-cQ8-b2 requires finite frequencies with 0 < fmin < fmax",
                section="anelastic_cQ8_b2_list",
            )
            return None
        if params.weight_method != "fixed-q50":
            self.add(
                "ERROR",
                "CFG-CQ8-B2-001",
                f"unsupported anelastic-cQ8-b2 weight_method: {params.weight_method}",
                section="anelastic_cQ8_b2_list",
            )
            return None
        if abs(params.fmin - 0.05) > 100.0 * MACHINE_EPS or abs(params.fmax - 20.0) > 100.0 * MACHINE_EPS:
            self.add(
                "ERROR",
                "CFG-CQ8-B2-001",
                "anelastic-cQ8-b2 fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz",
                section="anelastic_cQ8_b2_list",
            )
            return None
        return params

    def read_cq(self) -> QParams | None:
        assigns = self._q_assigns(
            "anelastic_cQ_list",
            "CFG-CQ-001",
            CQ_FIELDS,
            "response anelastic-cQ requires a valid &anelastic_cQ_list namelist",
        )
        if assigns is None:
            return None
        params = QParams(kind="cQ", qs0=[-1.0, -1.0], qp0=[-1.0, -1.0])
        try:
            self._fill_shared_q(params, assigns, pair=True)
            if raw := optional_assign(assigns, "n_mechanisms"):
                params.n_mechanisms = fortran_int(raw)
            if raw := optional_assign(assigns, "nnls_samples"):
                params.nnls_samples = fortran_int(raw)
            if raw := optional_assign(assigns, "coefficient_policy"):
                params.coefficient_policy = fortran_string(raw).strip()
            if raw := optional_assign(assigns, "nnls_objective"):
                params.nnls_objective = fortran_string(raw).strip()
            if raw := optional_assign(assigns, "nnls_tolerance"):
                params.nnls_tolerance = fortran_float(raw)
            if raw := optional_assign(assigns, "max_fit_error"):
                params.max_fit_error = fortran_float(raw)
        except InputError as error:
            self.add("ERROR", "CFG-CQ-001", str(error), section="anelastic_cQ_list")
            return None
        if not all(finite_positive(value) for value in params.qs0 + params.qp0):
            self.add(
                "ERROR",
                "CFG-CQ-001",
                "anelastic-cQ requires two finite, positive Qs0 and Qp0 values",
                section="anelastic_cQ_list",
                suggestion="Provide valid Q pairs, frequencies, policy, and n_mechanisms=4..8.",
            )
            return None
        if not (
            math.isfinite(params.fref)
            and math.isfinite(params.fmin)
            and math.isfinite(params.fmax)
            and params.fmin > 0.0
            and params.fmax > params.fmin
            and params.fmin <= params.fref <= params.fmax
        ):
            self.add(
                "ERROR",
                "CFG-CQ-001",
                "anelastic-cQ requires 0 < fmin <= fref <= fmax",
                section="anelastic_cQ_list",
            )
            return None
        if params.n_mechanisms < 4 or params.n_mechanisms > 8:
            self.add(
                "ERROR",
                "CFG-CQ-001",
                "anelastic-cQ n_mechanisms must be one of 4, 5, 6, 7, or 8",
                section="anelastic_cQ_list",
            )
            return None
        if params.nnls_samples < params.n_mechanisms:
            self.add(
                "ERROR",
                "CFG-CQ-001",
                "anelastic-cQ nnls_samples must be at least n_mechanisms",
                section="anelastic_cQ_list",
            )
            return None
        if params.nnls_objective != "relative-q":
            self.add(
                "ERROR",
                "CFG-CQ-001",
                "anelastic-cQ nnls_objective must be relative-q",
                section="anelastic_cQ_list",
            )
            return None
        if params.coefficient_policy == "fixed-q50":
            if params.n_mechanisms != 8:
                self.add(
                    "ERROR",
                    "CFG-CQ-001",
                    "anelastic-cQ fixed-q50 requires n_mechanisms=8",
                    section="anelastic_cQ_list",
                )
                return None
            if abs(params.fmin - 0.05) > 100.0 * MACHINE_EPS or abs(params.fmax - 20.0) > 100.0 * MACHINE_EPS:
                self.add(
                    "ERROR",
                    "CFG-CQ-001",
                    "anelastic-cQ fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz",
                    section="anelastic_cQ_list",
                )
                return None
        elif params.coefficient_policy not in {"nnls-shared", "nnls-block", "nnls-block-ps"}:
            self.add(
                "ERROR",
                "CFG-CQ-001",
                f"unsupported anelastic-cQ coefficient_policy: {params.coefficient_policy}",
                section="anelastic_cQ_list",
            )
            return None
        if not finite_positive(params.nnls_tolerance) or not finite_positive(params.max_fit_error):
            self.add(
                "ERROR",
                "CFG-CQ-001",
                "anelastic-cQ NNLS tolerances must be finite and positive",
                section="anelastic_cQ_list",
            )
            return None
        return params

    def read_fq8(self) -> QParams | None:
        assigns = self._q_assigns(
            "anelastic_fQ8_list",
            "CFG-FQ8-001",
            FQ8_FIELDS,
            "anelastic-fQ8 requires a valid &anelastic_fQ8_list (legacy c is unsupported)",
        )
        if assigns is None:
            return None
        params = QParams(
            kind="fQ8",
            qs0=[-1.0],
            qp0=[-1.0],
            fref=-1.0,
            coefficient_method="conventional-nnls",
            weight_policy="table-exact",
            coarse_grain=-1,
            gamma=-1.0,
            f_transition=-1.0,
        )
        try:
            if raw := optional_assign(assigns, "coefficient_method"):
                params.coefficient_method = fortran_string(raw).strip()
            if raw := optional_assign(assigns, "weight_policy"):
                params.weight_policy = fortran_string(raw).strip()
            if raw := optional_assign(assigns, "coarse_grain"):
                params.coarse_grain = fortran_int(raw)
            if raw := optional_assign(assigns, "qs0"):
                params.qs0 = [fortran_float(raw)]
            if raw := optional_assign(assigns, "qp0"):
                params.qp0 = [fortran_float(raw)]
            if raw := optional_assign(assigns, "gamma"):
                params.gamma = fortran_float(raw)
            if raw := optional_assign(assigns, "f_transition"):
                params.f_transition = fortran_float(raw)
            if raw := optional_assign(assigns, "fref"):
                params.fref = fortran_float(raw)
        except InputError as error:
            self.add("ERROR", "CFG-FQ8-001", str(error), section="anelastic_fQ8_list")
            return None
        suggestion = (
            "Use coarse_grain=2 with coefficient_method='withers-2015', "
            "or use coarse_grain=0 with coefficient_method='conventional-nnls'; "
            "also provide Qs0/Qp0 >= 15 and valid gamma/frequencies."
        )
        if params.coefficient_method not in {"withers-2015", "conventional-nnls"}:
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 coefficient_method must be conventional-nnls or withers-2015",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if params.weight_policy not in {"table-exact", "nonnegative-refit"}:
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 weight_policy must be table-exact or nonnegative-refit",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if params.coefficient_method != "withers-2015" and params.weight_policy != "table-exact":
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 nonnegative-refit applies only to coefficient_method=withers-2015",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if params.coarse_grain == -1:
            params.coarse_grain = 2 if params.coefficient_method == "withers-2015" else 0
        if params.coarse_grain not in (0, 2):
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 coarse_grain must be 0 or 2",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if params.coarse_grain == 0 and params.coefficient_method == "withers-2015":
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 coarse_grain=0 requires coefficient_method=conventional-nnls; "
                "raw withers-2015 strengths are for the 2x2x2 coarse layout",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if not (
            math.isfinite(params.qs0[0])
            and math.isfinite(params.qp0[0])
            and params.qs0[0] >= 15.0
            and params.qp0[0] >= 15.0
        ):
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 requires finite Qs0 and Qp0 >= 15",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if not (math.isfinite(params.gamma) and 0.0 <= params.gamma <= 0.9):
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 gamma must be finite and in [0,0.9]",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if not finite_positive(params.f_transition):
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 f_transition must be finite and positive",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        if not finite_positive(params.fref):
            self.add(
                "ERROR",
                "CFG-FQ8-001",
                "anelastic-fQ8 fref must be finite and positive",
                section="anelastic_fQ8_list",
                suggestion=suggestion,
            )
            return None
        return params

    def _fill_shared_q(self, params: QParams, assigns: dict[str, str], pair: bool = False) -> None:
        if raw := optional_assign(assigns, "weight_method"):
            params.weight_method = fortran_string(raw).strip()
        if raw := optional_assign(assigns, "qs0"):
            params.qs0 = list(require_floats(raw, 2 if pair else 1, "Qs0"))
        if raw := optional_assign(assigns, "qp0"):
            params.qp0 = list(require_floats(raw, 2 if pair else 1, "Qp0"))
        if raw := optional_assign(assigns, "fref"):
            params.fref = fortran_float(raw)
        if raw := optional_assign(assigns, "fmin"):
            params.fmin = fortran_float(raw)
        if raw := optional_assign(assigns, "fmax"):
            params.fmax = fortran_float(raw)

    def read_output(self) -> None:
        try:
            body = namelist_body(self.text, "output_list")
        except InputError as error:
            self.add("ERROR", "CFG-OUTPUT-001", str(error), section="output_list")
            return
        if body is None:
            self.add(
                "ERROR",
                "CFG-OUTPUT-001",
                "Cannot parse &output_list: namelist not found",
                section="output_list",
                suggestion="Fix namelist syntax or unknown fields.",
            )
            return
        try:
            assigns = parse_assignments(body)
        except InputError as error:
            self.add("ERROR", "CFG-OUTPUT-001", f"Cannot parse &output_list: {error}",
                     section="output_list")
            return
        unknown = unknown_fields(assigns, OUTPUT_FIELDS, "output_list")
        if unknown:
            self.add(
                "ERROR",
                "CFG-OUTPUT-001",
                f"Cannot parse &output_list: unknown field {unknown[0]}",
                section="output_list",
                field=unknown[0],
                suggestion="Fix namelist syntax or unknown fields.",
            )
            return
        try:
            if raw := optional_assign(assigns, "output_seismograms"):
                self.output.output_seismograms = fortran_bool(raw)
            if raw := optional_assign(assigns, "station_list"):
                self.output.station_list = fortran_string(raw).strip()
            if raw := optional_assign(assigns, "station_list_file"):
                self.output.station_list_file = fortran_string(raw)
            if raw := optional_assign(assigns, "station_number_in_list"):
                self.output.station_number_in_list = fortran_bool(raw)
        except InputError as error:
            self.add("ERROR", "CFG-OUTPUT-001", f"Cannot parse &output_list: {error}",
                     section="output_list")
            return
        if self.output.output_seismograms:
            self.validate_station_rows()

    def validate_station_rows(self) -> None:
        source_name = self.output.station_list.strip()
        external = source_name in {"extfile", "exfile"}
        if external:
            source = self.output.station_list_file
            if not source.strip():
                self.add(
                    "ERROR",
                    "CFG-STATION-001",
                    "station_list_file is required when station_list=extfile.",
                    section="output_list",
                    field="station_list_file",
                )
                return
            path = Path(source)
            if not path.is_file():
                self.add(
                    "ERROR",
                    "CFG-STATION-001",
                    f"Cannot open station list: {source}",
                    section="station_list",
                    field=source,
                )
                return
            lines = path.read_text(encoding="utf-8").splitlines()
            label = source
        else:
            lines = self.raw.splitlines()
            label = "input station list"

        begin = next(
            (i for i, line in enumerate(lines) if line.strip() == "!---begin:station_list---"),
            None,
        )
        if begin is None:
            if external:
                self.add(
                    "ERROR",
                    "CFG-STATION-001",
                    "!---begin:station_list--- was not found.",
                    section="station_list",
                    field=label,
                    suggestion="Add the required station-list markers.",
                )
            return

        expected = 4 if self.output.station_number_in_list else 3
        row = 0
        for line in lines[begin + 1 :]:
            stripped = line.strip()
            if stripped == "!---end:station_list---":
                return
            row += 1
            fields = stripped.split()
            if len(fields) != expected:
                kind = (
                    " fields: station_number x y z. Row: "
                    if self.output.station_number_in_list
                    else " fields: x y z. Row:                "
                )
                self.add(
                    "ERROR",
                    "CFG-STATION-002",
                    f"Station row {row} must contain exactly {expected}{kind}{line.rstrip()}",
                    section="station_list",
                    field=label,
                    suggestion="Make station_number_in_list match the station row format.",
                )
                return
        self.add(
            "ERROR",
            "CFG-STATION-003",
            "Station list is missing !---end:station_list---.",
            section="station_list",
            field=label,
        )

    def relaxation_dt_limit(self) -> float | None:
        if self.q is None:
            return None
        if self.q.kind == "Q4":
            return q4_relaxation_dt_limit(self.q)
        if self.q.kind == "Q8":
            return q8_relaxation_dt_limit(self.q)
        if self.q.kind == "cQ8-b2":
            return q8_relaxation_dt_limit(self.q)
        if self.q.kind == "cQ":
            return cq_relaxation_dt_limit(self.q)
        if self.q.kind == "fQ8":
            return fq8_relaxation_dt_limit(self.q)
        return None


def print_diagnostic(item: Diagnostic) -> None:
    print(f"\n{item.severity}  {item.code}")
    if item.section:
        print(f"  Section: {item.section}")
    if item.field:
        print(f"  Field: {item.field}")
    if item.block_id > 0:
        print(f"  Block: {item.block_id}")
    print(f"  {item.message}")
    if item.suggestion:
        print(f"  Suggested fix: {item.suggestion}")


def format_ints(values: tuple[int, ...] | list[int]) -> str:
    return " ".join(str(value) for value in values)


def format_owned(nqrs: tuple[int, int, int], dims: tuple[int, int, int]) -> str:
    owned = [nqrs[i] // dims[i] for i in range(3)]
    remainders = [nqrs[i] % dims[i] for i in range(3)]
    parts = []
    for i, axis in enumerate("qrs"):
        extra = f" ({remainders[i]} ranks own +1)" if remainders[i] else ""
        parts.append(f"{axis}:{owned[i]}{extra}")
    return ", ".join(parts)


def print_decomposition(check: Preflight) -> None:
    print("Resolved MPI decomposition:")
    print(f"  processes: {check.np}")
    if check.stencil.supported:
        print(
            f"  stencil: {check.stencil.operator_name}  "
            f"halo={check.stencil.halo_width}  "
            f"boundary={check.stencil.boundary_width}  "
            f"min global points={check.stencil.minimum_global_points}  "
            f"min owned={check.stencil.minimum_owned_points}"
        )
    if check.decomp is None:
        print("  no valid Cartesian topology for this np")
        return
    for block, dims, begin, end in zip(
        check.blocks, check.decomp.process_dims, check.decomp.rank_begin, check.decomp.rank_end
    ):
        nq, nr, ns = block.nqrs
        pq, pr, ps = dims
        print(
            f"  block {block.block_id}: grid {nq} {nr} {ns}  "
            f"topology {pq} {pr} {ps}  ranks {begin}-{end}"
        )
        print(f"           floor owned points/rank  {format_owned(block.nqrs, dims)}")
        print(
            f"           ranks on block: {end - begin + 1}  "
            f"points: {math.prod(block.nqrs)}"
        )
    if check.decomp.serial_shared_blocks:
        print("  serial two-block mode: rank 0 owns both blocks")


def print_time_parameters(check: Preflight) -> None:
    print("Resolved time parameters:")
    print(f"  problem name: {check.problem.name}")
    print(f"  response: {check.problem.response}")
    print(f"  CFL: {check.problem.cfl:.17g}")
    print(f"  mesh: {check.problem.type_of_mesh}  (nominal aqrs/bqrs/nqrs spacing)")
    print(f"  requested final time: {check.problem.t_final:.17g}")

    elastic_limit = math.inf
    limiting_block = 0
    for block in check.blocks:
        spacing = block.spacing
        block_limit = block.elastic_dt(check.problem.cfl)
        print(
            f"  block {block.block_id}: spacing {spacing[0]:.17g}  "
            f"{spacing[1]:.17g}  {spacing[2]:.17g}"
        )
        print(
            f"           h=min(spacing)={block.h:.17g}  "
            f"sqrt(Vs^2+Vp^2)={block.wave_speed_factor:.17g}  "
            f"elastic dt limit {block_limit:.17g}"
        )
        if block_limit < elastic_limit:
            elastic_limit = block_limit
            limiting_block = block.block_id

    relaxation_limit = check.relaxation_dt_limit()
    dt_limit = elastic_limit
    if relaxation_limit is not None:
        dt_limit = min(elastic_limit, relaxation_limit)

    if not math.isfinite(dt_limit) or dt_limit <= 0.0:
        print("  selected dt: unavailable")
        return

    nt = math.floor(check.problem.t_final / dt_limit) if dt_limit > 0.0 else 0
    covered = nt * dt_limit
    print(f"  elastic CFL limit: {elastic_limit:.17g} (block {limiting_block})")
    if relaxation_limit is not None:
        print(f"  relaxation limit: {relaxation_limit:.17g}")
        limiter = "relaxation" if relaxation_limit < elastic_limit else "elastic CFL"
        print(f"  limiting constraint: {limiter}")
    print(f"  selected dt: {dt_limit:.17g}")
    print(f"  number of time steps: {nt}")
    print(f"  time reached after nt steps: {covered:.17g}")
    print(f"  remainder to requested final time: {check.problem.t_final - covered:.17g}")
    if check.problem.type_of_mesh.lower() != "cartesian":
        print(
            "  note: this uses the same nominal-spacing formula as block_time_step(); "
            "it does not measure true cell sizes in a curvilinear mesh."
        )


def print_q_parameters(check: Preflight) -> None:
    q = check.q
    if q is None:
        return
    if q.kind == "cQ8-b2":
        print("anelastic-cQ8-b2 parameters:")
        print(f"  block 1: Qs0 = {q.qs0[0]:.4E}, Qp0 = {q.qp0[0]:.4E}")
        print(f"  block 2: Qs0 = {q.qs0[1]:.4E}, Qp0 = {q.qp0[1]:.4E}")
    elif q.kind == "cQ":
        print(f"anelastic-cQ mechanisms: {q.n_mechanisms}")
        print(f"  block 1: Qs0={q.qs0[0]:.4E}, Qp0={q.qp0[0]:.4E}")
        print(f"  block 2: Qs0={q.qs0[1]:.4E}, Qp0={q.qp0[1]:.4E}")


def report(check: Preflight) -> int:
    print(f"Input:      {check.path.resolve()}")
    print(f"Processes:  {check.np}")
    print()

    if check.issues:
        title = (
            f"Input preflight passed with warnings: {check.path}"
            if not check.has_errors()
            else f"Input preflight failed: {check.path}"
        )
        print(title)
        for item in check.issues:
            print_diagnostic(item)
        print(f"Summary: {check.error_count()} errors, {check.warning_count()} warnings.")
    else:
        print("Input preflight passed: 0 errors, 0 warnings.")

    if check.q is not None and check.q.kind in {"cQ8-b2", "cQ"}:
        print()
        print_q_parameters(check)

    if check.blocks:
        print()
        print_decomposition(check)

    geometry_ok = bool(check.blocks) and all(
        all(count >= 2 for count in block.nqrs)
        and all(upper > lower for lower, upper in zip(block.aqrs, block.bqrs))
        and all(value > 0.0 for value in block.rho_s_p)
        for block in check.blocks
    )
    if geometry_ok and check.problem.cfl > 0.0:
        print()
        print_time_parameters(check)

    return 2 if check.has_errors() else 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check a WaveQLab3D-Q .in file against src/input_preflight.f90, "
            "print the MPI decomposition for np ranks, and print dt/nt from "
            "the src/block_time_step CFL formula."
        )
    )
    parser.add_argument("input_file", type=Path, metavar="fname.in")
    parser.add_argument("np", type=int, help="number of MPI processes")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_file = args.input_file.expanduser()
    if not input_file.is_file():
        print(f"error: input file not found: {input_file}", file=sys.stderr)
        return 2
    if args.np < 1:
        print("error: np must be a positive integer", file=sys.stderr)
        return 2
    try:
        check = Preflight(input_file, args.np)
        check.run()
    except (OSError, UnicodeError, InputError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return report(check)


if __name__ == "__main__":
    raise SystemExit(main())
