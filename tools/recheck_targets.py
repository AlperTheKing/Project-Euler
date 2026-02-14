#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import subprocess
import urllib.request
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path

getcontext().prec = 80

SOLUTIONS_URL = "https://raw.githubusercontent.com/lucky-bai/projecteuler-solutions/master/Solutions.md"
LINE_RE = re.compile(r"^(\d+)\.\s*(.+?)\s*$")
INT_RE = re.compile(r"^[+-]?\d+$")
FRAC_RE = re.compile(r"^[+-]?\d+/[1-9]\d*$")
HEX_SPECIAL_RE = re.compile(r"^[0-9A-Fa-f]+(?:,[0-9A-Fa-f]+)*$")
NUM_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
TEXT_RE = re.compile(r"^[A-Za-z]+$")

NOISE_KEYWORDS = (
    "elapsed",
    "checkpoint",
    "pass",
    "fail",
    "validation",
    "test passed",
    "calculating",
    "build warning",
    "using",
    "scheduler",
    "probability mass",
    "states in canonical graph",
    "time:",
)


def fetch_expected() -> dict[int, str]:
    with urllib.request.urlopen(SOLUTIONS_URL, timeout=30) as r:
        text = r.read().decode("utf-8", errors="replace")
    out: dict[int, str] = {}
    for line in text.splitlines():
        m = LINE_RE.match(line.strip())
        if m:
            out[int(m.group(1))] = m.group(2).strip()
    return out


def is_atom(token: str) -> bool:
    return (
        FRAC_RE.fullmatch(token) is not None
        or NUM_RE.fullmatch(token) is not None
        or HEX_SPECIAL_RE.fullmatch(token) is not None
        or TEXT_RE.fullmatch(token) is not None
    )


def normalize_candidate(token: str) -> str:
    token = token.strip().strip("()[]{}")
    if token.startswith('"') and token.endswith('"') and len(token) >= 2:
        token = token[1:-1]
    return token.strip()


def parse_stdout(stdout: str) -> str | None:
    lines = [ln.strip() for ln in stdout.splitlines() if ln.strip()]
    for line in reversed(lines):
        low = line.lower()
        if any(k in low for k in NOISE_KEYWORDS) and "answer" not in low and "f(20!)" not in low:
            continue

        segments = [line]
        if "=" in line:
            segments.append(line.rsplit("=", 1)[1].strip())
        if ":" in line:
            segments.append(line.rsplit(":", 1)[1].strip())

        for segment in segments:
            cand = normalize_candidate(segment)
            if not cand:
                continue

            if is_atom(cand):
                return cand

            if "," in cand:
                parts = [normalize_candidate(p) for p in cand.split(",")]
                if all(parts) and all(is_atom(p) for p in parts):
                    return ",".join(parts)

            tokens = [normalize_candidate(t) for t in cand.split()]
            for tok in reversed(tokens):
                if tok and is_atom(tok):
                    return tok
    return None


def mantissa_decimals(value: str) -> int:
    head = re.split(r"[eE]", value, 1)[0]
    if "." not in head:
        return 0
    return len(head.split(".", 1)[1])


def norm_special(value: str) -> str:
    return ",".join(part.strip().upper() for part in value.split(","))


def compare_values(expected: str, actual: str) -> bool:
    e = expected.strip()
    a = actual.strip()

    if FRAC_RE.fullmatch(e):
        en, ed = e.split("/", 1)
        en_i = int(en)
        ed_i = int(ed)
        if FRAC_RE.fullmatch(a):
            an, ad = a.split("/", 1)
            return en_i * int(ad) == int(an) * ed_i
        try:
            da = Decimal(a)
        except InvalidOperation:
            return False
        return da == (Decimal(en_i) / Decimal(ed_i))

    if INT_RE.fullmatch(e) and INT_RE.fullmatch(a):
        return int(e) == int(a)

    if HEX_SPECIAL_RE.fullmatch(e) and HEX_SPECIAL_RE.fullmatch(a):
        return norm_special(e) == norm_special(a)

    try:
        de = Decimal(e)
        da = Decimal(a)
        d = mantissa_decimals(e)
        if d > 0:
            abs_tol = Decimal("0.5") * (Decimal(10) ** Decimal(-d))
        else:
            abs_tol = Decimal("1e-12")
        rel_tol = Decimal("1e-12")
        diff = abs(da - de)
        tol = max(abs_tol, rel_tol * max(abs(de), Decimal(1)))
        return diff <= tol
    except InvalidOperation:
        return e == a


def run_problem(repo_root: Path, pid: int, cc: str, timeout_sec: int) -> tuple[str, str, str | None]:
    src = repo_root / f"Euler{pid}.cpp"
    if not src.exists():
        return ("UNVERIFIED", "MISSING_SOURCE", None)

    binary = Path(f"/tmp/pe_target_{pid}")
    comp = subprocess.run(
        [cc, "-O2", "-std=c++20", "-pthread", str(src), "-o", str(binary)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if comp.returncode != 0:
        binary.unlink(missing_ok=True)
        return ("UNVERIFIED", "COMPILE_ERROR", None)

    try:
        run = subprocess.run(
            [str(binary)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout_sec,
        )
    except subprocess.TimeoutExpired:
        binary.unlink(missing_ok=True)
        return ("UNVERIFIED", "TIMEOUT", None)
    finally:
        binary.unlink(missing_ok=True)

    if run.returncode != 0:
        return ("UNVERIFIED", "RUNTIME_ERROR", None)

    parsed = parse_stdout(run.stdout)
    if parsed is None:
        return ("UNVERIFIED", "NO_PARSE", None)
    return ("OK", "PARSED", parsed)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--targets",
        default=str(Path(__file__).resolve().with_name("recheck_targets.json")),
    )
    parser.add_argument("--compiler", default="/opt/homebrew/bin/g++-15")
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    targets_path = Path(args.targets).resolve()
    data = json.loads(targets_path.read_text(encoding="utf-8"))
    expected = fetch_expected()

    mismatch_list: list[tuple[int, str, str]] = []
    unverified_list: list[tuple[int, str]] = []
    matched = 0
    checked = 0

    for section in data.get("order", []):
        ids = data.get(f"{section}_ids", [])
        print(f"[{section}]")
        for pid in ids:
            pid = int(pid)
            checked += 1
            exp = expected.get(pid)
            if exp is None:
                unverified_list.append((pid, "MISSING_EXPECTED"))
                print(f"{pid}\tUNVERIFIED\tMISSING_EXPECTED")
                continue
            status, reason, actual = run_problem(repo_root, pid, args.compiler, args.timeout)
            if status == "UNVERIFIED":
                unverified_list.append((pid, reason))
                print(f"{pid}\tUNVERIFIED\t{reason}")
                continue
            if compare_values(exp, actual or ""):
                matched += 1
                print(f"{pid}\tMATCH")
            else:
                mismatch_list.append((pid, exp, actual or ""))
                print(f"{pid}\tMISMATCH\texpected={exp}\tactual={actual}")

    print("MISMATCH_LIST")
    if mismatch_list:
        for pid, exp, act in mismatch_list:
            print(f"{pid}\texpected={exp}\tactual={act}")
    else:
        print("(none)")

    print("UNVERIFIED_LIST")
    if unverified_list:
        for pid, reason in unverified_list:
            print(f"{pid}\treason={reason}")
    else:
        print("(none)")

    print("SUMMARY")
    print(
        f"checked={checked}\tmatched={matched}\tmismatch={len(mismatch_list)}\tunverified={len(unverified_list)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
