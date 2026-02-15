#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

DEFAULT_START = 1
DEFAULT_END = 984
DEFAULT_THRESHOLD = 3.0
DEFAULT_TIMEOUT = 120
DEFAULT_STD = "c++20"
DEFAULT_OPT = "-O2"
DEFAULT_HOMEBREW_COMPILER = "/opt/homebrew/bin/g++-15"
REQUIRED_RUNTIME_PATHS = (
    "resources",
    "p424_kakuro200.txt",
    "I-expressions.txt",
)


def choose_default_compiler() -> str:
    homebrew = Path(DEFAULT_HOMEBREW_COMPILER)
    if homebrew.exists():
        return str(homebrew)
    fallback = shutil.which("g++")
    if fallback:
        return fallback
    return "g++"


def resolve_compiler(compiler: str) -> str:
    candidate = Path(compiler)
    if candidate.is_absolute():
        if not candidate.exists():
            raise FileNotFoundError(f"Compiler not found: {compiler}")
        return str(candidate)

    resolved = shutil.which(compiler)
    if resolved is None:
        raise FileNotFoundError(f"Compiler not found in PATH: {compiler}")
    return resolved


def parse_opt_flags(opt: str) -> list[str]:
    flags = shlex.split(opt)
    if not flags:
        raise ValueError("--opt must contain at least one compiler flag")
    return flags


def validate_range(start: int, end: int) -> None:
    if start < 1:
        raise ValueError("--start must be >= 1")
    if end < start:
        raise ValueError("--end must be >= --start")


def validate_sources(repo_root: Path, start: int, end: int) -> None:
    missing: list[int] = []
    for pid in range(start, end + 1):
        src = repo_root / f"Euler{pid}.cpp"
        if not src.exists():
            missing.append(pid)
    if missing:
        preview = ",".join(str(v) for v in missing[:20])
        suffix = "" if len(missing) <= 20 else f"...(+{len(missing) - 20} more)"
        raise FileNotFoundError(f"Missing source files in requested range: {preview}{suffix}")


def link_or_copy(src: Path, dst: Path) -> None:
    if dst.exists() or dst.is_symlink():
        return
    try:
        dst.symlink_to(src, target_is_directory=src.is_dir())
    except OSError:
        if src.is_dir():
            shutil.copytree(src, dst, symlinks=True)
        else:
            shutil.copy2(src, dst)


def prepare_runtime_tree(repo_root: Path, runtime_root: Path) -> None:
    for rel in REQUIRED_RUNTIME_PATHS:
        src = repo_root / rel
        if not src.exists():
            raise FileNotFoundError(f"Required runtime path not found: {src}")
        dst = runtime_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        link_or_copy(src, dst)


def compile_problem(
    compiler: str,
    opt_flags: list[str],
    std: str,
    src: Path,
    binary: Path,
) -> subprocess.CompletedProcess[str]:
    cmd = [compiler, *opt_flags, f"-std={std}", "-pthread", str(src), "-o", str(binary)]
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def run_problem(binary: Path, cwd: Path, timeout_sec: int) -> tuple[str, float, int | None]:
    start = time.perf_counter()
    try:
        proc = subprocess.run(
            [str(binary)],
            cwd=str(cwd),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout_sec,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return ("TIMEOUT", float(timeout_sec), None)

    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        return ("RUNTIME_ERROR", elapsed, proc.returncode)
    return ("OK", elapsed, 0)


def classify_records(records: list[dict[str, Any]], threshold: float) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    slow_list: list[dict[str, Any]] = []
    error_list: list[dict[str, Any]] = []
    ok_list: list[dict[str, Any]] = []

    for row in records:
        status = str(row["status"])
        runtime = row.get("runtime_s")

        if status in {"COMPILE_ERROR", "RUNTIME_ERROR"}:
            error_list.append(row)
            continue

        if status == "TIMEOUT":
            slow_list.append(row)
            continue

        if status == "OK" and isinstance(runtime, (int, float)):
            if runtime > threshold:
                slow_list.append(row)
            else:
                ok_list.append(row)

    slow_list.sort(key=lambda r: (-float(r["runtime_s"]), int(r["id"])))
    return slow_list, error_list, ok_list


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Sequential Project Euler benchmark runner")
    parser.add_argument("--start", type=int, default=DEFAULT_START)
    parser.add_argument("--end", type=int, default=DEFAULT_END)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    parser.add_argument("--compiler", default=choose_default_compiler())
    parser.add_argument("--std", default=DEFAULT_STD)
    parser.add_argument("--opt", default=DEFAULT_OPT)
    parser.add_argument("--json-out")
    return parser


def main() -> int:
    args = make_parser().parse_args()

    try:
        validate_range(args.start, args.end)
        if args.threshold < 0:
            raise ValueError("--threshold must be >= 0")
        if args.timeout <= 0:
            raise ValueError("--timeout must be > 0")

        repo_root = Path(__file__).resolve().parents[1]
        validate_sources(repo_root, args.start, args.end)

        compiler = resolve_compiler(args.compiler)
        opt_flags = parse_opt_flags(args.opt)
    except (ValueError, FileNotFoundError) as exc:
        print(f"error: {exc}")
        return 2

    records: list[dict[str, Any]] = []
    start_wall = time.perf_counter()

    with tempfile.TemporaryDirectory(prefix="pe_benchmark_") as tmpdir:
        temp_root = Path(tmpdir)
        runtime_root = temp_root / "runtime"
        runtime_root.mkdir(parents=True, exist_ok=True)
        prepare_runtime_tree(repo_root, runtime_root)
        bin_dir = temp_root / "bin"
        bin_dir.mkdir(parents=True, exist_ok=True)

        for pid in range(args.start, args.end + 1):
            src = repo_root / f"Euler{pid}.cpp"
            binary = bin_dir / f"Euler{pid}"

            comp = compile_problem(compiler, opt_flags, args.std, src, binary)
            if comp.returncode != 0:
                records.append(
                    {
                        "id": pid,
                        "source": str(src),
                        "status": "COMPILE_ERROR",
                        "runtime_s": None,
                        "compile_returncode": comp.returncode,
                        "runtime_returncode": None,
                    }
                )
                binary.unlink(missing_ok=True)
                continue

            status, runtime_s, runtime_rc = run_problem(binary, runtime_root, args.timeout)
            records.append(
                {
                    "id": pid,
                    "source": str(src),
                    "status": status,
                    "runtime_s": runtime_s,
                    "compile_returncode": 0,
                    "runtime_returncode": runtime_rc,
                }
            )
            binary.unlink(missing_ok=True)

    total_wall = time.perf_counter() - start_wall
    slow_list, error_list, ok_list = classify_records(records, args.threshold)

    ok_count = sum(1 for r in records if r["status"] == "OK")
    timeout_count = sum(1 for r in records if r["status"] == "TIMEOUT")
    compile_error_count = sum(1 for r in records if r["status"] == "COMPILE_ERROR")
    runtime_error_count = sum(1 for r in records if r["status"] == "RUNTIME_ERROR")

    print(f"RANGE Euler{args.start}..Euler{args.end}")
    print(f"COMPILER {compiler}")
    print(f"THRESHOLD {args.threshold:.3f}s")
    print(f"TIMEOUT {args.timeout}s")
    print("SLOW_LIST")
    if slow_list:
        for row in slow_list:
            print(
                f"Euler{row['id']} | runtime={float(row['runtime_s']):.6f}s | status={row['status']}"
            )
    else:
        print("(none)")

    print("ERROR_LIST")
    if error_list:
        for row in error_list:
            print(
                f"Euler{row['id']} | status={row['status']} | compile_rc={row['compile_returncode']} | runtime_rc={row['runtime_returncode']}"
            )
    else:
        print("(none)")

    print("SUMMARY")
    print(f"total={len(records)}")
    print(f"ok={ok_count}")
    print(f"timeout={timeout_count}")
    print(f"compile_error={compile_error_count}")
    print(f"runtime_error={runtime_error_count}")
    print(f"over_threshold={len(slow_list)}")
    print(f"ok_under_threshold={len(ok_list)}")
    print(f"wall_clock_total_s={total_wall:.6f}")

    if args.json_out:
        payload = {
            "config": {
                "start": args.start,
                "end": args.end,
                "threshold": args.threshold,
                "timeout": args.timeout,
                "compiler": compiler,
                "std": args.std,
                "opt_flags": opt_flags,
            },
            "summary": {
                "total": len(records),
                "ok": ok_count,
                "timeout": timeout_count,
                "compile_error": compile_error_count,
                "runtime_error": runtime_error_count,
                "over_threshold": len(slow_list),
                "ok_under_threshold": len(ok_list),
                "wall_clock_total_s": total_wall,
            },
            "slow_list": slow_list,
            "error_list": error_list,
            "ok_list": ok_list,
            "records": records,
        }
        out_path = Path(args.json_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"json_written={out_path.resolve()}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
