#!/usr/bin/env python3
"""Project Euler batch-mode orchestrator.

This tool maintains `batch_ledger.json`, selects batch candidates, tracks
blocked problems, and records push boundaries.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

LEDGER_VERSION = 1
RSS_URL = "https://projecteuler.net/rss2_euler.xml"
EULER_FILE_RE = re.compile(r"^Euler([0-9]+)\.cpp$")
PROBLEM_TITLE_RE = re.compile(r"<title>Problem\s+([0-9]+)</title>")


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def run_cmd(args: List[str], cwd: Path) -> str:
    proc = subprocess.run(args, cwd=str(cwd), check=True, capture_output=True, text=True)
    return proc.stdout


def detect_current_branch(repo_root: Path) -> str:
    return run_cmd(["git", "rev-parse", "--abbrev-ref", "HEAD"], repo_root).strip()


def detect_head_commit(repo_root: Path) -> str:
    return run_cmd(["git", "rev-parse", "HEAD"], repo_root).strip()


def detect_untracked(repo_root: Path) -> List[str]:
    out = run_cmd(["git", "ls-files", "--others", "--exclude-standard"], repo_root)
    values = [line.strip() for line in out.splitlines() if line.strip()]
    return sorted(set(values))


def discover_solved_ids(repo_root: Path) -> List[int]:
    solved: List[int] = []
    for entry in repo_root.iterdir():
        if not entry.is_file():
            continue
        match = EULER_FILE_RE.match(entry.name)
        if match:
            solved.append(int(match.group(1)))
    return sorted(set(solved))


def fetch_latest_problem_id() -> int:
    with urllib.request.urlopen(RSS_URL, timeout=30) as response:
        payload = response.read().decode("utf-8", errors="replace")
    ids = [int(m.group(1)) for m in PROBLEM_TITLE_RE.finditer(payload)]
    if not ids:
        raise RuntimeError("Could not parse latest problem id from rss feed")
    return max(ids)


def new_ledger(
    ceiling: int,
    batch_size: int,
    problem_budget_min: int,
    retry_slots: int,
    branch: str,
    completed_ids: List[int],
) -> Dict:
    return {
        "version": LEDGER_VERSION,
        "run_started_at_utc": utc_now_iso(),
        "ceiling_id": ceiling,
        "batch_size": batch_size,
        "problem_budget_min": problem_budget_min,
        "retry_slots_per_batch": retry_slots,
        "branch": branch,
        "rss_ceiling_snapshot": None,
        "baseline_untracked": [],
        "completed_ids": sorted(set(completed_ids)),
        "blocked": [],
        "batches": [],
    }


def load_ledger(path: Path) -> Dict:
    if not path.exists():
        raise FileNotFoundError(f"Ledger not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    validate_ledger(data)
    return data


def save_ledger(path: Path, ledger: Dict) -> None:
    normalize_ledger(ledger)
    with path.open("w", encoding="utf-8") as f:
        json.dump(ledger, f, indent=2, sort_keys=False)
        f.write("\n")


def normalize_ledger(ledger: Dict) -> None:
    ledger["completed_ids"] = sorted(set(int(v) for v in ledger.get("completed_ids", [])))

    blocked_map: Dict[int, Dict] = {}
    for row in ledger.get("blocked", []):
        pid = int(row["id"])
        blocked_map[pid] = {
            "id": pid,
            "attempts": int(row.get("attempts", 0)),
            "last_error": str(row.get("last_error", "")),
            "last_attempt_utc": str(row.get("last_attempt_utc", "")),
            "last_batch": int(row.get("last_batch", 0)),
        }
    ledger["blocked"] = sorted(blocked_map.values(), key=lambda x: (x["attempts"], x["id"]))

    normalized_batches = []
    for raw in ledger.get("batches", []):
        normalized_batches.append(
            {
                "batch_index": int(raw["batch_index"]),
                "solved_ids": sorted(set(int(v) for v in raw.get("solved_ids", []))),
                "skipped_ids": sorted(set(int(v) for v in raw.get("skipped_ids", []))),
                "push_commit": raw.get("push_commit"),
            }
        )
    ledger["batches"] = sorted(normalized_batches, key=lambda x: x["batch_index"])


def validate_ledger(ledger: Dict) -> None:
    required = [
        "run_started_at_utc",
        "ceiling_id",
        "batch_size",
        "problem_budget_min",
        "retry_slots_per_batch",
        "completed_ids",
        "blocked",
        "batches",
    ]
    for key in required:
        if key not in ledger:
            raise ValueError(f"Ledger missing key: {key}")


def ensure_batch(ledger: Dict, batch_index: int) -> Dict:
    for row in ledger["batches"]:
        if row["batch_index"] == batch_index:
            return row
    row = {
        "batch_index": batch_index,
        "solved_ids": [],
        "skipped_ids": [],
        "push_commit": None,
    }
    ledger["batches"].append(row)
    return row


def active_batch_index(ledger: Dict) -> int:
    if not ledger["batches"]:
        return 1
    tail = max(ledger["batches"], key=lambda x: x["batch_index"])
    if tail.get("push_commit"):
        return tail["batch_index"] + 1
    return tail["batch_index"]


def select_candidates(ledger: Dict, ceiling: int, batch_size: int, retry_slots: int) -> Tuple[int, List[int]]:
    batch_index = active_batch_index(ledger)
    batch = ensure_batch(ledger, batch_index)

    solved = set(int(v) for v in ledger["completed_ids"])
    for p in batch["solved_ids"]:
        solved.add(int(p))

    needed = max(0, batch_size - len(batch["solved_ids"]))
    if needed == 0:
        return batch_index, []

    blocked_rows = [b for b in ledger["blocked"] if int(b["id"]) not in solved and int(b["id"]) <= ceiling]
    blocked_rows.sort(key=lambda b: (int(b["attempts"]), int(b["id"])))

    blocked_ids = [int(r["id"]) for r in blocked_rows]
    retry_take = min(retry_slots, needed, len(blocked_ids))
    chosen = blocked_ids[:retry_take]

    blocked_set = set(blocked_ids)
    fresh = [pid for pid in range(1, ceiling + 1) if pid not in solved and pid not in blocked_set]

    remaining = needed - len(chosen)
    chosen.extend(fresh[:remaining])
    return batch_index, chosen


def mark_solved(ledger: Dict, pid: int, batch_index: int) -> None:
    if pid not in ledger["completed_ids"]:
        ledger["completed_ids"].append(pid)

    ledger["blocked"] = [b for b in ledger["blocked"] if int(b["id"]) != pid]

    batch = ensure_batch(ledger, batch_index)
    if pid not in batch["solved_ids"]:
        batch["solved_ids"].append(pid)
    if pid in batch["skipped_ids"]:
        batch["skipped_ids"].remove(pid)


def mark_blocked(ledger: Dict, pid: int, batch_index: int, error: str) -> None:
    if pid in ledger["completed_ids"]:
        return

    found = None
    for row in ledger["blocked"]:
        if int(row["id"]) == pid:
            found = row
            break

    now = utc_now_iso()
    if found is None:
        found = {
            "id": pid,
            "attempts": 1,
            "last_error": error,
            "last_attempt_utc": now,
            "last_batch": batch_index,
        }
        ledger["blocked"].append(found)
    else:
        found["attempts"] = int(found.get("attempts", 0)) + 1
        found["last_error"] = error
        found["last_attempt_utc"] = now
        found["last_batch"] = batch_index

    batch = ensure_batch(ledger, batch_index)
    if pid not in batch["skipped_ids"]:
        batch["skipped_ids"].append(pid)


def mark_batch_pushed(ledger: Dict, batch_index: int, push_commit: str) -> None:
    batch = ensure_batch(ledger, batch_index)
    batch["push_commit"] = push_commit


def cleanup_untracked(repo_root: Path, ledger: Dict) -> List[str]:
    baseline = set(str(v) for v in ledger.get("baseline_untracked", []))
    current = set(detect_untracked(repo_root))
    to_remove = sorted(current - baseline)

    removed: List[str] = []
    for rel in to_remove:
        path = repo_root / rel
        if not path.exists():
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
        removed.append(rel)
    return removed


def status_line(ledger: Dict) -> str:
    ceiling = int(ledger["ceiling_id"])
    done = len([x for x in ledger["completed_ids"] if int(x) <= ceiling])
    pending = max(0, ceiling - done)
    blocked = len([b for b in ledger["blocked"] if int(b["id"]) <= ceiling])
    return f"completed={done} pending={pending} blocked={blocked} batches={len(ledger['batches'])}"


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Project Euler batch runner")
    parser.add_argument("--ledger", default="batch_ledger.json")
    parser.add_argument("--ceiling", type=int, default=983)
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--problem-budget-min", type=int, default=45)
    parser.add_argument("--retry-slots", type=int, default=2)
    parser.add_argument("--branch", default="main")

    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("init", help="Create or refresh ledger")
    sub.add_parser("preflight", help="Snapshot rss ceiling + untracked baseline")
    sub.add_parser("status", help="Print summary")

    next_batch_parser = sub.add_parser("next-batch", help="Select candidates for active batch")
    next_batch_parser.add_argument("--json", action="store_true")

    solved_parser = sub.add_parser("mark-solved", help="Record solved problem")
    solved_parser.add_argument("--id", required=True, type=int)
    solved_parser.add_argument("--batch-index", required=True, type=int)

    blocked_parser = sub.add_parser("mark-blocked", help="Record blocked problem")
    blocked_parser.add_argument("--id", required=True, type=int)
    blocked_parser.add_argument("--batch-index", required=True, type=int)
    blocked_parser.add_argument("--error", required=True)

    pushed_parser = sub.add_parser("mark-pushed", help="Attach push commit to batch")
    pushed_parser.add_argument("--batch-index", required=True, type=int)
    pushed_parser.add_argument("--commit", required=False)

    sub.add_parser("cleanup-untracked", help="Remove newly-created untracked artifacts")
    return parser.parse_args(argv)


def cmd_init(repo_root: Path, args: argparse.Namespace) -> int:
    ledger_path = repo_root / args.ledger
    completed = discover_solved_ids(repo_root)

    if ledger_path.exists():
        ledger = load_ledger(ledger_path)
        ledger["ceiling_id"] = args.ceiling
        ledger["batch_size"] = args.batch_size
        ledger["problem_budget_min"] = args.problem_budget_min
        ledger["retry_slots_per_batch"] = args.retry_slots
        ledger["branch"] = args.branch
        merged = set(int(v) for v in ledger["completed_ids"])
        merged.update(completed)
        ledger["completed_ids"] = sorted(merged)
    else:
        ledger = new_ledger(
            ceiling=args.ceiling,
            batch_size=args.batch_size,
            problem_budget_min=args.problem_budget_min,
            retry_slots=args.retry_slots,
            branch=args.branch,
            completed_ids=completed,
        )

    save_ledger(ledger_path, ledger)
    print(f"initialized ledger: {ledger_path}")
    print(status_line(ledger))
    return 0


def cmd_preflight(repo_root: Path, args: argparse.Namespace) -> int:
    ledger_path = repo_root / args.ledger
    ledger = load_ledger(ledger_path)

    latest = fetch_latest_problem_id()
    ledger["rss_ceiling_snapshot"] = latest

    if latest < args.ceiling:
        raise RuntimeError(
            f"Requested ceiling {args.ceiling} exceeds latest published problem {latest}"
        )

    if latest != args.ceiling:
        print(
            f"warning: rss latest is {latest}, but run ceiling is fixed to {args.ceiling}",
            file=sys.stderr,
        )

    current_branch = detect_current_branch(repo_root)
    if current_branch != args.branch:
        raise RuntimeError(f"Current branch is {current_branch}, expected {args.branch}")

    if not ledger.get("baseline_untracked"):
        ledger["baseline_untracked"] = detect_untracked(repo_root)

    save_ledger(ledger_path, ledger)
    print(f"preflight ok: latest={latest} ceiling={args.ceiling} branch={current_branch}")
    return 0


def cmd_status(repo_root: Path, args: argparse.Namespace) -> int:
    ledger = load_ledger(repo_root / args.ledger)
    print(status_line(ledger))
    batch_idx, candidates = select_candidates(
        ledger,
        int(ledger["ceiling_id"]),
        int(ledger["batch_size"]),
        int(ledger["retry_slots_per_batch"]),
    )
    print(f"active_batch={batch_idx}")
    print("next_candidates=" + ",".join(str(v) for v in candidates))
    return 0


def cmd_next_batch(repo_root: Path, args: argparse.Namespace) -> int:
    ledger = load_ledger(repo_root / args.ledger)
    batch_idx, candidates = select_candidates(
        ledger,
        int(ledger["ceiling_id"]),
        int(ledger["batch_size"]),
        int(ledger["retry_slots_per_batch"]),
    )

    if args.json:
        print(json.dumps({"batch_index": batch_idx, "candidates": candidates}))
    else:
        print(f"batch_index={batch_idx}")
        print("candidates=" + ",".join(str(v) for v in candidates))
    return 0


def cmd_mark_solved(repo_root: Path, args: argparse.Namespace) -> int:
    ledger_path = repo_root / args.ledger
    ledger = load_ledger(ledger_path)
    mark_solved(ledger, args.id, args.batch_index)
    save_ledger(ledger_path, ledger)
    print(f"marked solved: id={args.id} batch={args.batch_index}")
    return 0


def cmd_mark_blocked(repo_root: Path, args: argparse.Namespace) -> int:
    ledger_path = repo_root / args.ledger
    ledger = load_ledger(ledger_path)
    mark_blocked(ledger, args.id, args.batch_index, args.error)
    save_ledger(ledger_path, ledger)
    print(f"marked blocked: id={args.id} batch={args.batch_index}")
    return 0


def cmd_mark_pushed(repo_root: Path, args: argparse.Namespace) -> int:
    ledger_path = repo_root / args.ledger
    ledger = load_ledger(ledger_path)
    commit = args.commit or detect_head_commit(repo_root)
    mark_batch_pushed(ledger, args.batch_index, commit)
    save_ledger(ledger_path, ledger)
    print(f"marked pushed: batch={args.batch_index} commit={commit}")
    return 0


def cmd_cleanup_untracked(repo_root: Path, args: argparse.Namespace) -> int:
    ledger_path = repo_root / args.ledger
    ledger = load_ledger(ledger_path)
    removed = cleanup_untracked(repo_root, ledger)
    if removed:
        print("removed:")
        for rel in removed:
            print(rel)
    else:
        print("removed: none")
    return 0


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    repo_root = Path.cwd()

    dispatch = {
        "init": cmd_init,
        "preflight": cmd_preflight,
        "status": cmd_status,
        "next-batch": cmd_next_batch,
        "mark-solved": cmd_mark_solved,
        "mark-blocked": cmd_mark_blocked,
        "mark-pushed": cmd_mark_pushed,
        "cleanup-untracked": cmd_cleanup_untracked,
    }

    func = dispatch[args.cmd]
    return func(repo_root, args)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
