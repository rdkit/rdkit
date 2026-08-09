#!/usr/bin/env bash

set -euo pipefail

readonly SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}")
readonly SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
readonly SOURCE_ROOT=$(readlink -f "$SCRIPT_DIR/../..")
readonly RUN_USER=${RDK_RUN_USER:-${SUDO_USER:-$(id -un)}}
readonly INPUT=${RDK_SMILES_INPUT:-/home/kevin/data/enamine_scrambled.cxsmiles}
readonly COUNT=${RDK_SMILES_COUNT:-135000}
readonly ROUNDS=${RDK_BENCHMARK_ROUNDS:-8}
readonly PERF_CPU=${RDK_MEASUREMENT_CPU:-27}
readonly RESERVED_CPUS=${RDK_RESERVED_CPUS:-11,27}
readonly ORDINARY_CPUS=${RDK_ORDINARY_CPUS:-0-10,12-26,28-31}
readonly LOCKED_KHZ=${RDK_LOCKED_KHZ:-4500000}
readonly MEMORY_COUNT=${RDK_RINGINFO_MEMORY_COUNT:-100000}
readonly CATCH_SAMPLES=${RDK_CATCH_BENCHMARK_SAMPLES:-3}
readonly CATCH_WARMUP_MS=${RDK_CATCH_BENCHMARK_WARMUP_MS:-100}
readonly PERF_LOCK=${RDK_PERF_LOCK:-/run/lock/rdkit-ringinfo-perf.lock}
readonly SERIES_ROOT=${RDK_SERIES_ROOT:-/tmp/rdkit-ringinfo-ab}
readonly RESULTS_DIR=${RDK_RESULTS_DIR:-$SERIES_ROOT/results}
readonly BEFORE_SOURCE=${RDK_BEFORE_SOURCE:-$SERIES_ROOT/source-before}
readonly AFTER_SOURCE=${RDK_AFTER_SOURCE:-$SOURCE_ROOT}
readonly BEFORE_BUILD=${RDK_BEFORE_BUILD:-$SERIES_ROOT/build-before}
readonly AFTER_BUILD=${RDK_AFTER_BUILD:-$SOURCE_ROOT/build}
readonly -a RING_DATASETS=(rings_2 rings_3 rings_4 rings_5 rings_6)

die() {
  echo "error: $*" >&2
  exit 1
}

usage() {
  cat <<EOF
Usage: $0 [all|measure|analyze]

Defaults:
  before build:    $BEFORE_BUILD
  after build:     $AFTER_BUILD
  rounds:          $ROUNDS interleaved pairs
  CPU:             $PERF_CPU (reserve SMT siblings $RESERVED_CPUS)
  locked clock:    $LOCKED_KHZ kHz
  results:         $RESULTS_DIR

The default 'all' action validates existing before/after binaries, obtains sudo
credentials, reserves the profiling core, locks every cpufreq policy, measures
all workloads, restores the host, and generates CSV summaries and SVG box
plots. It never configures, builds, or modifies either source tree.
EOF
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "$1 is required"
}

revision_at() {
  local source=$1
  git -c safe.directory="$source" -C "$source" rev-parse HEAD
}

require_prepared() {
  [[ -r $INPUT ]] || die "SMILES input is not readable: $INPUT"
  [[ $ROUNDS =~ ^[0-9]+$ && $ROUNDS -ge 8 && $((ROUNDS % 2)) -eq 0 ]] ||
    die "RDK_BENCHMARK_ROUNDS must be an even number of at least 8"
  local path
  for path in \
    "$BEFORE_BUILD/Code/Bench/bench" \
    "$BEFORE_BUILD/Code/Bench/ringinfo_memory_bench" \
    "$BEFORE_BUILD/Code/Bench/smiles_pipeline_bench" \
    "$AFTER_BUILD/Code/Bench/bench" \
    "$AFTER_BUILD/Code/Bench/ringinfo_memory_bench" \
    "$AFTER_BUILD/Code/Bench/smiles_pipeline_bench"; do
    [[ -x $path ]] || die "benchmark is not prepared: $path"
  done
}

run_scoped() {
  systemd-run --scope --quiet --slice=perf.slice --uid="$RUN_USER" \
    -p AllowedCPUs="$PERF_CPU" -p CPUWeight=10000 -- \
    "$@"
}

perf_value() {
  local file=$1 event=$2
  awk -F, -v wanted="$event:u" \
    '$3 == wanted {gsub(/[[:space:]]/, "", $1); print $1; found=1} END {exit !found}' \
    "$file"
}

run_micro() {
  local label=$1 round=$2 position=$3 source=$4 build=$5
  local stem
  stem=$(printf 'micro-round-%02d-position-%02d-%s' "$round" "$position" "$label")
  run_scoped env RDBASE="$source" \
    "$build/Code/Bench/bench" '[ringinfo]' \
    --benchmark-samples "$CATCH_SAMPLES" \
    --benchmark-warmup-time "$CATCH_WARMUP_MS" \
    --reporter XML >"$RESULTS_DIR/raw/$stem.xml"
}

run_memory() {
  local label=$1 round=$2 position=$3 source=$4 build=$5 dataset=$6
  local stem
  stem=$(printf 'memory-round-%02d-position-%02d-%s-%s' \
    "$round" "$position" "$label" "$dataset")
  run_scoped env RDBASE="$source" \
    "$build/Code/Bench/ringinfo_memory_bench" \
    "$source/Code/Bench/data/$dataset.smi" "$MEMORY_COUNT" \
    >"$RESULTS_DIR/raw/$stem.stdout"
}

run_e2e() {
  local label=$1 round=$2 position=$3 build=$4
  local stem stdout_file perf_file elapsed cycles instructions branches misses
  stem=$(printf 'e2e-round-%02d-position-%02d-%s' "$round" "$position" "$label")
  stdout_file=$RESULTS_DIR/raw/$stem.stdout
  perf_file=$RESULTS_DIR/raw/$stem.perf.csv
  run_scoped perf stat -x, -o "$perf_file" \
    -e cycles:u,instructions:u,branches:u,branch-misses:u -- \
    "$build/Code/Bench/smiles_pipeline_bench" "$INPUT" "$COUNT" \
    >"$stdout_file"

  grep -q "read=$COUNT" "$stdout_file" || die "wrong record count in $stdout_file"
  grep -q "sanitized=$COUNT" "$stdout_file" ||
    die "not all molecules sanitized in $stdout_file"
  grep -q 'failures=0' "$stdout_file" || die "failures reported in $stdout_file"
  elapsed=$(sed -n 's/.* elapsed_s=\([^ ]*\).*/\1/p' "$stdout_file")
  cycles=$(perf_value "$perf_file" cycles)
  instructions=$(perf_value "$perf_file" instructions)
  branches=$(perf_value "$perf_file" branches)
  misses=$(perf_value "$perf_file" branch-misses)
  printf '%s,%d,%d,%s,%s,%s,%s,%s\n' "$label" "$round" "$position" \
    "$elapsed" "$cycles" "$instructions" "$branches" "$misses" \
    >>"$RESULTS_DIR/smiles_135k.csv"
}

run_revision_workload() {
  local workload=$1 label=$2 round=$3 position=$4
  local source build
  if [[ $label == before ]]; then
    source=$BEFORE_SOURCE
    build=$BEFORE_BUILD
  else
    source=$AFTER_SOURCE
    build=$AFTER_BUILD
  fi
  case $workload in
    micro) run_micro "$label" "$round" "$position" "$source" "$build" ;;
    e2e) run_e2e "$label" "$round" "$position" "$build" ;;
    memory:*)
      run_memory "$label" "$round" "$position" "$source" "$build" \
        "${workload#memory:}"
      ;;
    *) die "unknown workload: $workload" ;;
  esac
}

measure_as_root() (
  ((EUID == 0)) || die "the internal measurement action must run as root"
  [[ -n ${SUDO_USER:-} || -n ${RDK_RUN_USER:-} ]] ||
    die "the original run user is unavailable"
  require_command flock
  require_command lscpu
  require_command perf
  require_command systemctl
  require_command systemd-run
  require_prepared
  [[ -d /sys/devices/system/cpu/cpu$PERF_CPU ]] ||
    die "measurement CPU $PERF_CPU does not exist"
  local siblings
  siblings=$(</sys/devices/system/cpu/cpu$PERF_CPU/topology/thread_siblings_list)
  [[ $siblings == "$RESERVED_CPUS" ]] ||
    die "CPU $PERF_CPU has SMT siblings $siblings, expected $RESERVED_CPUS"

  exec 9>"$PERF_LOCK"
  flock -x -w 600 9 || die "timed out waiting for exclusive measurement lock"

  readonly ORIGINAL_PARANOID=$(</proc/sys/kernel/perf_event_paranoid)
  readonly -a POLICIES=(/sys/devices/system/cpu/cpufreq/policy*)
  [[ -d ${POLICIES[0]} ]] || die "no cpufreq policies found"
  declare -a original_min=() original_max=()
  readonly -a CONTROLLED_UNITS=(system.slice user.slice init.scope)
  declare -A original_allowed_cpus=()
  local policy unit
  for policy in "${POLICIES[@]}"; do
    original_min+=("$(<"$policy/scaling_min_freq")")
    original_max+=("$(<"$policy/scaling_max_freq")")
    ((LOCKED_KHZ >= original_min[${#original_min[@]} - 1] &&
      LOCKED_KHZ <= original_max[${#original_max[@]} - 1])) ||
      die "requested clock $LOCKED_KHZ is outside the original bounds at $policy"
  done
  for unit in "${CONTROLLED_UNITS[@]}"; do
    original_allowed_cpus["$unit"]=$(
      systemctl show "$unit" -p AllowedCPUs --value
    )
  done

  restore_host() {
    local status=$?
    trap - EXIT INT TERM
    set +e
    local cleanup_failed=0
    local actual_allowed actual_max actual_min i saved_allowed unit
    for unit in "${CONTROLLED_UNITS[@]}"; do
      saved_allowed=${original_allowed_cpus["$unit"]-}
      if [[ -n $saved_allowed ]]; then
        systemctl set-property --runtime "$unit" \
          AllowedCPUs="$saved_allowed" || cleanup_failed=1
      else
        systemctl revert "$unit" || cleanup_failed=1
      fi
      actual_allowed=$(systemctl show "$unit" -p AllowedCPUs --value) ||
        cleanup_failed=1
      if [[ $actual_allowed != "$saved_allowed" ]]; then
        echo "error: failed to restore AllowedCPUs for $unit" >&2
        cleanup_failed=1
      fi
    done
    for i in "${!POLICIES[@]}"; do
      printf '%s\n' "${original_min[$i]}" \
        >"${POLICIES[$i]}/scaling_min_freq" || cleanup_failed=1
      printf '%s\n' "${original_max[$i]}" \
        >"${POLICIES[$i]}/scaling_max_freq" || cleanup_failed=1
      actual_min=$(<"${POLICIES[$i]}/scaling_min_freq")
      actual_max=$(<"${POLICIES[$i]}/scaling_max_freq")
      if [[ $actual_min != "${original_min[$i]}" ||
        $actual_max != "${original_max[$i]}" ]]; then
        echo "error: failed to restore clocks at ${POLICIES[$i]}" >&2
        cleanup_failed=1
      fi
    done
    sysctl -q -w "kernel.perf_event_paranoid=$ORIGINAL_PARANOID" ||
      cleanup_failed=1
    if [[ $(</proc/sys/kernel/perf_event_paranoid) != "$ORIGINAL_PARANOID" ]]; then
      echo "error: failed to restore kernel.perf_event_paranoid" >&2
      cleanup_failed=1
    fi
    chown -R "$RUN_USER" "$RESULTS_DIR" || true
    if ((cleanup_failed)); then
      echo "error: host cleanup was incomplete" >&2
      ((status == 0)) && status=1
    else
      echo "Host CPU scheduling, clocks, and perf permissions restored and verified."
    fi
    return "$status"
  }
  trap restore_host EXIT
  trap 'exit 130' INT
  trap 'exit 143' TERM

  mkdir -p "$RESULTS_DIR/raw"
  chown -R "$RUN_USER" "$RESULTS_DIR"
  cat >"$RESULTS_DIR/manifest.txt" <<EOF
before_source=$BEFORE_SOURCE
before_revision=$(revision_at "$BEFORE_SOURCE")
before_build=$BEFORE_BUILD
after_source=$AFTER_SOURCE
after_revision=$(revision_at "$AFTER_SOURCE")
after_build=$AFTER_BUILD
input=$INPUT
count=$COUNT
rounds=$ROUNDS
measurement_cpu=$PERF_CPU
reserved_cpus=$RESERVED_CPUS
locked_khz=$LOCKED_KHZ
catch_samples=$CATCH_SAMPLES
catch_warmup_ms=$CATCH_WARMUP_MS
memory_count=$MEMORY_COUNT
EOF
  printf '%s\n' \
    'revision,round,position,elapsed_s,cycles,instructions,branches,branch_misses' \
    >"$RESULTS_DIR/smiles_135k.csv"

  sysctl -q -w kernel.perf_event_paranoid=0
  for policy in "${POLICIES[@]}"; do
    # This order matches the host's proven lockClocks implementation.
    printf '%s\n' "$LOCKED_KHZ" >"$policy/scaling_min_freq"
    printf '%s\n' "$LOCKED_KHZ" >"$policy/scaling_max_freq"
  done
  local effective_locked_khz="" actual_min actual_max
  for policy in "${POLICIES[@]}"; do
    actual_min=$(<"$policy/scaling_min_freq")
    actual_max=$(<"$policy/scaling_max_freq")
    [[ $actual_min == "$actual_max" ]] ||
      die "clock is not fixed at $policy (min=$actual_min max=$actual_max)"
    if [[ -z $effective_locked_khz ]]; then
      effective_locked_khz=$actual_min
    elif [[ $actual_min != "$effective_locked_khz" ]]; then
      die "inconsistent fixed clock at $policy ($actual_min, expected $effective_locked_khz)"
    fi
  done
  echo "Requested $LOCKED_KHZ kHz; driver fixed all policies at $effective_locked_khz kHz."
  printf 'effective_locked_khz=%s\n' "$effective_locked_khz" \
    >>"$RESULTS_DIR/manifest.txt"

  {
    echo "timestamp=$(date --iso-8601=seconds)"
    echo "kernel=$(uname -srmo)"
    echo "model=$(lscpu | sed -n 's/^Model name:[[:space:]]*//p')"
    echo "measurement_cpu=$PERF_CPU"
    echo "thread_siblings=$siblings"
    echo "ordinary_cpus=$ORDINARY_CPUS"
    echo "perf_event_paranoid=$(</proc/sys/kernel/perf_event_paranoid)"
    for policy in "${POLICIES[@]}"; do
      echo "$(basename "$policy")_min=$(<"$policy/scaling_min_freq")"
      echo "$(basename "$policy")_max=$(<"$policy/scaling_max_freq")"
    done
  } >"$RESULTS_DIR/measurement_host.txt"

  systemctl set-property --runtime system.slice AllowedCPUs="$ORDINARY_CPUS"
  systemctl set-property --runtime user.slice AllowedCPUs="$ORDINARY_CPUS"
  systemctl set-property --runtime init.scope AllowedCPUs="$ORDINARY_CPUS"

  echo "Warming both revisions on CPU $PERF_CPU"
  run_scoped env RDBASE="$BEFORE_SOURCE" \
    "$BEFORE_BUILD/Code/Bench/bench" '[ringinfo]' \
    --benchmark-samples 1 --benchmark-warmup-time 0 >/dev/null
  run_scoped env RDBASE="$AFTER_SOURCE" \
    "$AFTER_BUILD/Code/Bench/bench" '[ringinfo]' \
    --benchmark-samples 1 --benchmark-warmup-time 0 >/dev/null
  run_scoped "$BEFORE_BUILD/Code/Bench/smiles_pipeline_bench" "$INPUT" 1000 \
    >/dev/null
  run_scoped "$AFTER_BUILD/Code/Bench/smiles_pipeline_bench" "$INPUT" 1000 \
    >/dev/null

  local round position first second workload label
  local -a workloads=(micro)
  for label in "${RING_DATASETS[@]}"; do
    workloads+=("memory:$label")
  done
  workloads+=(e2e)

  for ((round = 1; round <= ROUNDS; ++round)); do
    if ((round % 2)); then
      first=before
      second=after
    else
      first=after
      second=before
    fi
    echo "Measurement round $round/$ROUNDS ($first then $second)"
    for workload in "${workloads[@]}"; do
      position=1
      run_revision_workload "$workload" "$first" "$round" "$position"
      position=2
      run_revision_workload "$workload" "$second" "$round" "$position"
    done
  done
)

analyze() {
  require_command python3
  [[ -d $RESULTS_DIR/raw ]] || die "raw results are missing: $RESULTS_DIR/raw"
  mkdir -p "$RESULTS_DIR/plots" "$RESULTS_DIR/mplconfig"
  MPLCONFIGDIR=$RESULTS_DIR/mplconfig python3 - "$RESULTS_DIR" <<'PY'
import csv
import glob
import math
import os
import re
import statistics
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
except ImportError as exc:
    raise SystemExit(f"matplotlib is required to generate box plots: {exc}")

results = os.path.abspath(sys.argv[1])
raw = os.path.join(results, "raw")
plots = os.path.join(results, "plots")


def revision_from_name(path):
    name = os.path.basename(path)
    match = re.search(r"-(before|after)(?:-|\.)", name)
    if not match:
        raise ValueError(f"revision missing from filename: {name}")
    return match.group(1)


def round_position(path):
    name = os.path.basename(path)
    match = re.search(r"round-(\d+)-position-(\d+)", name)
    if not match:
        raise ValueError(f"round/position missing from filename: {name}")
    return int(match.group(1)), int(match.group(2))


micro_rows = []
for path in sorted(glob.glob(os.path.join(raw, "micro-*.xml"))):
    revision = revision_from_name(path)
    round_number, position = round_position(path)
    root = ET.parse(path).getroot()
    for result in root.iter("BenchmarkResults"):
        name = result.attrib["name"]
        mean = result.find("mean")
        if mean is None:
            raise ValueError(f"mean missing for {name} in {path}")
        dataset_match = re.search(r"(rings [2-6])$", name)
        if not dataset_match:
            raise ValueError(f"ring dataset missing from benchmark name: {name}")
        dataset = dataset_match.group(1).replace(" ", "_")
        operation = name[: dataset_match.start()].strip()
        suite = "perception" if operation.startswith("MolOps::") else "raw"
        micro_rows.append({
            "revision": revision,
            "round": round_number,
            "position": position,
            "suite": suite,
            "dataset": dataset,
            "operation": operation,
            "mean_ns": float(mean.attrib["value"]),
        })

with open(os.path.join(results, "ringinfo_micro.csv"), "w", newline="") as out:
    fields = ["revision", "round", "position", "suite", "dataset", "operation", "mean_ns"]
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    writer.writerows(micro_rows)


def parse_key_values(line):
    values = {}
    for key, quoted, plain in re.findall(r'(\w+)=(?:"([^"]*)"|(\S+))', line):
        values[key] = quoted if quoted else plain
    return values


memory_rows = []
for path in sorted(glob.glob(os.path.join(raw, "memory-*.stdout"))):
    revision = revision_from_name(path)
    round_number, position = round_position(path)
    dataset_match = re.search(r"-(rings_[2-6])\.stdout$", path)
    if not dataset_match:
        raise ValueError(f"dataset missing from memory filename: {path}")
    with open(path) as inp:
        values = parse_key_values(inp.read().strip())
    row = {
        "revision": revision,
        "round": round_number,
        "position": position,
        "dataset": dataset_match.group(1),
    }
    for key in (
        "count", "sizeof_ringinfo", "object_bytes", "baseline_rss_bytes",
        "object_rss_bytes", "populated_rss_bytes", "object_rss_delta_bytes",
        "population_rss_delta_bytes", "retained_delta_bytes",
        "object_bytes_per_ringinfo", "ring_bytes_per_ringinfo",
        "total_bytes_per_ringinfo", "rings", "atom_members",
    ):
        if key not in values:
            raise ValueError(f"{key} missing from {path}")
        row[key] = float(values[key])
    memory_rows.append(row)

with open(os.path.join(results, "ringinfo_memory.csv"), "w", newline="") as out:
    fields = list(memory_rows[0])
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    writer.writerows(memory_rows)


with open(os.path.join(results, "smiles_135k.csv"), newline="") as inp:
    e2e_rows = list(csv.DictReader(inp))
for row in e2e_rows:
    for key in ("round", "position"):
        row[key] = int(row[key])
    for key in ("elapsed_s", "cycles", "instructions", "branches", "branch_misses"):
        row[key] = float(row[key])


summary = []


def add_summary(suite, benchmark, unit, before_values, after_values):
    before_median = statistics.median(before_values)
    after_median = statistics.median(after_values)
    summary.append({
        "suite": suite,
        "benchmark": benchmark,
        "unit": unit,
        "before_median": before_median,
        "after_median": after_median,
        "change_percent": 100.0 * (after_median / before_median - 1.0),
        "before_min": min(before_values),
        "before_max": max(before_values),
        "after_min": min(after_values),
        "after_max": max(after_values),
        "samples_per_revision": len(before_values),
    })


micro_groups = defaultdict(lambda: {"before": [], "after": []})
for row in micro_rows:
    key = (row["suite"], row["dataset"], row["operation"])
    micro_groups[key][row["revision"]].append(row["mean_ns"])
for (suite, dataset, operation), values in sorted(micro_groups.items()):
    add_summary(f"ringinfo_{suite}", f"{dataset}: {operation}", "ns",
                values["before"], values["after"])

memory_groups = defaultdict(lambda: {"before": [], "after": []})
memory_object_groups = defaultdict(lambda: {"before": [], "after": []})
memory_ring_groups = defaultdict(lambda: {"before": [], "after": []})
for row in memory_rows:
    memory_groups[row["dataset"]][row["revision"]].append(
        row["total_bytes_per_ringinfo"])
    memory_object_groups[row["dataset"]][row["revision"]].append(
        row["object_bytes_per_ringinfo"])
    memory_ring_groups[row["dataset"]][row["revision"]].append(
        row["ring_bytes_per_ringinfo"])
for dataset, values in sorted(memory_groups.items()):
    add_summary("ringinfo_memory", dataset, "bytes_per_ringinfo",
                values["before"], values["after"])
    add_summary("ringinfo_memory_object", dataset, "bytes_per_ringinfo",
                memory_object_groups[dataset]["before"],
                memory_object_groups[dataset]["after"])
    add_summary("ringinfo_memory_ring_additions", dataset,
                "bytes_per_ringinfo",
                memory_ring_groups[dataset]["before"],
                memory_ring_groups[dataset]["after"])

for metric, unit in (
    ("elapsed_s", "seconds"),
    ("cycles", "cycles"),
    ("instructions", "instructions"),
    ("branches", "branches"),
    ("branch_misses", "branch_misses"),
):
    values = {
        revision: [row[metric] for row in e2e_rows if row["revision"] == revision]
        for revision in ("before", "after")
    }
    add_summary("smiles_135k", metric, unit, values["before"], values["after"])

with open(os.path.join(results, "summary.csv"), "w", newline="") as out:
    fields = list(summary[0])
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    writer.writerows(summary)


COLORS = {"before": "#6c757d", "after": "#255a91"}
DATASETS = [f"rings_{index}" for index in range(2, 7)]


for stale in glob.glob(os.path.join(plots, "ringinfo_raw_rings_*.svg")):
    os.remove(stale)

raw_operations = [
    "RingInfo construct",
    "RingInfo incremental addRing construction",
    "RingInfo numAtomRings",
    "RingInfo numBondRings",
    "RingInfo atomMembers traversal",
    "RingInfo bondMembers traversal",
    "RingInfo atomRings traversal",
    "RingInfo bondRings traversal",
    "RingInfo ring-size queries",
    "RingInfo ring-size vector materialization",
    "RingInfo same-ring queries",
    "RingInfo fused initialization and traversal",
]


def raw_values(operation, dataset, revision):
    return micro_groups[("raw", dataset, operation)][revision]


fig, axes = plt.subplots(4, 3, figsize=(16, 15))
for ax, operation in zip(axes.flat, raw_operations):
    positions = list(range(1, len(DATASETS) + 1))
    for offset, revision in ((-0.18, "before"), (0.18, "after")):
        values = [
            [sample / 1000.0 for sample in raw_values(operation, dataset, revision)]
            for dataset in DATASETS
        ]
        boxes = ax.boxplot(
            values,
            positions=[position + offset for position in positions],
            widths=0.30,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.3},
        )
        for box in boxes["boxes"]:
            box.set(facecolor=COLORS[revision], alpha=0.72)
        for position, samples in zip(positions, values):
            ax.scatter([position + offset] * len(samples), samples, s=12,
                       color=COLORS[revision], alpha=0.65, zorder=3)
    ax.set_title(operation.replace("RingInfo ", "", 1), fontsize=10)
    ax.set_xticks(positions, [dataset.replace("_", " ") for dataset in DATASETS])
    ax.set_ylabel("Actual time (µs)")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.22)
fig.legend(handles=[Patch(color=COLORS[key], label=key.title(), alpha=0.72)
                    for key in ("before", "after")],
           loc="upper center", ncol=2)
fig.suptitle("RingInfo raw operations: actual linear-scale time", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.975))
fig.savefig(os.path.join(plots, "ringinfo_raw_actual_time.svg"))
plt.close(fig)

reductions = {}
for operation in raw_operations:
    reductions[operation] = []
    for dataset in DATASETS:
        before = statistics.median(raw_values(operation, dataset, "before"))
        after = statistics.median(raw_values(operation, dataset, "after"))
        reductions[operation].append(100.0 * (1.0 - after / before))
all_reductions = [value for values in reductions.values() for value in values]
lower = min(0.0, min(all_reductions))
upper = max(0.0, max(all_reductions))
margin = max(3.0, 0.08 * (upper - lower))

fig, axes = plt.subplots(4, 3, figsize=(16, 15), sharey=True)
for ax, operation in zip(axes.flat, raw_operations):
    values = reductions[operation]
    colors = ["#2b8a3e" if value >= 0 else "#c92a2a" for value in values]
    bars = ax.bar(range(1, len(DATASETS) + 1), values, color=colors, alpha=0.82)
    ax.axhline(0, color="#333333", linewidth=0.9)
    ax.bar_label(bars, labels=[f"{value:+.1f}%" for value in values],
                 padding=3, fontsize=8)
    ax.set_title(operation.replace("RingInfo ", "", 1), fontsize=10)
    ax.set_xticks(range(1, len(DATASETS) + 1),
                  [dataset.replace("_", " ") for dataset in DATASETS])
    ax.set_ylabel("Median time reduction (%)")
    ax.set_ylim(lower - margin, upper + margin)
    ax.grid(axis="y", alpha=0.22)
fig.suptitle("RingInfo raw operations: percent reduction (positive is faster)",
             y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.975))
fig.savefig(os.path.join(plots, "ringinfo_raw_percent_reduction.svg"))
plt.close(fig)

perception_operations = ["MolOps::findSSSR", "MolOps::symmetrizeSSSR"]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
for ax, operation in zip(axes, perception_operations):
    positions = list(range(1, len(DATASETS) + 1))
    for offset, revision in ((-0.18, "before"), (0.18, "after")):
        values = [
            [sample / 1e6 for sample in
             micro_groups[("perception", dataset, operation)][revision]]
            for dataset in DATASETS
        ]
        boxes = ax.boxplot(values,
                           positions=[position + offset for position in positions],
                           widths=0.30, patch_artist=True, showfliers=False,
                           medianprops={"color": "black", "linewidth": 1.3})
        for box in boxes["boxes"]:
            box.set(facecolor=COLORS[revision], alpha=0.72)
        for position, samples in zip(positions, values):
            ax.scatter([position + offset] * len(samples), samples, s=16,
                       color=COLORS[revision], alpha=0.65, zorder=3)
    ax.set_title(operation.replace("MolOps::", ""))
    ax.set_xticks(positions, [dataset.replace("_", " ") for dataset in DATASETS])
    ax.set_ylabel("Actual time (ms)")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.22)
fig.legend(handles=[Patch(color=COLORS[key], label=key.title(), alpha=0.72)
                    for key in ("before", "after")],
           loc="upper center", ncol=2)
fig.suptitle("Ring perception: actual linear-scale time", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(plots, "ring_perception.svg"))
plt.close(fig)

perception_reductions = {}
for operation in perception_operations:
    perception_reductions[operation] = []
    for dataset in DATASETS:
        values = micro_groups[("perception", dataset, operation)]
        before = statistics.median(values["before"])
        after = statistics.median(values["after"])
        perception_reductions[operation].append(
            100.0 * (1.0 - after / before))
all_perception_reductions = [
    value
    for values in perception_reductions.values()
    for value in values
]
lower = min(0.0, min(all_perception_reductions))
upper = max(0.0, max(all_perception_reductions))
margin = max(0.5, 0.12 * (upper - lower))

fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharey=True)
for ax, operation in zip(axes, perception_operations):
    values = perception_reductions[operation]
    colors = ["#2b8a3e" if value >= 0 else "#c92a2a" for value in values]
    bars = ax.bar(range(1, len(DATASETS) + 1), values,
                  color=colors, alpha=0.82)
    ax.axhline(0, color="#333333", linewidth=0.9)
    ax.bar_label(bars, labels=[f"{value:+.1f}%" for value in values],
                 padding=3, fontsize=9)
    ax.set_title(operation.replace("MolOps::", ""))
    ax.set_xticks(range(1, len(DATASETS) + 1),
                  [dataset.replace("_", " ") for dataset in DATASETS])
    ax.set_ylabel("Median time reduction (%)")
    ax.set_ylim(lower - margin, upper + margin)
    ax.grid(axis="y", alpha=0.22)
fig.suptitle("Ring perception: percent reduction (positive is faster)",
             y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(plots, "ring_perception_percent_reduction.svg"))
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 5.8))
datasets = sorted(memory_groups)
positions = list(range(1, len(datasets) + 1))
for offset, revision in ((-0.19, "before"), (0.19, "after")):
    object_values = [
        statistics.median(memory_object_groups[dataset][revision]) / 1024.0
        for dataset in datasets
    ]
    ring_values = [
        statistics.median(memory_ring_groups[dataset][revision]) / 1024.0
        for dataset in datasets
    ]
    bar_positions = [position + offset for position in positions]
    ax.bar(bar_positions, object_values, width=0.34,
           color=COLORS[revision], alpha=0.95,
           label=f"{revision.title()}: object/default storage")
    ax.bar(bar_positions, ring_values, width=0.34, bottom=object_values,
           color=COLORS[revision], alpha=0.48, hatch="//",
           label=f"{revision.title()}: ring additions")
    for position, object_value, ring_value in zip(
            bar_positions, object_values, ring_values):
        total = object_value + ring_value
        ax.text(position, total, f"{total:.2f}", ha="center", va="bottom",
                fontsize=8, rotation=90)
ax.set_xticks(list(positions), [dataset.replace("_", " ") for dataset in datasets])
ax.set_ylabel("KiB per RingInfo")
ax.set_title("Total RingInfo memory: object/default storage + ring additions")
ax.grid(axis="y", alpha=0.22)
ax.legend(ncol=2, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(plots, "ringinfo_memory.svg"))
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.4))
for ax, metric, label, divisor in (
    (axes[0], "elapsed_s", "Elapsed time (seconds)", 1.0),
    (axes[1], "cycles", "User cycles (billions)", 1e9),
):
    values = [
        [row[metric] / divisor for row in e2e_rows if row["revision"] == revision]
        for revision in ("before", "after")
    ]
    boxes = ax.boxplot(values, tick_labels=["Before", "After"], patch_artist=True,
                       showfliers=False,
                       medianprops={"color": "black", "linewidth": 1.5})
    for box, revision in zip(boxes["boxes"], ("before", "after")):
        box.set(facecolor=COLORS[revision], alpha=0.72)
    for position, samples, revision in zip((1, 2), values, ("before", "after")):
        ax.scatter([position] * len(samples), samples, s=22,
                   color=COLORS[revision], alpha=0.72, zorder=3)
    ax.set_ylabel(label)
    ax.grid(axis="y", alpha=0.22)
fig.suptitle("SMILES parse + sanitize: 135,000 molecules")
fig.tight_layout()
fig.savefig(os.path.join(plots, "smiles_135k.svg"))
plt.close(fig)

print(f"Wrote summaries to {results}")
print(f"Wrote box plots to {plots}")
PY
  echo "Before/after summary: $RESULTS_DIR/summary.csv"
  echo "Box plots: $RESULTS_DIR/plots"
}

invoke_root_measurement() {
  require_command sudo
  sudo -v
  sudo env \
    RDK_RUN_USER="$RUN_USER" \
    RDK_BEFORE_SOURCE="$BEFORE_SOURCE" \
    RDK_AFTER_SOURCE="$AFTER_SOURCE" \
    RDK_BEFORE_BUILD="$BEFORE_BUILD" \
    RDK_AFTER_BUILD="$AFTER_BUILD" \
    RDK_SMILES_INPUT="$INPUT" \
    RDK_SMILES_COUNT="$COUNT" \
    RDK_BENCHMARK_ROUNDS="$ROUNDS" \
    RDK_MEASUREMENT_CPU="$PERF_CPU" \
    RDK_RESERVED_CPUS="$RESERVED_CPUS" \
    RDK_ORDINARY_CPUS="$ORDINARY_CPUS" \
    RDK_LOCKED_KHZ="$LOCKED_KHZ" \
    RDK_RINGINFO_MEMORY_COUNT="$MEMORY_COUNT" \
    RDK_CATCH_BENCHMARK_SAMPLES="$CATCH_SAMPLES" \
    RDK_CATCH_BENCHMARK_WARMUP_MS="$CATCH_WARMUP_MS" \
    RDK_PERF_LOCK="$PERF_LOCK" \
    RDK_SERIES_ROOT="$SERIES_ROOT" \
    RDK_RESULTS_DIR="$RESULTS_DIR" \
    "$SCRIPT_PATH" _measure-root
}

main() {
  case ${1:-all} in
    all)
      require_prepared
      invoke_root_measurement
      analyze
      ;;
    measure)
      require_prepared
      invoke_root_measurement
      ;;
    analyze) analyze ;;
    _measure-root) measure_as_root ;;
    -h|--help|help) usage ;;
    *) usage; die "unknown action: $1" ;;
  esac
}

main "${1:-all}"
