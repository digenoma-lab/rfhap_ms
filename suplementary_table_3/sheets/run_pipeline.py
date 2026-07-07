#!/usr/bin/env python3
"""Rebuild Supplementary Table 3 directly from RFhap/yak/hifiasm logs.

The pipeline performs four reproducible steps:
1. Select one complete RFhap Nextflow execution per dataset.
2. Parse task-level CPU/RSS and standalone yak/hifiasm logs.
3. Calculate component, stage, and workflow resource requirements.
4. Export auditable TSV files and a formula-driven XLSX workbook that opens
   directly in Apple Numbers.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

try:
    import xlsxwriter
except ImportError as exc:  # pragma: no cover - user-facing dependency check
    raise SystemExit(
        "Missing dependency 'XlsxWriter'. Run: python -m pip install -r requirements.txt"
    ) from exc


PROCESS_MAP = {
    "create_kmers_database": ("Parental databases", "K-mer database construction"),
    "sort_kmer_db": ("Parental databases", "K-mer database sorting"),
    "print_paths": ("Classification + assembly", "Input preparation"),
    "FastKM": ("Classification + assembly", "FASTKM feature extraction"),
    "trainRF": ("Classification + assembly", "Random Forest training"),
    "predictRF": ("Classification + assembly", "Random Forest prediction"),
    "setHaplotypes": ("Classification + assembly", "Haplotype assignment"),
    "seqtk": ("Classification + assembly", "Read partitioning"),
}

REQUIRED_PROCESSES = {"FastKM", "predictRF", "seqtk"}
STAGES = ("Parental databases", "Classification + assembly")
FAMILIES = ("RFhap workflow", "Hifiasm-trio")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_name("config.json"),
        help="JSON configuration file (default: config.json beside this script)",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove previous generated files from the configured output directory",
    )
    return parser.parse_args()


def read_config(path: Path) -> dict[str, Any]:
    path = path.resolve()
    config = json.loads(path.read_text(encoding="utf-8"))
    base = path.parent
    for key in ("rfhap_pipeline_info", "rfhap_assemblies", "hifiasm_logs", "output_dir"):
        config[key] = (base / config[key]).resolve()
    config["config_path"] = path
    return config


def clean_process_name(name: str) -> str:
    return re.sub(r" \([0-9]+\)$", "", name or "")


def parse_filename_time(path: Path) -> datetime:
    hit = re.search(r"(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})", path.name)
    if not hit:
        return datetime.min
    return datetime.strptime(hit.group(1), "%Y-%m-%d_%H-%M-%S")


def matching_report(trace_path: Path) -> Path:
    name = trace_path.name.replace("execution_trace_", "execution_report_").replace(".txt", ".html")
    return trace_path.with_name(name)


def inspect_trace(trace_path: Path, dataset: str) -> dict[str, Any]:
    try:
        with trace_path.open(newline="", encoding="utf-8", errors="replace") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except Exception as exc:
        return {
            "dataset": dataset,
            "trace_file": str(trace_path),
            "report_file": str(matching_report(trace_path)),
            "n_rows": 0,
            "n_good": 0,
            "n_bad": 1,
            "n_required": 0,
            "complete_candidate": False,
            "timestamp": parse_filename_time(trace_path).isoformat(),
            "error": str(exc),
        }

    processes = {clean_process_name(row.get("name", "")) for row in rows}
    statuses = [row.get("status", "") for row in rows]
    n_good = sum(status in {"COMPLETED", "CACHED"} for status in statuses)
    n_bad = sum(status in {"FAILED", "ABORTED"} for status in statuses)
    report = matching_report(trace_path)
    return {
        "dataset": dataset,
        "trace_file": str(trace_path.resolve()),
        "report_file": str(report.resolve()),
        "n_rows": len(rows),
        "n_good": n_good,
        "n_bad": n_bad,
        "n_required": len(REQUIRED_PROCESSES.intersection(processes)),
        "complete_candidate": bool(
            rows
            and n_bad == 0
            and REQUIRED_PROCESSES.issubset(processes)
            and report.exists()
        ),
        "timestamp": parse_filename_time(trace_path).isoformat(),
        "error": "",
    }


def select_traces(root: Path, datasets: list[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    candidates: list[dict[str, Any]] = []
    selected: list[dict[str, Any]] = []
    for dataset in datasets:
        dataset_dir = root / f"pipeline_info_{dataset}"
        trace_files = sorted(dataset_dir.glob("execution_trace_*.txt"))
        if not trace_files:
            raise FileNotFoundError(f"No execution traces found for {dataset}: {dataset_dir}")
        dataset_candidates = [inspect_trace(path, dataset) for path in trace_files]
        candidates.extend(dataset_candidates)
        ranked = sorted(
            dataset_candidates,
            key=lambda row: (
                bool(row["complete_candidate"]),
                int(row["n_required"]),
                int(row["n_good"]),
                int(row["n_rows"]),
                row["timestamp"],
                row["trace_file"],
            ),
            reverse=True,
        )
        best = dict(ranked[0])
        if not best["complete_candidate"]:
            raise RuntimeError(f"No complete RFhap trace candidate found for {dataset}")
        best["selected"] = True
        selected.append(best)

    selected_paths = {row["trace_file"] for row in selected}
    for row in candidates:
        row["selected"] = row["trace_file"] in selected_paths
    candidates.sort(key=lambda row: (datasets.index(row["dataset"]), not row["selected"], row["timestamp"]))
    return candidates, selected


def json_like_field(obj: str, key: str) -> str | None:
    match = re.search(r'"' + re.escape(key) + r'":"([^"\\]*(?:\\.[^"\\]*)*)"', obj, re.S)
    return html.unescape(match.group(1)) if match else None


def parse_nextflow_report(report_path: Path, dataset: str) -> list[dict[str, Any]]:
    source = report_path.read_text(encoding="utf-8", errors="replace")
    start_marker = 'window.data = { "trace":['
    start = source.find(start_marker)
    summary = source.find('"summary"', start)
    if start < 0 or summary < 0:
        raise RuntimeError(f"Could not find embedded Nextflow trace data in {report_path}")
    block = source[start + len(start_marker):summary]
    block = re.sub(r"\]\s*,\s*$", "", block)
    objects = re.split(r'\},\s*\{(?="task_id")', block)

    tasks: list[dict[str, Any]] = []
    for raw in objects:
        obj = raw if raw.startswith("{") else "{" + raw
        obj = obj if obj.endswith("}") else obj + "}"
        process = json_like_field(obj, "process")
        if process not in PROCESS_MAP:
            continue
        cpus = json_like_field(obj, "cpus")
        peak_rss = json_like_field(obj, "peak_rss")
        tasks.append({
            "dataset": dataset,
            "process": process,
            "major_stage": PROCESS_MAP[process][0],
            "substage": PROCESS_MAP[process][1],
            "cpus": float(cpus) if cpus not in (None, "", "-") else None,
            "peak_rss_gb": float(peak_rss) / (1024 ** 3) if peak_rss not in (None, "", "-") else None,
            "start": float(json_like_field(obj, "start") or 0),
            "complete": float(json_like_field(obj, "complete") or 0),
            "source_file": str(report_path.resolve()),
        })
    if not tasks:
        raise RuntimeError(f"No recognized RFhap tasks parsed from {report_path}")
    return tasks


def storage_gb(value: float, unit: str) -> float:
    multiplier = {
        "B": 1e-9,
        "KB": 1e-6,
        "MB": 1e-3,
        "GB": 1.0,
        "TB": 1e3,
        "PB": 1e6,
    }
    return value * multiplier[unit.upper()]


def dataset_from_path(path: Path, datasets: list[str]) -> str:
    hits = [dataset for dataset in datasets if dataset in path.parts or dataset in str(path)]
    if len(hits) != 1:
        raise RuntimeError(f"Could not assign one dataset to source path: {path}")
    return hits[0]


def parse_external_log(path: Path, datasets: list[str], family: str) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    metrics = re.findall(
        r"\[M::main\] Real time:\s*([0-9.]+) sec; CPU:\s*([0-9.]+) sec; Peak RSS:\s*([0-9.]+)\s*([KMGTPE]?B)",
        text,
    )
    commands = re.findall(r"\[M::main\] CMD:\s*(.+)", text)
    versions = re.findall(r"\[M::main\] Version:\s*(.+)", text)
    dataset = dataset_from_path(path, datasets)
    if not metrics:
        return {
            "dataset": dataset,
            "source_file": str(path.resolve()),
            "parsed_ok": False,
            "workflow_family": family,
            "error": "No [M::main] resource line",
        }

    wall_seconds, cpu_seconds, rss_value, rss_unit = metrics[-1]
    command = commands[-1].strip() if commands else ""
    thread_match = re.search(r"(?:^|\s)-t\s*([0-9]+)", command)

    path_text = str(path)
    if "logs_maternal_yak" in path_text:
        stage, substage = "Parental databases", "Yak maternal database"
    elif "logs_paternal_yak" in path_text:
        stage, substage = "Parental databases", "Yak paternal database"
    elif "logs_hifiasm_trio" in path_text:
        stage, substage = "Classification + assembly", "Integrated trio-binning + assembly"
    elif "logs_hifiasm_hapA" in path_text:
        stage, substage = "Classification + assembly", "Haplotype A assembly"
    elif "logs_hifiasm_hapB" in path_text:
        stage, substage = "Classification + assembly", "Haplotype B assembly"
    else:
        raise RuntimeError(f"Unmapped external log role: {path}")

    return {
        "dataset": dataset,
        "major_stage": stage,
        "workflow_family": family,
        "workflow": "RFhap + hifiasm-ONT" if family == "RFhap workflow" else "Hifiasm-trio",
        "substage": substage,
        "cpu_threads": float(thread_match.group(1)) if thread_match else None,
        "peak_rss_gb": storage_gb(float(rss_value), rss_unit),
        "wall_time_h": float(wall_seconds) / 3600,
        "cpu_hours": float(cpu_seconds) / 3600,
        "version": versions[-1].strip() if versions else "",
        "command": command,
        "source_file": str(path.resolve()),
        "parsed_ok": True,
        "error": "",
    }


def discover_external_logs(config: dict[str, Any], datasets: list[str]) -> list[dict[str, Any]]:
    sources: list[tuple[Path, str]] = []
    sources.extend((path, "Hifiasm-trio") for path in Path(config["hifiasm_logs"]).rglob("*.err"))
    sources.extend((path, "RFhap workflow") for path in Path(config["rfhap_assemblies"]).rglob("*.err"))
    logs = [parse_external_log(path, datasets, family) for path, family in sorted(sources)]
    failed = [row for row in logs if not row.get("parsed_ok")]
    if failed:
        details = "\n".join(f"- {row['source_file']}: {row['error']}" for row in failed)
        raise RuntimeError(f"External log parsing failed:\n{details}")
    return logs


def max_present(values: Iterable[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return max(present) if present else None


def peak_concurrent_cpus(tasks: Iterable[dict[str, Any]]) -> float | None:
    events: list[tuple[float, int, float]] = []
    for task in tasks:
        if task["cpus"] is None or not task["start"] or not task["complete"]:
            continue
        events.append((task["start"], 1, task["cpus"]))
        events.append((task["complete"], -1, -task["cpus"]))
    events.sort(key=lambda event: (event[0], event[1]))
    current = 0.0
    peak = 0.0
    for _, _, change in events:
        current += change
        peak = max(peak, current)
    return peak if events else None


def component_workflow(stage: str, family: str) -> str:
    if family == "Hifiasm-trio":
        return "Hifiasm-trio"
    return "RFhap" if stage == "Parental databases" else "RFhap + hifiasm-ONT"


def construct_metrics(
    datasets: list[str],
    tasks: list[dict[str, Any]],
    external_logs: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    components: list[dict[str, Any]] = []
    task_keys = sorted({(task["dataset"], task["major_stage"], task["process"]) for task in tasks})
    for dataset, stage, process in task_keys:
        group = [task for task in tasks if (task["dataset"], task["major_stage"], task["process"]) == (dataset, stage, process)]
        components.append({
            "dataset": dataset,
            "major_stage": stage,
            "workflow_family": "RFhap workflow",
            "workflow": component_workflow(stage, "RFhap workflow"),
            "substage": PROCESS_MAP[process][1],
            "cpu_threads": max_present(task["cpus"] for task in group),
            "peak_rss_gb": max_present(task["peak_rss_gb"] for task in group),
            "cpu_definition": "Maximum allocated CPUs per task",
            "ram_definition": "Maximum task RSS",
            "row_type": "component",
            "source_file": group[0]["source_file"],
        })

    for row in external_logs:
        family = row["workflow_family"]
        components.append({
            "dataset": row["dataset"],
            "major_stage": row["major_stage"],
            "workflow_family": family,
            "workflow": component_workflow(row["major_stage"], family),
            "substage": row["substage"],
            "cpu_threads": row["cpu_threads"],
            "peak_rss_gb": row["peak_rss_gb"],
            "cpu_definition": "Requested command threads",
            "ram_definition": "Peak RSS",
            "row_type": "component",
            "source_file": row["source_file"],
        })

    stage_totals: list[dict[str, Any]] = []
    for dataset in datasets:
        for stage in STAGES:
            task_group = [task for task in tasks if task["dataset"] == dataset and task["major_stage"] == stage]
            external_rfhap = [
                row for row in components
                if row["dataset"] == dataset
                and row["major_stage"] == stage
                and row["workflow_family"] == "RFhap workflow"
                and row["substage"] in {"Haplotype A assembly", "Haplotype B assembly"}
            ]
            cpu_candidates = [peak_concurrent_cpus(task_group)] + [row["cpu_threads"] for row in external_rfhap]
            rss_candidates = [task["peak_rss_gb"] for task in task_group] + [row["peak_rss_gb"] for row in external_rfhap]
            stage_totals.append({
                "dataset": dataset,
                "major_stage": stage,
                "workflow_family": "RFhap workflow",
                "workflow": component_workflow(stage, "RFhap workflow"),
                "substage": "TOTAL - " + stage,
                "cpu_threads": max_present(cpu_candidates),
                "peak_rss_gb": max_present(rss_candidates),
                "cpu_definition": "Peak allocated CPUs (RFhap) or maximum requested threads",
                "ram_definition": "Maximum reported task/command RSS",
                "row_type": "stage_total",
                "source_file": "Multiple files",
            })

            hif_components = [
                row for row in components
                if row["dataset"] == dataset
                and row["major_stage"] == stage
                and row["workflow_family"] == "Hifiasm-trio"
            ]
            stage_totals.append({
                "dataset": dataset,
                "major_stage": stage,
                "workflow_family": "Hifiasm-trio",
                "workflow": "Hifiasm-trio",
                "substage": "TOTAL - " + stage,
                "cpu_threads": max_present(row["cpu_threads"] for row in hif_components),
                "peak_rss_gb": max_present(row["peak_rss_gb"] for row in hif_components),
                "cpu_definition": "Maximum requested threads per standalone command",
                "ram_definition": "Maximum command peak RSS",
                "row_type": "stage_total",
                "source_file": "Multiple files",
            })

    end_totals: list[dict[str, Any]] = []
    for dataset in datasets:
        for family in FAMILIES:
            group = [row for row in stage_totals if row["dataset"] == dataset and row["workflow_family"] == family]
            end_totals.append({
                "dataset": dataset,
                "major_stage": "End-to-end",
                "workflow_family": family,
                "workflow": "RFhap + hifiasm-ONT" if family == "RFhap workflow" else "Hifiasm-trio",
                "substage": "TOTAL - end-to-end",
                "cpu_threads": max_present(row["cpu_threads"] for row in group),
                "peak_rss_gb": max_present(row["peak_rss_gb"] for row in group),
                "cpu_definition": "Maximum stage CPU/thread requirement",
                "ram_definition": "Maximum stage RSS",
                "row_type": "end_to_end_total",
                "source_file": "Stage totals",
            })

    stage_order = {"Parental databases": 0, "Classification + assembly": 1, "End-to-end": 2}
    family_order = {"RFhap workflow": 0, "Hifiasm-trio": 1}
    all_rows = components + stage_totals + end_totals
    all_rows.sort(key=lambda row: (
        datasets.index(row["dataset"]),
        stage_order[row["major_stage"]],
        family_order[row["workflow_family"]],
        row["row_type"] != "component",
        row["substage"],
    ))
    return components, stage_totals, end_totals


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def metric_lookup(
    rows: list[dict[str, Any]],
    dataset: str,
    family: str,
    stage: str,
    row_type: str,
    substage: str | None = None,
) -> dict[str, Any]:
    matches = [
        row for row in rows
        if row["dataset"] == dataset
        and row["workflow_family"] == family
        and row["major_stage"] == stage
        and row["row_type"] == row_type
        and (substage is None or row["substage"] == substage)
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one metric row, found {len(matches)}: "
            f"{dataset} | {family} | {stage} | {row_type} | {substage}"
        )
    return matches[0]


def build_workbook(
    path: Path,
    datasets: list[str],
    rows: list[dict[str, Any]],
) -> dict[str, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook = xlsxwriter.Workbook(path)
    workbook.set_properties({
        "title": "Supplementary Table 3 — CPU and memory requirements",
        "author": "RFhap computational comparison",
        "comments": "Generated reproducibly from Nextflow, yak, and hifiasm logs.",
    })

    navy, text, grid = "#18324B", "#202A33", "#C7D2D9"
    rf_light, rf_mid, rf_dark = "#E5F4F1", "#BFE3DC", "#238E83"
    hi_light, hi_mid, hi_dark = "#FCEAE0", "#F5CCB8", "#D9793D"
    gray = "#F3F6F8"

    title_fmt = workbook.add_format({"bold": True, "font_color": "#FFFFFF", "bg_color": navy, "font_size": 16, "valign": "vcenter"})
    subtitle_fmt = workbook.add_format({"font_color": "#52606D", "font_size": 10, "italic": True, "valign": "vcenter"})
    header_fmt = workbook.add_format({"bold": True, "font_color": "#FFFFFF", "bg_color": navy, "font_size": 10, "align": "center", "valign": "vcenter", "text_wrap": True})
    section_fmt = workbook.add_format({"bold": True, "font_color": "#FFFFFF", "bg_color": navy, "font_size": 11, "valign": "vcenter"})
    note_fmt = workbook.add_format({"font_color": "#52606D", "bg_color": gray, "font_size": 9, "italic": True, "text_wrap": True, "valign": "top"})
    calc_header_fmt = workbook.add_format({"bold": True, "font_color": text, "bg_color": gray, "border": 1, "border_color": grid, "text_wrap": True})
    calc_text_fmt = workbook.add_format({"font_color": text, "border": 1, "border_color": grid, "text_wrap": True})
    calc_int_fmt = workbook.add_format({"font_color": text, "border": 1, "border_color": grid, "num_format": "#,##0", "align": "right"})

    def row_format(fill: str, bold: bool = False, white: bool = False, top: int = 0) -> tuple[Any, Any, Any]:
        base = {
            "bg_color": fill,
            "font_color": "#FFFFFF" if white else text,
            "bold": bold,
            "font_size": 10,
            "valign": "vcenter",
            "text_wrap": True,
            "bottom": 1,
            "bottom_color": grid,
        }
        if top:
            base.update({"top": top, "top_color": navy})
        text_fmt = workbook.add_format(base)
        integer_fmt = workbook.add_format({**base, "num_format": "#,##0", "align": "right"})
        decimal_fmt = workbook.add_format({**base, "num_format": "0.00", "align": "right"})
        return text_fmt, integer_fmt, decimal_fmt

    rf_component = row_format(rf_light)
    rf_subtotal = row_format(rf_mid, bold=True, top=2)
    rf_total = row_format(rf_dark, bold=True, white=True, top=2)
    hi_component = row_format(hi_light)
    hi_subtotal = row_format(hi_mid, bold=True, top=2)
    hi_total = row_format(hi_dark, bold=True, white=True, top=2)

    formula_counts: dict[str, int] = {}
    for dataset in datasets:
        sheet = workbook.add_worksheet(dataset)
        sheet.hide_gridlines(2)
        sheet.freeze_panes(4, 0)
        sheet.set_landscape()
        sheet.fit_to_pages(1, 2)
        sheet.set_margins(0.25, 0.25, 0.35, 0.35)
        sheet.set_column("A:A", 23)
        sheet.set_column("B:B", 21)
        sheet.set_column("C:C", 34)
        sheet.set_column("D:D", 14)
        sheet.set_column("E:E", 15)
        sheet.set_column("F:F", 14)
        sheet.set_column("G:G", 31)

        sheet.merge_range("A1:G1", f"Supplementary Table 3 — CPU and memory requirements: {dataset}", title_fmt)
        sheet.set_row(0, 26)
        sheet.merge_range("A2:G2", "Process breakdown with formula-driven stage subtotals and workflow totals. RFhap is teal; hifiasm-trio is orange.", subtitle_fmt)
        sheet.set_row(1, 20)
        headers = ["Stage", "Workflow", "Process / subtotal / total", "CPUs / threads", "Peak RSS (GB)", "Row type", "Calculation basis"]
        for col, value in enumerate(headers):
            sheet.write(3, col, value, header_fmt)
        sheet.set_row(3, 25)

        rf_parent_components = [
            metric_lookup(rows, dataset, "RFhap workflow", "Parental databases", "component", "K-mer database construction"),
            metric_lookup(rows, dataset, "RFhap workflow", "Parental databases", "component", "K-mer database sorting"),
        ]
        rf_parent_total_row = metric_lookup(rows, dataset, "RFhap workflow", "Parental databases", "stage_total")
        rf_class_names = [
            "FASTKM feature extraction", "Haplotype A assembly", "Haplotype B assembly",
            "Haplotype assignment", "Input preparation", "Random Forest prediction",
            "Random Forest training", "Read partitioning",
        ]
        rf_class_components = [metric_lookup(rows, dataset, "RFhap workflow", "Classification + assembly", "component", name) for name in rf_class_names]
        rf_class_total_row = metric_lookup(rows, dataset, "RFhap workflow", "Classification + assembly", "stage_total")
        hif_parent_components = [
            metric_lookup(rows, dataset, "Hifiasm-trio", "Parental databases", "component", "Yak maternal database"),
            metric_lookup(rows, dataset, "Hifiasm-trio", "Parental databases", "component", "Yak paternal database"),
        ]
        trio_component = metric_lookup(rows, dataset, "Hifiasm-trio", "Classification + assembly", "component", "Integrated trio-binning + assembly")

        formulas = 0

        def write_component(excel_row: int, stage: str, workflow: str, row: dict[str, Any], fmts: tuple[Any, Any, Any]) -> None:
            values = [stage, workflow, row["substage"], row["cpu_threads"], round(float(row["peak_rss_gb"]), 2), "Component", "Parsed process allocation / RSS" if row["workflow_family"] == "RFhap workflow" else "Requested -t threads / command peak RSS"]
            for col, value in enumerate(values):
                fmt = fmts[1] if col == 3 else fmts[2] if col == 4 else fmts[0]
                sheet.write(excel_row - 1, col, value, fmt)

        write_component(5, "Parental databases", "RFhap", rf_parent_components[0], rf_component)
        write_component(6, "", "", rf_parent_components[1], rf_component)
        subtotal_values = ["", "", "SUBTOTAL — Parental databases", None, None, "Stage subtotal", "CPU: parsed peak concurrency; RAM: MAX of components"]
        for col, value in enumerate(subtotal_values):
            if col not in {3, 4}:
                sheet.write(6, col, value, rf_subtotal[0])
        sheet.write_formula("D7", "=$C$30", rf_subtotal[1], rf_parent_total_row["cpu_threads"]); formulas += 1
        sheet.write_formula("E7", "=MAX(E5:E6)", rf_subtotal[2], round(float(rf_parent_total_row["peak_rss_gb"]), 2)); formulas += 1

        for index, row in enumerate(rf_class_components):
            write_component(8 + index, "RFhap classification + assembly" if index == 0 else "", "RFhap" if index == 0 else "", row, rf_component)
        subtotal_values = ["", "", "SUBTOTAL — RFhap classification + assembly", None, None, "Stage subtotal", "CPU: parsed peak concurrency; RAM: MAX of components"]
        for col, value in enumerate(subtotal_values):
            if col not in {3, 4}:
                sheet.write(15, col, value, rf_subtotal[0])
        sheet.write_formula("D16", "=$C$31", rf_subtotal[1], rf_class_total_row["cpu_threads"]); formulas += 1
        sheet.write_formula("E16", "=MAX(E8:E15)", rf_subtotal[2], round(float(rf_class_total_row["peak_rss_gb"]), 2)); formulas += 1
        total_values = ["Workflow total", "RFhap", "TOTAL — maximum resource requirement", None, None, "Workflow total", "MAX of the two RFhap stage subtotals"]
        for col, value in enumerate(total_values):
            if col not in {3, 4}:
                sheet.write(16, col, value, rf_total[0])
        rf_end = metric_lookup(rows, dataset, "RFhap workflow", "End-to-end", "end_to_end_total")
        sheet.write_formula("D17", "=MAX(D7,D16)", rf_total[1], rf_end["cpu_threads"]); formulas += 1
        sheet.write_formula("E17", "=MAX(E7,E16)", rf_total[2], round(float(rf_end["peak_rss_gb"]), 2)); formulas += 1

        write_component(19, "Parental databases", "Hifiasm-trio", hif_parent_components[0], hi_component)
        write_component(20, "", "", hif_parent_components[1], hi_component)
        subtotal_values = ["", "", "SUBTOTAL — Parental databases", None, None, "Stage subtotal", "MAX of maternal and paternal commands"]
        for col, value in enumerate(subtotal_values):
            if col not in {3, 4}:
                sheet.write(20, col, value, hi_subtotal[0])
        hif_parent_total = metric_lookup(rows, dataset, "Hifiasm-trio", "Parental databases", "stage_total")
        sheet.write_formula("D21", "=MAX(D19:D20)", hi_subtotal[1], hif_parent_total["cpu_threads"]); formulas += 1
        sheet.write_formula("E21", "=MAX(E19:E20)", hi_subtotal[2], round(float(hif_parent_total["peak_rss_gb"]), 2)); formulas += 1
        write_component(22, "Hifiasm trio", "Hifiasm-trio", trio_component, hi_component)
        subtotal_values = ["", "", "SUBTOTAL — Hifiasm trio", None, None, "Stage subtotal", "MAX of stage commands"]
        for col, value in enumerate(subtotal_values):
            if col not in {3, 4}:
                sheet.write(22, col, value, hi_subtotal[0])
        hif_trio_total = metric_lookup(rows, dataset, "Hifiasm-trio", "Classification + assembly", "stage_total")
        sheet.write_formula("D23", "=MAX(D22:D22)", hi_subtotal[1], hif_trio_total["cpu_threads"]); formulas += 1
        sheet.write_formula("E23", "=MAX(E22:E22)", hi_subtotal[2], round(float(hif_trio_total["peak_rss_gb"]), 2)); formulas += 1
        total_values = ["Workflow total", "Hifiasm-trio", "TOTAL — maximum resource requirement", None, None, "Workflow total", "MAX of the two hifiasm-trio stage subtotals"]
        for col, value in enumerate(total_values):
            if col not in {3, 4}:
                sheet.write(23, col, value, hi_total[0])
        hif_end = metric_lookup(rows, dataset, "Hifiasm-trio", "End-to-end", "end_to_end_total")
        sheet.write_formula("D24", "=MAX(D21,D23)", hi_total[1], hif_end["cpu_threads"]); formulas += 1
        sheet.write_formula("E24", "=MAX(E21,E23)", hi_total[2], round(float(hif_end["peak_rss_gb"]), 2)); formulas += 1

        sheet.merge_range("A27:G27", "Calculation inputs and interpretation", section_fmt)
        calc_headers = ["Workflow", "Stage", "Parsed peak concurrent CPUs"]
        for col, value in enumerate(calc_headers):
            sheet.write(28, col, value, calc_header_fmt)
        sheet.merge_range("D29:G29", "Definition", calc_header_fmt)
        sheet.write_row(29, 0, ["RFhap", "Parental databases"], calc_text_fmt)
        sheet.write(29, 2, rf_parent_total_row["cpu_threads"], calc_int_fmt)
        sheet.merge_range("D30:G30", rf_parent_total_row["cpu_definition"], calc_text_fmt)
        sheet.write_row(30, 0, ["RFhap", "Classification + assembly"], calc_text_fmt)
        sheet.write(30, 2, rf_class_total_row["cpu_threads"], calc_int_fmt)
        sheet.merge_range("D31:G31", rf_class_total_row["cpu_definition"], calc_text_fmt)
        sheet.set_row(28, 28)
        sheet.set_row(29, 28)
        sheet.set_row(30, 28)

        note = (
            "Notes. RFhap component CPUs are Nextflow task allocations; RFhap stage CPUs are parsed peak concurrent allocations and therefore are not sums of the displayed component maxima. "
            "Yak/hifiasm CPUs are requested -t threads. Peak RSS is the maximum resident memory and is never summed. Workflow totals are maximum capacity requirements across stages."
        )
        sheet.merge_range("A35:G37", note, note_fmt)
        sheet.set_row(34, 20)
        sheet.set_row(35, 20)
        sheet.set_row(36, 20)
        sheet.print_area("A1:G37")
        formula_counts[dataset] = formulas

    workbook.close()
    return formula_counts


def validate_metrics(
    datasets: list[str],
    components: list[dict[str, Any]],
    stage_totals: list[dict[str, Any]],
    end_totals: list[dict[str, Any]],
    formula_counts: dict[str, int],
) -> dict[str, Any]:
    expected_component_names = {
        "K-mer database construction", "K-mer database sorting",
        "FASTKM feature extraction", "Haplotype A assembly", "Haplotype B assembly",
        "Haplotype assignment", "Input preparation", "Random Forest prediction",
        "Random Forest training", "Read partitioning", "Yak maternal database",
        "Yak paternal database", "Integrated trio-binning + assembly",
    }
    qc: dict[str, Any] = {"datasets": {}, "passed": True}
    for dataset in datasets:
        dataset_components = [row for row in components if row["dataset"] == dataset]
        names = {row["substage"] for row in dataset_components}
        missing = sorted(expected_component_names - names)
        duplicate_counts = Counter(row["substage"] for row in dataset_components)
        duplicates = sorted(name for name, count in duplicate_counts.items() if count > 1)
        dataset_stage = [row for row in stage_totals if row["dataset"] == dataset]
        dataset_end = [row for row in end_totals if row["dataset"] == dataset]
        passed = not missing and not duplicates and len(dataset_stage) == 4 and len(dataset_end) == 2 and formula_counts.get(dataset) == 12
        qc["datasets"][dataset] = {
            "component_rows": len(dataset_components),
            "stage_total_rows": len(dataset_stage),
            "end_to_end_rows": len(dataset_end),
            "formula_cells": formula_counts.get(dataset, 0),
            "missing_components": missing,
            "duplicate_components": duplicates,
            "passed": passed,
        }
        qc["passed"] = qc["passed"] and passed
    if not qc["passed"]:
        raise RuntimeError("QC validation failed; inspect outputs/qc_summary.json")
    return qc


def clean_outputs(output_dir: Path) -> None:
    if not output_dir.exists():
        return
    for path in output_dir.iterdir():
        if path.is_file() and path.suffix.lower() in {".tsv", ".json", ".xlsx"}:
            path.unlink()


def main() -> int:
    args = parse_args()
    config = read_config(args.config)
    datasets = list(config["datasets"])
    output_dir = Path(config["output_dir"])
    if args.clean:
        clean_outputs(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for key in ("rfhap_pipeline_info", "rfhap_assemblies", "hifiasm_logs"):
        if not Path(config[key]).exists():
            raise FileNotFoundError(f"Configured source directory does not exist: {config[key]}")

    trace_qc, selected = select_traces(Path(config["rfhap_pipeline_info"]), datasets)
    tasks: list[dict[str, Any]] = []
    for row in selected:
        tasks.extend(parse_nextflow_report(Path(row["report_file"]), row["dataset"]))

    external_logs = discover_external_logs(config, datasets)
    components, stage_totals, end_totals = construct_metrics(datasets, tasks, external_logs)
    all_rows = components + stage_totals + end_totals

    fields = [
        "dataset", "major_stage", "workflow_family", "workflow", "substage",
        "cpu_threads", "peak_rss_gb", "cpu_definition", "ram_definition",
        "row_type", "source_file",
    ]
    write_tsv(output_dir / "trace_selection_qc.tsv", trace_qc, list(trace_qc[0].keys()))
    write_tsv(output_dir / "selected_traces.tsv", selected, list(selected[0].keys()))
    write_tsv(output_dir / "parsed_external_logs.tsv", external_logs, sorted({key for row in external_logs for key in row.keys()}))
    write_tsv(output_dir / "Supplementary_Table_3_CPU_RAM_components.tsv", all_rows, fields)
    write_tsv(output_dir / "Supplementary_Table_3_CPU_RAM_totals.tsv", stage_totals + end_totals, fields)

    workbook_path = output_dir / config["workbook_name"]
    formula_counts = build_workbook(workbook_path, datasets, all_rows)
    qc = validate_metrics(datasets, components, stage_totals, end_totals, formula_counts)
    qc.update({
        "config": str(config["config_path"]),
        "workbook": str(workbook_path),
        "selected_traces": {row["dataset"]: row["trace_file"] for row in selected},
        "external_logs_parsed": len(external_logs),
        "component_rows_total": len(components),
        "stage_total_rows_total": len(stage_totals),
        "end_to_end_rows_total": len(end_totals),
    })
    (output_dir / "qc_summary.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")

    print(f"Workbook: {workbook_path}")
    print(f"QC: {output_dir / 'qc_summary.json'}")
    print("Status: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
