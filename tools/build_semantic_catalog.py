from __future__ import annotations

import csv
import hashlib
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "10_other_users"
MANIFEST = ARCHIVE / "manifest.tsv"
OUTPUT = ROOT / "11_semantic_catalog"
ANALYSIS_VIEW = OUTPUT / "by_analysis"
CELLCHAT_VIEW = OUTPUT / "03_cellchat"


def read_text(path: Path) -> str:
    data = path.read_bytes()
    for encoding in ("utf-8-sig", "gb18030", "utf-16", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            pass
    return data.decode("utf-8", errors="replace")


def code_only(text: str) -> str:
    return re.sub(r"(?m)^\s*#.*$", "", text)


def language_of(path: Path) -> str:
    return {
        ".r": "R",
        ".rmd": "R Markdown",
        ".py": "Python",
        ".sh": "Shell",
    }.get(path.suffix.lower(), path.suffix.lstrip(".") or "unknown")


def has_any(text: str, terms: tuple[str, ...]) -> bool:
    return any(term in text for term in terms)


def classify_role(name: str, suffix: str, code: str) -> str:
    key = f"{name}\n{code}".casefold()
    if suffix == ".rmd":
        return "report_template"
    if suffix == ".sh" or re.search(r"(^|[_-])(run|runner|submit|job)([_-]|\.)", name.casefold()):
        return "runner"
    if has_any(key, ("compare", "comparison", "contrast", "diffinteraction", "differential")):
        return "comparison"
    if has_any(key, ("plot", "visual", "heatmap", "dotplot", "violin", "bubble", "circlize")):
        return "visualization"
    if has_any(key, ("convert", "prepare", "preprocess", "split", "format", "extract", "reconstruct")):
        return "preprocessing"
    if re.search(r"(?m)^\s*[\w.]+\s*(?:<-|=)\s*function\s*\(", code):
        return "functions"
    return "analysis_pipeline"


def classify_general(category: str, name: str, code: str) -> str:
    key = f"{name}\n{code}".casefold()
    if category == "03_cellchat":
        return classify_cellchat(name, "", code)[0]
    if category == "04_cellphonedb":
        if has_any(key, ("ktplot", "heatmap", "dotplot", "plot_cpdb", "plot.")):
            return "cellphonedb_visualization"
        if has_any(key, ("compare", "contrast", "deg", "differential")):
            return "cellphonedb_comparison"
        if has_any(key, ("meta", "counts", "prepare", "convert")):
            return "cellphonedb_input_preparation"
        return "cellphonedb_analysis"
    if category == "05_trajectory_monocle":
        if has_any(key, ("monocle3", "learn_graph", "order_cells", "cell_data_set")):
            return "monocle3"
        if has_any(key, ("beam(", "branchedheatmap", "branch")):
            return "monocle2_branch_analysis"
        return "monocle2"
    if category == "06_trajectory_other":
        for term, label in (
            ("paga", "paga"),
            ("slingshot", "slingshot"),
            ("scvelo", "rna_velocity"),
            ("velocyto", "rna_velocity"),
            ("cytotrace", "cytotrace"),
            ("palantir", "palantir"),
        ):
            if term in key:
                return label
        return "other_trajectory"
    if category == "07_infercnv_copykat":
        if "copykat" in key:
            return "copykat"
        if "infercnv" in key:
            return "infercnv"
        return "cnv_supporting"
    if category == "08_scenic_regulon":
        if has_any(key, ("pyscenic", "ctx ", "aucell ", "grn ")):
            return "pyscenic"
        if has_any(key, ("rss", "regulonspecificityscore", "auc")) and has_any(key, ("plot", "heatmap")):
            return "regulon_visualization"
        return "scenic_r"
    if category == "09_qc_doublet":
        if "doubletfinder" in key:
            return "doubletfinder"
        if "scrublet" in key:
            return "scrublet"
        if "soupx" in key:
            return "soupx"
        return "qc"
    if category in {"10_scanpy", "11_clustering_annotation", "15_core_scrna"}:
        if has_any(key, ("doublet", "scrublet", "percent.mt", "pct_counts_mt", "filter_cells")):
            return "qc_doublet"
        if has_any(key, ("subcluster", "recluster", "subset(")):
            return "subclustering"
        if has_any(key, ("singler", "celltypist", "azimuth", "annotation", "annotate")):
            return "celltype_annotation"
        if has_any(key, ("findallmarkers", "findmarkers", "rank_genes_groups", "deg", "differential")):
            return "markers_deg"
        if has_any(key, ("harmony", "integratedata", "scanorama", "bbknn", "combat", "scvi")):
            return "integration_batch_correction"
        if has_any(key, ("findclusters", "leiden", "louvain", "neighbors", "runumap", "sc.tl.umap")):
            return "clustering_embedding"
        if has_any(key, ("convert", "h5ad", "loom", "sceasy", "seuratdisk")):
            return "format_conversion"
        return "scrna_core"
    if category == "12_scores_enrichment":
        if has_any(key, ("gsea", "fgsea")):
            return "gsea"
        if has_any(key, ("gsva", "ssgsea")):
            return "gsva_ssgsea"
        if has_any(key, ("addmodulescore", "score_genes", "module_score", "geneset_score")):
            return "module_score"
        if has_any(key, ("enrichgo", "enrichkegg", "clusterprofiler", "enrichment")):
            return "go_kegg_enrichment"
        return "scores_other"
    if category == "13_spatial_xenium":
        if "rctd" in key:
            return "rctd_deconvolution"
        if has_any(key, ("cellchat", "cell_chat")):
            return "spatial_cellchat"
        if has_any(key, ("segmentation", "nucleus", "boundary")):
            return "segmentation_qc"
        if has_any(key, ("neighbor", "niche", "microenvironment", "milo")):
            return "spatial_neighborhood"
        if has_any(key, ("annotation", "celltypist", "singler", "cell_type")):
            return "spatial_annotation"
        if has_any(key, ("cluster", "leiden", "umap")):
            return "spatial_clustering"
        return "spatial_preprocessing"
    return "other"


def comparison_group_count(code: str) -> int:
    patterns = (
        r"(?:group|groups|comparison|compare\.list)\s*=\s*c\(([^)]*)\)",
        r"(?:group|groups|comparison|compare\.list)\s*<-\s*c\(([^)]*)\)",
        r"object\.list\s*(?:<-|=)\s*list\(([^)]*)\)",
    )
    largest = 0
    for pattern in patterns:
        for match in re.finditer(pattern, code, flags=re.I | re.S):
            largest = max(largest, match.group(1).count(",") + 1)
    return largest


def classify_cellchat(name: str, source_path: str, text: str) -> tuple[str, str, str, dict[str, bool]]:
    suffix = Path(name).suffix.casefold()
    code = code_only(text)
    key = f"{name}\n{source_path}\n{code}".casefold()
    has_create = bool(re.search(r"\bcreatecellchat\s*\(", code, flags=re.I))
    has_compute = bool(re.search(r"\bcomputecommunprob(?:pathway)?\s*\(", code, flags=re.I))
    has_compare = bool(
        re.search(
            r"\b(?:mergecellchat|compareinteractions|netvisual_diffinteraction)\s*\(|"
            r"ranknet\s*\([^)]*mode\s*=\s*['\"]comparison",
            code,
            flags=re.I | re.S,
        )
    )
    has_plot = has_any(
        key,
        ("netvisual_", "netanalysis_", "ranknet(", "ggplot(", "heatmap", "dotplot", "bubble", "circlize"),
    )
    is_runner = suffix == ".sh" or (
        suffix == ".py" and has_any(key, ("subprocess", "rscript", "os.system", "run_cellchat"))
    )
    is_report = suffix == ".rmd"
    group_count = comparison_group_count(code)
    is_multi = has_any(
        key, ("多组", "multi_group", "muti_group", "multigroup", "three_group", "3group", "3groups")
    ) or group_count >= 3
    is_two = has_any(key, ("two_group", "2group", "pairwise")) or group_count == 2
    is_spatial = has_any(
        key,
        ("xenium", "visium", "spatial", "rctd", "computeregiondistance", "spatial.factors", "scale.factors"),
    )

    if is_report:
        scenario = "07_report_template"
    elif is_runner:
        scenario = "06_batch_runner"
    elif has_compare and has_create:
        if is_multi:
            scenario = "03_multi_group_comparison"
        elif is_two:
            scenario = "02_two_group_comparison"
        else:
            scenario = "09_group_comparison_unspecified"
    elif has_compare:
        scenario = "04_comparison_postprocessing"
    elif has_create or has_compute:
        scenario = "01_single_dataset"
    elif has_plot:
        scenario = "05_visualization"
    elif has_any(key, ("cellchat", "cell_chat", "convert", "split", "prepare", "print_data")):
        scenario = "08_preprocessing_support"
    else:
        scenario = "99_unclear"

    if has_compare:
        if is_multi:
            comparison_scope = "three_or_more_groups"
        elif is_two:
            comparison_scope = "two_groups"
        else:
            comparison_scope = "group_count_unspecified"
    elif has_create or has_compute:
        comparison_scope = "single_dataset"
    else:
        comparison_scope = "not_applicable"

    modality = "spatial_source" if is_spatial else "scrna_or_unspecified"
    flags = {
        "has_create_cellchat": has_create,
        "has_compute_communication": has_compute,
        "has_comparison": has_compare,
        "has_visualization": has_plot,
    }
    return scenario, comparison_scope, modality, flags


def safe_stem(name: str, limit: int = 72) -> str:
    stem = Path(name).stem
    stem = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", stem)
    stem = re.sub(r"\s+", "_", stem)
    stem = re.sub(r"_+", "_", stem).strip(" ._")
    return (stem or "script")[:limit].rstrip(" ._")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

    copied = [row for row in manifest_rows if row["status"] == "copied"]
    source_counts = Counter(row["sha256"] for row in manifest_rows if row.get("target_file"))
    name_counts = Counter(row["basename"].casefold() for row in copied)

    OUTPUT.mkdir(exist_ok=True)
    if ANALYSIS_VIEW.exists():
        shutil.rmtree(ANALYSIS_VIEW)
    ANALYSIS_VIEW.mkdir(parents=True)
    readme_path = CELLCHAT_VIEW / "README.md"
    cellchat_readme = readme_path.read_bytes() if readme_path.exists() else None
    if CELLCHAT_VIEW.exists():
        shutil.rmtree(CELLCHAT_VIEW)
    CELLCHAT_VIEW.mkdir(parents=True)
    if cellchat_readme is not None:
        readme_path.write_bytes(cellchat_readme)

    catalog_rows: list[dict[str, object]] = []
    cellchat_rows: list[dict[str, object]] = []
    duplicate_rows: list[dict[str, object]] = []
    checksum_errors: list[str] = []

    sorted_rows = sorted(copied, key=lambda row: (row["category"], row["basename"].casefold(), row["owner"], row["sha256"]))
    variant_numbers: defaultdict[str, int] = defaultdict(int)

    for row in sorted_rows:
        archived_path = ROOT / Path(row["target_file"])
        if not archived_path.exists():
            raise FileNotFoundError(archived_path)

        digest = hashlib.sha256(archived_path.read_bytes()).hexdigest()
        if digest != row["sha256"]:
            checksum_errors.append(row["target_file"])

        text = read_text(archived_path)
        code = code_only(text)
        suffix = archived_path.suffix.casefold()
        normalized_name = row["basename"].casefold()
        role = classify_role(row["basename"], suffix, code)
        scenario = classify_general(row["category"], row["basename"], code)
        modality = "spatial" if has_any(
            f"{row['source_path']}\n{text}".casefold(), ("xenium", "visium", "spatial", "rctd")
        ) else "scrna_or_unspecified"
        semantic_name = f"{scenario}__{row['owner']}__{safe_stem(row['basename'])}__{row['sha256'][:12]}{archived_path.suffix}"
        semantic_relative = Path("11_semantic_catalog") / "by_analysis" / row["category"] / scenario / semantic_name
        semantic_path = ROOT / semantic_relative
        semantic_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(archived_path, semantic_path)

        catalog = {
            "category": row["category"],
            "analysis_scenario": scenario,
            "workflow_role": role,
            "modality": modality,
            "owner": row["owner"],
            "original_name": row["basename"],
            "same_name_variant_count": name_counts[normalized_name],
            "source_count_same_content": source_counts[row["sha256"]],
            "semantic_name": semantic_name,
            "semantic_file": semantic_relative.as_posix(),
            "language": language_of(archived_path),
            "reuse_level": row["reuse_level"],
            "line_count": row["line_count"],
            "size_bytes": row["size_bytes"],
            "sha256": row["sha256"],
            "archived_file": row["target_file"],
            "source_path": row["source_path"],
        }
        catalog_rows.append(catalog)

        if name_counts[normalized_name] > 1:
            variant_numbers[normalized_name] += 1
            duplicate_rows.append(
                {
                    "normalized_name": normalized_name,
                    "original_name": row["basename"],
                    "variant_no": variant_numbers[normalized_name],
                    "variant_count": name_counts[normalized_name],
                    "category": row["category"],
                    "analysis_scenario": scenario,
                    "workflow_role": role,
                    "owner": row["owner"],
                    "semantic_name": semantic_name,
                    "sha256": row["sha256"],
                    "archived_file": row["target_file"],
                    "source_path": row["source_path"],
                }
            )

        if row["category"] == "03_cellchat":
            cell_scenario, comparison_scope, cell_modality, flags = classify_cellchat(
                row["basename"], row["source_path"], text
            )
            cell_name = (
                f"{cell_modality}__{row['owner']}__{safe_stem(row['basename'])}__"
                f"{row['sha256'][:12]}{archived_path.suffix}"
            )
            semantic_relative = Path("11_semantic_catalog") / "03_cellchat" / cell_scenario / cell_modality / cell_name
            semantic_path = ROOT / semantic_relative
            semantic_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(archived_path, semantic_path)
            cellchat_rows.append(
                {
                    "scenario": cell_scenario,
                    "comparison_scope": comparison_scope,
                    "modality": cell_modality,
                    "workflow_role": role,
                    "owner": row["owner"],
                    "original_name": row["basename"],
                    "same_name_variant_count": name_counts[normalized_name],
                    "source_count_same_content": source_counts[row["sha256"]],
                    "semantic_file": semantic_relative.as_posix(),
                    "archived_file": row["target_file"],
                    "source_path": row["source_path"],
                    "reuse_level": row["reuse_level"],
                    "line_count": row["line_count"],
                    "size_bytes": row["size_bytes"],
                    "sha256": row["sha256"],
                    **flags,
                }
            )

    if checksum_errors:
        raise RuntimeError(f"Checksum mismatch for {len(checksum_errors)} archived files")

    catalog_fields = [
        "category", "analysis_scenario", "workflow_role", "modality", "owner", "original_name",
        "same_name_variant_count", "source_count_same_content", "semantic_name", "semantic_file", "language", "reuse_level",
        "line_count", "size_bytes", "sha256", "archived_file", "source_path",
    ]
    duplicate_fields = [
        "normalized_name", "original_name", "variant_no", "variant_count", "category", "analysis_scenario",
        "workflow_role", "owner", "semantic_name", "sha256", "archived_file", "source_path",
    ]
    cellchat_fields = [
        "scenario", "comparison_scope", "modality", "workflow_role", "owner", "original_name",
        "same_name_variant_count", "source_count_same_content", "semantic_file", "archived_file", "source_path",
        "reuse_level", "line_count", "size_bytes", "sha256", "has_create_cellchat",
        "has_compute_communication", "has_comparison", "has_visualization",
    ]
    write_tsv(OUTPUT / "code_catalog.tsv", catalog_rows, catalog_fields)
    write_tsv(OUTPUT / "duplicate_original_names.tsv", duplicate_rows, duplicate_fields)
    write_tsv(CELLCHAT_VIEW / "catalog.tsv", cellchat_rows, cellchat_fields)

    summary_counter = Counter((row["scenario"], row["comparison_scope"], row["modality"]) for row in cellchat_rows)
    summary_rows = [
        {"scenario": key[0], "comparison_scope": key[1], "modality": key[2], "file_count": count}
        for key, count in sorted(summary_counter.items())
    ]
    write_tsv(CELLCHAT_VIEW / "summary.tsv", summary_rows, ["scenario", "comparison_scope", "modality", "file_count"])

    print(f"Cataloged {len(catalog_rows)} unique archived files")
    print(f"Built by-analysis semantic view with {len(catalog_rows)} byte-identical copies")
    print(f"Indexed {len(duplicate_rows)} files from repeated original-name groups")
    print(f"Built semantic CellChat view with {len(cellchat_rows)} files")
    print("Verified all archived SHA-256 values")


if __name__ == "__main__":
    main()
