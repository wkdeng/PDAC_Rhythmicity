#!/usr/bin/env python3
"""Standalone ChIP background sensitivity; never writes notebooks or legacy results."""
import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import false_discovery_control, fisher_exact, hypergeom


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, records):
    records = list(records)
    if not records:
        return
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def clean_id(value):
    return value.split(".")[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    data = root / "data"
    out = data / "chip_background_sensitivity"
    log = data / "logs" / "chip_background_sensitivity"
    out.mkdir(parents=True, exist_ok=True)
    log.mkdir(parents=True, exist_ok=True)
    inputs = []

    def read(relative):
        path = data / relative
        inputs.append(path)
        return rows(path)

    groups = ["GTEx_Normal", "TCGA-PAAD_Tumor"]
    tfs = ["ARNTL", "CLOCK", "CRY1", "DBP", "NPAS2", "NR1D1", "RORA", "RORC", "TEF"]
    tested, rhythmic = {}, {}
    for group in groups:
        all_rows = read(f"Merged/regression/all_genes_regression_periodic_12h_{group}.txt")
        tested[group] = {clean_id(r["ENSG"]) for r in all_rows
                         if math.isfinite(float(r["p"])) and 0 <= float(r["p"]) <= 1}
        assert len(tested[group]) == len(all_rows), "Invalid or duplicate tested genes"
        sig = read(f"Merged/regression/sig_genes_regression_periodic_12h_{group}.txt")
        rhythmic[group] = {clean_id(r["ENSG"]) for r in sig}
        assert rhythmic[group] <= tested[group]
    original_targets = {"GTEx_Normal": rhythmic["GTEx_Normal"]}
    legacy_identifiers = {"GTEx_Normal": rhythmic["GTEx_Normal"]}
    selected_sets = [rhythmic["GTEx_Normal"]]
    legacy_selected_sets = [rhythmic["GTEx_Normal"]]
    for group in ["TCGA-PAAD_Tumor", "CPTAC-3_Tumor"]:
        gene_rows = read(f"enrichment/cancer_specific_genes/method_A_presence/{group}_cancer_specific_genes.txt")
        id_key = next(k for k in ["ENSG", "ENSG_clean", "ensg", "gene_id"] if k in gene_rows[0])
        ids = {clean_id(r[id_key]) for r in gene_rows}
        selected_sets.append(ids)
        legacy_key = next(k for k in ["ENSG", "Gene", "gene_id", "ENSG_clean"] if k in gene_rows[0])
        legacy_ids = {clean_id(r[legacy_key]) for r in gene_rows}
        legacy_selected_sets.append(legacy_ids)
        if group in groups:
            original_targets[group] = ids
            legacy_identifiers[group] = legacy_ids
    selected_union = set.union(*selected_sets)
    shared_tested = set.intersection(*tested.values())

    # Preserve the exact cached annotation and promoter convention of the source
    # for the primary background-only comparison. Check BED conversion separately.
    cached = read("reference/gencode_gene_coordinates.tsv")
    coordinates = {clean_id(r["gene_id"]): r for r in cached}
    candidate_ids = selected_union | set.union(*tested.values())
    eligible_coordinates = {g: coordinates[g] for g in candidate_ids if g in coordinates
                            and coordinates[g]["chrom"].startswith("chr")
                            and "_" not in coordinates[g]["chrom"]
                            and coordinates[g]["strand"] in ["+", "-"]}
    gtf_path = data / "Ref/gencode.v49.basic.annotation.gtf"
    inputs.append(gtf_path)
    gtf_seen = {}
    with gtf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = re.search(r'gene_id "([^"]+)"', fields[8])
            if not match:
                continue
            gene = clean_id(match.group(1))
            if gene in eligible_coordinates and "_" not in fields[0]:
                gtf_seen[gene] = (fields[0], int(fields[3]) if fields[6] == "+" else int(fields[4]), fields[6])
    annotation_mismatch = [g for g, r in eligible_coordinates.items()
                           if gtf_seen.get(g) != (r["chrom"], int(r["tss"]), r["strand"])]
    assert not annotation_mismatch, f"Cache/GTF annotation mismatch: {annotation_mismatch[:10]}"
    genes = sorted(eligible_coordinates)
    gene_index = {g: i for i, g in enumerate(genes)}
    promoter_by_chrom = defaultdict(list)
    for i, gene in enumerate(genes):
        r = eligible_coordinates[gene]
        tss = int(r["tss"])
        promoter_by_chrom[r["chrom"]].append((i, max(1, tss - 2000), tss + 2000))
    counts = np.zeros((len(genes), len(tfs)), dtype=np.int64)
    bed_counts = np.zeros_like(counts)
    peak_qc = []
    for j, tf in enumerate(tfs):
        path = data / "ChIP_Atlas" / f"Oth.ALL.05.{tf}.AllCell.bed"
        inputs.append(path)
        by_chrom = defaultdict(lambda: ([], []))
        n_raw = 0
        with path.open() as handle:
            for line in handle:
                if line.startswith(("track", "#")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                n_raw += 1
                if int(f[4]) < 200:
                    continue
                start, end = int(f[1]), int(f[2])
                assert start < end
                by_chrom[f[0]][0].append(start)
                by_chrom[f[0]][1].append(end)
        for chrom, promoters in promoter_by_chrom.items():
            starts = np.sort(np.asarray(by_chrom[chrom][0], dtype=np.int64))
            ends = np.sort(np.asarray(by_chrom[chrom][1], dtype=np.int64))
            p = np.asarray(promoters, dtype=np.int64)
            # Number of peaks with start < promoter end AND end > promoter start.
            counts[p[:, 0], j] = np.searchsorted(starts, p[:, 2], side="left") - np.searchsorted(ends, p[:, 1], side="right")
            # Correct 1-based inclusive interval -> zero-based half-open BED.
            bed_counts[p[:, 0], j] = np.searchsorted(starts, p[:, 2], side="left") - np.searchsorted(ends, p[:, 1] - 1, side="right")
        assert np.all(counts[:, j] >= 0)
        peak_qc.append(dict(tf=tf, raw_peak_rows=n_raw,
                            retained_peak_rows=sum(len(x[0]) for x in by_chrom.values()),
                            genes_with_peak=int((counts[:, j] > 0).sum()),
                            bed_conversion_binding_changes=int(((counts[:, j] > 0) != (bed_counts[:, j] > 0)).sum())))
        print(tf, peak_qc[-1], flush=True)

    # Require exact reproduction of the saved matrix before interpreting changes.
    old_matrix = read("enrichment/chip_atlas_enrichment/peak_gene_matrix.txt")
    id_col = list(old_matrix[0])[0]
    mismatches = []
    for r in old_matrix:
        gene = clean_id(r[id_col])
        for j, tf in enumerate(tfs):
            if int(float(r[tf])) != counts[gene_index[gene], j]:
                mismatches.append((gene, tf))
    assert not mismatches, f"Legacy binding count mismatch: {mismatches[:10]}"
    old_universe = {clean_id(r[id_col]) for r in old_matrix}
    annotated = set(genes)
    universes = {
        "legacy_peak_positive": old_universe,
        "selected_pool_with_zeros": selected_union & annotated,
        "shared_tested": shared_tested & annotated,
    }
    # Audit coverage: the legacy loader chose Gene symbols before ENSG_clean
    # for cancer-specific tables. Keep that bug only as a diagnostic control.
    mapped_selected_positive = {g for g in universes["selected_pool_with_zeros"] if counts[gene_index[g]].sum() > 0}
    legacy_selected_union = set.union(*legacy_selected_sets)
    legacy_annotated = legacy_selected_union & annotated
    legacy_remapped_positive = {g for g in legacy_annotated if counts[gene_index[g]].sum() > 0}
    assert legacy_remapped_positive == old_universe
    mapping_coverage = dict(selected_pool_size=len(selected_union), selected_annotated=len(universes["selected_pool_with_zeros"]),
                            legacy_positive=len(old_universe), remapped_positive=len(mapped_selected_positive),
                            legacy_mixed_id_pool=len(legacy_selected_union), legacy_mixed_id_annotated=len(legacy_annotated),
                            additional_positive=sorted(mapped_selected_positive-old_universe),
                            missing_legacy=sorted(old_universe-mapped_selected_positive))
    (log / "mapping_coverage.json").write_text(json.dumps(mapping_coverage, indent=2)+"\n")
    print("Mapping coverage:", {k:len(v) if isinstance(v,list) else v for k,v in mapping_coverage.items()},flush=True)
    assert not mapping_coverage["missing_legacy"]
    configs = [
        ("legacy_reproduction", "legacy_peak_positive", original_targets, counts),
        ("shared_tested_legacy_identifier_control", "shared_tested", legacy_identifiers, counts),
        ("selected_pool_with_zeros", "selected_pool_with_zeros", original_targets, counts),
        ("shared_tested_original_targets", "shared_tested", original_targets, counts),
        ("shared_tested_all_rhythmic", "shared_tested", rhythmic, counts),
        ("shared_tested_all_rhythmic_bed_coordinates", "shared_tested", rhythmic, bed_counts),
    ]
    results, gene_qc = [], []
    for scenario, universe_name, targets, matrix in configs:
        universe = universes[universe_name]
        bg_idx = [gene_index[g] for g in sorted(universe)]
        M = len(bg_idx)
        for group in groups:
            target = targets[group] & universe
            target_idx = [gene_index[g] for g in sorted(target)]
            N = len(target)
            assert 0 < N < M
            gene_qc.append(dict(scenario=scenario, group=group, full_target_genes=len(targets[group]),
                                target_in_universe=N, excluded_target_genes=len(targets[group] - universe),
                                background_genes=M,
                                zero_binding_background=int((matrix[bg_idx].sum(axis=1) == 0).sum()),
                                zero_binding_target=int((matrix[target_idx].sum(axis=1) == 0).sum())))
            for j, tf in enumerate(tfs + ["ANY_CLOCK_TF"]):
                bound = matrix[:, j] > 0 if j < len(tfs) else matrix.sum(axis=1) > 0
                n, k = int(bound[bg_idx].sum()), int(bound[target_idx].sum())
                table = [[k, N - k], [n - k, M - N - n + k]]
                assert min(v for row in table for v in row) >= 0
                p = float(hypergeom.sf(k - 1, M, n, N))
                fisher = fisher_exact(table, alternative="greater")
                assert np.isclose(p, fisher.pvalue, rtol=1e-10, atol=1e-300)
                results.append(dict(scenario=scenario, group=group, tf=tf, M=M, n=n, N=N, k=k,
                                    target_bound=k, target_unbound=N-k, remainder_bound=n-k,
                                    remainder_unbound=M-N-n+k, target_pct=100*k/N,
                                    background_pct=100*n/M, fold_enrichment=(k/N)/(n/M) if n else float("nan"),
                                    odds_ratio=float(fisher.statistic), pvalue=p,
                                    q_bh_18_tfs=None, q_bh_9_tfs=None, q_bh_2_any=None))
    # Each scenario is a distinct sensitivity analysis. Main FDR family: all 18
    # pre-specified cohort x TF tests; also retain per-cohort BH for legacy comparability.
    for scenario, _, _, _ in configs:
        subset = [r for r in results if r["scenario"] == scenario and r["tf"] != "ANY_CLOCK_TF"]
        for r, q in zip(subset, false_discovery_control([r["pvalue"] for r in subset])):
            r["q_bh_18_tfs"] = float(q)
        for group in groups:
            subset9 = [r for r in subset if r["group"] == group]
            for r, q in zip(subset9, false_discovery_control([r["pvalue"] for r in subset9])):
                r["q_bh_9_tfs"] = float(q)
        combined = [r for r in results if r["scenario"] == scenario and r["tf"] == "ANY_CLOCK_TF"]
        for r, q in zip(combined, false_discovery_control([r["pvalue"] for r in combined])):
            r["q_bh_2_any"] = float(q)
    saved = read("enrichment/chip_atlas_enrichment/enrichment_summary.txt")
    for r in results:
        if r["scenario"] == "legacy_reproduction" and r["tf"] in tfs:
            ref = next(x for x in saved if x["Group"] == r["group"] and x["TF"] == r["tf"] and x["Type"] == "ChIP")
            assert abs(r["fold_enrichment"] - float(ref["Fold_enrichment"])) < 1e-12
            assert abs(r["pvalue"] - float(ref["pvalue"])) < 1e-12

    write_table(out / "enrichment_results.tsv", results)
    write_table(out / "gene_set_qc.tsv", gene_qc)
    write_table(out / "peak_mapping_qc.tsv", peak_qc)
    write_table(out / "eligible_gene_membership.tsv", [dict(ENSG=g,
        tested_gtex=g in tested[groups[0]], tested_tcga=g in tested[groups[1]],
        valid_promoter=g in annotated, shared_background=g in universes["shared_tested"],
        rhythmic_gtex=g in rhythmic[groups[0]], rhythmic_tcga=g in rhythmic[groups[1]],
        tcga_original_target=g in original_targets[groups[1]]) for g in sorted(candidate_ids)])
    write_table(out / "assessed_promoter_peak_counts.tsv", [dict(ENSG=g, **{tf:int(counts[i,j]) for j,tf in enumerate(tfs)}) for i,g in enumerate(genes)])
    input_manifest = []
    for path in sorted(set(inputs)):
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1048576), b""):
                digest.update(chunk)
        input_manifest.append(dict(path=str(path.relative_to(root)),bytes=path.stat().st_size,sha256=digest.hexdigest()))
    write_table(out / "input_manifest.tsv", input_manifest)
    info = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                tested_gtex=len(tested[groups[0]]), tested_tcga=len(tested[groups[1]]),
                shared_tested=len(shared_tested), shared_annotated=len(universes["shared_tested"]),
                selected_annotated=len(universes["selected_pool_with_zeros"]),
                legacy_matrix_genes=len(old_universe), cached_gtf_mismatches=len(annotation_mismatch),
                legacy_peak_count_mismatches=len(mismatches),
                score_threshold=200,promoter_upstream_bp=2000,promoter_downstream_bp=2000,
                primary_coordinates="Exact legacy promoter coordinates; BED conversion sensitivity reported separately",
                tests="One-sided enrichment; hypergeometric and Fisher agree; BH across 18 TF/cohort tests per scenario")
    (log / "validation.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2), flush=True)
    for r in results:
        if r["scenario"] == "shared_tested_all_rhythmic":
            print(r["group"],r["tf"],f"FE={r['fold_enrichment']:.4f}",f"p={r['pvalue']:.5g}",
                  f"q={r['q_bh_18_tfs'] if r['tf'] != 'ANY_CLOCK_TF' else r['q_bh_2_any']}",flush=True)


if __name__ == "__main__":
    main()
