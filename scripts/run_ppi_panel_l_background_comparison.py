#!/usr/bin/env python3
"""Nine-clock PPI test with exact panel L targets/background and legacy controls."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, fisher_exact, hypergeom


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    data = root / "data"
    out = data / "ppi_panel_l_background_comparison"
    out.mkdir(parents=True, exist_ok=True)
    notebooks = {p: digest(p) for p in list(root.glob("*.ipynb")) + list((root / "notebook_publication").glob("*.ipynb"))}
    nine = ["ARNTL", "CLOCK", "CRY1", "DBP", "NPAS2", "NR1D1", "RORA", "RORC", "TEF"]
    ten = ["ARNTL", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "CIART", "NPAS2"]
    partners = {gene: set() for gene in set(nine + ten)}
    background = set()
    bio = data / "PPI/BIOGRID-ORGANISM-Homo_sapiens-5.0.253.tab3.txt"
    physical_rows = 0
    with bio.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Experimental System Type"] != "physical":
                continue
            physical_rows += 1
            a, b = (row[f"Official Symbol Interactor {side}"].upper() for side in ["A", "B"])
            background.update([a, b])
            if a in partners:
                partners[a].add(b)
            if b in partners:
                partners[b].add(a)
    for gene in partners:
        partners[gene].discard(gene)
    coords_path = data / "reference/gencode_gene_coordinates.tsv"
    coords = {r["gene_id"].split(".")[0]: r["gene_name"] for r in read(coords_path)}
    paths = {
        "GTEx_Normal": data / "Merged/regression/sig_genes_regression_periodic_12h_GTEx_Normal.txt",
        "TCGA-PAAD_Tumor": data / "enrichment/cancer_specific_genes/method_A_presence/TCGA-PAAD_Tumor_cancer_specific_genes.txt",
        "CPTAC-3_Tumor": data / "enrichment/cancer_specific_genes/method_A_presence/CPTAC-3_Tumor_cancer_specific_genes.txt",
    }
    targets, mapping = {}, []
    for group, path in paths.items():
        rows = read(path)
        # Preserve panel L's identifier selection and conversion, including direct symbol input.
        key = next(k for k in ["ENSG", "Gene", "gene_id", "ENSG_clean"] if k in rows[0])
        symbols = set()
        for row in rows:
            gene = row[key].strip('"').split(".")[0]
            symbol = coords.get(gene, "") if gene.startswith("ENSG") else gene
            symbol = symbol.upper()
            if symbol:
                symbols.add(symbol)
            mapping.append(dict(group=group, source_column=key, input_identifier=gene,
                                symbol=symbol, in_background=symbol in background))
        targets[group] = symbols & background
    endpoints = {gene: partners[gene] for gene in nine}
    endpoints["POOLED_NINE_CCG"] = set.union(*(partners[g] for g in nine))
    legacy_pool = set.union(*(partners[g] for g in ten))

    def test(group, label, partner_set):
        target = targets[group]
        M, n, N, k = len(background), len(partner_set), len(target), len(target & partner_set)
        assert partner_set <= background
        table = [[k, N-k], [n-k, M-N-n+k]]
        assert min(v for row in table for v in row) >= 0
        p = float(hypergeom.sf(k-1, M, n, N))
        assert np.isclose(p, fisher_exact(table, alternative="greater").pvalue, rtol=1e-10, atol=1e-300)
        return dict(group=group, clock_gene=label, M=M, n=n, N=N, k=k,
                    target_partner=k, target_nonpartner=N-k, remainder_partner=n-k,
                    remainder_nonpartner=M-N-n+k, fold_enrichment=(k/N)/(n/M), pvalue=p)

    results = [test(group, label, values) for group in list(paths)[:2] for label, values in endpoints.items()]
    for row, q in zip(results, false_discovery_control([r["pvalue"] for r in results])):
        row["q_bh_20_tests"] = float(q)
    controls = [test(group, "PANEL_L_TEN_CCG", legacy_pool) for group in paths]
    legacy_path = data / "enrichment/ppi_enrichment/ppi_enrichment_combined.txt"
    legacy = {r["group"]: r for r in read(legacy_path)}
    for row in controls:
        old = legacy[row["group"]]
        for key in ["M", "n", "N", "k"]:
            assert row[key] == int(old[key]), (row["group"], key, row[key], old[key])
        assert np.isclose(row["fold_enrichment"], float(old["fold_enrichment"]), rtol=1e-12)
        # Legacy used 1-cdf; sf avoids cancellation for small p-values.
        assert np.isclose(row["pvalue"], float(old["pvalue"]), rtol=1e-6, atol=1e-15)
    assert all(digest(p) == h for p, h in notebooks.items())
    write(out / "enrichment_results.tsv", results)
    write(out / "panel_l_reproduction_controls.tsv", controls)
    write(out / "target_identifier_mapping.tsv", mapping)
    write(out / "background_symbols.tsv", [dict(symbol=g) for g in sorted(background)])
    write(out / "partner_symbols.tsv", [dict(clock_gene=k, symbol=g) for k, v in endpoints.items() for g in sorted(v)])
    write(out / "target_overlap_symbols.tsv", [dict(group=group, clock_gene=k, symbol=g)
          for group in list(paths)[:2] for k, v in endpoints.items() for g in sorted(targets[group] & v)])
    inputs = [bio, coords_path, legacy_path, *paths.values()]
    write(out / "input_manifest.tsv", [dict(path=str(p.relative_to(root)), sha256=digest(p)) for p in inputs])
    validation = dict(notebooks_checked=len(notebooks), notebooks_changed=0, physical_rows=physical_rows,
                      background_symbols=len(background), panel_l_controls_reproduced=3,
                      tests=len(results), fisher_hypergeometric_agreement=True)
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    report = """# PPI comparison with panel L targets and background

This standalone analysis uses the same nine clock genes as the recent ChIP/PPI comparison plus their unique pooled direct partners. TCGA uses panel L's method_A_presence cancer-specific rhythmic genes; GTEx uses its rhythmic-gene reference. The background is exactly panel L's 22,695 uppercase official-symbol entries from all physical interactions in BioGRID 5.0.253. As required to match panel L, no protein-coding, shared-tested-gene, or both-interactors-human restriction is applied. Identifier selection/conversion also matches panel L. This is a legacy-comparable sensitivity analysis, not a revision of notebooks or publication figures.

Self-interactions are excluded from each seed's partner set; partners are deduplicated by symbol, and a pooled partner counts once. Fold enrichment is (target partner fraction)/(full-background partner fraction). Upper-tail hypergeometric p-values are verified against one-sided Fisher exact tests. BH q-values cover 20 tests (two cohorts × nine individual genes plus one pooled endpoint). Panel L's plotted stars use raw p-values; both raw p and q are supplied here.

The original ten-clock-gene pool is also rerun for GTEx, TCGA and CPTAC. All stored panel L M/n/N/k counts, fold enrichments and p-values reproduce. That control isolates the effect of changing clock-seed membership in this run. See panel_l_reproduction_controls.tsv. The broader legacy background includes symbols from interactions with nonhuman partners and genes not tested in the expression analysis; matching it does not establish it as the preferred inferential background.

| Clock gene | GTEx k/N | Fold | raw p | BH q | TCGA k/N | Fold | raw p | BH q |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
"""
    for label in endpoints:
        pair = [next(r for r in results if r["group"] == group and r["clock_gene"] == label) for group in list(paths)[:2]]
        report += "| " + label + " | " + " | ".join(f"{r['k']}/{r['N']} | {r['fold_enrichment']:.3f} | {r['pvalue']:.3g} | {r['q_bh_20_tests']:.3g}" for r in pair) + " |\n"
    (out / "README.md").write_text(report)
    print(report, flush=True)
    print(json.dumps(validation), flush=True)


if __name__ == "__main__":
    main()
