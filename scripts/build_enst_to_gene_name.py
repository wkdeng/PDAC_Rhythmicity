"""Parse gencode.v49.primary.gtf and write a transcript_id -> gene_name lookup TSV.

One-time helper. Output:
  data/Ref/enst_to_gene_name.tsv

Columns:
  transcript_id (no version)  gene_id (no version)  gene_name
"""
from __future__ import annotations
import re
from pathlib import Path
import time

GTF      = Path('data/Ref/gencode.v49.primary.gtf')
OUT_TSV  = Path('data/Ref/enst_to_gene_name.tsv')

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GN_RE = re.compile(r'gene_name "([^"]+)"')
GI_RE = re.compile(r'gene_id "([^"]+)"')


def strip_version(s: str) -> str:
    return s.split('.', 1)[0]


def main() -> None:
    t0 = time.time()
    print(f'Parsing {GTF} ...')
    seen: dict[str, tuple[str, str]] = {}
    n_lines = 0
    n_tx = 0
    with GTF.open() as fh:
        for line in fh:
            n_lines += 1
            if line.startswith('#'):
                continue
            # quick filter: 3rd tab-delimited field is feature type
            try:
                second_tab = line.index('\t', line.index('\t') + 1)
                third_tab = line.index('\t', second_tab + 1)
            except ValueError:
                continue
            if line[second_tab + 1:third_tab] != 'transcript':
                continue
            n_tx += 1
            tx_m = TX_RE.search(line)
            gn_m = GN_RE.search(line)
            gi_m = GI_RE.search(line)
            if not (tx_m and gn_m and gi_m):
                continue
            tx = strip_version(tx_m.group(1))
            gn = gn_m.group(1)
            gi = strip_version(gi_m.group(1))
            if tx not in seen:
                seen[tx] = (gi, gn)
            if n_tx % 50000 == 0:
                print(f'  {n_tx:>8,} transcript records scanned ({n_lines:,} lines)')
    elapsed = time.time() - t0
    print(f'Parsed {n_tx:,} transcript records / {n_lines:,} lines in {elapsed:.1f}s')
    print(f'Unique transcript IDs: {len(seen):,}')

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open('w') as out:
        out.write('transcript_id\tgene_id\tgene_name\n')
        for tx, (gi, gn) in sorted(seen.items()):
            out.write(f'{tx}\t{gi}\t{gn}\n')
    print(f'Wrote -> {OUT_TSV}')


if __name__ == '__main__':
    main()
