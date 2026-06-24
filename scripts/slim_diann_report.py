#!/usr/bin/env python3
"""Slim a large DIA-NN report TSV down to the columns needed for §1 counting.

Purpose: some deposited "original" analyses ship a full precursor-level DIA-NN
report that is far too large to download per rebuild (e.g. the ProCan original
``ProCan-DepMapSanger_DIANN_output.tsv`` on PRIDE is ~237 GB). This utility
streams that TSV in chunks and writes a column-slimmed parquet carrying only the
columns ``count_report`` (in scripts/rebuild.py) needs, so the reanalysis-figure
baseline can be RECOMPUTED from a small public file instead of a hard-coded
constant.

Run it where the big file is local/fast (e.g. an HPC node), then move the small
output parquet to the quantmsdiann-benchmarks FTP, e.g.:
    cell-lines/PXD030304/original_v1_7/quant_tables/diann_report.parquet

Usage:
    python scripts/slim_diann_report.py INPUT.tsv OUTPUT.parquet
    # INPUT may also be a .parquet (re-slimmed) or a gzipped .tsv.gz

It also prints the §1 counts (prot_global, prot_2pep, prec_global) of the slim
output so the recomputed baseline can be verified immediately (ProCan expected
prot_2pep ~ 6,692).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Columns count_report needs (scripts/rebuild.py:_NEEDED_COLS). Decoy/Proteotypic
# are kept when present; everything else is dropped.
WANTED = [
    'Run', 'Precursor.Id', 'Protein.Group', 'Stripped.Sequence',
    'Q.Value', 'PG.Q.Value', 'Lib.Q.Value', 'Lib.PG.Q.Value',
    'Decoy', 'Proteotypic',
]
STR_COLS = {'Run', 'Precursor.Id', 'Protein.Group', 'Stripped.Sequence'}
FLOAT_COLS = {'Q.Value', 'PG.Q.Value', 'Lib.Q.Value', 'Lib.PG.Q.Value'}
CHUNK = 1_000_000
Q = 0.01


def _present_columns(path: Path) -> list[str]:
    if path.suffix == '.parquet':
        have = set(pq.ParquetFile(path).schema_arrow.names)
    else:
        have = set(pd.read_csv(path, sep='\t', nrows=0).columns)
    cols = [c for c in WANTED if c in have]
    missing = [c for c in ('Protein.Group', 'Stripped.Sequence', 'Lib.PG.Q.Value', 'Lib.Q.Value') if c not in have]
    if missing:
        raise SystemExit(f'ERROR: {path} is missing required column(s): {missing}')
    return cols


def _iter_chunks(path: Path, cols: list[str]):
    if path.suffix == '.parquet':
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(batch_size=CHUNK, columns=cols):
            yield batch.to_pandas()
    else:
        for chunk in pd.read_csv(path, sep='\t', usecols=cols, dtype=str, chunksize=CHUNK):
            yield chunk


def slim(input_path: Path, output_path: Path) -> None:
    cols = _present_columns(input_path)
    print(f'Keeping columns: {cols}', flush=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer: pq.ParquetWriter | None = None
    rows = 0
    try:
        for i, chunk in enumerate(_iter_chunks(input_path, cols)):
            for c in cols:
                if c in FLOAT_COLS:
                    chunk[c] = pd.to_numeric(chunk[c], errors='coerce').astype('float32')
                elif c in STR_COLS:
                    chunk[c] = chunk[c].astype('string')
            table = pa.Table.from_pandas(chunk[cols], preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output_path, table.schema, compression='zstd')
            writer.write_table(table)
            rows += len(chunk)
            if i % 20 == 0:
                print(f'  ... {rows:,} rows', flush=True)
    finally:
        if writer is not None:
            writer.close()
    size_mb = output_path.stat().st_size / 1e6
    print(f'wrote {output_path} ({rows:,} rows, {size_mb:,.1f} MB)', flush=True)


def verify(output_path: Path) -> None:
    """Print the §1 counts of the slim parquet (mirrors count_report)."""
    have = set(pq.ParquetFile(output_path).schema_arrow.names)
    df = pq.read_table(output_path, columns=list(have)).to_pandas()
    if 'Decoy' in df.columns:
        df = df[df['Decoy'].astype(str).isin(('0', '0.0', 'False', 'false'))]
    prot_global = int(df.loc[df['Lib.PG.Q.Value'] <= Q, 'Protein.Group'].nunique())
    prec_global = int(df.loc[df['Lib.Q.Value'] <= Q, 'Precursor.Id'].nunique()) if 'Precursor.Id' in df else -1
    dq = df[df['Lib.Q.Value'] <= Q]
    dp = dq[dq['Lib.PG.Q.Value'] <= Q]
    prot_2pep = int((dp.groupby('Protein.Group')['Stripped.Sequence'].nunique() >= 2).sum())
    gps = df.loc[df['Lib.PG.Q.Value'] <= Q, 'Protein.Group'].dropna().unique()
    multi = 100.0 * sum(';' in str(g) for g in gps) / max(len(gps), 1)
    print('--- §1 counts of slim parquet (Lib rule) ---')
    print(f'  prot_global (Lib.PG.Q.Value<=0.01): {prot_global:,}')
    print(f'  prot_2pep   (>=2 stripped peptides): {prot_2pep:,}')
    print(f'  prec_global (Lib.Q.Value<=0.01):     {prec_global:,}')
    print(f'  multi-accession protein groups:      {multi:.2f}%  ({"OK" if multi <= 1.2 else "OVER 1.2% CAP — flag"})')


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2
    input_path, output_path = Path(argv[1]), Path(argv[2])
    if not input_path.exists():
        raise SystemExit(f'ERROR: input not found: {input_path}')
    slim(input_path, output_path)
    verify(output_path)
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
