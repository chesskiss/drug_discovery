# _alias — gene symbol → Entrez resolver

Shared helper. Gene symbols get renamed by HGNC over time (`ACPP` → `ACP3`, `ATP5E` → `ATP5F1E`,
`FGFR1OP` → `CEP43`), but **Entrez ids are stable**. A dataset keyed by current symbols will
silently miss genes whose symbol differs from the target's vintage.

This maps *any* symbol — current or historical synonym — to its Entrez id, so alignment can join
on Entrez instead of fragile symbols.

## Build

```bash
uv run --project data/bio_context _alias/build_alias_map.py
```

Source: NCBI [`Homo_sapiens.gene_info.gz`](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)
(~5 MB), columns `GeneID`, `Symbol`, `Synonyms` (pipe-separated). Both the primary symbol and
every synonym are exploded into rows.

Output: `data/alias_to_entrez.csv` — `alias_symbol`, `entrez`. ~268k rows / 262k unique symbols
covering 194k genes.

## Caveat

~9.7k rows have a symbol that maps to **more than one** gene (an old symbol reused, or ambiguous
shorthand). Consumers should drop ambiguous symbols rather than guess —
`progeny/align_progeny_to_tdc.py` does this via `drop_duplicates("alias_symbol", keep=False)`,
so an ambiguous symbol is treated as unmatched instead of being resolved to an arbitrary gene.

## Impact

On PROGENy: raw symbol matching hits 1,246/1,295 (96.2%); adding this alias step recovers the
remaining **49 genes → 100%**.
