# KEGG (gene-to-pathway membership + graph degree)

**KEGG PATHWAY** (e.g. https://www.kegg.jp/pathway/map01100, "Metabolic pathways") is normally
presented as a diagram, but the underlying data is available as plain tables through the
[KEGG REST API](https://rest.kegg.jp) — no graph/image parsing required. This script pulls
human gene → pathway membership across **all 372 KEGG pathways** (not just the metabolic
overview map), via three REST endpoints:

- `link/pathway/hsa` — every human gene linked to every KEGG pathway it belongs to.
- `list/pathway/hsa` — pathway id → pathway name.
- `list/hsa` — gene id → primary gene symbol.

It then downloads **each pathway's own KGML graph** (372 individual files, not the combined
`map01100` overview) and adds a `degree` column: the number of reaction/relation edges that
gene's node touches within that specific pathway's graph.

**What `degree` is, and isn't:** node position/distance in the KEGG diagram is a manual layout
choice (verified by inspecting `graphics` elements — only `x`/`y`/`width`/`height`/color, no
strength field) and carries no biological meaning. Graph *adjacency* does, however: `degree`
counts real edges — `reaction` steps (substrate → enzyme → product) in metabolic pathways, or
typed `relation` edges (activation, inhibition, phosphorylation, binding, etc.) in signaling
pathways. It is a topological signal ("how wired-in is this gene to this pathway's graph"),
**not** a statistically-derived effect size like PROGENy's `weight` — KEGG has no equivalent of
that. Caveats:
- A KGML node can bundle several isozyme genes together; they all get the same `degree`.
- `degree` doesn't distinguish edge direction or sign (`activation` vs `inhibition` count the same).
- 36 of the 372 pathways (e.g. Ribosome, DNA replication, ABC transporters) are genuinely
  edgeless complex-membership diagrams in KEGG, not reaction/signaling networks — every gene in
  those gets `degree=0` correctly, not due to a fetch failure.

- Homepage: https://www.kegg.jp/kegg/pathway.html
- REST API docs: https://www.kegg.jp/kegg/rest/keggapi.html

## KEGG is not just metabolism

KEGG's own category hierarchy (BRITE `br08901`) splits our 372 human pathways as:

| Category | Count | Example |
|---|---:|---|
| Human Diseases | 96 | EGFR tyrosine kinase inhibitor resistance |
| Metabolism | 95 | Glycolysis / Gluconeogenesis |
| Organismal Systems | 86 | PPAR signaling pathway |
| Environmental Information Processing | 37 | ABC transporters |
| Genetic Information Processing | 33 | Aminoacyl-tRNA biosynthesis |
| Cellular Processes | 25 | Cell cycle |

Metabolism is only ~25% of what's here — the rest is signaling, disease mechanisms, transport,
and cellular machinery. Relevant for drug synergy, since most drug targets (receptors, kinases)
live outside the metabolic category.

## What an edge means, concretely

- **Metabolic pathways:** an edge is a `<reaction>` element — a real chemical conversion,
  substrate → enzyme → product (e.g. reaction `R00351`/`R00352` converts compounds `C00024` +
  `C00036` → `C00158`).
- **Signaling pathways:** an edge is a `<relation>` element with a labeled subtype —
  `activation`, `inhibition`, `phosphorylation`, `binding/association`, `expression`,
  `indirect effect`, etc.
- Either way: an edge means *"this gene's enzyme/protein is a documented participant in a
  specific reaction or interaction step drawn in this pathway's diagram."* It says nothing about
  magnitude or direction of effect on the pathway's overall output, only that a documented
  mechanistic link exists.

## `degree` is per-(gene, pathway), not a pathway-wide or cross-pathway total

`degree` is computed separately inside each pathway's own graph. A gene in 10 pathways gets 10
independent rows with 10 independent `degree` values — it is **not** summed across all pathways
a gene belongs to, and it is **not** a distance to some single "pathway node" (no such node
exists; the pathway *is* the whole graph). It is literally: for this one gene, in this one
pathway's graph, how many reaction/relation edges touch its node.

## Degree is not the same as importance or "bottleneck-ness"

Higher degree loosely tracks the "centrality-lethality" pattern seen in biological networks
(highly-connected hub nodes are disproportionately essential), but it is a purely structural
count — it says nothing about flux, kinetics, or redundancy. In particular, **rate-limiting /
committed-step enzymes** (e.g. PFK1 in glycolysis, HMG-CoA reductase in cholesterol synthesis)
are classic low-degree bottlenecks: each catalyzes essentially one reaction, yet because it's the
sole, often allosterically-regulated committed step, small changes in its activity control the
whole pathway's flux disproportionately (the subject of Metabolic Control Analysis, Kacser &
Burns). In graph terms this is closer to a **bottleneck / articulation point** or **betweenness
centrality** — a node with few direct edges that nonetheless sits on the only path connecting two
otherwise-separate parts of the network. `degree` cannot detect this: a degree-1 gene could be a
critical bridge, or a totally redundant dead end, and `degree` looks identical either way.
Capturing that would require a different graph metric (betweenness centrality or explicit
articulation-point detection per pathway) computed from the same KGML data — not yet done here.

## Download

```bash
uv run --project data/bio_context data/bio_context/kegg/download_kegg.py
```

Output: `data/kegg_human_gene_pathways.csv`, one row per (gene, pathway) pair, columns:

| Column | Meaning |
|---|---|
| `gene_id` | KEGG `hsa:<entrez>` id |
| `gene_symbol` | primary gene symbol |
| `pathway_id` | KEGG `path:hsa#####` id |
| `pathway_name` | human-readable pathway name |
| `category` | BRITE `br08901` top-level category (Metabolism, Human Diseases, ...) |
| `is_enzyme` | gene has an EC/enzyme link (`link/enzyme/hsa`) — a drug-target hint, not a filter |
| `degree` | # reaction/relation edges touching that gene's node in this pathway's graph |
| `betweenness` | networkx betweenness centrality (0–1) of that gene's node in this pathway's graph |
| `is_articulation` | removing that gene's node disconnects this pathway's graph (hard bottleneck) |
| `importance` | combined score: `max(1 − 1/(1+degree), betweenness)`, overridden to `1.0` if `is_articulation` |

Fetching the 372 per-pathway KGML files takes a couple of minutes (one request per pathway,
~0.1 s delay between requests to be polite to the free public API), plus two small extra REST
calls (`link/enzyme/hsa` for `is_enzyme`, `get/br:br08901/json` for `category`).

### The graph metrics (`degree` / `betweenness` / `is_articulation` / `importance`)

- **`degree`** is local connectivity — how many edges touch the node. Loosely tracks importance
  (centrality-lethality), but misses low-degree bottlenecks.
- **`betweenness`** catches those bottlenecks: a gene with few edges can still lie on the only path
  bridging two halves of a pathway (e.g. rate-limiting enzymes like PFK1 / HMG-CoA reductase —
  Metabolic Control Analysis, Kacser & Burns). networkx normalizes it to 0–1, comparable across
  pathways of different sizes.
- **`is_articulation`** is the discrete, provable version of a bottleneck: remove this node and the
  pathway graph splits into disconnected pieces. Treated as the strongest signal, so `importance`
  is forced to `1.0` there.
- **Graph construction:** metabolic pathways connect enzymes only *through* shared compound nodes
  (enzyme → compound → enzyme), so the per-pathway graph includes compound/ortholog/group nodes
  for correct connectivity (`map` link-nodes excluded to avoid artificial shortcuts); metrics are
  reported only for `gene` nodes. Isozyme-bundled nodes share their metrics across all their genes.
- **Still not a `weight`.** All four are topology, not a statistically-derived per-gene effect size
  like PROGENy's `weight`. They rank structural centrality within a pathway, not magnitude/direction
  of effect on its output.

## Contents (actual downloaded file)

- **File:** `data/kegg_human_gene_pathways.csv` — ~5 MB, 39,545 rows (KEGG updates the underlying
  annotations over time, so exact counts drift slightly between downloads).
- **9,416 unique genes**, **372 unique pathways**, zero null `category`.
- `category` row split: Human Diseases 13,264 · Organismal Systems 8,527 · Environmental Info
  Processing 6,262 · Metabolism 5,238 · Cellular Processes 3,642 · Genetic Info Processing 2,612.
- `is_enzyme`: 2,491 unique genes flagged (genes that are both EC-linked *and* in ≥1 pathway;
  KEGG lists 3,503 EC-linked human genes total, but the rest are in no pathway).
- `importance`: median 0.5, max 1.0; `betweenness` median 0.0, max 0.78; 11,369 articulation rows.
- Just the metabolic-pathways map (`map01100`/`hsa01100`) alone is much smaller: 1,587 genes,
  37 KB — use `rest.kegg.jp/link/hsa/path:hsa01100` directly if you only want that one map.

## Alignment to the TDC gene axis

`align_kegg_to_tdc.py` places KEGG's per-(gene, pathway) `importance` onto the 23,808-dim TDC
NCI-60 gene axis, so a cell line's expression vector `x` (`CellLine[0]`) compresses to a 372-dim
pathway-activity vector via `z = W @ x`.

```bash
# prerequisites (build once, shared by all datasets)
uv run --project data/bio_context tdc_gene_index/build_tdc_gene_index.py
uv run --project data/bio_context _alias/build_alias_map.py

uv run --project data/bio_context kegg/align_kegg_to_tdc.py
```

Matching strategy — the **mirror image of PROGENy's**. KEGG's `gene_id` is `hsa:<entrez>`, so
Entrez is already embedded and the primary join is a direct number match; symbols are only a
fallback (and KEGG's `gene_symbol` is corrupted for ~625 genes — see below — so we never rely on
it when Entrez is present):
1. direct Entrez: `hsa:<entrez>` from `gene_id` matches an Entrez in [`tdc_gene_index`](../tdc_gene_index) → `matched_via=direct`
2. else KEGG `gene_symbol` == a TDC symbol → `matched_via=symbol`
3. else KEGG symbol → Entrez via the NCBI [`_alias`](../_alias) map → Entrez match → `matched_via=alias`
4. else dropped (logged)

**Result: 8,424/9,416 genes matched (89.5%)** — 8,419 by direct Entrez, 5 by symbol, 0 net by
alias. 372/372 pathways retained. The 992 unmatched are genes whose KEGG Entrez id isn't in
CellMiner's 23,808-gene panel: non-coding RNAs (rRNA/tRNA), readthrough/fusion transcripts,
mitochondrially-encoded genes, and newer/uncharacterized ids. At the (gene, pathway) row grain,
37,880/39,545 rows (95.8%) are matched.

Outputs:
- `data/kegg_tdc_weights.npz` — `W` `[372 x 23808]` float32 (27,662 nonzero), `pathways`
  (`path:hsa#####` ids), `pathway_names`, plus `gene_symbol`/`entrez` per column. Note: `W` holds
  `importance`, which is legitimately `0.0` for degree-0 / edgeless-pathway members, so ~10.2k
  matched cells land as 0 (populated 37,867, nonzero 27,662) — a zero here means "in the pathway
  but structurally peripheral", not "absent".
- `data/kegg_tdc_alignment.csv` — audit trail, one row per source (gene, pathway) pair. Carries
  **two** symbol columns on purpose: `kegg_gene_symbol` (KEGG's raw, possibly-corrupted label,
  kept for traceability) and `gene_symbol` (TDC's clean label, looked up via `tdc_idx`), plus
  `entrez`, `matched_via`, and the structural metrics (`category`, `is_enzyme`, `degree`,
  `betweenness`, `is_articulation`, `importance`).

The script warns loudly if two genes collide on the same `(pathway, tdc_idx)` with *conflicting*
`importance` (would silently lose data). For KEGG the only collisions are 26 rows where a paralog
rescued by symbol lands on a gene already matched directly (e.g. CBS, CRYAA, CCL4L1) — all agree
on `importance` (KEGG bundles them into one KGML node), so harmless.

Verified end-to-end: `X @ W.T` over all 59 cell lines yields a finite `[59 x 372]` activity matrix
with real cross-cell-line variation (Metabolic pathways varies most).

## Other available forms (not used by this script)

- **KGML** (`rest.kegg.jp/get/hsa01100/kgml`, 2.7 MB) — the literal reaction graph (nodes +
  edges) behind the map01100 diagram, in XML. Needed only if reaction/mechanistic structure
  (not just membership) is required; needs an XML/KGML parser.
- The organism-agnostic reference graph (`map01100`) is not servable as KGML directly — only
  per-organism instances (e.g. `hsa01100`) are.
