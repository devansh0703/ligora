# BM25 Search Index

Ligora ships a local BM25 text search over searchable entities. The index is
**built from live sources** — it is a cache, not a data source. Querying
never invents results: every indexed document was fetched from a provider
endpoint listed in `../DATA_SOURCES.md`.

## What is indexed

| Corpus | Document | Fields | Fetched from |
|---|---|---|---|
| `structures` | PDB entry | id, title, citation title, journal, year, method, resolution | RCSB Data API + Search API |
| `chem` | chemical component (CCD) / PubChem compound | id, name, formula, SMILES, classification | RCSB CCD + PubChem PUG REST |
| `workspace` | user's saved artifacts | structure id, title, ligand name, notes | local workspace session files |

## Refresh policy

- `POST/GET search_index_refresh` triggers re-fetching the entities you
  already touched in this workspace (session structures, resolved ligands)
  plus any PDB IDs / compound names the user explicitly asks to add.
- Nothing is prefetched "for the sake of it": a desktop app must not hammer
  public APIs. The index grows with your work.
- Search ranks with BM25 (Okapi, k1=1.5, b=0.75) across fields with weights
  (id > name/title > other fields). Purely local scoring math over real
  documents — no network calls at query time.

## API

- `search({query, scope?, limit?})` — ranked hits with source and fetch date.
- `search_index_refresh({pdb_ids?, chem_names?})` — re-fetch/rebuild.

The UI exposes this in the **Search** tab.
