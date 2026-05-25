# 🧬 KEGGtangled

[![PyPI version](https://badge.fury.io/py/keggtangled.svg)](https://pypi.org/project/keggtangled/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Tie all your KEGG data together. Untangle it your way.**

`keggtangled` is a Python library that lets you **download, cache, and interrogate** the entire metabolic network of any organism in the [KEGG](https://www.genome.jp/kegg/) database.  
It builds a rich, in‑memory knowledge graph connecting **genes → KOs → reactions → compounds**, so you can answer biological questions with clean, Pythonic code.

---

## 📦 Installation

```bash
pip install keggtangled
```

`keggtangled` requires Python 3.8 or later and depends on `biopython`.  
Optionally, install `tqdm` for progress bars during bulk downloads:

```bash
pip install tqdm
```

---

## 🚀 Quick start

```python
import keggtangled as kt

# Create an organism – first run fetches and caches all KO‑gene links and reactions
mta = kt.Organism("mta", cache_dir="kegg_cache")

# Load a pathway (downloads + parses KGML, caches everything)
glycolysis = mta.load_pathway("mta00010")
print(glycolysis)
# Pathway(mta00010, 44 genes, 23 reactions)

# From gene → KOs → reactions → compounds
gene = "Moth_0033"
kos = mta.get_kos_for_gene(gene)            # frozenset of ko:K...
rxns = mta.get_reactions_for_gene(gene)     # frozenset of reaction IDs

# Get a specific reaction and its readable formula
rxn = mta.get_reaction("R00200")
print(rxn.formula_per_pathway["mta00010"]["formula_read"])
# e.g., "ATP + D-Glucose --> ADP + D-Glucose 6-phosphate"

# Work with a compound
atp = mta.get_compound("C00074")
print(f"{atp.name} – {atp.formula} – mass {atp.mass}")
# ATP – C10H16N5O13P3 – mass 507.181

# Which genes are linked to ATP?
print(atp.get_genes())    # frozenset of locus tags
```

---

## 🧠 What can you do with KEGGtangled?

- **Fetch & cache all KEGG data** for any organism code (e.g., `eco`, `hsa`, `mta`) with a single line.  
  No more manual downloads or API wrangling.

- **Trace metabolic connections** effortlessly:
  - Gene → KOs → Reactions → Compounds (and back)
  - Built‑in walkers: `get_reactions_for_gene()`, `get_compounds_for_reaction()`, `get_genes_for_reaction()`

- **Inspect pathways in detail**:
  - Readable reaction formulas (substrate / product names)
  - Reversible vs. irreversible arrows (`<=>` / `-->`)
  - Gene lists per pathway, pathways per gene, etc.

- **Save & resume your work**:
  - `organism.save("my_eco.pkl")` snapshots the fully‑loaded state
  - `Organism.load("my_eco.pkl")` restores it instantly – no network calls, no recomputation

- **Iterate over everything**:
  ```python
  for comp in mta.compounds:
      print(comp.id, comp.name)
  for rxn in mta.all_reactions:
      print(rxn.reaction_id)
  ```

- **Get a quick overview** with `mta.summary()` or just `print(mta)`.

---

## 🔍 Key API at a glance

| Method / Property | Description |
|------------------|-------------|
| `Organism(org_code, cache_dir)` | Create or load a KEGG organism |
| `.load_pathway(pathway_id)` | Load a pathway (genes, reactions, compounds, formulas) |
| `.load_all_cached_pathways()` | Load every pathway already cached on disk |
| `.get_compound("C00022")` | Retrieve a compound (fetches details if needed) |
| `.get_reaction("R00200")` | Retrieve a reaction |
| `.get_genes_for_ko(ko)` | All genes annotated with a KO |
| `.get_reactions_for_ko(ko)` | All reactions linked to a KO |
| `.get_kos_for_gene(locus_tag)` | KOs linked to a gene |
| `.get_reactions_for_gene(locus_tag)` | Reactions linked to a gene (via KOs) |
| `.get_compounds_for_reaction(rxn_id)` | Compound IDs involved in a reaction |
| `.get_genes_for_reaction(rxn_id)` | Genes linked to a reaction |
| `.get_pathways_for_gene(locus_tag)` | Pathways containing a gene |
| `.get_pathways_for_reaction(rxn_id)` | Pathways containing a reaction |
| `.save(filename)` | Pickle the entire organism with version check |
| `Organism.load(filename)` | Load a previously saved organism |
| `.compounds` / `.all_reactions` / `.all_pathways` | Iterate over loaded data |
| `.summary()` | Dictionary with counts of pathways, reactions, compounds, genes, KOs |

---

## 📂 Caching & performance

`keggtangled` stores everything it downloads in the `cache_dir` folder:
- KO‑gene and KO‑reaction mappings (JSON)
- Pathway flat files and parsed KGML (both raw XML and pickled parsed objects)
- Compound details (name, formula, mass)

After the first run, subsequent `Organism()` creations will load from cache in seconds.  
When you save an organism with `.save()`, **all in‑memory computed relations are preserved**, so loading is nearly instant.

---

## 🧪 Example: discover all compounds linked to a gene

```python
# Find every compound that 'Moth_0033' touches
reactions = mta.get_reactions_for_gene("Moth_0033")
all_compounds = set()
for rxn in reactions:
    all_compounds.update(mta.get_compounds_for_reaction(rxn))

for cid in sorted(all_compounds):
    comp = mta.get_compound(cid)
    print(f"{comp.id:10s} {comp.name or 'unknown'}")
```

---

## 📚 Requirements & dependencies

- Python ≥ 3.8
- [Biopython](https://biopython.org/) (for KEGG API access and KGML parsing)
- Optional: [`tqdm`](https://github.com/tqdm/tqdm) – shows progress bars during large fetches.

---

## 🤝 Contributing

Bug reports, feature requests, and pull requests are welcome!  
Please check the [issues page](https://github.com/emarquezz/keggtangled/issues) first.

---

## 📄 License

`keggtangled` is released under the MIT License.  
KEGG data is provided by Kanehisa Laboratories – please see their [terms of use](https://www.kegg.jp/kegg/legal.html).

