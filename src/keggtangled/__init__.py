#!/usr/bin/env python
# coding: utf-8
"""
keggtangled v0.6.0 – Fully linked KEGG metabolism library.

Now every object "knows" about its neighbours:
- Compound ↔ reactions, pathways, genes
- Reaction ↔ compounds, pathways, genes
- Pathway ↔ compounds, reactions, genes
- Gene ↔ reactions, pathways, compounds (via Organism)

All relationships are stored as frozensets of string IDs,
ensuring easy serialisation and no circular pickling issues.
"""

__version__ = "0.6.0"

import re
import os
import json
import io
import hashlib
import logging
import pickle
from collections import defaultdict
from typing import Dict, Set, FrozenSet, Tuple, Optional, Iterator, Union, List

from Bio.KEGG.REST import kegg_link, kegg_get
from Bio.KEGG.KGML.KGML_parser import read as kgml_read

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# ----------------------------------------------------------------------
# Colored logging
# ----------------------------------------------------------------------
class ColoredFormatter(logging.Formatter):
    COLORS = {
        'DEBUG': '\033[36m',
        'INFO': '\033[32m',
        'WARNING': '\033[33m',
        'ERROR': '\033[31m',
        'CRITICAL': '\033[1;31m'
    }
    RESET = '\033[0m'

    def format(self, record):
        levelname = record.levelname
        if levelname in self.COLORS:
            record.levelname = f"{self.COLORS[levelname]}{levelname}{self.RESET}"
            record.msg = f"{self.COLORS.get(levelname, '')}{record.msg}{self.RESET}"
        return super().format(record)

logger = logging.getLogger("keggtangled")
handler = logging.StreamHandler()
handler.setFormatter(ColoredFormatter('%(levelname)s: %(message)s'))
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# ----------------------------------------------------------------------
# Regex
# ----------------------------------------------------------------------
KO_RE = re.compile(r'\[KO:(K\d+)\]')

# ----------------------------------------------------------------------
# Compound (now fully linked)
# ----------------------------------------------------------------------
class Compound:
    __slots__ = ('id', 'organism', 'name', 'formula', 'mass', 'bigg_id',
                 '_reaction_ids', '_pathway_ids', '_gene_ids')

    def __init__(self, compound_id: str, organism: 'Organism') -> None:
        self.id = compound_id                # bare ID, e.g. 'C00011'
        self.organism = organism
        self.name: Optional[str] = None
        self.formula: Optional[str] = None
        self.mass: Optional[str] = None
        self.bigg_id: Optional[str] = None   # user can fill this in later
        self._reaction_ids: Optional[FrozenSet[str]] = None
        self._pathway_ids: Optional[FrozenSet[str]] = None
        self._gene_ids: Optional[FrozenSet[str]] = None

    def __repr__(self) -> str:
        return f"Compound({self.id}, {self.name})"

    @property
    def reactions(self) -> FrozenSet[str]:
        """Reaction IDs that involve this compound."""
        if self._reaction_ids is None:
            # Should be filled during Organism.finalize()
            return frozenset()
        return self._reaction_ids

    @property
    def pathways(self) -> FrozenSet[str]:
        """Pathway IDs that contain this compound (via reactions)."""
        if self._pathway_ids is None:
            return frozenset()
        return self._pathway_ids

    @property
    def genes(self) -> FrozenSet[str]:
        """Gene IDs associated with reactions involving this compound."""
        if self._gene_ids is None:
            return frozenset()
        return self._gene_ids

    # Legacy methods for compatibility (still useful for direct lookups)
    def get_kos(self) -> FrozenSet[str]:
        kos = set()
        for rxn_id in self.reactions:
            kos.update(self.organism._reaction_to_kos.get(rxn_id, frozenset()))
        return frozenset(kos)

    def get_genes(self) -> FrozenSet[str]:   # kept as method to not conflict with property 'genes'
        return self.genes                    # now identical to property


# ----------------------------------------------------------------------
# Reaction (now holds substrate/product IDs, pathways, etc.)
# ----------------------------------------------------------------------
class Reaction:
    __slots__ = ('reaction_id', 'organism', 'ko_to_genes',
                 'formula_per_pathway',
                 '_substrate_ids', '_product_ids', '_pathway_ids',
                 '_kos_cache', '_genes_cache')

    def __init__(self, reaction_id: str, organism: 'Organism') -> None:
        self.reaction_id = reaction_id
        self.organism = organism
        self.ko_to_genes: Dict[str, Set[str]] = {}
        self.formula_per_pathway: Dict[str, dict] = {}
        self._substrate_ids: Set[str] = set()     # bare compound IDs
        self._product_ids: Set[str] = set()
        self._pathway_ids: Optional[FrozenSet[str]] = None
        self._kos_cache: Optional[Tuple[str, ...]] = None
        self._genes_cache: Optional[FrozenSet[str]] = None

        kos = organism._reaction_to_kos.get(reaction_id, frozenset())
        for ko in kos:
            genes = organism.get_genes_for_ko(ko)
            if genes:
                self.ko_to_genes[ko] = genes

    def __repr__(self) -> str:
        return (f"Reaction({self.reaction_id}, {self.organism.org_code}) – "
                f"{len(self.ko_to_genes)} KOs mapped to genes, "
                f"{len(self.formula_per_pathway)} pathway formulas")

    @property
    def substrates(self) -> FrozenSet[str]:
        """Bare compound IDs of substrates."""
        return frozenset(self._substrate_ids)

    @property
    def products(self) -> FrozenSet[str]:
        """Bare compound IDs of products."""
        return frozenset(self._product_ids)

    @property
    def compounds(self) -> FrozenSet[str]:
        """All compound IDs (substrates + products)."""
        return self.substrates | self.products

    @property
    def pathways(self) -> FrozenSet[str]:
        """Pathway IDs this reaction belongs to."""
        if self._pathway_ids is None:
            return frozenset()
        return self._pathway_ids

    def get_genes(self) -> FrozenSet[str]:
        if self._genes_cache is None:
            all_genes = set()
            for genes in self.ko_to_genes.values():
                all_genes.update(genes)
            self._genes_cache = frozenset(all_genes)
        return self._genes_cache

    @property
    def genes(self) -> FrozenSet[str]:
        """Gene IDs associated with this reaction."""
        return self.get_genes()

    def get_kos(self) -> Tuple[str, ...]:
        if self._kos_cache is None:
            self._kos_cache = tuple(self.ko_to_genes.keys())
        return self._kos_cache


# ----------------------------------------------------------------------
# Pathway (now knows its compounds)
# ----------------------------------------------------------------------
class Pathway:
    __slots__ = ('id', 'title', 'description', 'dblinks',
                 '_gene_ids', '_reaction_ids', '_compound_ids',
                 'organism')

    def __init__(self, pathway_id: str, gene_kos: Optional[Dict[str, Set[str]]] = None,
                 organism: Optional['Organism'] = None) -> None:
        self.id = pathway_id
        self.title: Optional[str] = None
        self.description: Optional[str] = None
        self.dblinks: Dict[str, str] = {}
        if gene_kos is None:
            gene_kos = {}
        self._gene_ids: FrozenSet[str] = frozenset(gene_kos.keys())
        self._reaction_ids: Union[Set[str], FrozenSet[str]] = set()
        self._compound_ids: Optional[FrozenSet[str]] = None
        self.organism: Optional['Organism'] = organism

    def __setstate__(self, state):
        """Handle unpickling of old Pathway objects that lack the 'organism' slot."""
        if isinstance(state, dict):
            for slot in self.__slots__:
                setattr(self, slot, state.get(slot, None))
        else:
            expected = len(self.__slots__)
            if len(state) < expected:
                state = state + (None,) * (expected - len(state))
            for slot, value in zip(self.__slots__, state):
                setattr(self, slot, value)
        if not hasattr(self, 'organism'):
            self.organism = None

    @property
    def gene_ids(self) -> FrozenSet[str]:
        return self._gene_ids

    @gene_ids.setter
    def gene_ids(self, value):
        self._gene_ids = frozenset(value)

    @property
    def reaction_ids(self) -> FrozenSet[str]:
        return frozenset(self._reaction_ids) if not isinstance(self._reaction_ids, frozenset) else self._reaction_ids

    def add_reactions(self, reaction_ids: Union[str, Set[str]]) -> None:
        if isinstance(reaction_ids, str):
            self._reaction_ids.add(reaction_ids)
        else:
            self._reaction_ids.update(reaction_ids)

    @property
    def compounds(self) -> FrozenSet[str]:
        """Compound IDs present in this pathway (union of all reaction compounds)."""
        if self._compound_ids is None:
            if self.organism is None:
                return frozenset()
            # compute lazily
            cids = set()
            for rxn_id in self.reaction_ids:
                rxn = self.organism.reactions.get(rxn_id)
                if rxn:
                    cids.update(rxn.substrates | rxn.products)
            self._compound_ids = frozenset(cids)
        return self._compound_ids

    @compounds.setter
    def compounds(self, value):
        """Allow setting (used by Organism.finalize)."""
        self._compound_ids = frozenset(value)

    def __repr__(self) -> str:
        return f"Pathway({self.id}, {len(self.gene_ids)} genes, {len(self.reaction_ids)} reactions)"

    # ------------------------------------------------------------------
    # DataFrame helper (updated with clean compound IDs)
    # ------------------------------------------------------------------
    def get_reactions_df(self, other_pathways: bool = True, include_genes: bool = True):
        if self.organism is None:
            raise RuntimeError(
                "This Pathway object is not linked to an Organism. "
                "Reload the pathway via organism.load_pathway() instead."
            )
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("pandas is required for get_reactions_df(). Install with `pip install pandas`.")

        rows = []
        for rxn_id in sorted(self.reaction_ids):
            rxn = self.organism.reactions.get(rxn_id)
            if rxn is None:
                continue

            genes = rxn.genes if include_genes else frozenset()
            gene_str = ", ".join(sorted(genes)) if genes else ""

            if not other_pathways:
                info = rxn.formula_per_pathway.get(self.id)
                if info is None:
                    continue
                row = {
                    'Reaction': rxn_id,
                    'Type': info.get('type', ''),
                    'Substrates_kegg': info.get('substrates', []),
                    'Substrates_read': info.get('substrates_read', []),
                    'Products_kegg': info.get('products', []),
                    'Products_read': info.get('products_read', []),
                    'Formula_read': info.get('formula_read', ''),
                    'Formula_kegg': info.get('formula_kegg', ''),
                }
                if include_genes:
                    row['Genes'] = gene_str
                rows.append(row)
            else:
                for pw_id, info in rxn.formula_per_pathway.items():
                    row = {
                        'Reaction': rxn_id,
                        'Pathway': pw_id,
                        'Type': info.get('type', ''),
                        'Substrates_kegg': info.get('substrates', []),
                        'Substrates_read': info.get('substrates_read', []),
                        'Products_kegg': info.get('products', []),
                        'Products_read': info.get('products_read', []),
                        'Formula_read': info.get('formula_read', ''),
                        'Formula_kegg': info.get('formula_kegg', ''),
                    }
                    if include_genes:
                        row['Genes'] = gene_str
                    rows.append(row)

        df = pd.DataFrame(rows)

        if other_pathways and not df.empty:
            list_cols = ['Substrates_kegg', 'Substrates_read', 'Products_kegg', 'Products_read']
            for col in list_cols:
                if col in df.columns:
                    df[col] = df[col].apply(tuple)

            group_cols = [c for c in df.columns if c != 'Pathway']
            df = (
                df.groupby(group_cols, as_index=False)['Pathway']
                  .agg(lambda pws: ', '.join(sorted(pws)))
            )

            for col in list_cols:
                if col in df.columns:
                    df[col] = df[col].apply(list)

            first_cols = ['Reaction', 'Pathway']
            other_cols = [c for c in df.columns if c not in first_cols]
            df = df[first_cols + other_cols]

        return df


# ----------------------------------------------------------------------
# KGML fetcher (unchanged)
# ----------------------------------------------------------------------
def get_pathway_kgml(pathway_id: str, cache_dir: str = "kegg_cache"):
    os.makedirs(cache_dir, exist_ok=True)
    full_id = f"path:{pathway_id}" if not pathway_id.startswith("path:") else pathway_id
    raw_cache_file = os.path.join(cache_dir, f"{pathway_id}.kgml")
    pkl_cache_file = os.path.join(cache_dir, f"{pathway_id}.kgml.pkl")

    if os.path.exists(pkl_cache_file):
        try:
            with open(pkl_cache_file, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            logger.warning(f"Could not load cached KGML pickle for {pathway_id}, re‑parsing: {e}")

    if not os.path.exists(raw_cache_file):
        raw_kgml = kegg_get(full_id, "kgml").read()
        with open(raw_cache_file, 'w', encoding='utf-8') as f:
            f.write(raw_kgml)
    else:
        with open(raw_cache_file, 'r', encoding='utf-8') as f:
            raw_kgml = f.read()

    kgml = kgml_read(io.StringIO(raw_kgml))

    try:
        with open(pkl_cache_file, 'wb') as f:
            pickle.dump(kgml, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        logger.warning(f"Could not save parsed KGML pickle: {e}")

    return kgml


# ----------------------------------------------------------------------
# Organism (now builds full linkage during finalize)
# ----------------------------------------------------------------------
class Organism:
    def __init__(self, org_code: str, batch_size: int = 10, cache_dir: str = "kegg_cache"):
        self.org_code = org_code
        self.batch_size = batch_size
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        self.pathways: Dict[str, Pathway] = {}
        self.reactions: Dict[str, Reaction] = {}
        self._pathway_reaction_map: Dict[str, Set[str]] = defaultdict(set)
        self._gene_pathway_map: Dict[str, Set[str]] = defaultdict(set)

        self._ko_to_genes: Dict[str, FrozenSet[str]] = {}
        self._gene_to_kos: Dict[str, FrozenSet[str]] = {}
        self._ko_to_reactions: Dict[str, FrozenSet[str]] = {}
        self._reaction_to_kos: Dict[str, FrozenSet[str]] = {}

        self._compounds: Dict[str, Compound] = {}
        self._compound_cache_file = os.path.join(cache_dir, f"{org_code}_compounds.json")

        self._is_finalized: bool = False

        self._load_all_ko_genes()
        self._prefetch_all_ko_reactions()
        self._load_compounds()
        self._freeze_mappings()

    # ------------------------------------------------------------------
    # Caching helper (unchanged)
    # ------------------------------------------------------------------
    def _cache_get(self, key: str, subdir: str = "", fetcher_func=None, *args, **kwargs) -> str:
        cache_subdir = os.path.join(self.cache_dir, subdir)
        os.makedirs(cache_subdir, exist_ok=True)
        key_str = str(key)
        key_hash = hashlib.md5(key_str.encode()).hexdigest()
        cache_file = os.path.join(cache_subdir, f"{key_hash}.txt")
        if os.path.exists(cache_file):
            with open(cache_file, 'r', encoding='utf-8') as f:
                return f.read()
        raw_data = fetcher_func(*args, **kwargs)
        with open(cache_file, 'w', encoding='utf-8') as f:
            f.write(raw_data)
        return raw_data

    # ------------------------------------------------------------------
    # KO–gene mapping (unchanged)
    # ------------------------------------------------------------------
    def _load_all_ko_genes(self) -> None:
        json_file = os.path.join(self.cache_dir, f"{self.org_code}_ko_genes.json")
        if os.path.exists(json_file):
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self._ko_to_genes = {ko: frozenset(genes) for ko, genes in data.items()}
        else:
            raw = kegg_link("ko", self.org_code).read().strip()
            ko_to_genes = defaultdict(set)
            if raw:
                for line in raw.splitlines():
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    gene_part = parts[0]
                    ko_part = parts[1]
                    locus = gene_part.split(":", 1)[-1]
                    ko_to_genes[ko_part].add(locus)
            self._ko_to_genes = {ko: frozenset(genes) for ko, genes in ko_to_genes.items()}
            json_data = {ko: list(genes) for ko, genes in self._ko_to_genes.items()}
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(json_data, f, indent=2)

        self._gene_to_kos = defaultdict(set)
        for ko, genes in self._ko_to_genes.items():
            for gene in genes:
                self._gene_to_kos[gene].add(ko)
        self._gene_to_kos = {gene: frozenset(kos) for gene, kos in self._gene_to_kos.items()}

    # ------------------------------------------------------------------
    # KO–reaction mapping (unchanged)
    # ------------------------------------------------------------------
    def _prefetch_all_ko_reactions(self) -> None:
        all_kos = list(self._ko_to_genes.keys())
        if not all_kos:
            return

        ko_reactions_file = os.path.join(self.cache_dir, f"{self.org_code}_ko_reactions.json")
        reaction_kos_file  = os.path.join(self.cache_dir, f"{self.org_code}_reaction_kos.json")

        if os.path.exists(ko_reactions_file) and os.path.exists(reaction_kos_file):
            with open(ko_reactions_file, 'r') as f:
                data = json.load(f)
                self._ko_to_reactions = {ko: frozenset(rxn_list) for ko, rxn_list in data.items()}
            with open(reaction_kos_file, 'r') as f:
                data = json.load(f)
                self._reaction_to_kos = {rn: frozenset(ko_list) for rn, ko_list in data.items()}
            return

        ko_to_reactions: Dict[str, Set[str]] = defaultdict(set)
        batches = list(range(0, len(all_kos), self.batch_size))
        batch_iter = tqdm(batches, desc="Fetching KO‑reaction links", unit="batch") if tqdm else batches

        for i in batch_iter:
            batch = all_kos[i:i + self.batch_size]
            ko_query = "+".join(batch)
            url_key = f"link/rn/{ko_query}"
            raw = self._cache_get(
                url_key,
                subdir="ko_reactions_batches",
                fetcher_func=lambda: kegg_link("rn", ko_query).read().strip()
            )
            if not raw:
                for ko in batch:
                    self._get_reactions_for_ko_fallback(ko, ko_to_reactions)
                continue

            for line in raw.splitlines():
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                ko_id = parts[0]
                rn_id = parts[1].split(":")[1]
                ko_to_reactions[ko_id].add(rn_id)

        self._ko_to_reactions = {ko: frozenset(rxns) for ko, rxns in ko_to_reactions.items()}

        reaction_to_kos: Dict[str, Set[str]] = defaultdict(set)
        for ko, rxn_set in self._ko_to_reactions.items():
            for rn in rxn_set:
                reaction_to_kos[rn].add(ko)
        self._reaction_to_kos = {rn: frozenset(kos) for rn, kos in reaction_to_kos.items()}

        with open(ko_reactions_file, 'w') as f:
            json.dump({ko: list(rxns) for ko, rxns in self._ko_to_reactions.items()}, f, indent=2)
        with open(reaction_kos_file, 'w') as f:
            json.dump({rn: list(kos) for rn, kos in self._reaction_to_kos.items()}, f, indent=2)

    def _get_reactions_for_ko_fallback(self, ko: str, ko_to_reactions_dict: Dict[str, Set[str]]) -> None:
        raw = kegg_link("rn", ko).read().strip()
        for line in raw.splitlines():
            if line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    rn_id = parts[1].split(":")[1]
                    ko_to_reactions_dict[ko].add(rn_id)

    # ------------------------------------------------------------------
    # Internal freezing (unchanged)
    # ------------------------------------------------------------------
    def _freeze_mappings(self) -> None:
        self._ko_to_genes = {k: frozenset(v) for k, v in self._ko_to_genes.items()}
        self._gene_to_kos = {k: frozenset(v) for k, v in self._gene_to_kos.items()}
        self._ko_to_reactions = {k: frozenset(v) for k, v in self._ko_to_reactions.items()}
        self._reaction_to_kos = {k: frozenset(v) for k, v in self._reaction_to_kos.items()}

    # ------------------------------------------------------------------
    # Compound management (enhanced with bigg_id)
    # ------------------------------------------------------------------
    def _load_compounds(self) -> None:
        if os.path.exists(self._compound_cache_file):
            with open(self._compound_cache_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            for cid, info in data.items():
                comp = Compound(cid, self)
                comp.name = info.get("name")
                comp.formula = info.get("formula")
                comp.mass = info.get("mass")
                comp.bigg_id = info.get("bigg_id")
                self._compounds[cid] = comp

    def _save_compounds(self) -> None:
        data = {}
        for cid, comp in self._compounds.items():
            if comp.name:
                data[cid] = {
                    "name": comp.name,
                    "formula": comp.formula,
                    "mass": comp.mass,
                    "bigg_id": comp.bigg_id
                }
        with open(self._compound_cache_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

    def get_compound(self, compound_id: str, fetch_if_missing: bool = True) -> Compound:
        # Ensure bare ID (strip cpd: if accidentally given)
        cid = compound_id.split(':')[-1] if ':' in compound_id else compound_id
        if cid not in self._compounds:
            self._compounds[cid] = Compound(cid, self)
            if fetch_if_missing:
                self._fetch_compound_details(cid)
        return self._compounds[cid]

    def _fetch_compound_details(self, compound_id: str) -> None:
        comp = self._compounds.get(compound_id)
        if not comp:
            return
        try:
            full_id = f"cpd:{compound_id}" if not compound_id.startswith("cpd:") else compound_id
            raw = kegg_get(full_id).read()
            for line in raw.splitlines():
                if line.startswith("NAME") and comp.name is None:
                    parts = line.split(maxsplit=1)
                    if len(parts) > 1:
                        comp.name = parts[1].strip()
                elif line.startswith("FORMULA"):
                    parts = line.split(maxsplit=1)
                    if len(parts) > 1:
                        comp.formula = parts[1].strip()
                elif line.startswith("MASS"):
                    parts = line.split()
                    if len(parts) >= 2:
                        comp.mass = parts[1]
            self._save_compounds()
        except Exception as e:
            logger.warning(f"Could not fetch details for compound {compound_id}: {e}")

    def _parse_compound_names_from_text(self, flat_text: str) -> None:
        in_compound = False
        for line in flat_text.splitlines():
            if line.startswith("COMPOUND"):
                in_compound = True
                parts = line.split()
                if len(parts) >= 3:
                    cid = parts[1]
                    name = ' '.join(parts[2:])
                    comp = self.get_compound(cid, fetch_if_missing=False)
                    comp.name = name
                continue
            if in_compound:
                if line.startswith(" "):
                    parts = line.split()
                    if len(parts) >= 2:
                        cid = parts[0]
                        name = ' '.join(parts[1:])
                        comp = self.get_compound(cid, fetch_if_missing=False)
                        comp.name = name
                else:
                    in_compound = False
        self._save_compounds()

    # ------------------------------------------------------------------
    # Public lookups
    # ------------------------------------------------------------------
    def get_genes_for_ko(self, ko: str) -> FrozenSet[str]:
        return self._ko_to_genes.get(ko, frozenset())

    def get_reactions_for_ko(self, ko: str) -> FrozenSet[str]:
        return self._ko_to_reactions.get(ko, frozenset())

    def get_kos_for_gene(self, locus_tag: str) -> FrozenSet[str]:
        return self._gene_to_kos.get(locus_tag, frozenset())

    def get_reactions_for_gene(self, locus_tag: str) -> FrozenSet[str]:
        reactions = set()
        for ko in self.get_kos_for_gene(locus_tag):
            reactions.update(self.get_reactions_for_ko(ko))
        return frozenset(reactions)

    def get_compounds_for_reaction(self, reaction_id: str) -> FrozenSet[str]:
        rxn = self.reactions.get(reaction_id)
        if not rxn:
            return frozenset()
        return rxn.substrates | rxn.products

    def get_genes_for_reaction(self, reaction_id: str) -> FrozenSet[str]:
        rxn = self.reactions.get(reaction_id)
        if rxn:
            return rxn.genes
        return frozenset()

    # ------------------------------------------------------------------
    # Pathway loading (modified to strip cpd: prefix and collect compound IDs)
    # ------------------------------------------------------------------
    def _parse_gene_kos_from_text(self, flat_text: str) -> Dict[str, Set[str]]:
        gene_kos = {}
        in_gene_section = False
        for line in flat_text.splitlines():
            if line.startswith("GENE"):
                in_gene_section = True
                parts = line.split()
                if len(parts) > 1:
                    gene_id = parts[1]
                    kos = KO_RE.findall(line)
                    if kos:
                        gene_kos[gene_id] = {f"ko:{ko}" for ko in kos}
                continue
            if in_gene_section:
                if line.startswith(" "):
                    parts = line.split()
                    if parts:
                        gene_id = parts[0]
                        kos = KO_RE.findall(line)
                        if kos:
                            gene_kos[gene_id] = {f"ko:{ko}" for ko in kos}
                else:
                    in_gene_section = False
        return gene_kos

    def load_pathway(self, pathway_id: str) -> Pathway:
        if pathway_id in self.pathways:
            return self.pathways[pathway_id]

        cache_file = os.path.join(self.cache_dir, "pathways", f"{pathway_id}.txt")
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)

        if os.path.exists(cache_file):
            with open(cache_file, 'r', encoding='utf-8') as f:
                flat_text = f.read()
        else:
            flat_text = kegg_get(pathway_id).read()
            with open(cache_file, 'w', encoding='utf-8') as f:
                f.write(flat_text)

        gene_kos = self._parse_gene_kos_from_text(flat_text)
        pw = Pathway(pathway_id, gene_kos, organism=self)
        self.pathways[pathway_id] = pw

        # ---------- Parse additional metadata ----------
        in_desc = False
        in_dblinks = False
        desc_lines = []
        for line in flat_text.splitlines():
            if line.startswith("NAME") and pw.title is None:
                title = line[5:].strip()
                pw.title = title
            if line.startswith("DESCRIPTION"):
                in_desc = True
                rest = line[12:].strip()
                if rest:
                    desc_lines.append(rest)
                continue
            if in_desc:
                if line.startswith(" "):
                    desc_lines.append(line.strip())
                else:
                    in_desc = False
            if line.startswith("DBLINKS"):
                in_dblinks = True
                rest = line[8:].strip()
                if rest:
                    db_parts = rest.split(": ", 1)
                    if len(db_parts) == 2:
                        pw.dblinks[db_parts[0]] = db_parts[1]
                continue
            if in_dblinks:
                if line.startswith(" "):
                    rest = line.strip()
                    db_parts = rest.split(": ", 1)
                    if len(db_parts) == 2:
                        pw.dblinks[db_parts[0]] = db_parts[1]
                else:
                    in_dblinks = False

        if desc_lines:
            pw.description = " ".join(desc_lines)
        # ------------------------------------------------

        self._parse_compound_names_from_text(flat_text)

        all_kos = set().union(*gene_kos.values()) if gene_kos else set()
        for ko in all_kos:
            rxn_ids = self._ko_to_reactions.get(ko, frozenset())
            for rn_id in rxn_ids:
                pw.reaction_ids.add(rn_id)
                if rn_id not in self.reactions:
                    self.reactions[rn_id] = Reaction(rn_id, self)

        try:
            kgml = get_pathway_kgml(pathway_id, self.cache_dir)
        except Exception as e:
            logger.warning(f"Could not fetch/parse KGML for {pathway_id}: {e}")
            kgml = None

        if kgml is not None:
            pw.title = kgml.title

            for kgml_rxn in kgml.reactions:
                rxn_id = kgml_rxn.name.split(':')[-1]

                if rxn_id not in self.reactions:
                    self.reactions[rxn_id] = Reaction(rxn_id, self)
                rxn_obj = self.reactions[rxn_id]

                # --- Store bare compound IDs (strip 'cpd:' prefix) ---
                substrates_kegg = [s.name.split(':')[-1] for s in kgml_rxn.substrates]
                products_kegg   = [p.name.split(':')[-1] for p in kgml_rxn.products]

                # Update reaction's compound sets
                rxn_obj._substrate_ids.update(substrates_kegg)
                rxn_obj._product_ids.update(products_kegg)

                subs_short = substrates_kegg   # already bare
                prod_short = products_kegg
                arrow = ' <=> ' if kgml_rxn.type == 'reversible' else ' --> '

                formula_kegg = ' + '.join(subs_short) + arrow + ' + '.join(prod_short) if (subs_short or prod_short) else ''

                substrates_read = []
                for cid in substrates_kegg:
                    comp = self.get_compound(cid, fetch_if_missing=False)
                    substrates_read.append(comp.name if comp.name else cid)

                products_read = []
                for cid in products_kegg:
                    comp = self.get_compound(cid, fetch_if_missing=False)
                    products_read.append(comp.name if comp.name else cid)

                formula_read = ' + '.join(substrates_read) + arrow + ' + '.join(products_read)

                rxn_obj.formula_per_pathway[pathway_id] = {
                    'type': kgml_rxn.type,
                    'substrates': substrates_kegg,        # bare IDs
                    'products': products_kegg,
                    'substrates_read': substrates_read,
                    'products_read': products_read,
                    'formula_kegg': formula_kegg,
                    'formula_read': formula_read
                }

                pw.reaction_ids.add(rxn_id)

        # Freeze reaction_ids
        pw.reaction_ids = frozenset(pw.reaction_ids) if not isinstance(pw.reaction_ids, frozenset) else pw.reaction_ids

        for rn_id in pw.reaction_ids:
            self._pathway_reaction_map[rn_id].add(pathway_id)
        for locus in pw.gene_ids:
            self._gene_pathway_map[locus].add(pathway_id)

        return pw

    # ------------------------------------------------------------------
    # Load all cached pathways (unchanged except now calls finalize)
    # ------------------------------------------------------------------
    def load_all_cached_pathways(self) -> list:
        pathways_dir = os.path.join(self.cache_dir, "pathways")
        if not os.path.isdir(pathways_dir):
            logger.warning(f"No pathways directory found at {pathways_dir}")
            return []

        files = [f for f in os.listdir(pathways_dir)
                 if f.startswith(self.org_code) and f.endswith('.txt')]
        if tqdm:
            files = tqdm(files, desc="Loading pathways", unit="pathway")

        loaded = []
        for fname in files:
            pw_id = fname[:-4]
            pw = self.load_pathway(pw_id)
            loaded.append(pw)

        logger.info(f"Loaded {len(loaded)} cached pathways for {self.org_code}")
        self.finalize()
        return loaded

    # ------------------------------------------------------------------
    # Finalize – builds all cross‑references
    # ------------------------------------------------------------------
    def finalize(self) -> None:
        if self._is_finalized:
            return

        logger.info("Building compound–reaction–gene–pathway links...")

        # 1. Freeze reaction substrate/product sets
        for rxn in self.reactions.values():
            rxn._substrate_ids = frozenset(rxn._substrate_ids)
            rxn._product_ids = frozenset(rxn._product_ids)

        # 2. Build compound → reactions index
        compound_to_reactions = defaultdict(set)
        for rxn_id, rxn in self.reactions.items():
            for cid in rxn.substrates | rxn.products:
                compound_to_reactions[cid].add(rxn_id)

        # 3. Assign reaction_ids, pathway_ids, gene_ids to each compound
        for cid, comp in self._compounds.items():
            rxn_set = compound_to_reactions.get(cid, set())
            comp._reaction_ids = frozenset(rxn_set)
            # pathway_ids = union of all pathways for those reactions
            pw_set = set()
            gene_set = set()
            for rxn_id in rxn_set:
                pw_set.update(self._pathway_reaction_map.get(rxn_id, set()))
                gene_set.update(self.reactions[rxn_id].genes)
            comp._pathway_ids = frozenset(pw_set)
            comp._gene_ids = frozenset(gene_set)

        # 4. Assign pathway_ids to each reaction
        for rxn_id, rxn in self.reactions.items():
            rxn._pathway_ids = frozenset(self._pathway_reaction_map.get(rxn_id, set()))

        # 5. Assign compound_ids to each pathway
        for pw in self.pathways.values():
            cids = set()
            for rxn_id in pw.reaction_ids:
                rxn = self.reactions.get(rxn_id)
                if rxn:
                    cids.update(rxn.substrates | rxn.products)
            pw._compound_ids = frozenset(cids)

        self._is_finalized = True
        logger.debug("Linkage finalization complete.")

    # ------------------------------------------------------------------
    # Save / Load (unchanged)
    # ------------------------------------------------------------------
    def save(self, filepath: str) -> None:
        data = {
            'version': __version__,
            'organism': self
        }
        with open(filepath, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        logger.info(f"Organism saved to {filepath} (keggtangled v{__version__})")

    @staticmethod
    def load(filepath: str) -> 'Organism':
        with open(filepath, 'rb') as f:
            data = pickle.load(f)

        saved_version = data.get('version', '0.0.0')
        if saved_version != __version__:
            logger.warning(
                f"Loading organism saved with keggtangled v{saved_version} "
                f"(current is v{__version__}). Compatibility not guaranteed."
            )

        obj = data['organism']
        if not isinstance(obj, Organism):
            raise TypeError(f"File {filepath} does not contain a keggtangled Organism")
        # After unpickling, ensure linkages are frozen (they should be, but just in case)
        if not obj._is_finalized:
            obj.finalize()
        logger.info(f"Organism loaded from {filepath}")
        return obj

    # ------------------------------------------------------------------
    # Summary (unchanged)
    # ------------------------------------------------------------------
    def summary(self) -> dict:
        return {
            "organism": self.org_code,
            "pathways": len(self.pathways),
            "reactions": len(self.reactions),
            "compounds": len(self._compounds),
            "genes": len(self._gene_to_kos),
            "KOs": len(self._ko_to_genes),
            "finalized": self._is_finalized
        }

    def __str__(self) -> str:
        return (
            f"Organism({self.org_code}) – "
            f"{len(self.pathways)} pathways, {len(self.reactions)} reactions, "
            f"{len(self._compounds)} compounds, {len(self._gene_to_kos)} genes"
        )

    # ------------------------------------------------------------------
    # Convenience iterators
    # ------------------------------------------------------------------
    @property
    def compounds(self) -> Iterator[Compound]:
        return iter(self._compounds.values())

    @property
    def all_reactions(self) -> Iterator[Reaction]:
        return iter(self.reactions.values())

    @property
    def all_pathways(self) -> Iterator[Pathway]:
        return iter(self.pathways.values())

    # ------------------------------------------------------------------
    # Convenience methods
    # ------------------------------------------------------------------
    def get_pathway(self, pathway_id: str) -> Pathway:
        return self.load_pathway(pathway_id)

    def get_reaction(self, reaction_id: str) -> Reaction:
        if reaction_id not in self.reactions:
            self.reactions[reaction_id] = Reaction(reaction_id, self)
        return self.reactions[reaction_id]

    def get_pathways_for_reaction(self, reaction_id: str) -> FrozenSet[str]:
        return frozenset(self._pathway_reaction_map.get(reaction_id, set()))

    def get_pathways_for_gene(self, locus_tag: str) -> FrozenSet[str]:
        return frozenset(self._gene_pathway_map.get(locus_tag, set()))

    def get_genes_for_pathway(self, pathway_id: str) -> FrozenSet[str]:
        pw = self.pathways.get(pathway_id)
        return pw.gene_ids if pw else frozenset()

    def get_reactions_for_pathway(self, pathway_id: str) -> FrozenSet[str]:
        pw = self.pathways.get(pathway_id)
        return pw.reaction_ids if pw else frozenset()