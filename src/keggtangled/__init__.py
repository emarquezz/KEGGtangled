#!/usr/bin/env python
# coding: utf-8
__version__ = "0.4.0"

import re
import os
import json
import io
import hashlib
import logging
import pickle
from collections import defaultdict
from typing import Dict, Set, FrozenSet, Tuple, Optional, Iterator, Union

from Bio.KEGG.REST import kegg_link, kegg_get
from Bio.KEGG.KGML.KGML_parser import read as kgml_read

# Optional progress bar
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# ----------------------------------------------------------------------
# Colored logging (ANSI, no external dependency)
# ----------------------------------------------------------------------
class ColoredFormatter(logging.Formatter):
    COLORS = {
        'DEBUG': '\033[36m',     # cyan
        'INFO': '\033[32m',      # green
        'WARNING': '\033[33m',   # yellow
        'ERROR': '\033[31m',     # red
        'CRITICAL': '\033[1;31m' # bold red
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
logger.setLevel(logging.INFO)   # set to DEBUG for more details

# ----------------------------------------------------------------------
# Pre‑compiled regular expressions
# ----------------------------------------------------------------------
KO_RE = re.compile(r'\[KO:(K\d+)\]')

# ----------------------------------------------------------------------
# Compound class (optimized with __slots__ and manual caching)
# ----------------------------------------------------------------------
class Compound:
    """Represents a KEGG compound (e.g., C00022). Memory‑efficient with __slots__."""

    __slots__ = ('id', 'organism', 'name', 'formula', 'mass', 'reactions',
                 '_kos', '_genes')

    def __init__(self, compound_id: str, organism: 'Organism') -> None:
        self.id: str = compound_id
        self.organism: 'Organism' = organism
        self.name: Optional[str] = None
        self.formula: Optional[str] = None
        self.mass: Optional[str] = None
        self.reactions: Union[Set[str], FrozenSet[str]] = set()   # mutable until finalize()
        self._kos: Optional[FrozenSet[str]] = None
        self._genes: Optional[FrozenSet[str]] = None

    def __repr__(self) -> str:
        return f"Compound({self.id}, {self.name})"

    def get_kos(self) -> FrozenSet[str]:
        """Return KOs linked to this compound via any reaction (cached)."""
        if self._kos is None:
            kos = set()
            for rxn_id in self.reactions:
                kos.update(self.organism._reaction_to_kos.get(rxn_id, frozenset()))
            self._kos = frozenset(kos)
        return self._kos

    def get_genes(self) -> FrozenSet[str]:
        """Return genes linked to this compound through KOs and reactions (cached)."""
        if self._genes is None:
            genes = set()
            for ko in self.get_kos():
                genes.update(self.organism.get_genes_for_ko(ko))
            self._genes = frozenset(genes)
        return self._genes


# ----------------------------------------------------------------------
# Reaction class (optimized with __slots__ and manual caching)
# ----------------------------------------------------------------------
class Reaction:
    """Represents a single KEGG reaction with per‑pathway formulas."""

    __slots__ = ('reaction_id', 'organism', 'ko_to_genes',
                 'formula_per_pathway', '_kos_cache', '_genes_cache')

    def __init__(self, reaction_id: str, organism: 'Organism') -> None:
        self.reaction_id: str = reaction_id
        self.organism: 'Organism' = organism
        self.ko_to_genes: Dict[str, Set[str]] = {}
        self.formula_per_pathway: Dict[str, dict] = {}
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

    def get_genes(self) -> FrozenSet[str]:
        """Return all genes associated with this reaction (cached)."""
        if self._genes_cache is None:
            all_genes = set()
            for genes in self.ko_to_genes.values():
                all_genes.update(genes)
            self._genes_cache = frozenset(all_genes)
        return self._genes_cache

    def get_kos(self) -> Tuple[str, ...]:
        """Return KOs associated with this reaction (cached, immutable)."""
        if self._kos_cache is None:
            self._kos_cache = tuple(self.ko_to_genes.keys())
        return self._kos_cache


# ----------------------------------------------------------------------
# Pathway class (optimized with __slots__, gene_ids is a frozenset)
# ----------------------------------------------------------------------
class Pathway:
    """Represents a KEGG pathway with immutable gene_ids and reaction_ids."""

    __slots__ = ('id', 'gene_ids', 'reaction_ids')

    def __init__(self, pathway_id: str, gene_kos: Optional[Dict[str, Set[str]]] = None) -> None:
        self.id: str = pathway_id
        if gene_kos is None:
            gene_kos = {}
        self.gene_ids: FrozenSet[str] = frozenset(gene_kos.keys())
        self.reaction_ids: Union[Set[str], FrozenSet[str]] = set()  # frozen after loading

    def add_reactions(self, reaction_ids: Union[str, Set[str]]) -> None:
        """Add reaction IDs to the pathway (used during loading)."""
        if isinstance(reaction_ids, str):
            self.reaction_ids.add(reaction_ids)
        else:
            self.reaction_ids.update(reaction_ids)

    def __repr__(self) -> str:
        return f"Pathway({self.id}, {len(self.gene_ids)} genes, {len(self.reaction_ids)} reactions)"


# ----------------------------------------------------------------------
# KGML fetcher (cached, with pickle of parsed object)
# ----------------------------------------------------------------------
def get_pathway_kgml(pathway_id: str, cache_dir: str = "kegg_cache") -> 'Pathway':
    """Return a parsed KGML pathway, caching both raw XML and parsed object."""
    os.makedirs(cache_dir, exist_ok=True)
    full_id = f"path:{pathway_id}" if not pathway_id.startswith("path:") else pathway_id
    raw_cache_file = os.path.join(cache_dir, f"{pathway_id}.kgml")
    pkl_cache_file = os.path.join(cache_dir, f"{pathway_id}.kgml.pkl")

    # 1. Try parsed pickle cache first
    if os.path.exists(pkl_cache_file):
        try:
            with open(pkl_cache_file, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            logger.warning(f"Could not load cached KGML pickle for {pathway_id}, re‑parsing: {e}")

    # 2. Raw XML cache or fetch
    if not os.path.exists(raw_cache_file):
        raw_kgml = kegg_get(full_id, "kgml").read()
        with open(raw_cache_file, 'w', encoding='utf-8') as f:
            f.write(raw_kgml)
    else:
        with open(raw_cache_file, 'r', encoding='utf-8') as f:
            raw_kgml = f.read()

    kgml = kgml_read(io.StringIO(raw_kgml))

    # 3. Save parsed pickle for next time
    try:
        with open(pkl_cache_file, 'wb') as f:
            pickle.dump(kgml, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        logger.warning(f"Could not save parsed KGML pickle: {e}")

    return kgml


# ----------------------------------------------------------------------
# Organism class (fully featured)
# ----------------------------------------------------------------------
class Organism:
    """Central class for a KEGG organism – holds all pathways, reactions, compounds."""

    def __init__(self, org_code: str, batch_size: int = 10, cache_dir: str = "kegg_cache"):
        self.org_code: str = org_code
        self.batch_size: int = batch_size
        self.cache_dir: str = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        self.pathways: Dict[str, Pathway] = {}
        self.reactions: Dict[str, Reaction] = {}
        self._pathway_reaction_map: Dict[str, Set[str]] = defaultdict(set)
        self._gene_pathway_map: Dict[str, Set[str]] = defaultdict(set)

        # Core mappings (frozen after loading)
        self._ko_to_reactions: Dict[str, FrozenSet[str]] = {}
        self._ko_to_genes: Dict[str, FrozenSet[str]] = {}
        self._reaction_to_kos: Dict[str, FrozenSet[str]] = {}
        self._gene_to_kos: Dict[str, FrozenSet[str]] = {}

        self._compounds: Dict[str, Compound] = {}
        self._compound_cache_file = os.path.join(cache_dir, f"{org_code}_compounds.json")

        self._is_finalized: bool = False

        self._load_all_ko_genes()
        self._prefetch_all_ko_reactions()
        self._load_compounds()
        self._freeze_mappings()

    # ------------------------------------------------------------------
    # Caching helper
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
    # KO–gene mapping
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
    # KO–reaction mapping (with progress bar)
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
    # Internal freezing of mappings
    # ------------------------------------------------------------------
    def _freeze_mappings(self) -> None:
        """Convert all internal mapping sets to frozensets (call once after init)."""
        self._ko_to_genes = {k: frozenset(v) for k, v in self._ko_to_genes.items()}
        self._gene_to_kos = {k: frozenset(v) for k, v in self._gene_to_kos.items()}
        self._ko_to_reactions = {k: frozenset(v) for k, v in self._ko_to_reactions.items()}
        self._reaction_to_kos = {k: frozenset(v) for k, v in self._reaction_to_kos.items()}

    # ------------------------------------------------------------------
    # Compound management
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
                self._compounds[cid] = comp

    def _save_compounds(self) -> None:
        data = {}
        for cid, comp in self._compounds.items():
            if comp.name:
                data[cid] = {"name": comp.name, "formula": comp.formula, "mass": comp.mass}
        with open(self._compound_cache_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

    def get_compound(self, compound_id: str, fetch_if_missing: bool = True) -> Compound:
        """Return Compound object, creating and fetching if needed."""
        if compound_id not in self._compounds:
            self._compounds[compound_id] = Compound(compound_id, self)
            if fetch_if_missing:
                self._fetch_compound_details(compound_id)
        return self._compounds[compound_id]

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
    # Public lookups (return immutable views)
    # ------------------------------------------------------------------
    def get_genes_for_ko(self, ko: str) -> FrozenSet[str]:
        return self._ko_to_genes.get(ko, frozenset())

    def get_reactions_for_ko(self, ko: str) -> FrozenSet[str]:
        return self._ko_to_reactions.get(ko, frozenset())

    def get_kos_for_gene(self, locus_tag: str) -> FrozenSet[str]:
        return self._gene_to_kos.get(locus_tag, frozenset())

    # ------------------------------------------------------------------
    # Gene ↔ reaction ↔ compound walkers
    # ------------------------------------------------------------------
    def get_reactions_for_gene(self, locus_tag: str) -> FrozenSet[str]:
        """All reactions linked to a gene (via KOs)."""
        reactions = set()
        for ko in self.get_kos_for_gene(locus_tag):
            reactions.update(self.get_reactions_for_ko(ko))
        return frozenset(reactions)

    def get_compounds_for_reaction(self, reaction_id: str) -> FrozenSet[str]:
        """Compound IDs participating in a given reaction."""
        rxn = self.reactions.get(reaction_id)
        if not rxn:
            return frozenset()
        compounds = set()
        for pw_data in rxn.formula_per_pathway.values():
            for s in pw_data.get('substrates', []):
                compounds.add(s.split(':')[-1])
            for p in pw_data.get('products', []):
                compounds.add(p.split(':')[-1])
        return frozenset(compounds)

    def get_genes_for_reaction(self, reaction_id: str) -> FrozenSet[str]:
        """Genes linked to a reaction (via KOs)."""
        rxn = self.reactions.get(reaction_id)
        if rxn:
            return rxn.get_genes()
        return frozenset()

    # ------------------------------------------------------------------
    # Pathway loading
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
        """Load a single pathway, including compounds and KGML formulas."""
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
        pw = Pathway(pathway_id, gene_kos)
        self.pathways[pathway_id] = pw

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
            for kgml_rxn in kgml.reactions:
                rxn_id = kgml_rxn.name.split(':')[-1]

                if rxn_id not in self.reactions:
                    self.reactions[rxn_id] = Reaction(rxn_id, self)
                rxn_obj = self.reactions[rxn_id]

                substrates_kegg = [s.name for s in kgml_rxn.substrates]
                products_kegg   = [p.name for p in kgml_rxn.products]

                subs_short = [s.split(':')[-1] for s in substrates_kegg]
                prod_short = [p.split(':')[-1] for p in products_kegg]
                arrow = ' <=> ' if kgml_rxn.type == 'reversible' else ' --> '

                formula_kegg = ' + '.join(subs_short) + arrow + ' + '.join(prod_short) if (subs_short or prod_short) else ''

                substrates_read = []
                for s in substrates_kegg:
                    cid = s.split(':')[-1]
                    comp = self.get_compound(cid, fetch_if_missing=False)
                    substrates_read.append(comp.name if comp.name else cid)
                    comp.reactions.add(rxn_id)

                products_read = []
                for p in products_kegg:
                    cid = p.split(':')[-1]
                    comp = self.get_compound(cid, fetch_if_missing=False)
                    products_read.append(comp.name if comp.name else cid)
                    comp.reactions.add(rxn_id)

                formula_read = ' + '.join(substrates_read) + arrow + ' + '.join(products_read)

                rxn_obj.formula_per_pathway[pathway_id] = {
                    'type': kgml_rxn.type,
                    'substrates': substrates_kegg,
                    'products': products_kegg,
                    'substrates_read': substrates_read,
                    'products_read': products_read,
                    'formula_kegg': formula_kegg,
                    'formula_read': formula_read
                }

                pw.reaction_ids.add(rxn_id)

        pw.reaction_ids = frozenset(pw.reaction_ids)

        for rn_id in pw.reaction_ids:
            self._pathway_reaction_map[rn_id].add(pathway_id)
        for locus in pw.gene_ids:
            self._gene_pathway_map[locus].add(pathway_id)

        return pw

    # ------------------------------------------------------------------
    # Load all cached pathways with progress bar
    # ------------------------------------------------------------------
    def load_all_cached_pathways(self) -> list:
        """Load every cached pathway, showing a progress bar if tqdm is available."""
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
        self.finalize()   # automatically freeze compounds after bulk load
        return loaded

    # ------------------------------------------------------------------
    # Finalize: freeze mutable sets
    # ------------------------------------------------------------------
    def finalize(self) -> None:
        """Freeze compound reactions and any remaining mutable structures."""
        if self._is_finalized:
            return
        for comp in self._compounds.values():
            if isinstance(comp.reactions, set):
                comp.reactions = frozenset(comp.reactions)
        self._is_finalized = True
        logger.debug("Organism finalized – compound reactions frozen.")

    # ------------------------------------------------------------------
    # Save / Load with version check
    # ------------------------------------------------------------------
    def save(self, filepath: str) -> None:
        """Serialize the whole Organism to a pickle file, preserving cached state."""
        data = {
            'version': __version__,
            'organism': self
        }
        with open(filepath, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        logger.info(f"Organism saved to {filepath} (keggtangled v{__version__})")

    @staticmethod
    def load(filepath: str) -> 'Organism':
        """Load a previously saved Organism from a pickle file."""
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
        logger.info(f"Organism loaded from {filepath}")
        return obj

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    def summary(self) -> dict:
        """Return a dictionary with key statistics about the organism."""
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
        """Iterator over all Compound objects."""
        return iter(self._compounds.values())

    @property
    def all_reactions(self) -> Iterator[Reaction]:
        """Iterator over all Reaction objects."""
        return iter(self.reactions.values())

    @property
    def all_pathways(self) -> Iterator[Pathway]:
        """Iterator over all Pathway objects."""
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