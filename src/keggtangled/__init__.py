#!/usr/bin/env python
# coding: utf-8
__version__ = "0.1.0"

import re
import os
import json
import io
import hashlib
from Bio.KEGG.REST import kegg_link, kegg_get
from Bio.KEGG.KGML.KGML_parser import read as kgml_read

# ----------------------------------------------------------------------
# Pre‑compiled regular expressions
# ----------------------------------------------------------------------
KO_RE = re.compile(r'\[KO:(K\d+)\]')

# ----------------------------------------------------------------------
# Compound class
# ----------------------------------------------------------------------
class Compound:
    """
    Represents a KEGG compound (e.g., C00022).
    """
    def __init__(self, compound_id, organism):
        self.id = compound_id           # e.g., 'C00022'
        self.organism = organism        # Organism instance (needed for KOs / genes)
        self.name = None                # human‑readable name
        self.formula = None             # molecular formula
        self.mass = None                # molecular weight (string from KEGG)
        self.reactions = set()          # reaction IDs this compound participates in

    def __repr__(self):
        return f"Compound({self.id}, {self.name})"

    def get_kos(self):
        """Return the set of KOs linked to this compound via any reaction."""
        kos = set()
        for rxn_id in self.reactions:
            kos.update(self.organism._reaction_to_kos.get(rxn_id, set()))
        return kos

    def get_genes(self):
        """Return the set of genes linked to this compound (through KOs and reactions)."""
        genes = set()
        for ko in self.get_kos():
            genes.update(self.organism.get_genes_for_ko(ko))
        return genes


# ----------------------------------------------------------------------
# Reaction class
# ----------------------------------------------------------------------
class Reaction:
    def __init__(self, reaction_id, organism):
        self.reaction_id = reaction_id
        self.organism = organism
        self.ko_to_genes = {}

        # pathway_id -> {type, substrates, products, substrates_read, products_read,
        #                formula_kegg, formula_read}
        self.formula_per_pathway = {}

        kos = organism._reaction_to_kos.get(reaction_id, set())
        for ko in kos:
            genes = organism.get_genes_for_ko(ko)
            if genes:
                self.ko_to_genes[ko] = genes

    def __repr__(self):
        return (f"Reaction({self.reaction_id}, {self.organism.org_code}) – "
                f"{len(self.ko_to_genes)} KOs mapped to genes, "
                f"{len(self.formula_per_pathway)} pathway formulas")

    def get_genes(self):
        all_genes = set()
        for genes in self.ko_to_genes.values():
            all_genes.update(genes)
        return all_genes

    def get_kos(self):
        return list(self.ko_to_genes.keys())


# ----------------------------------------------------------------------
# Pathway class
# ----------------------------------------------------------------------
class Pathway:
    def __init__(self, pathway_id, gene_kos=None):
        self.id = pathway_id
        if gene_kos is None:
            gene_kos = {}
        self.gene_ids = set(gene_kos.keys())
        self.reaction_ids = set()

    def add_reactions(self, reaction_ids):
        self.reaction_ids.update(reaction_ids)

    def __repr__(self):
        return f"Pathway({self.id}, {len(self.gene_ids)} genes, {len(self.reaction_ids)} reactions)"


# ----------------------------------------------------------------------
# KGML fetcher (cached)
# ----------------------------------------------------------------------
def get_pathway_kgml(pathway_id, cache_dir="kegg_cache"):
    """
    Return a parsed KGML pathway object (Bio.KEGG.KGML.Pathway).
    The raw XML is cached on disk.
    """
    os.makedirs(cache_dir, exist_ok=True)
    full_id = f"path:{pathway_id}" if not pathway_id.startswith("path:") else pathway_id
    cache_file = os.path.join(cache_dir, f"{pathway_id}.kgml")

    if not os.path.exists(cache_file):
        raw_kgml = kegg_get(full_id, "kgml").read()
        with open(cache_file, 'w', encoding='utf-8') as f:
            f.write(raw_kgml)
    else:
        with open(cache_file, 'r', encoding='utf-8') as f:
            raw_kgml = f.read()

    return kgml_read(io.StringIO(raw_kgml))


# ----------------------------------------------------------------------
# Organism class (fully integrated with Compound, Reaction, Pathway)
# ----------------------------------------------------------------------
class Organism:
    def __init__(self, org_code, batch_size=10, cache_dir="kegg_cache"):
        self.org_code = org_code
        self.batch_size = batch_size
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        self.pathways = {}
        self.reactions = {}
        self._pathway_reaction_map = {}
        self._gene_pathway_map = {}
        self._ko_to_reactions = {}
        self._ko_to_genes = {}
        self._reaction_to_kos = {}
        self._gene_to_kos = {}                     # reverse mapping: locus_tag -> set of KOs

        # Compound objects dictionary
        self._compounds = {}                       # "C00022" -> Compound instance
        self._compound_cache_file = os.path.join(cache_dir, f"{org_code}_compounds.json")

        # Bulk pre‑fetching
        self._load_all_ko_genes()                  # builds _ko_to_genes and _gene_to_kos
        self._prefetch_all_ko_reactions()          # builds _ko_to_reactions and _reaction_to_kos
        self._load_compounds()

    # ------------------------------------------------------------------
    # Caching helper (unchanged)
    # ------------------------------------------------------------------
    def _cache_get(self, key, subdir="", fetcher_func=None, *args, **kwargs):
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
    # KO–gene mapping (JSON) + reverse gene→KO mapping
    # ------------------------------------------------------------------
    def _load_all_ko_genes(self):
        json_file = os.path.join(self.cache_dir, f"{self.org_code}_ko_genes.json")
        if os.path.exists(json_file):
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
                self._ko_to_genes = {ko: set(genes) for ko, genes in data.items()}
        else:
            raw = kegg_link("ko", self.org_code).read().strip()
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
                    self._ko_to_genes.setdefault(ko_part, set()).add(locus)

                json_data = {ko: list(genes) for ko, genes in self._ko_to_genes.items()}
                with open(json_file, 'w', encoding='utf-8') as f:
                    json.dump(json_data, f, indent=2)

        # Build reverse mapping (gene -> KOs)
        self._gene_to_kos = {}
        for ko, genes in self._ko_to_genes.items():
            for gene in genes:
                self._gene_to_kos.setdefault(gene, set()).add(ko)

    # ------------------------------------------------------------------
    # KO–reaction mapping (JSON)
    # ------------------------------------------------------------------
    def _prefetch_all_ko_reactions(self):
        all_kos = list(self._ko_to_genes.keys())
        if not all_kos:
            return

        ko_reactions_file = os.path.join(self.cache_dir, f"{self.org_code}_ko_reactions.json")
        reaction_kos_file  = os.path.join(self.cache_dir, f"{self.org_code}_reaction_kos.json")

        if os.path.exists(ko_reactions_file) and os.path.exists(reaction_kos_file):
            with open(ko_reactions_file, 'r') as f:
                data = json.load(f)
                self._ko_to_reactions = {ko: set(rxn_list) for ko, rxn_list in data.items()}
            with open(reaction_kos_file, 'r') as f:
                data = json.load(f)
                self._reaction_to_kos = {rn: set(ko_list) for rn, ko_list in data.items()}
            return

        for i in range(0, len(all_kos), self.batch_size):
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
                    self._get_reactions_for_ko_fallback(ko)
                continue

            for line in raw.splitlines():
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                ko_id = parts[0]
                rn_id = parts[1].split(":")[1]
                self._ko_to_reactions.setdefault(ko_id, set()).add(rn_id)

        self._reaction_to_kos = {}
        for ko, rxn_set in self._ko_to_reactions.items():
            for rn in rxn_set:
                self._reaction_to_kos.setdefault(rn, set()).add(ko)

        with open(ko_reactions_file, 'w') as f:
            json.dump({ko: list(rxn_set) for ko, rxn_set in self._ko_to_reactions.items()}, f, indent=2)
        with open(reaction_kos_file, 'w') as f:
            json.dump({rn: list(ko_set) for rn, ko_set in self._reaction_to_kos.items()}, f, indent=2)

    def _get_reactions_for_ko_fallback(self, ko):
        raw = kegg_link("rn", ko).read().strip()
        rxn_set = set()
        for line in raw.splitlines():
            if line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    rn_id = parts[1].split(":")[1]
                    rxn_set.add(rn_id)
        self._ko_to_reactions[ko] = rxn_set
        for rn in rxn_set:
            self._reaction_to_kos.setdefault(rn, set()).add(ko)

    # ------------------------------------------------------------------
    # Compound management
    # ------------------------------------------------------------------
    def _load_compounds(self):
        """Load cached compound info (name, formula, mass) from JSON."""
        if os.path.exists(self._compound_cache_file):
            with open(self._compound_cache_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            for cid, info in data.items():
                comp = Compound(cid, self)
                comp.name = info.get("name")
                comp.formula = info.get("formula")
                comp.mass = info.get("mass")
                self._compounds[cid] = comp
        # (reaction sets will be populated when pathways are loaded)

    def _save_compounds(self):
        """Save compound names, formulas, masses to JSON (reaction sets are transient)."""
        data = {}
        for cid, comp in self._compounds.items():
            if comp.name:   # only save if at least a name is known
                data[cid] = {"name": comp.name, "formula": comp.formula, "mass": comp.mass}
        with open(self._compound_cache_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

    def get_compound(self, compound_id, fetch_if_missing=True):
        """
        Return a Compound object. If not yet loaded, create it.
        Optionally fetch full details from KEGG (name, formula, mass) via flat file.
        """
        if compound_id not in self._compounds:
            self._compounds[compound_id] = Compound(compound_id, self)
            if fetch_if_missing:
                self._fetch_compound_details(compound_id)
        return self._compounds[compound_id]

    def _fetch_compound_details(self, compound_id):
        """Fetch compound flat file from KEGG and populate name, formula, mass."""
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
            print(f"Warning: could not fetch details for compound {compound_id}: {e}")

    def _parse_compound_names_from_text(self, flat_text):
        """
        Extract compound IDs and names from the COMPOUND section of a pathway flat file.
        Updates Compound objects in self._compounds and saves the cache.
        """
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
                if line.startswith(" "):   # continuation line
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
    def get_genes_for_ko(self, ko):
        return self._ko_to_genes.get(ko, set())

    def get_reactions_for_ko(self, ko):
        return self._ko_to_reactions.get(ko, set())

    def get_kos_for_gene(self, locus_tag):
        """Return the set of KOs associated with a gene locus tag."""
        return self._gene_to_kos.get(locus_tag, set())

    # ------------------------------------------------------------------
    # Pathway loading (with compound integration)
    # ------------------------------------------------------------------
    def _parse_gene_kos_from_text(self, flat_text):
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

    def load_pathway(self, pathway_id):
        """Load pathway, extract genes, reactions, and compounds; attach formulas."""
        if pathway_id in self.pathways:
            return self.pathways[pathway_id]

        # 1. Load / fetch the flat file
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

        # 2. Update compound names from the flat file and save
        self._parse_compound_names_from_text(flat_text)

        # 3. Add reactions from KO annotations (from flat file)
        all_kos = set().union(*gene_kos.values()) if gene_kos else set()
        for ko in all_kos:
            rxn_ids = self._ko_to_reactions.get(ko, set())
            for rn_id in rxn_ids:
                pw.reaction_ids.add(rn_id)
                if rn_id not in self.reactions:
                    self.reactions[rn_id] = Reaction(rn_id, self)

        # 4. Parse KGML and attach per‑pathway formulas, link compounds to reactions
        try:
            kgml = get_pathway_kgml(pathway_id, self.cache_dir)
        except Exception as e:
            print(f"Warning: could not fetch/parse KGML for {pathway_id}: {e}")
            kgml = None

        if kgml is not None:
            for kgml_rxn in kgml.reactions:
                rxn_id = kgml_rxn.name.split(':')[-1]

                if rxn_id not in self.reactions:
                    self.reactions[rxn_id] = Reaction(rxn_id, self)
                rxn_obj = self.reactions[rxn_id]

                # Substrates / products as 'cpd:C00022'
                substrates_kegg = [s.name for s in kgml_rxn.substrates]
                products_kegg   = [p.name for p in kgml_rxn.products]

                # Short IDs for formulas
                subs_short = [s.split(':')[-1] for s in substrates_kegg]
                prod_short = [p.split(':')[-1] for p in products_kegg]
                arrow = ' <=> ' if kgml_rxn.type == 'reversible' else ' --> '

                # KEGG‑ID formula
                formula_kegg = ' + '.join(subs_short) + arrow + ' + '.join(prod_short) if (subs_short or prod_short) else ''

                # Readable names (using Compound objects)
                substrates_read = []
                for s in substrates_kegg:
                    cid = s.split(':')[-1]
                    comp = self.get_compound(cid, fetch_if_missing=False)
                    substrates_read.append(comp.name if comp.name else cid)
                    # Link compound to this reaction
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

        # 5. Cross‑reference maps
        for rn_id in pw.reaction_ids:
            self._pathway_reaction_map.setdefault(rn_id, set()).add(pathway_id)
        for locus in pw.gene_ids:
            self._gene_pathway_map.setdefault(locus, set()).add(pathway_id)

        return pw

    # ------------------------------------------------------------------
    # Convenience methods
    # ------------------------------------------------------------------
    def get_pathway(self, pathway_id):
        return self.load_pathway(pathway_id)

    def get_reaction(self, reaction_id):
        if reaction_id not in self.reactions:
            self.reactions[reaction_id] = Reaction(reaction_id, self)
        return self.reactions[reaction_id]

    def get_pathways_for_reaction(self, reaction_id):
        return self._pathway_reaction_map.get(reaction_id, set())

    def get_pathways_for_gene(self, locus_tag):
        return self._gene_pathway_map.get(locus_tag, set())

    def get_genes_for_pathway(self, pathway_id):
        pw = self.pathways.get(pathway_id)
        return pw.gene_ids if pw else set()

    def get_reactions_for_pathway(self, pathway_id):
        pw = self.pathways.get(pathway_id)
        return pw.reaction_ids if pw else set()