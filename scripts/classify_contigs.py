#!/usr/bin/env python

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import taxopy

OUTPUT_CSV = "dataset_challenge_classification.csv"
DATASET_FASTA = "dataset_challenge.fasta"
MMSEQS2_OUTPUT = Path("mmseqs2_output")
QUERY_H = MMSEQS2_OUTPUT / "query_h"
ORFS_TAXONOMY = MMSEQS2_OUTPUT / "orfs_taxonomy.tsv"
ORFS_TAXONOMY_ALN = MMSEQS2_OUTPUT / "orfs_taxonomy_aln.tsv"


@dataclass
class ContigTaxonomicAssignment:
    taxa: list[taxopy.Taxon]
    weights: list[float]
    identities: list[float]

    # Create a `taxonomic_assignment_dict` dictionary that will store the
    # contig-level taxonomic assignments. The taxonomic assignment of each
    # contig corresponds to the majority vote of its ORFs. The support of each
    # taxon is calculated as the sum of the weights of the ORFs that are
    # assigned to it, divided by the sum of the weights of all ORFs. Read more
    # about the majority vote algorithm here:
    # https://apcamargo.github.io/taxopy/guide/#the-weights-parameter
    def assign_taxonomy(self, taxdb):
        assignment_dict = {}
        if len(self.taxa) > 1:
            aggregated_taxon = taxopy.find_majority_vote(
                self.taxa,
                taxdb,
                weights=self.weights,
            )
        else:
            # If there is only one ORF with a hit, use its taxonomic assignment.
            aggregated_taxon = self.taxa[0]
        # If identity is low and contig is assigned to species, shift it one
        # rank up
        avg_ident = sum(self.identities) / len(self.identities)
        if aggregated_taxon.rank == "species" and avg_ident < 0.8:
            aggregated_taxon = taxopy.Taxon(aggregated_taxon.taxid_lineage[1], taxdb)
        # Process each rank in the lineage
        for rank, taxid in reversed(aggregated_taxon.ranked_taxid_lineage):
            taxon = taxopy.Taxon(taxid, taxdb)
            if taxon.taxid == 1:
                continue
            support = get_support_taxon(taxon, self.taxa, self.weights)
            assignment_dict[rank] = (taxon.name, support)
        return assignment_dict


taxdb = taxopy.TaxDb(
    names_dmp="ictv_taxdump/names.dmp",
    nodes_dmp="ictv_taxdump/nodes.dmp",
    keep_files=True,
)


# Function to calculate the support of a taxon within a list of taxonomic
# lineages, each with a corresponding weight.
def get_support_taxon(taxon, taxon_list, weight_list):
    support = 0
    for t, w in zip(taxon_list, weight_list):
        if t.rank_taxid_dictionary.get(taxon.rank) == taxon.taxid:
            support += w / sum(weight_list)
    return support


# Create a dictionary that maps MMseqs2 interval ids to contig headers
header_dict = {}
with open(QUERY_H) as fi:
    c = 0
    for line in fi:
        header_dict[c] = line.strip("\n\x00")
        c += 1

# Store the ORF-level taxonomic assignments
orf_taxon_dict = {}
with open(ORFS_TAXONOMY) as fi:
    for qheader, taxid, rank, taxname in csv.reader(fi, delimiter="\t"):
        taxid = int(taxid)
        if taxid:
            orf_taxon_dict[qheader] = taxopy.Taxon(taxid, taxdb)

# Read the ORF alignments and store the weights and alignment identities of each
# ORF. The weight is the negative log of the E-value of the alignment to the
# best hit in the reference database.
weight_dict = {}
identity_dict = {}
with open(ORFS_TAXONOMY_ALN) as fi:
    for qheader, target, fident, bits, evalue in csv.reader(fi, delimiter="\t"):
        evalue, fident = float(evalue), float(fident)
        if qheader not in orf_taxon_dict:
            continue
        if evalue == 0:
            weight = 1000
        else:
            weight = -math.log(evalue)
        weight_dict[qheader] = weight
        identity_dict[qheader] = fident

# Create a `contig_taxonomy_dict` dictionary that will hold, for each contig,
# a dataclass containing three lists: the taxonomic assignments of its ORFs, the
# weights of these assignments, and the identities of the alignments of the ORFs
# to their best hits.
contig_taxonomy_dict = {}
for qheader in orf_taxon_dict:
    contig = header_dict[int(qheader.split("_", 1)[0])]
    if contig not in contig_taxonomy_dict:
        contig_taxonomy_dict[contig] = ContigTaxonomicAssignment([], [], [])
    contig_taxonomy_dict[contig].taxa.append(orf_taxon_dict[qheader])
    contig_taxonomy_dict[contig].weights.append(weight_dict[qheader])
    contig_taxonomy_dict[contig].identities.append(identity_dict[qheader])


# Write the output to a CSV file
ranks = [
    "realm",
    "subrealm",
    "kingdom",
    "subkingdom",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "genus",
    "subgenus",
    "species",
]
with open(DATASET_FASTA) as fi, open(OUTPUT_CSV, "w") as fo:
    fo.write(
        "SequenceID,Realm (-viria),Realm_score,Subrealm (-vira),Subrealm_score,Kingdom (-virae),Kingdom_score,Subkingdom (-virites),Subkingdom_score,Phylum (-viricota),Phylum_score,Subphylum (-viricotina),Subphylum_score,Class (-viricetes),Class_score,Subclass (-viricetidae),Subclass_score,Order (-virales),Order_score,Suborder (-virineae),Suborder_score,Family (-viridae),Family_score,Subfamily (-virinae),Subfamily_score,Genus (-virus),Genus_score,Subgenus (-virus),Subgenus_score,Species (binomial),Species_score\n"
    )
    for line in fi:
        if not line.startswith(">"):
            continue
        contig = line.strip()[1:]
        if contig not in contig_taxonomy_dict:
            row = [contig] + ["NA"] * 30
        else:
            row = [contig]
            contig_taxonomy = contig_taxonomy_dict[contig].assign_taxonomy(taxdb)
            for r in ranks:
                taxon_name, support = contig_taxonomy.get(r, ("NA", "NA"))
                if support != "NA":
                    support = f"{support:.4f}"
                row.extend([taxon_name, support])
        fo.write(",".join(row) + "\n")
