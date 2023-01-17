import ensembl_rest
import pandas as pd

from collections import defaultdict

###############################################################
################ Query parameter settings #####################
###############################################################
to_find = [
    ("EMX1", "GTGATGGGAGCCCTTCTTCT"),
    ("AGBL1", "TGTTGGCTCAAACACCAGAT"),
    ("AGBL1", "AATGAATGGCTACTCTAACC"),
    ("AGBL1", "ATCTGCGCATCTGGTAAAGC"),
    ("AGBL1", "ATGAGATGCAAGTCCACATC"),
    ("AGBL1", "ttttGTAACCTATGAAtttt"),
    ("AGBL1", "AGGAACATAGATGAAGGGTT"),
    ("AGBL1", "GTTCCAGAAAACTGttttCT"),
    ("AGBL1", "GTTTCAGACATGTAATTCAA"),
    ("AGBL1", "AATCCGAGTCTGCATTATGC"),
    ("AGBL1", "TTCATCAGGTTCTAttttCA"),
    ("AGBL1", "TAAAGAAAAACCCATAAAGG"),
    ("AGBL1", "GATGATGAAACTTAACATGC"),
    ("AGBL1", "AGAAATAAAATTCTAATTTA"),
    ("AGBL1", "TCCAACCTGGTAATATTCTC"),
    ("AGBL1", "TCCCTAGAGTAAAACAAAAC"),
]
upstream_flanking = 600
downstream_flanking = 600

upstream_context = 8  # 8bp upstream of the guide sequence
downstream_context = 5  # 8bp downstream of the guide sequence

###############################################################
###############################################################
###############################################################


if __name__ == "__main__":
    # https://www.geeksforgeeks.org/reverse-complement-of-dna-strand-using-python/
    def verify(sequence):
        """This code verifies if a sequence is a DNA or RNA"""

        # set the input sequence
        seq = set(sequence)

        # confirm if its elements is equal to
        # the set of valid DNA bases
        # Use a union method to ensure the
        # sequence is verified if does not
        # contain all the bases
        if seq == {"A", "T", "C", "G"}.union(seq):
            return "DNA"
        elif seq == {"A", "U", "C", "G"}.union(seq):
            return "RNA"
        else:
            return "Invalid sequence"

    def rev_comp_if(seq):
        comp = []
        if verify(seq) == "DNA":
            for base in seq:
                if base == ("A" or "a"):
                    comp.append("T")
                elif base == ("G" or "g"):
                    comp.append("C")
                elif base == ("T" or "t"):
                    comp.append("A")
                elif base == ("C" or "c"):
                    comp.append("G")
        elif verify(seq) == "RNA":
            for base in seq:
                if base == ("U" or "u"):
                    comp.append("A")
                elif base == ("G" or "g"):
                    comp.append("C")
                elif base == ("A" or "a"):
                    comp.append("U")
                elif base == ("C" or "c"):
                    comp.append("G")
        else:
            return "Invalid Sequence"

        # reverse the sequence
        comp_rev = comp[::-1]

        # convert list to string
        comp_rev = "".join(comp_rev)
        return comp_rev

    def constant_factory(value: list):
        if value is None:
            return {"Gene": None, "Sequence": None, "Direction": 0}
        return {
            "Gene": value[0],
            "Sequence": value[1],
            "Direction": value[2],
        }

    query_results = defaultdict(constant_factory)

    for gene, seq in to_find:
        print(gene, seq)

        response = ensembl_rest.symbol_lookup(species="homo_sapiens", symbol=gene)

        fwd = ensembl_rest.sequence_region(
            region=f"{response['seq_region_name']}:{response['start'] - upstream_flanking}..{response['end'] + downstream_flanking}:1",
            species="homo_sapiens",
            coord_system_version="GRCh38",
        )
        # ::-1 does not mean reverse complement
        # rev = ensembl_rest.sequence_region(
        #     region=f"{response['seq_region_name']}:{response['start']}..{response['end']}:-1",
        #     species="homo_sapiens",
        #     coord_system_version="GRCh38",
        # )

        fwd_seq = fwd["seq"].lower()
        rev_seq = rev_comp_if(fwd["seq"]).lower()

        # debug
        with open(f"{gene}_fwd.fa", "w") as f:
            f.write(f">{gene}\n{fwd_seq}")
        with open(f"{gene}_rev.fa", "w") as f:
            f.write(f">{gene}\n{rev_seq}")

        if seq.lower() in fwd_seq:
            query_results[seq] = constant_factory(
                [
                    gene,
                    fwd_seq[
                        fwd_seq.index(seq.lower())
                        - upstream_context : fwd_seq.index(seq.lower())
                        + len(seq)
                        + downstream_context
                    ],
                    1,
                ]
            )

        elif seq.lower() in rev_seq:
            query_results[seq] = constant_factory(
                [
                    gene,
                    rev_seq[
                        rev_seq.index(seq.lower())
                        - upstream_context : rev_seq.index(seq.lower())
                        + len(seq)
                        + downstream_context
                    ],
                    -1,
                ]
            )
        else:
            query_results[seq] = {"gene": gene, "sequence": None, "Direction": 0}

        # print(f"fwd: {fwd}")
        # print(f"rev: {rev}")

        pd.DataFrame.from_dict(query_results, orient="index").to_csv(
            "query_results.csv"
        )
