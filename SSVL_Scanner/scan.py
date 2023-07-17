import csv
from collections import defaultdict
import re

import pandas as pd

from aa_table import codonlst
from aa_table import codontab
from dna_table import dna_nt
from dna_table import complement_nt
from restriction_enzyme import enzymetab

# SOME GLOBAL VARIABLES
UPSTREAM_FLANKING_SEQUENCE = (
    "ggatcctctggcagcgagacaccaggaacaagcgagtcagcaacaccagagagcagtggcggcagcagcggcggcagcagc"
)
DOWNSTREAM_FLANKING_SEQUENCE = "tctggcggctcaaaaagaaccgcc"
SSVL_TARGET = "./SFV_RT.txt"

INIT_POS = 0
CODON_LENGTH = 3


def get_rev_strand(dna: str):
    return "".join(
        map(str, [complement_nt[dna[i]] for i in range(len(dna) - 1, -1, -1)])
    )


def translate_seq(dna, offset=0):
    """
    No need to translate the sequence at this time (20220830)



    :param dna:
    :param offset: amino acid offset, multiples of 3
    :return:
    """
    """https://python.plainenglish.io/bioinformatics-in-python-dna-toolkit-part-4-translation-codon-usage-7179d679135"""
    """Translates a DNA sequence (fwd strand) into an amino acid sequence """

    return codontab[dna[offset : offset + 3]]


def seq_handler(dna: str) -> list:
    # TODO: Only handling 'N' this time
    return [dna.replace("N", nt) for nt in dna_nt] if "N" in dna else [dna]


def SSVL_generator(
    lib: list,
    logs: list,
):
    with open(SSVL_TARGET, mode="r") as target, open("log.csv", mode="w") as log:
        sequence = target.read()

        fieldnames = [
            "FWD_SRC",
            "FWD_DST",
            "FWD_RESULT",
            "REV_SRC",
            "REV_DST",
            "REV_RESULT",
        ]
        writer = csv.DictWriter(
            log, fieldnames=fieldnames
        )  # For logging variant generation activities
        writer.writeheader()

        fwd = sequence
        rev = get_rev_strand(fwd)

        assert len(fwd) == len(
            rev
        ), f"fwd and rev sequence have different length, fwd_len: {len(fwd)}, rev_len: {len(rev)}"

        for offset in range(INIT_POS, len(fwd) - 2, CODON_LENGTH):
            # Fetch AA code

            # Fwd & Rev
            for codon in codonlst:
                fwd_var = "".join([fwd[0:offset], codon, fwd[offset + CODON_LENGTH :]])
                rev_var = "".join([rev[0:offset], codon, rev[offset + CODON_LENGTH :]])
                lib.append(fwd_var)
                lib.append(rev_var)

                logs.append(
                    {
                        "FWD_SRC": fwd[offset : offset + CODON_LENGTH],
                        "FWD_DST": codon,
                        "FWD_RESULT": fwd_var,
                        "REV_SRC": rev[offset : offset + CODON_LENGTH],
                        "REV_DST": codon,
                        "REV_RESULT": rev_var,
                    }
                )

        writer.writerows(logs)

        return None


def main():
    # Storage initialization
    lib = []
    logs = []
    search_result = defaultdict(int)

    # SSVL library generation
    SSVL_generator(lib, logs)
    print(f"All codon variants generated.\n" f"LIBRARY_SIZE: {len(lib) // 2}")

    # Finding restriction sites
    # TODO: Refactoring
    for enz, site in enzymetab.items():
        # Parsing restriction sites for regex operation...
        substrings = seq_handler(site)  # TODO: How to deal with 'N'?
        for candidate in lib:
            # Find all possible cutting sites
            for subs in substrings:
                # Can be multiple sequences (ex. 'N')
                search_result[enz] += candidate.count(subs)

    pd.DataFrame(
        [(k, v) for k, v in search_result.items()], columns=["Enzyme", "Occurrences"]
    ).to_csv("scan_result.csv")


if __name__ == "__main__":
    main()
