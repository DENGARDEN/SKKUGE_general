import random

SEQ_PRIMERS = {
    "WPRE_seqrR1": "ggcattaaagcagcgtatccac",
    "EF1a_seqF1": "tgtactggctccgcc",
    "3a_seqF4": "ACAGGTCAACATTGTGAAGAAGAC",
    "3a_seqF3": "GCTCCCAGATCCTGAAGGAG",
    "3a_seqF2": "ACAAAGGTGAAGTACGTGACAGA",
    "3a_seqF1": "CCAGCTGCCTGGAGAG",
}

with open(f"num", "w") as f:
    rvals = random.sample(range(21, 33, 1), 6)
    for i in rvals:
        for k, v in SEQ_PRIMERS.items():
            f.write(f"{i}\t{k}\t{v}\n")
