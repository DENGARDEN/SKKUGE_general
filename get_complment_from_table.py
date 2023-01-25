import pandas as pd

INPUT_CSV = "PBS RTT conversion.csv"
TARGET_COLUMN = "PBS RTT"

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

if __name__ == "__main__":
    df = pd.read_csv(INPUT_CSV)
    df["REV_COMP"] = df[TARGET_COLUMN].apply(rev_comp_if)
    df.to_csv(f"{INPUT_CSV}_REVCOMP.csv", index=False)
    

    

    
