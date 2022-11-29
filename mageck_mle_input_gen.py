import pandas as pd
from tqdm import tqdm

import os
import pathlib

DATE = "20221112"
FILE = "/home/dengarden/Dropbox/KRAS_drug_resistance/KRAS NGS/Read Analysis/mageck/src/20221112_MRTX_merged.csv"
DESIGN_FILE = "/home/dengarden/Dropbox/KRAS_drug_resistance/KRAS NGS/Read Analysis/mageck/src/design.txt"
TARGET_DIR = "/home/dengarden/Dropbox/KRAS_drug_resistance/KRAS NGS/Read Analysis/mageck/src/MRTX+no_read_filtered+20221112"
READ_FILTER_CUTOFF = 0
FRAGMENTS = ["FRAG1", "FRAG2", "FRAG3"]
GENE_COLUMN_NAME = "index"


def design_reader(path):
    with open(path, "r") as f:
        return [s.split() for s in f.readlines() if s[0] not in "#\n"]


if __name__ == "__main__":
    experiments = design_reader(DESIGN_FILE)

    df = pd.read_csv(FILE, sep=",")
    df[["Fragment", "Pos", "Position", "SRC_AA", "DST_AA"]] = df["index"].str.split(
        "_", expand=True
    )
    df["AA_NOTATION"] = (
        df["SRC_AA"].astype(str) + df["Position"].astype(str) + df["DST_AA"].astype(str)
    )

    # Column order rearrangement
    cols = df.columns.tolist()
    cols = cols[0:2] + cols[-6:] + cols[2:-6]

    df = df[cols].drop(["Pos"], axis=1)

    for exp in tqdm(experiments):

        sample_info = "_".join(exp[-1].split("_")[:-1])
        exp.insert(0, GENE_COLUMN_NAME)
        for frag in FRAGMENTS:
            # Writing read count table
            df_in_focus = df.loc[df["Fragment"] == frag][exp]

            df_table = pd.concat(
                [
                    pd.DataFrame(
                        [0 for i in range(df_in_focus.shape[0])], columns=["sgRNA"]
                    ),
                    df_in_focus.reset_index(drop=True),
                ],
                axis=1,
            )

            colnames = list(df_table.columns)
            colnames.remove(GENE_COLUMN_NAME)
            dtypes = ["int" for i in range(len(colnames))]
            dtype_dict = {colnames[i]: dtypes[i] for i in range(len(colnames))}

            df_table = df_table.astype(dtype_dict)
            df_table = df_table.loc[df_table[df_table.columns[2]] >= READ_FILTER_CUTOFF]

            pathlib.Path(TARGET_DIR).mkdir(parents=True, exist_ok=True)

            df_table.to_csv(
                os.path.join(TARGET_DIR, f"{DATE}+{sample_info}+table+{frag}.txt"),
                sep="\t",
                index=False,
            )

            # Writing a design file
            samples = colnames[1:]
            with open(
                os.path.join(TARGET_DIR, f"{DATE}+{sample_info}+design+{frag}.txt"), "w"
            ) as design_file:
                # First line
                design_file.write("Samples\t")
                exp_groups = "\t".join(samples)
                design_file.write(f"{exp_groups}\n")

                for idx, sample in enumerate(samples):
                    design_file.write(f"{sample}")
                    for i in range(len(samples)):
                        if i == idx or i == 0:
                            design_file.write("\t1")
                        else:
                            design_file.write("\t0")

                    design_file.write("\n")
