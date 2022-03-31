# requires Python 3.5+
import os
import sys
from subprocess import run

import pandas as pd

from extractor import do

# OS specific running
# if platform == "linux" or platform == "linux2":
#     # linux
#     run()
# elif platform == "darwin":
#     # OS X
#     completed = run()
# elif platform == "win32":
#     # Windows...

EXTRACTOR_ONLY = False

NGS_file_info = pd.read_excel("NGS_input.xlsx", engine="openpyxl")
leading_name = "HKK_220314"
format_name = ".fastq.gz"
combined_format = ".extendedFrags.fastq"
# TBD
src1 = str()
src2 = str()
output_name = str()
# iterate through items

# execute FLASH
if not EXTRACTOR_ONLY:
    for idx, rows in NGS_file_info.iterrows():
        src1 = f"{leading_name}__{rows['#num']}_1{format_name}"
        src2 = f"{leading_name}__{rows['#num']}_2{format_name}"
        output_name = f"{rows['#num']}"
        run(
            [
                "./flash",
                f"{src1}",
                f"{src2}",
                "-M",
                "400",
                "-m",
                "10",
                "-O",
                "-o",
                f"{output_name}",
            ]
        )

# extract with barcodes
# hyperparameter
START = 20
# BASE_DIR = os.path.dirname(sys.executable)

for i in range(7):
    src_file_name = os.path.join(f"{START + i}{combined_format}")
    dest_file_name = os.path.join(f"{START + i}{combined_format}")
    do(src_file_name, dest_file_name, START + i)
# anaylze indel

