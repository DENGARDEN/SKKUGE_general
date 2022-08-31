import subprocess
import shlex

example_cmd = "/flash HKK_220314__10_1.fastq.gz HKK_220314__10_2.fastq.gz -M 400 -m 10 -O -o 0314_10"
DATE = "220805"
FILENAMEPOS = 2
INPUT_FILES = """
HKK_220805_1_1.fastq.gz
HKK_220805_1_2.fastq.gz
HKK_220805_2_1.fastq.gz
HKK_220805_2_2.fastq.gz
HKK_220805_3_1.fastq.gz
HKK_220805_3_2.fastq.gz
HKK_220805_4_1.fastq.gz
HKK_220805_4_2.fastq.gz
HKK_220805_5_1.fastq.gz
HKK_220805_5_2.fastq.gz
HKK_220805_6_1.fastq.gz
HKK_220805_6_2.fastq.gz
HKK_220805_7_1.fastq.gz
HKK_220805_7_2.fastq.gz
HKK_220805_8_1.fastq.gz
HKK_220805_8_2.fastq.gz
HKK_220805_9_1.fastq.gz
HKK_220805_9_2.fastq.gz
HKK_220805_10_1.fastq.gz
HKK_220805_10_2.fastq.gz
HKK_220805_11_1.fastq.gz
HKK_220805_11_2.fastq.gz
HKK_220805_12_1.fastq.gz
HKK_220805_12_2.fastq.gz
"""

# https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list
def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)

if __name__ == "__main__":
    files = INPUT_FILES.split('\n')[1:] # Additional '\n' in the head

    waiting_queue = []
    for fwd, rev in grouped(files,2):
        waiting_queue.append((fwd, rev))

        print(f"Target loaded successfully: {waiting_queue[-1]}")

    # Generate child processes for FLASH
    for fwd, rev in waiting_queue:
        filename = fwd.split('_')[FILENAMEPOS]


        cmd_input = shlex.split(
            f"./flash {fwd} {rev} -M 400 -m 10 -O -o {DATE}_{filename}")

        process = subprocess.Popen(cmd_input)

        

    
