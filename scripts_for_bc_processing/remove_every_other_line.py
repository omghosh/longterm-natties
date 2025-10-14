# read in an inp file and delete every other line in the file
# save the new file as a new file 

import os
import pandas as pd

def remove_every_other_line(file_path, output_file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    with open(output_file_path, 'w') as f:
        for i, line in enumerate(lines):
            if i % 2 == 0:
                f.write(line)


remove_every_other_line('/Users/olivia/Desktop/PetrovLab/natty/longterm/whole_genome/wgs_merged_samples.inp', '/Users/olivia/Desktop/PetrovLab/natty/longterm/whole_genome/wgs_merged_samples.inp')
