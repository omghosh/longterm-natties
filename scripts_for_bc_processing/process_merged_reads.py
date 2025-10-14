# process merged reads to extract first step indices and barcodes
# written Nov 20, 2024 
# Olivia M. Ghosh 

# Usage: python3 process_merged_reads.py <sample_name>

# sample_name: name of the sample
# primer_reference_table: path to the primer reference table

# Output: a csv file with the extracted indices and barcodes

# Importing required libraries
import regex
import numpy as np
import pandas as pd
import sys
import os
from concurrent.futures import ThreadPoolExecutor

## inputs
sample_name = sys.argv[1]
seq_lane_name = sample_name.split('-')[0]+ '-' + sample_name.split('-')[1]

path_to_sample = f'merged_fastqs/{sample_name}.assembled.fastq'

primer_ref = "primer_table_natties.csv"


primer_reference_table = pd.read_csv(primer_ref)
primer_reference_table['index'] = primer_reference_table['index'].astype(str).str.strip().str.upper()


direc= 'regex_output_April2025'
if not os.path.exists(f'{direc}/'):
    os.makedirs(f'{direc}/')

############
# Functions #
############

def FourLineFastq(handle):
    """
    Reads 4 lines from the file and returns 1, 2, and 4 with newlines stripped
    The only check for fastq format is that line 3 starts with '+'
    """
    while True:
        line = handle.readline()    
        if not line:
            # end of file
            break       
        title = line.rstrip()
        seq = handle.readline().rstrip()
        jnk_line = handle.readline()
        if jnk_line[0] != "+":
            print(title, seq, jnk_line)
            raise ValueError("Looks like this isn’t a strictly 4-line fastq file")
        qual = handle.readline().rstrip()
        yield (title, seq, qual)

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a':'T', 'c':'G', 'g':'C', 't':'A', 'N':'N'}
    return ''.join([complement[base] for base in reversed(seq)])


count = 0
qual_fail = 0
reg_fail = 0
reg_fail_indices = 0
bc_info = {'bc_start': 82, 'bc_len': 20, 'regex_group': 2}
extracted_data = []

# Define regexes for barcode extraction
regexes = [
            regex.compile('\D*?(TGCAGAT)(\D{20})(GACATG)\D*'), # normal barcode
            regex.compile('\D*?(TGCAGAT){e<=1}(\D{20})(GACATG){e<=1}\D*'), # normal barcode with single mutation
            regex.compile('\D*?(TGTGCG){e<=1}(\D{20})(GACATG){e<=1}\D*'), #truncated upstream by 5 
            regex.compile('\D*?(TGCAGAT){e<=1}(\D{28})(GGAGGC){e<=1}\D*'), #truncated downstream by 5
            regex.compile('\D*?(TGCAGAT){e<=1}(\D{28})(CATGGAG){e<=1}\D*'), #truncated downstream by 2
            regex.compile('\D*?(TGCAGAT){e<=1}(\D{12})(GACATG){e<=1}\D*'), # 12 bp barcode
            regex.compile('\D*?(TGCAGAT){e<=1}(\D{14})(GACATG){e<=1}\D*'), # 14 bp barcode
            regex.compile('\D*?(TGCAGAT){e<=1}(\D{18,22})(GACATG){e<=1}\D*'), # 18-22 bp barcode
            regex.compile('\D*?(TGCAGAT){e<=1}(\D{26})(GACATG){e<=1}\D*'), # 26 bp barcode
        ]

strict_regexes = [regex.compile('.*?(GTTATGTGCGCAGAT)(.*?)(GACATGGAGGCCCAG).*')]
all_bc_regexes = strict_regexes + regexes

# Define regex for index extraction
index_pattern = r"^(.*?)GACTATTCTGAT{e<=3}.*CTTGACGTGC{e<=3}(.*?)$"


# Read the fastq file and extract barcodes

file_name =  path_to_sample #f'{sample_name}.assembled.fastq'
with open(file_name, "rt") as handle:
    for title, seq, qual in FourLineFastq(handle):
        count += 1
        quality_score = np.mean([ord(c) - 33 for c in qual[bc_info['bc_start']:bc_info['bc_start'] + bc_info['bc_len']]])
        if quality_score < 30:
            qual_fail += 1
            continue

        # Compile regex and search
        match = regex.match(index_pattern, seq)

        if match:
            index_1 = match.group(1)  # Extract the first index
            index_2 = match.group(2)  # Extract the second index
            index_2_reversed = reverse_complement(index_2)  # Reverse complement the second index

        else:
            reg_fail_indices += 1
            continue
        # Extract barcode
        barcode = None
        for regex_pattern in all_bc_regexes:
            match = regex_pattern.match(seq)
            if match:
                barcode = match.group(bc_info['regex_group'])
                break
        if not barcode:
            reg_fail += 1
            continue
        extracted_data.append((sample_name, index_1, index_2_reversed, barcode))

df = pd.DataFrame(extracted_data, columns=['sample_name', 'index_1', 'index_2', 'barcode'])

# Print summary
print(f'File contained {count} reads \n')
print(f'Quality score failed for {qual_fail} reads \n')
print(f'Regex for barcodes failed for {reg_fail} reads \n')
print(f'Index regex failed for {reg_fail_indices} reads \n')


# Precompile regex patterns with fuzzy matching
compiled_patterns = {
    index: regex.compile(f"{index}{{e<=2}}")
    for index in primer_reference_table['index']
}

def process_row(row, column_name):
    """Processes a row and matches the specified index column."""
    index_value = row[column_name]
    for index, pattern in compiled_patterns.items():
        if pattern.fullmatch(index_value):
            match = {
                "matched_index": index,
                "confidence": 1 - sum(pattern.match(index_value).fuzzy_counts) / len(index),
                "name": primer_reference_table[primer_reference_table['index'] == index]['name'].values[0],
            }
            return match
    return {"matched_index": None, "confidence": 0, "name": None}

def process_all_rows(df, column_name):
    """Processes all rows for a specific index column in parallel."""
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(lambda row: process_row(row, column_name), df.to_dict('records')))

    # Extract results into separate columns
    matched_column = f"{column_name}_match"
    confidence_column = f"{column_name}_confidence"
    name_column = f"{column_name}_matched_name"

    df[matched_column] = [res['matched_index'] for res in results]
    df[confidence_column] = [res['confidence'] for res in results]
    df[name_column] = [res['name'] for res in results]

# Process both index_1 and index_2
process_all_rows(df, 'index_1')
process_all_rows(df, 'index_2')

df['index_1_index_2'] = df['index_1_matched_name'] + '/' + df['index_2_matched_name']
df.to_csv(f'{direc}/{sample_name}_demultiplexed.csv', index=False)