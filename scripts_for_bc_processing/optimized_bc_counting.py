import numpy as np
import pandas as pd
from Levenshtein import distance as levenshtein_distance
import sys
from rapidfuzz import process, fuzz
import os

direc= 'bc_counting_output'
if not os.path.exists(f'{direc}/'):
    os.makedirs(f'{direc}/')


# Load Data
sample = sys.argv[1]
sample_sheet_reference_file = sys.argv[2]
bc_reference_file = sys.argv[3]

path_to_data = '/scratch/groups/dpetrov/fit_natty_2024/regex_output/merged_fastqs/18032-175/'
print(f'{path_to_data}{sample}_demultiplexed.csv')

data_df = pd.read_csv(f'{path_to_data}{sample}_demultiplexed.csv')
data_df['admera_ID'] = data_df['sample_name'].str.split('/').str[-1]
sample_reference_df = pd.read_csv(sample_sheet_reference_file)
barcode_reference_file = pd.read_csv(bc_reference_file)

# Diagnostic Information
diagnostic_info = {
    'num_index1_nan': data_df['index_1_match'].isna().sum(),
    'num_index2_nan': data_df['index_2_match'].isna().sum(),
    'frac_index1_nan': data_df['index_1_match'].isna().mean(),
    'frac_index2_nan': data_df['index_2_match'].isna().mean(),
}

# pd.DataFrame(diagnostic_info, index=[0]).to_csv(f'{sample}_diagnostic_info.csv')

# Clean Data
cleaned_df = data_df.dropna(subset=['index_1_match', 'index_2_match']).copy()
cleaned_df['index_1_num'] = pd.to_numeric(cleaned_df['index_1_matched_name'].str.extract(r'^NF(\d+)$')[0], errors='coerce')
cleaned_df['index_2_num'] = pd.to_numeric(cleaned_df['index_2_matched_name'].str.extract(r'^NR(\d+)$')[0], errors='coerce')

cleaned_df['row_status'] = np.where(
    (cleaned_df['index_1_num'] == cleaned_df['index_2_num']),
    'Valid Match',
    'Mismatched Numbers'
)
row_status_counts = cleaned_df['row_status'].value_counts().to_dict()
diagnostic_info.update({
    'valid_matches': row_status_counts.get('Valid Match', 0),
    'mismatched_numbers': row_status_counts.get('Mismatched Numbers', 0),
    'invalid_patterns': row_status_counts.get('Invalid Pattern', 0)
})
diagnostic_info.update({
    'valid_matches_fraction': row_status_counts.get('Valid Match', 0) / len(cleaned_df),
    'mismatched_numbers_fraction': row_status_counts.get('Mismatched Numbers', 0) / len(cleaned_df),
    'invalid_patterns_fraction': row_status_counts.get('Invalid Pattern', 0) / len(cleaned_df)
})
print(diagnostic_info)

# Save updated diagnostic info
pd.DataFrame(diagnostic_info, index=[0]).to_csv(f'{direc}/{sample}_diagnostic_info.csv', mode='a', header=False)
0,72445,112521,0.012742423978587733,0.019791431962104634,5516649,217,0,0.999960666073818,3.933392618200261e-05,0.0

working_df = cleaned_df[cleaned_df['row_status'] == 'Valid Match']

# Prepare Reference Data
sky_bc_list = ['CTCGTGCTGGGTATGCCGGAC', 'TCAAAAGTACCGTATTAACAC', 'TACTTACAAGCAAAACCTATC', 
               'CGCGAGTACCTTTTATCGTTC', 'TCAGTTTACGCTGGTACGGTC', 'ATGTCCCTCGGATGTTCTCAC', 
               'TTTCTTTGATGTTTGTAGCTC', 'GCCTGGTAACAGTGGTAACTC', 'CTATGTATGCGTTAGTCGAGC', 
               'TTGTTTATATTCCCACGAGTC', 'AAAATCTTTAATGGCCAGTAC', 'CCTAAGCCGGTCGAACGTAGC', 
               'AGCTAATGATCTTAAATTTAC', 'CATTACTGAATCATCTTATCC']

reference_sequences = np.concatenate([barcode_reference_file['consensus barcode'], sky_bc_list])
reference_labels = np.concatenate([barcode_reference_file['Strain_id'], [f'sky_bc_{i+1}' for i in range(len(sky_bc_list))]])

reference_df = pd.DataFrame({'Barcode': reference_sequences, 'Barcode_ID': reference_labels})

# Define function to find the best match using rapidfuzz
def find_best_match_rapidfuzz(query_seq, reference_seqs, threshold=90):
    """
    Find the best match for a query sequence using RapidFuzz.
    Args:
        query_seq (str): The query barcode sequence.
        reference_seqs (array-like): Array of reference barcode sequences.
        threshold (int): Minimum similarity score to consider a match.
    Returns:
        tuple: (best_match, best_score) or (None, None) if no match found.
    """
    if pd.isna(query_seq):
        return None, None

    # Use rapidfuzz's process.extractOne for fast similarity search
    result = process.extractOne(query_seq, reference_seqs, scorer=fuzz.ratio, score_cutoff=threshold)
    if result:
        best_match, best_score = result[0], result[1]
        return best_match, best_score
    return None, None

# Updated barcode mapping function
def map_barcodes_with_rapidfuzz(input_df, reference_df, threshold=90):
    """
    Map barcodes in the input DataFrame to reference barcodes using RapidFuzz.
    Args:
        input_df (DataFrame): The input DataFrame containing barcodes to map.
        reference_df (DataFrame): Reference DataFrame with barcodes and IDs.
        threshold (int): Minimum similarity score to consider a match.
    Returns:
        DataFrame: Updated input DataFrame with mapped barcodes and counts.
    """
    # Prepare reference sequences and mapping
    reference_seqs = reference_df['Barcode'].values
    barcode_to_id = reference_df.set_index('Barcode')['Barcode_ID'].to_dict()

    # Map each barcode using rapidfuzz
    corrected_barcodes = []
    barcode_ids = []

    for e, barcode in enumerate(input_df['barcode']):
        best_match, best_score = find_best_match_rapidfuzz(barcode, reference_seqs, threshold=threshold)
        if best_match is not None:
            corrected_barcodes.append(best_match)
            barcode_ids.append(barcode_to_id[best_match])
        else:
            corrected_barcodes.append(None)
            barcode_ids.append(None)

        # Optional progress update for large datasets
        if e % 1000 == 0:
            print(f'Processing barcode {e} out of {len(input_df)}')

    # Add results to the DataFrame
    input_df['Corrected_Barcode'] = corrected_barcodes
    input_df['Barcode_ID'] = barcode_ids

    # Count occurrences of corrected barcodes
    count_df = input_df.groupby(['Corrected_Barcode']).size().reset_index(name='Count')

    # Merge counts back into the input DataFrame
    input_df = input_df.merge(count_df, on='Corrected_Barcode', how='left')

    return input_df

output_df = map_barcodes_with_rapidfuzz(working_df, reference_df)



# Count Unique Combinations
counts_df = (
    output_df
    .groupby(['admera_ID', 'index_1_index_2', 'Barcode_ID'], as_index=False)
    .size()
    .rename(columns={'size': 'count'})  # Rename the size column to 'count'
)

# Final Merge
sample_reference_df.rename(columns={'1st_step_primers': 'index_1_index_2'}, inplace=True)
final_df = pd.merge(sample_reference_df, counts_df, on=['admera_ID', 'index_1_index_2'], how='inner')
# Save Final Counts
final_df.to_csv(f'{direc}/{sample}_final_counts.csv', index=False)

