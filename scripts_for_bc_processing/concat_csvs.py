import pandas as pd
import glob
import sys
# path_to_files = sys.argv[1]

# Define the path to your CSV files (adjust the path as needed)
path_to_files = "/scratch/groups/dpetrov/Natty_LongTerm/bc_counting_output_July2025/*_final_counts.csv"

# Get a list of all CSV files
csv_files = glob.glob(path_to_files)

# Combine all CSV files into one DataFrame
dataframes = [pd.read_csv(file) for file in csv_files]
combined_df = pd.concat(dataframes, ignore_index=True)

# Save the combined DataFrame to a new CSV file
output_file = "all_bc_counts_July2025.csv"
combined_df.to_csv(output_file, index=False)

print(f"Combined file saved to {output_file}")
