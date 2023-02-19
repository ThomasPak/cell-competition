import sys
import csv

# Collect all input file names from command line arguments
input_file_names = sys.argv[1:]

# Create a list to hold all the rows from all input files
merged_rows = []

# Iterate over each input file
for input_file_index, input_file_name in enumerate(input_file_names):
    # Open the input file and create a CSV reader object
    with open(input_file_name, 'r') as input_file:
        csv_reader = csv.reader(input_file)

        # Iterate over each row in the CSV file and append it to the merged rows list
        for row_index, row in enumerate(csv_reader):
            if input_file_index == 0 or row_index > 0:
                merged_rows.append(row)

# Create a CSV writer object for stdout
csv_writer = csv.writer(sys.stdout)

# Write each merged row to stdout using the CSV writer object
for merged_row in merged_rows:
    csv_writer.writerow(merged_row)
