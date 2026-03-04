import csv
import argparse
import os
import sys

def convert_csv_to_tsv(input_path, output_path):
    try:
        with open(input_path, 'r', encoding='utf-8') as f_in:
            # Use excel dialect to handle quotes and special characters correctly
            reader = csv.reader(f_in)
            
            with open(output_path, 'w', encoding='utf-8', newline='') as f_out:
                writer = csv.writer(f_out, delimiter='\t')
                
                rows_converted = 0
                for row in reader:
                    writer.writerow(row)
                    rows_converted += 1
                    
        print(f"✅ Success: {rows_converted} rows written to '{output_path}'")
        
    except FileNotFoundError:
        print(f"❌ Error: The file '{input_path}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a CSV file to a TSV file.")
    
    # Arguments
    parser.add_argument("input", help="Path to the source CSV file")
    parser.add_argument("-o", "--output", help="Path to the output TSV file (optional)")

    args = parser.parse_args()

    # Logic to auto-generate output name if not provided
    if not args.output:
        base = os.path.splitext(args.input)[0]
        args.output = f"{base}.tsv"

    convert_csv_to_tsv(args.input, args.output)
