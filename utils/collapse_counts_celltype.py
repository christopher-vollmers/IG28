#!/usr/bin/env python3
"""
Aggregates isoform count data from the Isofroms.filtered.clean.quant file into a new quantification file per cell type.

Inputs:
  --IG28_output_dir   Path to the output directory of IG28.

Example usage:
  python isoforms_counts_cell_type.py --IG28_output_dir IG28_output
"""

import argparse
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregates isoform count data from the Isofroms.filtered.clean.quant file into a new quantification file per cell type."
    )
    parser.add_argument(
        "--IG28_output_dir", "-i", required=True,
        help="Path to the output directory of IG28."
    )
    args = parser.parse_args()
    return args

def get_variables(IG28_output_dir):
    """
    Get the variables from the IG28 output directory.
    """
    quant_file_genes = f"{IG28_output_dir}/Combined/Isoforms.filtered.clean.genes.quant"
    output_quant_genes = f"{IG28_output_dir}/Combined/Isoforms.filtered.clean.genes.celltype.quant"
    output_rpm_genes = f"{IG28_output_dir}/Combined/Isoforms.filtered.clean.genes.celltype.rpms"
    fastainfo_file = f"{IG28_output_dir}/Combined/tmp/fastaFiles.info"
    return quant_file_genes, output_quant_genes, output_rpm_genes, fastainfo_file

def build_cell_type_mapping(fastainfo_path):
    """
    Parse the fastafile.info file and return a dict: { cell_id -> cell_type }.
    Uses column index 1 (cell ID) and index 2 (cell type).
    """
    cell_type_map = {}
    with open(fastainfo_path, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            cell_id = fields[1]
            cell_type = fields[2]
            cell_reads = int(fields[3])
            if cell_reads == 0:
                print(f"Warning: Cell {cell_id} has 0 reads")
            cell_type_map[cell_id] = (cell_type, cell_reads)
    print(f"Cell-cell type map built.")
    return cell_type_map


def aggregate_counts(quants, cell_type_map):
    """
    Read the counts Isoforms.filtered.clean.quant file, group cell columns by cell type, and sum counts.
    Returns two lists of strings: one for the quantification data and one for the RPM data.
    """
    output_list_quants = []
    output_list_rpms = []
    # Read the first line to get the header
    with open(quants, "r") as counts_file:
        header = counts_file.readline().strip().split("\t")
        row_string = ""
        # Define a dictionary and a list to store the sample_celltype and the index of the column in which 
        # they will be positioned in the aggregated count output file
        new_column_map = defaultdict(list)
        new_column_order = []
        cell_type_reads = defaultdict(int)
        header_cols = []
        for col_idx, col_name in enumerate(header):
            fields = col_name.strip().split("_")
            # Header fields can't be split, so continue
            if len(fields) == 1:
                row_string += col_name + "\t"
                header_cols.append(col_idx)
                continue
            cell_type, cell_reads = cell_type_map.get(col_name, ("Unknown", 0))
            new_col_name = f"{cell_type}"
            if new_col_name not in new_column_map:
                new_column_order.append(new_col_name)
                row_string += new_col_name + "\t"
                cell_type_reads[new_col_name] = 0
            new_column_map[new_col_name].append(col_idx)
            cell_type_reads[new_col_name] += int(cell_reads)

        row_string += "\n"
        output_list_quants.append(row_string)
        output_list_rpms.append(row_string)
        
        for line in counts_file:
            row_quants = ""
            row_rpms = ""
            fields = line.strip().split("\t")
            for col_idx in header_cols:
                row_quants += fields[col_idx] + "\t"
                row_rpms += fields[col_idx] + "\t"
            for new_col_name in new_column_order:
                group_sum = 0
                for col_idx in new_column_map[new_col_name]:
                    group_sum += int(fields[col_idx])

                row_quants += str(group_sum) + "\t"
                row_rpms += f"{float((group_sum / cell_type_reads[new_col_name]) * 1000000):.3f}" + "\t"
            row_quants += "\n"
            row_rpms += "\n"
            output_list_quants.append(row_quants)
            output_list_rpms.append(row_rpms)
        counts_file.close()
        print(f"Aggregated counts")
        return output_list_quants, output_list_rpms

def write_output(output_path, output_list):
    with open(output_path, "w") as f:
        for line in output_list:
            f.write(line)
    print(f"Wrote output file to {output_path}")

def main():
    args = parse_args()
    quant_file_genes, output_quant_genes, output_rpm_genes, fastainfo_file = get_variables(args.IG28_output_dir)
    celltype_map = build_cell_type_mapping(fastainfo_file)
    output_list_quants, output_list_rpms = aggregate_counts(quant_file_genes, celltype_map)
    write_output(output_quant_genes, output_list_quants)
    write_output(output_rpm_genes, output_list_rpms)

 
if __name__ == "__main__":
    main()
