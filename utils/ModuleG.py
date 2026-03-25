import sys
import os

def gene_cell_counts(core, tsv):
    print(f"Reading CORE file {core}")
    if os.path.exists(core) and os.path.getsize(core)>0:
        print(f"CORE file {core} exists and is not empty")
    else:
        print(f"CORE file {core} does not exist or is empty")
        sys.exit(1)

    countDict={}
    cells=set()
    genes=set()
    with open(core, 'r') as core_file:
        for line in core_file:
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue
            if parts[1] != "Assigned":
                continue
            gene = str(parts[3].split(',')[0])
            genes.add(gene)
            read_name = parts[0].split('|')
            cell = read_name[0]
            cells.add(cell)
            if gene not in countDict:
                countDict[gene] = {}
            if cell not in countDict[gene]:
                countDict[gene][cell] = 0
            countDict[gene][cell] += 1

    cells=sorted(list(cells))
    genes=sorted(list(genes))

    tsv_file = open(tsv, 'w')
    tsv_file.write(f'gene\t' + '\t'.join(cells) + '\n')
    for gene in genes:
        row=[]
        for cell in cells:
            if cell in countDict[gene]:
                count=countDict[gene][cell]
            else:
                count=0
            row.append(str(count))
        tsv_file.write(f'{gene}\t' + '\t'.join(row) + '\n')

    # Close the TSV file
    tsv_file.close()
    # Print the time taken
    print(f'Cell-gene counts matrix written to {tsv}')
