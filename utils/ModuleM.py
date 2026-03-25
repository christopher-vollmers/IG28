import re
import sys
import os
import pysam


def modifyMetaData(bamFile,fastaInfo,metaData):
    """
    Load readname->cell_type_or_cluster mapping from a two-column text file.
    Returns a dictionary of cell identifiers and cell types/clusters.
    """

    outInfo=open(f'{fastaInfo}.tmp','w')
    mapping = {}
    print('reading meta data')
    with open(metaData) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            cell, cluster = parts[0], parts[1]
            cluster = re.sub(r'[^a-zA-Z0-9_]', '', cluster)
            mapping[cell] = cluster

    print('rewriting meta data file', fastaInfo)
    with open(fastaInfo) as f:
        for line in f:
#            print(line)
            a = line.strip().split('\t')
            cell=a[1]
#            print(cell)
            if cell in mapping:
#                print(cell,mapping[cell])
                a[2]=mapping[cell]
            outInfo.write('\t'.join(a)+'\n')

    print('rewriting bam file')
    with pysam.AlignmentFile(bamFile, "rb", ignore_truncation=True) as bam:
        outbam = pysam.AlignmentFile(f'{bamFile}.tmp', mode="wb", header=bam.header)
        for aln in bam:
            read_name = aln.query_name
            cell,tissue,name,*rest=read_name.split('|')
            if cell in mapping:
                tissue = mapping[cell]
                read_name = f'{cell}|{tissue}|{name}'
                aln.query_name = read_name
                outbam.write(aln)
            else:
                outbam.write(aln)

        outbam.close()
    os.system(f'mv {fastaInfo}.tmp {fastaInfo}')
    os.system(f'mv {bamFile}.tmp {bamFile}')
