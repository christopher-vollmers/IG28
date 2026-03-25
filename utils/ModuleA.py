import sys
import argparse
import re
import os
import numpy as np
import mappy as mp
import gzip
import time
import pysam




def parseFile(fasta,out,outInfo,reads):

    tissue='None'
    if len(fasta.split('\t'))==2:
        fasta,id = fasta.split('\t')
    if len(fasta.split('\t'))>=3:
        fasta,id,tissue,*rest = fasta.split('\t')

    use_file=False
    if os.path.exists(fasta) and os.path.getsize(fasta)>0:
        use_file=True
    if use_file:
        reads=True
        if fasta.endswith('fasta') or fasta.endswith('fastq'):
            readCounter=0
            for name,seq,qual in mp.fastx_read(fasta):
                readCounter+=1
                out.write(f'>{id}|{tissue}|{name}\n{seq}\n')
            outInfo.write(f'{fasta}\t{id}\t{tissue}\t{str(readCounter)}\n')
        elif fasta.endswith('bam'):
            if id.startswith('scRNA'):
                state,base,tag = id.split(';')
                cellDict={}
                with pysam.AlignmentFile(fasta, "rb", check_sq=False) as sam:
                    for aln in sam.fetch(until_eof=True):
                        name = aln.query_name
                        seq = aln.query_sequence
                        barcode = aln.get_tag(tag)
                        out.write(f'>{base}_{barcode}|{tissue}|{name}\n{seq}\n')
                        if f'{base}_{barcode}' not in cellDict:
                            cellDict[f'{base}_{barcode}']=0
                        cellDict[f'{base}_{barcode}']+=1
                for cell,readCounter in cellDict.items():
                    outInfo.write(f'{fasta}\t{cell}\t{tissue}\t{str(readCounter)}\n')
            else:
                readCounter=0
                with pysam.AlignmentFile(fasta, "rb", check_sq=False) as sam:
                    for aln in sam.fetch(until_eof=True):
                        name = aln.query_name
                        seq = aln.query_sequence
                        readCounter+=1
                        out.write(f'>{id}|{tissue}|{name}\n{seq}\n')
                outInfo.write(f'{fasta}\t{id}\t{tissue}\t{str(readCounter)}\n')

    else:
        print(f'\t{fasta} does not exist or is an empty file')

    return reads
