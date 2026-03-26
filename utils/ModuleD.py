#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import argparse
import re
import os
import numpy as np
import mappy as mp
import gzip
import time
import pysam
import multiprocessing as mupro


#def call_variants(infile,genome,root,out_temp):
#    os.system(f'samtools addreplacerg -O bam -r "@RG\tID:1\tSM:sample_name" -o {out_temp}{root}.fixed.bam {infile} > /dev/null 2>&1')
#    os.system(f'samtools index {out_temp}{root}.fixed.bam > /dev/null 2>&1')
#    os.system(f'bcftools mpileup -Ou -f {genome} {out_temp}{root}.fixed.bam 2> /dev/null | bcftools call -mv -Oz -o {out_temp}{root}.vcf.gz > /dev/null 2>&1')
#    os.system(f'bcftools index {out_temp}{root}.vcf.gz > /dev/null 2>&1')
#    os.system(f'whatshap phase -o {out_temp}{root}.phased.vcf.gz --reference {genome} {out_temp}{root}.vcf.gz {out_temp}{root}.fixed.bam > /dev/null 2>&1')
#    os.system(f'bcftools index {out_temp}{root}.phased.vcf.gz > /dev/null 2>&1')
#    os.system(f'whatshap haplotag --sample sample_name --reference {genome} -o {out_temp}{root}.tagged.bam {out_temp}{root}.phased.vcf.gz {out_temp}{root}.fixed.bam > /dev/null 2>&1')
#    os.system(f'samtools index {out_temp}{root}.tagged.bam > /dev/null 2>&1')


def split_files_into_tissues(temp_path,tissueList,bamFile,fastaInfo):
    info_dict = {}
    bam_dict = {}

    with pysam.AlignmentFile(bamFile, "rb", ignore_truncation=True) as bam:

        for tissue in tissueList:
            infoSubFile = f'{temp_path}/{tissue}/tmp/fastaFiles.info'
            bamSubFile = f'{temp_path}/{tissue}/tmp/mm2Alignments.sorted.bam'
            info_dict[tissue] = open(infoSubFile, 'w')
            bam_dict[tissue] = pysam.AlignmentFile(bamSubFile, mode="wb", header=bam.header)

        for line in open(fastaInfo):
            columns = line.strip().split('\t')
            tissue = columns[2]
            info_dict[tissue].write(line)

        for aln in bam:
            read_name = aln.query_name
            tissue = read_name.split('|')[1]
            bam_dict[tissue].write(aln)

def myround(x, base=10):
    '''Rounds down to the nearest base'''
    return int(base * round(float(x) / base))

def prepare_locus(chrom,start,end,left_bounds,right_bounds):
    if chrom not in left_bounds:
        left_bounds[chrom] = {'5': [], '3': []}
    if chrom not in right_bounds:
        right_bounds[chrom] = {'5': [], '3': []}

    left_bounds_sub = {}
    right_bounds_sub = {}
    left_bounds_sub =  {'5': [], '3': []}
    right_bounds_sub = {'5': [], '3': []}
    for side in ['5','3']:
        for pos in left_bounds[chrom][side]:
            if start<pos<end:
                left_bounds_sub[side].append(pos)
        for pos in right_bounds[chrom][side]:
            if start<pos<end:
                 right_bounds_sub[side].append(pos)
    return left_bounds_sub,right_bounds_sub



def clean_psl(psl_file, clean_psl_file,primary):

    minimum_intron_size=10

    used=set()
    out = open(clean_psl_file, 'w')
    for line in open(psl_file):
        a = line.strip().split('\t')
        start = int(a[15])
        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]

        name=a[9]
        if primary and name in used:
            continue

        blockstarts_clean = []
        readstarts_clean = []
        blocksizes_clean = []
        size_gap = []

        for x in range(len(blocksizes)):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            blockend = blockstart + blocksize
            size_gap.append(blocksize)

            try:
                next_blockstart = int(blockstarts[x + 1])
                gap = next_blockstart - blockend
                size_gap.append(gap)
            except IndexError:
                pass

        new_size_gap = []
        block = 0
        for index in range(len(size_gap)):
            if index % 2 == 0:
                block += size_gap[index]
            if index % 2 == 1:
                if size_gap[index] < minimum_intron_size:
                    block += size_gap[index]
                else:
                    new_size_gap.append(block)
                    new_size_gap.append(size_gap[index])
                    block = 0

        new_size_gap.append(block)

        current_position = start
        current_readposition = int(a[11])
        for index in range(0, len(new_size_gap), 1):
            inter = new_size_gap[index]
            if index % 2 == 0:
                blockstarts_clean.append(str(current_position))
                blocksizes_clean.append(str(inter))
                readstarts_clean.append(str(current_readposition))

                current_position += inter
                current_readposition += inter

            if index % 2 == 1:
                current_position += inter

        a[17] = str(len(blockstarts_clean))

        blockstarts_clean.append('')
        blocksizes_clean.append('')
        readstarts_clean.append('')

        a[18] = (',').join(blocksizes_clean)
        a[19] = (',').join(readstarts_clean)
        a[20] = (',').join(blockstarts_clean)



        new_line = ('\t').join(a)
        out.write(new_line + '\n')
        used.add(name)



def get_parsed_files(out_tmp):
    roots=[]
    for file1 in os.listdir(out_tmp):
        if os.path.isfile(out_tmp+'/'+file1):
            if file1.endswith('.bam'):
#              if 'chr1~' in file1: # for testing only
                root=file1.split('.bam')[0]
                roots.append(root)
    return roots

def getCSaroundSS(CS,genomePosition,start,end):
    '''
    Goes through the CIGAR string by pairing numbers with the letters.
    Counts the number of matches/mismatches and indels.
    Returns the total number of all, matches/mismatches (M), and indels.
    '''

    status=''
    p = re.compile(r'([\\=\\+\-\\*\\~])')
    splitCstr = [i+j for i,j in zip(p.split(CS)[1::2], p.split(CS)[2::2])]
    record=[]
    spliceIndex=False
    for i in range(len(splitCstr)):
        sub=splitCstr[i]
        status=sub[0]
        entry=sub[1:]
        if status == '=':
            for item in entry:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '+':
            for item in entry:
                record.append((status,item))
        if status == '-':
            for item in entry:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '*':
            for item in entry[::2]:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '~':
            intronLength=int(entry[2:-2])
            genomePosition+=intronLength
            record.append('|'+entry+'|')
            if start<=genomePosition<=end:
                  spliceIndex=len(record)
    bases='nnnn'
    left=''
    right=''

    if spliceIndex:
        for index in range(max(spliceIndex-10,0),min(spliceIndex+10,len(record))):
            item=record[index]
            if '|' in item:
               bases=item[1:3]+item[-3:-1]
               left = record[index-5:index]
               right = record[index+1:index+6]
    return bases,left,right

def scan_for_best_bin(entry, dist_range, density_dict, chrom, side, reverse):
    best_extra_list, peak_center = 0, 0
    best_dirs = {'+': 0, '-': 0}
    best_names = []
    for x in dist_range:
        extra_list_exp = 0
        dirs = {'+': 0, '-': 0}
        names=[]
        cov_set = {}
        for y in dist_range:
            if entry + x + y in density_dict:
                count=density_dict[entry + x + y][0]
                extra_list_exp+=count
                for item in density_dict[entry + x + y][1]:
                    dirs[item[1]] += 1
                    names.append(item[0])

        if extra_list_exp > best_extra_list:
            best_extra_list = extra_list_exp
            best_names=names
            peak_center = entry + x
            best_dirs = dirs

    Left_to_Right = best_dirs['+']
    Right_to_Left = best_dirs['-']
    Type = ''

    if Left_to_Right < Right_to_Left:
        Type = '3' if reverse else '5'
    elif Left_to_Right > Right_to_Left:
        Type = '5' if reverse else '3'


    return peak_center, Type, best_names


def determine_cov(cov_area, chrom, reverse, peak_center, histo_cov):
    cov = [0]
    cov_area2 = []
    for covered_position, count in cov_area.items():
        if count > 1:
            cov_area2.append(covered_position)

    cov_area = sorted(cov_area2, reverse=reverse)
    counter = 0
    for base_f in cov_area:
        count = False
        if reverse and base_f < peak_center:
            count = True
        elif not reverse and base_f > peak_center:
            count = True
        if count:
            if counter <= 3:
                counter += 1
                base_f = myround(base_f)
                if base_f in histo_cov[chrom]:
                    cov.append(histo_cov[chrom][base_f])
            else:
                break
    cov = max(cov)
    return cov, cov_area





def find_peaks(density_dict, reverse, side, peak_areas, chrom,start,end,splice_site_width,minimum_read_count,junctions,inFile):
    dist_range=[0]
    toWrite=[]
    for shift in range(1,splice_site_width+1):
        dist_range.append(shift)
        dist_range.append(-shift)

    entry_list = []
    for entry,density in density_dict.items():
        if density[0]>=minimum_read_count:
            entry_list.append([entry, density])
    CSdict={}
    peakDict={}
    peakList=[]
    for entry, density in sorted(entry_list, key=lambda x: x[1][0], reverse=True):
        if entry in peak_areas[side]:
            continue
        density_sub={}
        peak_sub={}
        called=False
        for x in dist_range:
            if entry + x in peak_areas[side]:
                    called=True
            else:
                for y in dist_range:
                    if entry + x + y in density_dict:
                        density_sub[entry + x + y] = density_dict[entry + x + y]

        if not called:
            peak_center, Type, best_names = scan_for_best_bin(entry, dist_range, density_sub, chrom, side, reverse)
            if not Type:
                continue

            peakList.append((peak_center, Type))
            for name in best_names:
                if name not in peakDict:
                    peakDict[name]=[]
                peakDict[name].append((peak_center,Type))

    with pysam.AlignmentFile(inFile, "rb") as sam:
        for aln in sam.fetch():
            if aln.is_mapped:
                begin = aln.reference_start            # 0-based leftmost coordinate
                cs=aln.get_tag("cs")
                name=aln.query_name
                if name in peakDict:
                    keys=peakDict[name]
                    for key in keys:
                        if key not in CSdict:
                            CSdict[key]=[]
                        CSdict[key].append((cs,begin))

    for peak_center, Type in peakList:
        CSs=CSdict[(peak_center,Type)]
        called=False
        for x in dist_range:
            for y in dist_range:
                if peak_center + x + y in peak_areas[side]:
                    called=True
        if called:
            continue
        passed = characterize_splicing_event([chrom, peak_center - splice_site_width, peak_center + splice_site_width], junctions, CSs)
        if passed:
            start=str(peak_center - splice_site_width)
            end=str(peak_center + splice_site_width)
            toWrite.append([chrom,start,end,Type,side,'1'])
            for base in range(int(start), int(end) + 1):
                peak_areas[side][base] = 1

    return peak_areas, toWrite



def rewrite_cigar(cigartuples, max_indel=10):
    MATCH, INS, DEL = 0, 1, 2
    out = []
    pending = 0

    for op, length in (cigartuples or []):
        if op == MATCH or (op == DEL and length < max_indel):
            pending += length
        elif op == INS and length < max_indel:
            pass  # drop short insertions
        else:
            if pending:
                out.append((MATCH, pending))
                pending = 0
            out.append((op, length))

    if pending:
        out.append((MATCH, pending))

    return out


def collect_reads(reads, peak_areas,side):
    histo_bases = {}

    lineCounter=0
    with pysam.AlignmentFile(reads, "rb") as sam:
        for aln in sam.fetch():
            if aln.is_mapped:
                lineCounter+=1
                name = aln.query_name
                tid= aln.reference_id
                chrom = sam.get_reference_name(tid)
                begin = aln.reference_start            # 0-based leftmost coordinate
                span = aln.reference_end
                blocks = aln.get_blocks()
                aligned_bases = aln.query_alignment_length
                nm = aln.get_tag("NM") if aln.has_tag("NM") else 0
                accuracy=(aligned_bases-nm)/aligned_bases
                aln.cigartuples = rewrite_cigar(aln.cigartuples, max_indel=10)
                clean_blocks=aln.get_blocks()
                dirn='+'
                if aln.is_reverse:
                    dirn='-'
                cov_set = set()
                low_bounds, up_bounds = [], []
                for blockstart,blockend in clean_blocks:

                    if blockstart != begin:
                        up_bounds.append(blockstart)
                    if blockend != span:
                        low_bounds.append(blockend)

                if accuracy<0.9:
                    continue
                if side =='l':
                    bounds=low_bounds
                elif side=='r':
                    bounds=up_bounds

                for bound in bounds:
                    if bound not in peak_areas[side]:
                        if bound not in histo_bases:
                            histo_bases[bound] = [0,[]]
                        histo_bases[bound][0]+=1
                        if histo_bases[bound][0]<=500:
                            histo_bases[bound][1].append([name, dirn])
                        else:
                            m = np.random.randint(histo_bases[bound][0])          # 0…t‑1
                            if m < 500:
                                histo_bases[bound][1][m]=[name, dirn]

    return histo_bases, lineCounter


def parse_genome(input_file, left_bounds, right_bounds, white_list_polyA,chromSet):
    polyAWhiteList=[]
    chrom_list = set()
    gene_dict = {}
    if input_file.endswith('.gtf.gz'):
        print('\t\tgtf file ends on .gz and will be treated as gzipped')
        input=gzip.open(input_file,'rt')
    elif input_file.endswith('.gtf'):
        input=open(input_file,'r')

    for line in input:
        pAwl=False
        for element in white_list_polyA:
            if element in line:
                pAwl=True
        a = line.strip().split('\t')
        if len(a) <= 7:
            continue
        if a[0] in chromSet:
            if a[2] == 'exon':
                testKey = a[8].split('transcript_id "')[1].split('"')[0]
                if not gene_dict.get(testKey):
                    gene_dict[testKey] = []
                gene_dict[testKey].append((a[0], a[3], a[4], a[6], pAwl))

    for transcript_id in gene_dict:
        transcript_data = gene_dict[transcript_id]

        chrom = transcript_data[0][0]
        direction=transcript_data[0][3]
        chrom_list.add(chrom)
        if chrom not in right_bounds:
            left_bounds[chrom] = {'5': [], '3': []}
            right_bounds[chrom] = {'5': [], '3': []}

        start = sorted(transcript_data, key=lambda x: int(x[1]))[0][1]
        end = sorted(transcript_data, key=lambda x: int(x[2]), reverse=True)[0][2]
        direction=transcript_data[0][3]
        pAwl=transcript_data[0][4]
        if pAwl:
            if direction=='+':
                polyAWhiteList.append((chrom,direction,end,transcript_id))
            elif direction=='-':
                polyAWhiteList.append((chrom,direction,start,transcript_id))

        for entry in transcript_data:
            if entry[1] != start:
                if entry[3] == '+':
                    right_bounds[chrom]['3'].append(int(entry[1]) - 1)
                elif entry[3] == '-':
                    right_bounds[chrom]['5'].append(int(entry[1]) - 1)
            if entry[2] != end:
                if entry[3] == '+':
                    left_bounds[chrom]['5'].append(int(entry[2]))
                if entry[3] == '-':
                    left_bounds[chrom]['3'].append(int(entry[2]))
    return chrom_list, left_bounds, right_bounds, polyAWhiteList


def make_genome_bins(bounds, side, chrom, peak_areas,splice_site_width):
    toWrite=[]
    for type1 in ['5', '3']:
        covered = {}
        position_list = sorted(bounds[type1], key=int)
        for index1 in range(0, len(position_list)):
            if index1 in covered:
                continue
            sub_list = [position_list[index1]]
            for index2 in range(index1, len(position_list)):
                if position_list[index2] - max(sub_list) <= splice_site_width:
                    sub_list.append(position_list[index2])
                    covered[index2] = 1
                else:
                    break
            single = False
            if len(sub_list) > 1:
                splice_dists = []
                for splice_pos in range(0, len(sub_list) - 1):
                    splice_dists.append(int(sub_list[splice_pos + 1]) - int(sub_list[splice_pos]))
                if min(splice_dists) > 3:
                    for x in range(0, len(sub_list), 1):
                        if x != 0:
                            start = int(sub_list[x] - ((sub_list[x] - sub_list[x - 1]) / 2))
                        else:
                            start = int(sub_list[x]) - splice_site_width
                        if x != len(sub_list) - 1:
                            end = int(sub_list[x] + ((sub_list[x + 1] - sub_list[x]) / 2))
                        else:
                            end = int(sub_list[x]) + splice_site_width

                        toWrite.append([chrom,str(start),str(end),type1,side,'A'])

                        for base in range(start, end + 1):
                            peak_areas[side][base] = 1

                else:
                    single = True
            else:
                single = True
            if single:
                start = min(sub_list) - splice_site_width
                end = max(sub_list) + splice_site_width
                toWrite.append([chrom,str(start),str(end),type1,side,'A'])
                for base in range(start, end + 1):
                    peak_areas[side][base] = 1

    return peak_areas,toWrite

def error_handler(error):
    # This function is called in the main process with the exception object
    print(f"Error: {error}")


def get_loci(infile,out_tmp,minCount,numThreads):
    print('\t\tparsing sam and splitting into loci')

    tissue_bam_fh = pysam.AlignmentFile(infile, "rb")
    tissue_chromosomes = tissue_bam_fh.references
    mupro.set_start_method('fork', force=True)
    pool = mupro.Pool(min(numThreads,5))
    for tissue_chrom in tissue_chromosomes:
#        print(f'submitted {tissue_chrom} for splitting', end='\r')
        pool.apply_async(process_chromosome,[tissue_chrom,infile,out_tmp,minCount], error_callback=error_handler)
#        process_chromosome(tissue_chrom,infile,out_tmp,minCount)
    pool.close()
    pool.join()

def process_chromosome(chrom,infile,out_tmp,minCount):
    new={}
    reads={}
    previous_start={}
    previous_end={}
    tmpBam={}

    for dir in ['+','-']:
       reads[dir]=0
       new[dir]=False
       previous_start[dir]=-1
       previous_end[dir]=-1
       tmpBam[dir]=f'{out_tmp}/{chrom}.{dir}.tmp.bam'

    roots=set()
    total_psl=0
    minCount=int(minCount)
    fh={}

    with pysam.AlignmentFile(infile,ignore_truncation=True) as bam:
        for aln in bam.fetch(chrom):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.mapping_quality<30:
                continue  # skip unmapped reads for overlap logic
            previous_chrom=chrom
            start = aln.reference_start            # 0-based leftmost coordinate
            end = aln.reference_end
            direction='+'
            if aln.is_reverse:
                direction='-'

                # reference_start + aligned length
            if start>previous_end[direction]:
                new[direction]=True
            if not new[direction]:
                previous_end[direction]=max(end,previous_end[direction])
                fh[direction].write(aln)
                reads[direction]+=1
                total_psl+=1
            else:
                if reads[direction]:
                    root=previous_chrom+'~'+direction+'~'+str(previous_start[direction])+'~'+str(previous_end[direction])
                    roots.add(root)
                    fh[direction].close()
                    if reads[direction]>=minCount:
                        os.system(f'mv {tmpBam[direction]} {out_tmp}/{root}.bam')
                    reads[direction]=0

                fh[direction] = pysam.AlignmentFile(tmpBam[direction], mode="wb", header=bam.header)
                fh[direction].write(aln)
                reads[direction]+=1
                total_psl+=1
                previous_end[direction]=end
                previous_start[direction]=start
                new[direction]=False

        for dir in ['+','-']:
            if reads[dir]:
                root=previous_chrom+'~'+dir+'~'+str(previous_start[dir])+'~'+str(previous_end[dir])
                roots.add(root)
                fh[dir].close()
                if reads[dir]>=minCount:
                    os.system(f'mv {tmpBam[dir]} {out_tmp}/{root}.bam')
    for dir in ['+','-']:
        if os.path.isfile(f'{tmpBam[dir]}'):
            os.system(f'rm {tmpBam[dir]}')
    print(f'\t\tsplit {total_psl:,} sam entries into {len(roots):,} loci for {chrom}',' '*50, end='\r')



def get_chromosomes(infile,out_tmp,fastaList):
    chrom_list=set()
    outDict={}
    reads=[]
    new=False
    previous_chrom=''
    previous_start=0
    previous_end=0
    roots=set()
    total_psl=0
    print('\t\tparsing psl and splitting into loci')
    for line in open(infile):
        a=line.strip().split('\t')
        chrom=a[13]
        chrom_list.add(chrom)
        start=int(a[15])
        end=int(a[16])
        if chrom!=previous_chrom or start>previous_end:
            new=True

        if not new:
            previous_end=max(end,previous_end)
            reads.append(line)
        else:
            if reads:
                root=previous_chrom+'~'+str(previous_start)+'~'+str(previous_end)
                print('\t\tdefining locus',(' ').join(root.split('~')),' '*20,end='\r')

                fh=open(out_tmp+'/'+root+'.psl', 'w')
                roots.add(root)

                for read in reads:
                    fh.write(read)
                    total_psl+=1
                    name=read.strip().split('\t')[9]
                    outDict[name]=root
                fh.close()
                reads=[]
            reads.append(line)   ### bug in v4.0.0. Left the first read in each locus unused
            previous_chrom=chrom
            previous_end=end
            previous_start=start
            new=False
    if reads:
        root=previous_chrom+'~'+str(previous_start)+'~'+str(previous_end)
        fh=open(out_tmp+'/'+root+'.psl', 'w')
        roots.add(root)
        for line in reads:
            total_psl+=1
            fh.write(line)
            name=line.strip().split('\t')[9]
            outDict[name]=root
        fh.close()

    print('\t\tsplit', total_psl, 'psl entries into loci',' '*20)



def characterize_splicing_event(a,junctions,CSs):
    passed = False
    chromosome = a[0]
    splice_left = int(a[1])
    splice_right = int(a[2])

    basecontext={}
    basecontext['allowed']=0
    basecontext['all']=0
    leftCS={}
    leftCS['Total']=0
    rightCS={}
    rightCS['Total']=0
    for status in ['*','+','-','=','|']:
        leftCS[status]=0
        rightCS[status]=0
    for cs,genomePosition in CSs:
        bases,left,right = getCSaroundSS(cs,genomePosition,splice_left,splice_right)
        basecontext['all'] += 1
#        if bases == 'gtag' or bases == 'ctac':
        if bases in junctions:
            basecontext['allowed'] += 1
        else:
            if bases not in basecontext:
                basecontext[bases]=0
            basecontext[bases]+=1
        if left:
            for thing in left:
                leftstatus = thing[0]
                leftCS[leftstatus]+=1
                leftCS['Total']+=1
        if right:
            for thing in right:
                rightstatus = thing[0]
                rightCS[rightstatus]+=1
                rightCS['Total']+=1

    basepass = False
    if basecontext['allowed']/basecontext['all']>0.85:
        basepass = True
        leftacc=leftCS['=']/leftCS['Total']
        rightacc=rightCS['=']/rightCS['Total']
    if basepass:
        if leftacc>0.85 and rightacc>0.85:
            passed = True
    CIGARs, CigarDict = 0, 0
    return passed



def find_ends(starts, ends, identity, count_dict, upstream_buffer,downstream_buffer,minimum_feature_count):
    start_peaks, end_peaks, extend = {}, {}, True
    if '+' in identity:
        count_dict['+'].add(identity)
    else:
        count_dict['-'].add(identity)


    start_count = {}
    end_count = {}
    for position in starts:
        if position not in start_count:
            start_count[position] = 0
        start_count[position] += 1

    for position in ends:
        if position not in end_count:
            end_count[position] = 0
        end_count[position] += 1

    show=False
    for position in sorted(starts, key = lambda x: start_count[x],reverse = True):
        if position - upstream_buffer not in start_peaks:
            window_count = 0
            for i in np.arange(0,10):
                window_position = position + i
                if window_position in start_count:
                    window_count += start_count[window_position]
            if window_count >= minimum_feature_count:
                original_bin=[]
                for shift in np.arange(-upstream_buffer, downstream_buffer):
                    start_peaks[position + shift] = position
                    original_bin.append(position + shift)
                if extend:
                    bins=[]
                    for i in np.arange(min(original_bin),max(original_bin)):
                        bin=0
                        for shift in np.arange(0,10):
                            if i+shift in start_count:
                                 bin+=start_count[i+shift]
                        bins.append(bin)

                    best_bin=max(bins)
                    extended=True
                    adjacent=position - upstream_buffer
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent - i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in start_count:
                                window_count += start_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in start_peaks:
                                    start_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['start_left'].add(identity)

                    extended = True
                    adjacent = position + downstream_buffer - 1
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent + i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in start_count:
                                window_count += start_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in start_peaks:
                                    start_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                     extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['start_right'].add(identity)

    for position in sorted(ends, key = lambda x: end_count[x],reverse=True):
        if position + upstream_buffer -1 not in end_peaks:
            window_count = 0
            for i in np.arange(0,10):
                window_position = position - i
                if window_position in end_count:
                    window_count += end_count[window_position]
            if window_count >= minimum_feature_count:
                original_bin=[]
                for shift in np.arange(-downstream_buffer, upstream_buffer):
                    end_peaks[position + shift] = position
                    original_bin.append(position + shift)
                if extend:
                    bins=[]
                    for i in np.arange(min(original_bin),max(original_bin)):
                        bin=0
                        for shift in np.arange(0,10):
                            if i+shift in end_count:
                                 bin+=end_count[i+shift]
                        bins.append(bin)

                    best_bin=max(bins)
                    extended = True
                    adjacent = position - downstream_buffer
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent - i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in end_count:
                                window_count += end_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in end_peaks:
                                    end_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent=adjacent_pos
                        if extended:
                            count_dict['end_left'].add(identity)

                    extended=True
                    adjacent=position + upstream_buffer - 1
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent + i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in end_count:
                                window_count += end_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in end_peaks:
                                    original_bin.append(element)
                                    end_peaks[element] = position
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['end_right'].add(identity)

    return start_peaks, end_peaks, count_dict


def sort_reads_into_splice_junctions(splice_dict, infile, multi_exon_only,overhangs):
    start_end_dict, start_end_dict_mono = {}, {}
    fail_counter=0
    pass_counter=0
    with pysam.AlignmentFile(infile, "r", ignore_truncation=True) as sam:
        for aln in sam.fetch(until_eof=True):
            cigar = aln.cigartuples
            aln.cigartuples = rewrite_cigar(cigar, max_indel=10)
            clean_blocks=aln.get_blocks()
            exon_number=len(clean_blocks)
            if multi_exon_only=='0' or exon_number>1:
                name = aln.query_name
                tid= aln.reference_id
                read_chromosome = sam.get_reference_name(tid)
                start = aln.reference_start            # 0-based leftmost coordinate
                end = aln.reference_end
                blocks = aln.get_blocks()

                read_direction='+'
                if aln.is_reverse:
                    read_direction='-'

                left_clip=right_clip=0
                op, length = cigar[0]
                if op in (4, 5):  # soft or hard clip
                    left_clip = length

                op, length = cigar[-1]
                if op in (4, 5):
                    right_clip = length

                if read_direction == '+':
                    left_extra, right_extra = left_clip,right_clip
                if read_direction == '-':
                    right_extra, left_extra = left_clip,right_clip

                failed = False
                identity = read_chromosome + '_'

                if overhangs[0] <= left_extra <= overhangs[1] and overhangs[2] <= right_extra <= overhangs[3]:
                    failed=False
                else:
                    failed=True

                if not failed:
                    for x in range(0, len(clean_blocks) - 1):
                        blockstart = clean_blocks[x][0]
                        blockend = clean_blocks[x][1]
                        left_splice = blockend
                        right_splice = clean_blocks[x+1][0]

                        if left_splice not in splice_dict or right_splice not in splice_dict:
                            failed = True
                            fail_counter+=1
                            break

                        left_splice_site = splice_dict[left_splice]
                        right_splice_site = splice_dict[right_splice]
                        identity += str(left_splice_site) + '-' + str(right_splice_site) + '~'

                if not failed:

                        pass_counter+=1
#                    if identity.split('_')[1] != '':
                        if identity not in start_end_dict:
                            start_end_dict[identity] = [[],0]
                        start_end_dict[identity][1]+=1
                        if start_end_dict[identity][1]<=10000:
                            start_end_dict[identity][0].append((start, end,
                                                                read_direction))
                        m = np.random.randint(start_end_dict[identity][1])          # 0…t‑1
                        if m < 10000:
                            start_end_dict[identity][0][m]=(start, end,
                                                                read_direction)



#                    else:
#                        if identity not in start_end_dict_mono:
#                            start_end_dict_mono[identity] = [[],0]
#                        start_end_dict_mono[identity][1]+=1
#                        if start_end_dict_mono[identity][1]<=10000:
#                            start_end_dict_mono[identity][0].append((start, end,
#                                                                read_direction))
#                        m = np.random.randint(start_end_dict_mono[identity][1])          # 0…t‑1
#                        if m < 10000:
#                            start_end_dict_mono[identity][0][m]=(start, end,
#                                                                read_direction)

#    for identity in start_end_dict_mono:
#        iso_counter = 0
#        new_identity = identity + 'Mono' + str(iso_counter)
#        for direction in ['+','-']:
#            previous_end = 0
#            positions = start_end_dict_mono[identity][0]
#            for position in sorted(positions):
#                start = position[0]
#                end = position[1]
#                read_direction=position[2]
#                if direction==read_direction:
#                    if start > previous_end:
#                        iso_counter += 1
#                        new_identity = identity + 'Mono' + str(iso_counter)
#                        if new_identity not in start_end_dict:
#                            start_end_dict[new_identity] = [[],0]
#                        start_end_dict[new_identity][1]+=1
#                        start_end_dict[new_identity][0].append(position)
#                        previous_end = end
#                    else:
#                        if new_identity not in start_end_dict:
#                            start_end_dict[new_identity] = [[],0]
#                        start_end_dict[new_identity][1]+=1
#                        start_end_dict[new_identity][0].append(position)
#                        previous_end = max(end, previous_end)

    return start_end_dict


def define_start_end_sites(start_end_dict,upstream_buffer,downstream_buffer,minimum_feature_count,splice_dict,infile,overhangs,multi_exon_only):
    left_extras, right_extras = {}, {}
    isoform_counter, isoform_dict, IsoDict = 0, {}, {}

    number_of_isoforms = len(start_end_dict)
    counter = 0
    count_dict={}
    count_dict['start_left']=set()
    count_dict['start_right']=set()
    count_dict['end_left']=set()
    count_dict['end_right']=set()
    count_dict['+']=set()
    count_dict['-']=set()
    start_dicts={}
    end_dicts={}


    for identity in sorted(start_end_dict):
        counter += 1
        starts, ends = [], []
        ends = []
        dir_dict={}
        positions,count = start_end_dict[identity]
        for position in positions:
            starts.append(int(position[0]))
            ends.append(int(position[1]))

        start_dict, end_dict, count_dict = find_ends(starts, ends, identity, count_dict,upstream_buffer,downstream_buffer,minimum_feature_count)
        start_dicts[identity]=start_dict
        end_dicts[identity]=end_dict
        left_extras, right_extras = {}, {}
        combinations=set()
        for start, end, read_direction in positions:
            if int(start) in start_dict and int(end) in end_dict:
                left = start_dict[int(start)]
                right = end_dict[int(end)]
                if (left,right) not in dir_dict:
                    dir_dict[(left,right)]={}
                    dir_dict[(left,right)]['+']=0
                    dir_dict[(left,right)]['-']=0
                dir_dict[(left,right)][read_direction]+=1
                combinations.add((left,right))

        for combination in combinations:
            if dir_dict[combination]['+']>=dir_dict[combination]['-']:
                read_direction='+'
            else:
                read_direction='-'

            left,right=combination

            new_identity = identity + '_' + read_direction + '_' + str(left) + '_' + str(right)
            IsoDict[new_identity]=[[],0,{}]

    with pysam.AlignmentFile(infile, "r", ignore_truncation=True) as sam:
        for aln in sam.fetch(until_eof=True):
            cigar = aln.cigartuples
            aln.cigartuples = rewrite_cigar(cigar, max_indel=10)
            clean_blocks=aln.get_blocks()
            exon_number=len(clean_blocks)
            if multi_exon_only=='0' or exon_number>1:
                name = aln.query_name
                tid= aln.reference_id
                read_chromosome = sam.get_reference_name(tid)
                start = aln.reference_start            # 0-based leftmost coordinate
                end = aln.reference_end
                blocks = aln.get_blocks()

                read_direction='+'
                opposite_direction='-'
                if aln.is_reverse:
                    read_direction='-'
                    opposite_direction='+'
                aligned_bases = aln.query_alignment_length

                failed = False
                identity = read_chromosome + '_'

                for x in range(0, len(clean_blocks) - 1):
                    blockstart = clean_blocks[x][0]
                    blockend = clean_blocks[x][1]
                    left_splice = blockend
                    right_splice = clean_blocks[x+1][0]

                    if left_splice not in splice_dict or right_splice not in splice_dict:
                        failed = True
                        break

                    left_splice_site = splice_dict[left_splice]
                    right_splice_site = splice_dict[right_splice]
                    identity += str(left_splice_site) + '-' + str(right_splice_site) + '~'

                if not failed and identity in start_dicts:
                    sample=name.split('|')[0]
                    read=(name,aligned_bases,tuple(blocks))
                    start_dict=start_dicts[identity]
                    end_dict=end_dicts[identity]
                    if int(start) in start_dict and int(end) in end_dict:
                        left = start_dict[int(start)]
                        right = end_dict[int(end)]
                        new_identity = identity + '_' + read_direction + '_' + str(left) + '_' + str(right)
#                        new_identity2 = identity + '_' + opposite_direction + '_' + str(left) + '_' + str(right)
                        if new_identity in IsoDict:
                            if sample not in IsoDict[new_identity][2]:
                                IsoDict[new_identity][2][sample]=0
                            IsoDict[new_identity][2][sample]+=1

                            IsoDict[new_identity][1]+=1
                            if IsoDict[new_identity][1]<=10000:
                                IsoDict[new_identity][0].append(read)
                            else:
                                m = np.random.randint(IsoDict[new_identity][1])          # 0…t‑1
                                if m < 10000:
                                    IsoDict[new_identity][0][m] = read
#                        elif new_identity2 in IsoDict:
#                            print('opposite')

    return IsoDict

def revComp(seq):
    '''Returns the reverse complement of a seq'''
    bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    return ''.join([bases[x] for x in list(seq)])[::-1]


def parse_seqDict(seqDict,minimum_isoform_length,minimum_reads,overhangs,multi_exon_only):

    psl_dict = {}
    psl_info = {}
    isoform_list = []
    first_alignment = set()
    doneSet=set()
    filtered_isoforms=[]

    for isoform in seqDict:
        reads,abundance,quantDict=seqDict[isoform]
        chromosome,spliceSites,read_direction,readStart,readEnd = isoform.split('_')
        lengths=[]
        consensus_blocks=[]
        block_positions={}
        passed_filters=False
        if abundance >= minimum_reads:
            for name,length,blocks in reads:
                lengths.append(length)

                counter=0
                for block in blocks:
                    if counter not in block_positions:
                        block_positions[counter]=[]
                        block_positions[counter+1]=[]
                    block_positions[counter].append(block[0])
                    block_positions[counter+1].append(block[1])

                    counter+=2

            del reads

            consensus_blocks.append(np.percentile(block_positions[0], 10, method='lower'))

            if counter>2:
                for pos in range(1,counter-1):
                    consensus_blocks.append(np.percentile(block_positions[pos], 50,  method='lower'))

            consensus_blocks.append(np.percentile(block_positions[counter-1], 90, method='higher'))

            del block_positions

            readstart = consensus_blocks[0]
            readend = consensus_blocks[-1]

            readlength = np.median(lengths)
            exon_number=len(consensus_blocks)/2

            name = isoform
            if readlength >=minimum_isoform_length:
                if multi_exon_only=='0' or exon_number>1:
                    if name not in first_alignment:
                        isoform_list.append(name)
                        psl_dict[name] = [consensus_blocks, read_direction,readstart,readend,abundance,quantDict]
                else:
                    filtered_isoforms.append(name)
            else:
                filtered_isoforms.append(name)
        else:
            filtered_isoforms.append(isoform)

    return psl_dict, set(isoform_list),filtered_isoforms

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for name,seq,qual in mp.fastx_read(inFile):
        readDict[name]=seq
    return readDict


def filter_isoforms(count, isoform_names, chromosome, psl_dict, minimum_ratio):
    filtered_isoforms = []
    discard_isoforms = []
    for isoform in sorted(isoform_names):
        coverage_list = []
        info = psl_dict[isoform]
        direction = info[1]
        start = int(info[2])
        end = int(info[3])
        number = info[4]
        for base in np.arange(round(start, -1), round(end, -1), 10):
            coverage_list.append(count[chromosome + '_' + direction][base])
        max_coverage = max(coverage_list)

        if (number / max_coverage) >= minimum_ratio:
            filtered_isoforms.append(isoform)
        else:
            discard_isoforms.append(isoform)
    return filtered_isoforms,discard_isoforms



def look_for_contained_isoforms(isoform_list, chromosome, psl_dict, chr_sequence,polyAWhiteList,target_start,sw,Acutoff,internal_ratio,downstream_buffer):
    filtered_isoforms = []
    discard_isoforms = []
    bases_dict={}
    internal_buffer = 20
    coverageDict={}
    coverageDict['+']={}
    coverageDict['-']={}
    for isoform in isoform_list:
        info = psl_dict[isoform]
        status = set()
        coordinates = info[0]
        start,end = coordinates[0],coordinates[-1]
        minimum = min(coordinates[0] + internal_buffer, coordinates[1])
        maximum = max(coordinates[-1] - internal_buffer, coordinates[-2])
        direction = info[1]
        isoform_abundance=info[4]
        bases_dict[isoform]=[direction,isoform_abundance,minimum,maximum,start,end,coordinates]
        for i in range(0,len(coordinates)-1,2):
            for base in range(coordinates[i],coordinates[i+1],1):
                rounded_base=myround(base,5)
                if rounded_base not in coverageDict[direction]:
                    coverageDict[direction][rounded_base]=set()
                coverageDict[direction][rounded_base].add(isoform)

#        bases_dict[isoform][6] = {base for i in range(0, len(coordinates)-1, 2) for base in range(coordinates[i], coordinates[i+1])}


    for isoform in isoform_list:
        compatible_isoforms=set(isoform_list)
        decision=True
        direction,abundance,minimum,maximum,start,end,coordinates=bases_dict[isoform]
        if direction == '+':
            Acontent = chr_sequence[end-target_start:end-target_start + 15].upper().count('A') / 15
            polyArange = set(np.arange(end + 3, end + 23, 1))
            polyApos=end
        elif direction == '-':
            Acontent = chr_sequence[start - 15 - target_start:start-target_start].upper().count('T') / 15
            polyArange = set(np.arange(start - 23, start - 3, 1))
            polyApos=start

        if Acontent > Acutoff:
            overhangDict={}
            for base in polyArange:
                rounded_base=myround(base,5)
                if rounded_base in coverageDict[direction]:
                    for overhanging_isoform in coverageDict[direction][rounded_base]:
                        if overhanging_isoform!=isoform:
                            if overhanging_isoform not in overhangDict:
                                overhangDict[overhanging_isoform]=0
                            overhangDict[overhanging_isoform]+=1

            for overhanging_isoform,overhang_bases in overhangDict.items():
                if overhang_bases>=2:
                    if polyApos not in polyAWhiteList[direction]:
                        decision=False
                        break

        if decision==True:
            compare_bases=set()
            for i in range(0,len(coordinates)-1,2):
                for base in range(coordinates[i],coordinates[i+1],1):
                    if minimum <= base <= maximum:
                        compare_bases.add(base)
                        compatible_isoforms=compatible_isoforms&coverageDict[direction][myround(base,5)]

            for isoform2 in compatible_isoforms:
                if isoform2 != isoform:
                    direction2,abundance2,minimum2,maximum2,start2,end2,coordinates2=bases_dict[isoform2]
                    if direction == direction2:
                        if abundance / abundance2 < internal_ratio:
                            if start2<=minimum<=end2 and start2<=maximum<=end2:
                                compare_bases2=set()
                                for i in range(0,len(coordinates2)-1,2):
                                    for base in range(coordinates2[i],coordinates2[i+1],1):
                                        if minimum <= base <= maximum:
                                            compare_bases2.add(base)

                                if compare_bases == compare_bases2:
                                    decision=False
                                    break

        if  decision==True:
            filtered_isoforms.append(isoform)
        else:
            discard_isoforms.append(isoform)

    return filtered_isoforms,discard_isoforms

def readWhiteList(polyA,root):

    chromosome,dir,start,end=root.split('~')
    start,end=int(start),int(end)
    WhiteList={}
    WhiteList['+']=set()
    WhiteList['-']=set()

    for line in open(polyA):
        a=line.strip().split('\t')
        if chromosome == a[0]:
            if start<int(a[1])<end or start<int(a[2])<end:
                for pos in np.arange(int(a[1]),int(a[2]),1):
                    WhiteList[a[5]].add(pos)
    return WhiteList

def get_count(isoform_list, chromosome, psl_dict):
    count = {}
    for isoform in isoform_list:
        info = psl_dict[isoform]
        coordinates = info[0]
        direction = info[1]
        number = info[4]
        isoStart, isoEnd = coordinates[0], coordinates[-1]
        for base in np.arange(round(isoStart, -1), round(isoEnd, -1), 10):
            if (chromosome + '_' + direction) not in count:
                count[chromosome + '_' + direction] = {}
            if base not in count[chromosome + '_' + direction]:
                count[chromosome + '_' + direction][base] = 0
            count[chromosome + '_' + direction][base] += number
    return count


#def split_consensus(sequences):
