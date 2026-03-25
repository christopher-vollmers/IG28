#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import argparse
import re
import os
import numpy as np
from os.path import isfile
import multiprocessing as mp
import mappy
import time
import pysam

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/utils/'
sys.path.append(os.path.abspath(PATH))

import ModuleD as ModuleD

parser = argparse.ArgumentParser()

parser.add_argument('--path', '-p', type=str, action='store')
parser.add_argument('--roots', '-r', type=str, action='store')
parser.add_argument('--genome_sequence', '-G', type=str)
parser.add_argument('--genome_file', '-g', type=str, action='store')
parser.add_argument('--splice_site_width', '-w', type=int, action='store')
parser.add_argument('--minimum_feature_count', '-m', type=int, action='store')
parser.add_argument('--white_list_polyA', '-W', type=str, action='store')
parser.add_argument('--numThreads', '-n', type=str, action='store')
parser.add_argument('--junctions', '-j', type=str, action='store')
parser.add_argument('--upstream_buffer', '-u', type=str, action='store')
parser.add_argument('--downstream_buffer', '-d', type=str, action='store')
parser.add_argument('--minimum_read_count', '-M', type=int, action='store')
parser.add_argument('--internal_ratio','-ir', type=float, action='store')
parser.add_argument('--minimum_ratio', '-mr', type=float)
parser.add_argument('--overhangs', '-O', type=str)
parser.add_argument('--Acutoff', '-A', type=str)
parser.add_argument('--minimum_isoform_length', '-I', type=str)
parser.add_argument('--multi_exon_only', '-meo', type=str)




if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()
out_path = args.path + '/'  # path where you want your output files to go
root_file = args.roots
genome_file = args.genome_file
genome=args.genome_sequence
splice_site_width = int(args.splice_site_width)
minimum_feature_count = int(args.minimum_feature_count)
minimum_read_count = int(args.minimum_read_count)
white_list_polyA=args.white_list_polyA.split(',')
threads=int(args.numThreads)
junctions=args.junctions.split(',')
upstream_buffer=int(args.upstream_buffer)
downstream_buffer=int(args.downstream_buffer)
internal_ratio=float(args.internal_ratio)
minimum_ratio=float(args.minimum_ratio)
overhangs=np.array(args.overhangs.split(','), dtype=int)
Acutoff=float(args.Acutoff)
minimum_isoform_length=int(args.minimum_isoform_length)
multi_exon_only=args.multi_exon_only


polyAWhiteListFile=out_path + '/polyAWhiteList.bed'

def process_locus(out_tmp,root,chrom,left_bounds_chrom, right_bounds_chrom,start,end,splice_site_width,minimum_read_count,junctions,verbose,minimum_feature_count,locus_sequence,polyA,overhangs,minimum_isoform_length,internal_ratio,multi_exon_only,Acutoff,downstream_buffer,out_path_coord):
    infile=out_tmp+'/'+root+'.bam'
    pysam.index(infile)
#    SpliceDefineConsensus.call_variants(infile,genome,root,out_tmp+'/')
#    infile=out_tmp+'/'+root+'.tagged.bam'
#    print(infile)
#    verbose=True
    peak_areas={}
    peak_areas['l'] = {}
    peak_areas['r'] = {}
    filter_reasons=[]
    WhiteList = ModuleD.readWhiteList(polyA,root)

    if verbose:
        print(f'\t\t\t\tcollecting annotated splice sites in locus {chrom} {start} {end}')
    peak_areas, toWrite_A_l = ModuleD.make_genome_bins(left_bounds_chrom, 'l', chrom, peak_areas,splice_site_width)
    peak_areas, toWrite_A_r = ModuleD.make_genome_bins(right_bounds_chrom, 'r', chrom, peak_areas,splice_site_width)
    if verbose:
        print(f'\t\t\t\tfinding high-confidence unannotated splice sites in locus {chrom} {start} {end}')


    histo_left_bases,lineCounter = ModuleD.collect_reads(infile, peak_areas,'l')
    peak_areas, toWrite_N_l = ModuleD.find_peaks(histo_left_bases, True, 'l', peak_areas, chrom,start,end,splice_site_width,minimum_feature_count,junctions,infile)
    del histo_left_bases

    histo_right_bases,lineCounter = ModuleD.collect_reads(infile, peak_areas,'r')
    peak_areas, toWrite_N_r = ModuleD.find_peaks(histo_right_bases, False, 'r', peak_areas, chrom,start,end,splice_site_width,minimum_feature_count,junctions,infile)
    del histo_right_bases

    if verbose:
        print(f'\t\t\tprocessing locus {chrom} {start} {end} covering {end-start:,} genomic bases and {lineCounter:,} sequencing reads',' '*20)

    peakCounter={}
    peakCounter['l']=0
    peakCounter['r']=0
    spliceDict={}
    for toWrite in [toWrite_A_l,toWrite_A_r,toWrite_N_l,toWrite_N_r]:
       for sChrom,sStart,sEnd,type1,side,prop in toWrite:
            peakCounter[side]+=1
            peaks=str(peakCounter[side])
            splice_left = int(sStart)
            splice_right = int(sEnd)
#            print(splice_left,type1 + side + peaks,prop)
            for base in np.arange(splice_left, splice_right + 1):
                spliceDict[base] = type1 + side + peaks

#    obj_counter=1
#    for object in [histo_left_bases, histo_right_bases, peak_areas, toWrite_A_l, toWrite_A_r, toWrite_N_l, toWrite_N_r]:
#        print(f'{obj_counter} {asizeof.asizeof(object)/ (1024**2):,}')
#        obj_counter+=1
    del peak_areas, toWrite_A_l, toWrite_A_r, toWrite_N_l, toWrite_N_r

    if verbose:
        print(f'\t\t\t\tsorting read alignments into splice junctions for locus {chrom} {start} {end}')
    start_end_dict = ModuleD.sort_reads_into_splice_junctions(spliceDict, infile, multi_exon_only,overhangs)

    if verbose:
        print(f'\t\t\t\tfinding TSS and polyA sites for {len(start_end_dict):,} splice junctions chains for locus {chrom} {start} {end}')
    seqDict=ModuleD.define_start_end_sites(start_end_dict,upstream_buffer,downstream_buffer,minimum_feature_count,spliceDict,infile,overhangs,multi_exon_only)

    del start_end_dict

    if verbose:
        print(f'\t\t\t\tstarted filtering {len(seqDict)} isoforms for locus {chrom} {start} {end}')
    psl_dict, isoform_list,filter_reasons,filtered_isoforms = ModuleD.parse_seqDict(seqDict,minimum_isoform_length,minimum_read_count,overhangs,multi_exon_only,filter_reasons)

    del seqDict

    if verbose:
        print(f'\t\t\t\t\tfiltered down to {len(isoform_list)} isoforms with length,abs abundance, and readlength for locus {chrom} {start} {end}')
    count = ModuleD.get_count(isoform_list,chrom, psl_dict)

    if verbose:
        print(f'\t\t\t\t\tfiltering isoforms for relative read coverage starting with', len(isoform_list), f'isoforms for locus {chrom} {start} {end}')
    isoform_list,filter_reasons,filtered_isoforms = ModuleD.filter_isoforms(count, isoform_list, chrom, psl_dict, minimum_ratio,filter_reasons)
    for filtered_isoform in filtered_isoforms:
        psl_dict.pop(filtered_isoform)

    if verbose:
        print(f'\t\t\t\t\tfinding fully contained isoforms in', len(isoform_list), f'remaining isoforms for locus {chrom} {start} {end}')
    isoform_list,filter_reasons,filtered_isoforms = ModuleD.look_for_contained_isoforms(isoform_list, chrom, psl_dict, locus_sequence,WhiteList,start,splice_site_width,Acutoff,internal_ratio,downstream_buffer,filter_reasons)
    for filtered_isoform in filtered_isoforms:
        psl_dict.pop(filtered_isoform)

    if verbose:
        print(f'\t\t\t\tended on {len(isoform_list)} filtered isoforms for locus {chrom} {start} {end}')

#    for object in [psl_dict, count, locus_sequence, WhiteList]:
#        print(f'{obj_counter} {asizeof.asizeof(object)/ (1024**2):,}')
#        obj_counter+=1


    IsoData={}
    counter=0
    if len(isoform_list)>0:
        out_coords=open(f'{out_path_coord}/{root}.coords','w')
        out_counts=open(f'{out_path_coord}/{root}.counts','w')
        for isoform in isoform_list:
            counter+=1
            info = psl_dict[isoform]
            coordinates = info[0]
            dir=info[1]
            abundance=info[4]
            quantDict=info[5]

            spaced_root = "/".join(root.split("~"))
            nameString=f'IG28_Locus={spaced_root}_IsoformNumber={counter}_Coverage={abundance}'
            parsed_root = "\t".join(root.split("~"))
            parsed_consensus = "\t".join(map(str,coordinates))
            out_coords.write(f'{nameString}\t{parsed_root}\t{dir}\t{parsed_consensus}\n')
            out_counts.write(f'{nameString}\t')
            quantList=[]
            for sample,count in quantDict.items():
               out_counts.write(f'{sample}~{count}\t')
            out_counts.write('\n')
        out_coords.close()
        out_counts.close()

def main():

    chromSet=set()
    for root_line in open(root_file):
        root,root_sequence=root_line.strip().split('\t')
        chrom,dir,start,end=root.split('~')
        chromSet.add(chrom)

    out_tmp = out_path + '/tmp_SS'
    left_bounds, right_bounds = {}, {}
    if genome_file!='None':
        if genome_file.endswith('.gtf.gz') or genome_file.endswith('.gtf'):
            print('\t\tparsing annotated splice sites')
            chrom_list, left_bounds, right_bounds, polyAWhiteList = ModuleD.parse_genome(genome_file, left_bounds, right_bounds, white_list_polyA,chromSet)
        else:
            print('\t\tgenome annotation file does not end with .gtf or .gtf.gz. File will be ignored. Splice sites will be entirely read derived and no polyA sites will be white-listed')
            chrom_list=set()
            polyAWhiteList=[]
    else:
        print('\t\tNo genome annotation provided, so splice sites will be entirely read derived and no polyA sites will be white-listed')
        chrom_list=set()
        polyAWhiteList=[]

    outPolyA=open(out_path+'/polyAWhiteList.bed','w')

    if '0' not in white_list_polyA:
        print('\t\t'+str(len(polyAWhiteList)), 'poly(A) sites whitelisted')
        for chrom,direction,end,transcript_id in polyAWhiteList:
            polyA=int(end)
            polyAstart=polyA-20
            polyAend=polyA+20
            outPolyA.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,str(polyAstart),str(polyAend),transcript_id,'0',direction))
    outPolyA.close()


    numThreads=threads


    SampleTotalDict={}
    sampleSet=set()
    for line in open(out_path+'/fastaFiles.info'):
        a=line.strip().split('\t')
        sample=a[1]
        sampleSet.add(sample)
        count=int(a[3])
        if sample not in SampleTotalDict:
            SampleTotalDict[sample]=0
        SampleTotalDict[sample]+=count

    sampleList=sorted(list(sampleSet))
    out_path_coord = f'{out_path}/tmp_coords/'
    if  os.path.isdir(out_path_coord):
       os.system(f'rm -r {out_path_coord}')
    os.system(f'mkdir {out_path_coord}')
    out=open(out_path+'/Isoform_Coordinates.tsv','a')
    out_filter=open(out_path+'/filter_reasons.txt','a')
    out_r2i=open(out_path+'/reads2isoforms.txt','a')
    out_i2s=open(out_path+'/Isoforms.filtered.clean.quant','a')
    out_i2sn=open(out_path+'/Isoforms.filtered.clean.rpms','a')
    out_i2l=open(out_path+'/isoform2locus.txt','a')



    results={}
    print(f'\t\tstarting multithreaded processing of loci in batch')
    mp.set_start_method('fork', force=True)
    pool = mp.Pool(numThreads)
    locus_counter=0
    for root_line in open(root_file):
        root,root_sequence=root_line.strip().split('\t')
        locus_counter+=1
        chrom,dir,start,end=root.split('~')
        start=int(start)
        end=int(end)
        left_bounds_sub,right_bounds_sub=ModuleD.prepare_locus(chrom,start,end,left_bounds,right_bounds)
        pool.apply_async(process_locus,[out_tmp,root,chrom,left_bounds_sub, right_bounds_sub,start,end,splice_site_width,minimum_read_count,junctions,False,minimum_feature_count,root_sequence,polyAWhiteListFile,overhangs,minimum_isoform_length,internal_ratio,multi_exon_only,Acutoff,downstream_buffer,out_path_coord])
#        process_locus(out_tmp,root,chrom,left_bounds_sub, right_bounds_sub,start,end,splice_site_width,minimum_read_count,junctions,True,minimum_feature_count,root_sequence,polyAWhiteListFile,overhangs,minimum_isoform_length,internal_ratio,multi_exon_only,Acutoff,downstream_buffer,out_path_coord)
        print(f'\t\tloading locus {locus_counter:,} -- {root}',' '*20,end='\r')


    print(f'\n\t\t{len(results):,} loci submitted for mutlithreading')
    pool.close()
    pool.join()

    for file in os.listdir(out_path_coord):
        if file.endswith('.coords'):
            for line in open(f'{out_path_coord}/{file}'):
                out.write(line)

    for file in os.listdir(out_path_coord):
        if file.endswith('.counts'):
            for line in open(f'{out_path_coord}/{file}'):
                quantDict={}
                a=line.strip().split('\t')
                quants=a[1:]
                for entry in quants:
                    sample,count = entry.split('~')
                    quantDict[sample]=int(count)
                sampleDict={}
                for sample in sampleList:
                    if sample in quantDict:
                        sampleDict[sample]=quantDict[sample]
                    else:
                        sampleDict[sample]=0

                sampleCounts=[]
                sampleRPMS=[]
                for sample in sampleList:
                    sampleCount=sampleDict[sample]
                    sampleCounts.append(str(sampleCount))
                    sampleRPMS.append(str(round(sampleCount/SampleTotalDict[sample]*1000000,3)))
                out_i2s.write(f'{a[0]}\t'+'\t'.join(sampleCounts)+'\n')
                out_i2sn.write(f'{a[0]}\t'+'\t'.join(sampleRPMS)+'\n')



    out.close()
    out_r2i.close()
    out_i2s.close()
    out_i2l.close()

main()
