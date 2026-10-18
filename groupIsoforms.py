import sys
import argparse
import numpy
import gzip

def argParser():
    parser = argparse.ArgumentParser(
        description='Groups isoforms into loci and matches them with genes if annotation gtf is given',
        add_help=True,
        prefix_chars='-',
    )

    parser.add_argument('-p', '--temp_path', type=str, action='store')
    parser.add_argument('-g', '--genome_annotation', type=str, action='store')


    return vars(parser.parse_args())



args = argParser()
temp_path = args['temp_path']
infile = f'{temp_path}/Isoform_Coordinates.sorted.tsv'
gtffile = f'{temp_path}/Isoforms.filtered.clean.gtf'
outfile = f'{temp_path}/Isoforms.filtered.clean.genes'
rpmsfile = f'{temp_path}/Isoforms.filtered.clean.rpms'
quantfile = f'{temp_path}/Isoforms.filtered.clean.quant'
rpmsfile_genes = f'{temp_path}/Isoforms.filtered.clean.genes.rpms'
quantfile_genes = f'{temp_path}/Isoforms.filtered.clean.genes.quant'
genome_annotation = args['genome_annotation']

out=open(outfile,'w')

def read_annotation(genome_annotation):
    geneDict={}
    geneDict['+']={}
    geneDict['-']={}
    coordDict={}
    coordDict['+']={}
    coordDict['-']={}
    if genome_annotation == 'None':
        return coordDict

    else:
        if genome_annotation.endswith('.gtf.gz'):
            print('\t\t\tgtf file ends on .gz and will be treated as gzipped')
            input=gzip.open(genome_annotation,'rt')
        elif genome_annotation.endswith('.gtf'):
            input=open(genome_annotation,'r')

        with input as f:
            for line in f:
                if line[0]!='#':
                    a=line.strip().split('\t')
                    chromosome=a[0]
                    element=a[2]
                    left=int(a[3])-1
                    right=int(a[4])
                    direction=a[6]
                    gene_name=a[8].split('gene_id "')[1].split('"')[0]
                    if element == 'exon':
                        if 'gene_name' in a[8]:
                            gene_name+='~'+a[8].split('gene_name "')[1].split('"')[0]
                        else:
                            gene_name+='~'+gene_name
                        if gene_name not in geneDict[direction]:
                            geneDict[direction][gene_name]=[chromosome,set()]
                        geneDict[direction][gene_name][1].add((left,right))


        return geneDict


def group_isoforms(infile,geneDict):
    outDict={}
    isoforms=[]
    new=False
    previous_chrom=''
    previous_start=0
    previous_end=0
    roots=set()
    locus=0
    i2g={}
    covered=set()
    for direction in ['+','-']:
        print('\t\t\tgrouping isoforms  on',direction,'strand')
        outDict={}
        isoforms=[]
        new=False
        new_chrom=False
        previous_chrom=''
        previous_start=0
        previous_end=0
        roots=set()
        isoform_count=0
        coordDict={}
        for line in open(infile):
            a=line.strip().split('\t')
            isodirection=a[5]
            if direction==isodirection:
                isoform_count+=1

        current=0
        for line in open(infile):
            a=line.strip().split('\t')
            isodirection=a[5]
            if direction==isodirection:
                current+=1
                print('\t\t\t'+str(current),'of',isoform_count,str(round((current/isoform_count)*100,2))+'%'+' '*40,end='\r')

                chrom=a[1]
                start=int(a[6])
                end=int(a[-1])

                if chrom!=previous_chrom:
                    new=True
                    new_chrom=True
                else:
                    if start>previous_end:
                        new=True
                    else:
                        previous_end=max(end,previous_end)
                        isoforms.append(line)
                if new:
                    if isoforms:
                        locus,i2g = match_isoforms(isoforms,previous_chrom,previous_start,previous_end,direction,locus,i2g,coordDict)
                    isoforms=[]
                    isoforms.append(line)
                    previous_chrom=chrom
                    previous_end=end
                    previous_start=start
                    if new_chrom:
                        coordDict={}
                        coordDict[direction]={}
                        for gene,coords in geneDict[direction].items():
                            chromosome=coords[0]
                            if chromosome==previous_chrom:
                                for exon in coords[1]:
                                    start=exon[0]
                                    end=exon[1]
                                    if chromosome not in coordDict[direction]:
                                        coordDict[direction][chromosome]={}
                                    for i in range(start,end,2):
                                        if i not in coordDict[direction][chromosome]:
                                            coordDict[direction][chromosome][i]=set()
                                        coordDict[direction][chromosome][i].add(gene)

                    new=False
                    new_chrom=False
        if isoforms:
            locus,i2g = match_isoforms(isoforms,previous_chrom,previous_start,previous_end,direction,locus,i2g,coordDict)
        print('\n')
    return i2g

def match_isoforms(isoforms,previous_chrom,previous_start,previous_end,direction,locus,i2g,coordDict):
    genes={}
    covered=set()
    for line in isoforms:
        a=line.strip().split('\t')
        coords=numpy.array(a[6:],dtype=int)
        for i in range(0,len(coords)-1,2):
            blockstart = coords[i]
            blockend = coords[i+1]
            for i in range(blockstart,blockend,1):
                covered.add(i)


    for i in covered:
        if previous_chrom in coordDict[direction]:
            if i in coordDict[direction][previous_chrom]:
                for gene in coordDict[direction][previous_chrom][i]:
                    if gene not in genes:
                        genes[gene]=0
                    genes[gene]+=1

    geneList=[]
    geneSet=set()
    for gene,count in genes.items():
        geneList.append((count,gene))
        geneSet.add(gene)
    coords=[0]

    if geneList:
        best=sorted(geneList,reverse=True)[0][1]
    else:
        best=''
    locus+=1
    LocusName='Locus'+str(locus)
    if geneSet:
        GeneOverlaps=(',').join(geneSet)
    else:
        GeneOverlaps=''
    for line in isoforms:
        name=line.strip().split('\t')[0]
        if best:
            i2g[name]=best
        else:
            i2g[name]=LocusName+'~'+LocusName
        out.write(f'{name}\t{LocusName}\t{previous_chrom}\t{previous_start}\t{previous_end}\t{best}\t{GeneOverlaps}\n')
    return locus,i2g

def psl_to_gtf(psl_file,gtf_file,i2g):
    '''
    converts a psl file to a gtf file. It also adjusts the start position by +1 to but leave the end position of features alone. \
    This way it accomodates format definitions and coordinate systems.
    '''
    out=[]
    doneDict=set()
    for line in open(psl_file):
        a = line.strip().split('\t')
        direction, name, chromosome, start, end=a[5], a[0], a[1], int(a[6]), int(a[-1])
        coords=numpy.array(a[6:],dtype=int)
        gene_id,gene_name=i2g[name].split('~')
        out_tmp=[]
        out_tmp.append(f'{chromosome}\tIG28\ttranscript\t{start+1}\t{end}\t.\t{direction}\t.\ttranscript_id "{name}"; gene_id "{gene_id}"; gene_name "{gene_name}"\n')

        for i in range(0,len(coords)-1,2):
            blockstart = coords[i]
            blockend = coords[i+1]
            out_tmp.append(f'{chromosome}\tIG28\texon\t{blockstart+1}\t{blockend}\t.\t{direction}\t.\ttranscript_id "{name}"; gene_id "{gene_id}"; gene_name "{gene_name}"\n')
        out.append(out_tmp)
    gtf_handle=open(gtf_file,'w')
    for transcript in out:
        for feature in transcript:
            gtf_handle.write(feature)
    gtf_handle.close()

def add_gene(inFile,outFile,i2g):
    out=open(outFile,'w')
    with open(inFile) as f:
        line=f.readline()
        out.write('gene\t'+line)
        for line in f:
            a=line.strip().split('\t')
            name = a[0]
            gene=i2g[name]
            out.write(gene+'\t'+line)
    out.close()

print('\t\tparsing genome annotation')
geneDict=read_annotation(genome_annotation)
print('\t\tprocessing isoforms')
i2g=group_isoforms(infile,geneDict)
print('\tconverting psl to gtf')
psl_to_gtf(infile,gtffile,i2g)
print(f'\tadding genes to quantFiles')
add_gene(quantfile, quantfile_genes,i2g)
add_gene(rpmsfile,rpmsfile_genes,i2g)
