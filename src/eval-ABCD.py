from sys import argv
from math import log
from numpy import corrcoef

mode=argv[1]
count_file=argv[2]
truth_file=argv[3]
length_file="Gene-lengths.txt"

gene_length_table={}
for fl in open(length_file):
	fli=fl.split()
	gene_length_table[fli[0]]=int(fli[1])

dexseq_exon_table={}
if 'dexseq' in count_file:
	for fl in open("hg38-DEXseq.GTF"):
		if "exonic_part_number" in fl:
			fli=fl.split("\t")
			chro=fli[0]
			start=int(fli[3])
			stop=int(fli[4])
			geneid=fl.split('gene_id "')[1].split('"')[0]
			exonid=fl.split('exonic_part_number "')[1].split('"')[0]
			exonkey="%s:%s"%(geneid,exonid)
			dexseq_exon_table[exonkey]=(( chro, start, stop ))

count_table={}
for fl in open(count_file):
	if "start" in fl.lower() or "bam" in fl.lower() or "gene" in fl.lower(): continue # header line
	if fl[0]=="_": continue # stat lines
	fli=fl.split()

	if "exons" in mode:
		if 'dexseq' in count_file:
			chro, start, stop = dexseq_exon_table[fli[0]]
			key="%s:%d:%d"%(chro, start, stop)
			count=int(fli[1])
		else:
			chro=fli[0]
			if 'chrM' in chro or 'NT_' in chro or 'NW_' in chro:continue
			key=":".join(fli[0:3])
			count=int(fli[3])
	else:
		key=fli[0]
		if key not in gene_length_table:continue
		count=float(fli[1])
	count_table[key]=count

truth_L2RPKM_table={}
truth_count_table={}
repeated_genes={}
rt_pcr_genes = {}
rt_pcr_nas = 0
missed_in_len_tab = {} 

for fl in open(truth_file):
	if 'A1_value' in fl or "Fragment" in fl:continue
	fli=fl.split()
	if 'RTPCR' in truth_file:
		geneid=fli[0]
		if "NA" == geneid:
			rt_pcr_nas +=1
		else:
			rt_pcr_genes[geneid] = 1 
			exp_a1 = sum([ log(float(fli[ii]))/log(2) for ii in (2,4,6,8) ])/4
			exp_b1 = sum([ log(float(fli[ii]))/log(2) for ii in (10,12,14,16) ])/4
			exp_c1 = sum([ log(float(fli[ii]))/log(2) for ii in (18,20,22,24) ])/4
			exp_d1 = sum([ log(float(fli[ii]))/log(2) for ii in (26,28,30,32) ])/4
			if "UHRR" in mode: truth_L2RPKM_table[geneid]=exp_a1
			elif "HBRR" in mode: truth_L2RPKM_table[geneid]=exp_b1
			elif "MIXC" in mode: truth_L2RPKM_table[geneid]=exp_c1
			elif "MIXD" in mode: truth_L2RPKM_table[geneid]=exp_d1
			else:raise( "ERROR : NO MODE: "+ mode )
		if geneid not in gene_length_table: missed_in_len_tab[geneid]=1
		if geneid not in repeated_genes: repeated_genes[geneid]=1
		else: repeated_genes[geneid]+=1
	else:
		if len(fli) < 3 and "genes" in mode:
			geneid=fli[0]
			if geneid in gene_length_table: truth_count_table[geneid]=int(fli[1])
		if len(fli) > 3 and "exons" in mode:
			truth_count_table[":".join(fli[1:4])]=int(fli[4])

if 'RTPCR' in truth_file:
	for gene, rep in repeated_genes.items():
		if rep>1 and gene in truth_L2RPKM_table: del truth_L2RPKM_table[gene]

if 'RTPCR' in truth_file:
	used_keys = [ky for ky in truth_L2RPKM_table.keys() if ky in gene_length_table]
else:
	if 'genes' in mode:
		used_keys = list(set(list(count_table.keys())+list(truth_count_table.keys())))
	else:
		used_keys = list(count_table.keys())

log2_offset=0.5+0.5 

if 'RTPCR' not in truth_file:
	all_reads=sum(truth_count_table.values())
	for key in used_keys:
		count=truth_count_table[key] if key in truth_count_table else 0
		feature_length=gene_length_table[key] if "genes" in mode else (int(key.split(":")[2]) - int(key.split(":")[1])+1)
		log2RPKM=log((log2_offset+count)*1000/feature_length*1000000/(log2_offset+all_reads)) / log(2)
		truth_L2RPKM_table[key]=log2RPKM

pipeline_L2RPKM_table = {}
all_reads=sum(count_table.values())

for key in used_keys:
	count=count_table[key] if key in count_table else 0
	feature_length=gene_length_table[key] if "genes" in mode else (int(key.split(":")[2]) - int(key.split(":")[1])+1)
	pipeline_L2RPKM_table[key]=log((log2_offset+count)*1000/feature_length*1000000/(log2_offset+all_reads))/log(2)

truth_L2RPKM_vec=[ truth_L2RPKM_table[ky] for ky in used_keys ]
pipeline_L2RPKM_vec=[ pipeline_L2RPKM_table[ky] for ky in used_keys ]

if False:
	for ky in used_keys:
		tcnt = truth_count_table[ky] if ky  in truth_count_table else 0
		ccnt = count_table[ky] if ky in count_table else 0
		print("DEBUG_CNTS_TH_FC\t%s\t%d\t%d\t%d"%( ky, gene_length_table[ky], tcnt, ccnt))

	for ky in used_keys:
		print("DEBUG_L2RPKM_TH_FC\t%s\t%.3f\t%.3f"%( ky, truth_L2RPKM_table[ky], pipeline_L2RPKM_table[ky] ))

corr=corrcoef(truth_L2RPKM_vec, pipeline_L2RPKM_vec)[0][1]
print("Pearson correlation between the known truth and %s (log2 RPKM) is %.7f ; using RT-PCR data with %d NAs and %d meaningful genes. Overlapping %d genes"%(count_file, corr, rt_pcr_nas, len(rt_pcr_genes), len(used_keys)))
