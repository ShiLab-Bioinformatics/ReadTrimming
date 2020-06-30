from sys import argv

rawbam=argv[1]
tmbam=argv[2]
pureadapter=argv[3]

puretab = {}
for fl in open(pureadapter):
	if fl[0]=='#':continue

	fli=fl.strip().split("\t")
	(rname, is_R2_str, flag_str) = fli[0].split("/")
	is_R2 = is_R2_str=="True"
	flag = int(flag_str)
	purerange = [ (int(xx.split('-')[0]), int(xx.split('-')[1])) for xx in fli[1:] ]
	puretab[ (rname, is_R2) ] = ( flag, purerange )

def calc_uniq(Rl, Rc, Tl, Tp, Pa):
	Rgaps=[]
	Rgap_len=0
	tmpi = 0
	Rcur = 0
	for ch in Rc:
		if ch.isdigit():
			tmpi = tmpi*10+int(ch)
		else:
			if ch=='S':
				Rgaps.append( [ Rcur, Rcur+tmpi ] )
				Rgap_len += tmpi
			if ch in ('S','M','I'): Rcur += tmpi
			tmpi = 0

	Tgaps=[]
	Tgap_len=0
	if Tp>0:
		Tgap_len += Tp
		Tgaps=[ [0,Tp] ]
	if Tp + Tl < Rl:
		Tgaps.append( [ Tp+Tl, Rl ] )
		Tgap_len += Rl - (Tp+Tl)

	common = 0
	common_are_adaptor = 0
	for Rgap in Rgaps:
		for Tgap in Tgaps:
			common_end = min(Rgap[1], Tgap[1])
			common_start = max( Rgap[0], Tgap[0] )
			overlap = common_end - common_start

			if overlap>0:
				common+= overlap
				for pure_gap in Pa:
					coco_start = max(pure_gap[0], common_start )
					coco_end = min(pure_gap[1], common_end )
					coco_overlap = coco_end-coco_start
					if coco_overlap>0:
						common_are_adaptor += coco_overlap

	ada_unique = sum( [ aa[1]-aa[0] for aa in Pa ]) - common_are_adaptor
	return ( Rgap_len - common, common, Tgap_len - common, common_are_adaptor, ada_unique ) 




trimbam_tab = {}

def reverse(seq):
	ret=''
	for b in seq:
		if b == 'N': ret='N'+ret
		if b == 'A': ret='T'+ret
		if b == 'C': ret='G'+ret
		if b == 'G': ret='C'+ret
		if b == 'T': ret='A'+ret
	return ret

Softclipped_bases_unique = 0
Common_Clipped_Trimmed = 0
Trimmed_bases_unique = 0
for fl in open(tmbam):
	fli=fl.split("\t")
	if fli[0][0]=='@': continue
	rname=fli[0]
	flags=int(fli[1])
	cigar=fli[5]

	if flags & 2048 >0: continue
	assert(flags & 256 == 0)
	if (flags & 1)>0: is_R2 = (flags & 128)>0
	else:
		# DO NOT ALLOW single-end reads
		assert(0)

	seq = fli[9]
	trimbam_tab[(rname, is_R2)] = ( flags, cigar, seq )

seq="!!!!!!!!"

Total_Raw_Bases = 0
rawreads = 0
ambiguous_adapters = 0
rawreads_mapped = 0
unmatched_reads = 0
for_debug=0
Common_CliTri_also_Adapter = 0
Adapter_but_no_CommonCliTri = 0

for fl in open(rawbam):
	fli=fl.split("\t")
	rname=fli[0]
	rawflags=int(fli[1])
	is_R2 = (rawflags & 128)>0
	rawcigar = fli[5]
	rawseq = fli[9]
	rawreads +=1
	Total_Raw_Bases += len(rawseq)

	if (rname, is_R2) in trimbam_tab:
		trimflags, trimcigar, trimseq = trimbam_tab[(rname, is_R2)] 
	else:
		trimflags, trimcigar, trimseq = (1+(128 if is_R2 else 64), "0M","")

	if (trimflags & 16) != (rawflags & 16):
		trimseq=reverse(trimseq)

	matches = 0
	inraw_pos = -1

	if trimseq:
		for test_start in range( len(rawseq) - len(trimseq) + 1):
			if rawseq[ test_start:test_start+len(trimseq) ] == trimseq:
				matches+=1
				inraw_pos = test_start

		if 0 == matches:
			if rawflags & 4 == 0: print(rname)
			assert(rawflags & 4 == 4)
			unmatched_reads +=1
			trimseq=reverse(trimseq)
			for test_start in range( len(rawseq) - len(trimseq) + 1):
				if rawseq[ test_start:test_start+len(trimseq) ] == trimseq:
					matches+=1
					inraw_pos = test_start
	else:
		matches = 1
		inraw_pos = 0

	if(matches<1):
		print( "## UMR="+rname )
	#	print( trimseq )
	#	print( rawseq )
	else:
		#print (rname, is_R2)
		puflag, puranges = puretab[ (rname, is_R2) ]
		is_ambiguous_adapter = 0
		if puranges and puranges[0][0]>99999:
			is_ambiguous_adapter = 1
			ambiguous_adapters +=1
	
		if (puflag & 16) != ( rawflags & 16 ):
			pu2 = []
			for ra in puranges:
				pu2 = [ [ len(rawseq) - ra[1], len(rawseq) - ra[0]] ] + pu2
			puranges = pu2


		if for_debug:
			print( "TS="+trimseq )
			print( "RS="+rawseq )
			print( "RCigar="+rawcigar )
			print( "PA="+(", ".join( [ "%d-%d"%(a[0],a[1]) for a in puranges ] )) )
			print( "RL=%d, TL=%d, TPos=%d"%(  len(rawseq), len(trimseq) , inraw_pos ) )

		if rawflags &4 ==0: rawreads_mapped+=1
		if matches==1 and rawflags & 4 ==0 and not is_ambiguous_adapter:
			raw_uniq, common, truniq, common_in_adapter, adapter_uniq = calc_uniq( len(rawseq), rawcigar, len(trimseq), inraw_pos, puranges )
			if for_debug: print("RQ, CO, TQ=%d, %d, %d  ;  CO_ADA, AU=%d, %d"%(raw_uniq, common, truniq, common_in_adapter, adapter_uniq))
			Softclipped_bases_unique += raw_uniq
			Common_Clipped_Trimmed += common
			Trimmed_bases_unique += truniq
			Common_CliTri_also_Adapter += common_in_adapter
			Adapter_but_no_CommonCliTri += adapter_uniq

		if for_debug:
			print()
			print()

#print("%d\t%d\t%d\t%d\t%d\t%d\tUnmatch_recovered=%d\tTrimmedReads=%d\tRawReads=%d\tRawMapped=%d\tAmbiguousAdapterReads=%d" %(Softclipped_bases_unique, Common_Clipped_Trimmed, Trimmed_bases_unique, Common_Clipped_Trimmed - Common_CliTri_also_Adapter, Common_CliTri_also_Adapter, Adapter_but_no_CommonCliTri, unmatched_reads, len(trimbam_tab), rawreads , rawreads_mapped, ambiguous_adapters))

Trimmed_base_proportion=(Common_Clipped_Trimmed+Trimmed_bases_unique)*1./Total_Raw_Bases
if Common_Clipped_Trimmed+Trimmed_bases_unique>0:
	Softclipped_in_trimmed=Common_Clipped_Trimmed*1./(Common_Clipped_Trimmed+Trimmed_bases_unique)
else:
	Softclipped_in_trimmed=-9999
if Common_Clipped_Trimmed>0:
	Adapter_in_common_CliTri=Common_CliTri_also_Adapter*1./Common_Clipped_Trimmed
else:
	Adapter_in_common_CliTri=-9999

print( "%.5f,%.5f,%.5f;"%( Trimmed_base_proportion, Softclipped_in_trimmed, Adapter_in_common_CliTri))
#print( "%.5f,%.5f,%.5f::%d::%d;"%( Trimmed_base_proportion, Softclipped_in_trimmed, Adapter_in_common_CliTri , Common_CliTri_also_Adapter, Common_Clipped_Trimmed))
