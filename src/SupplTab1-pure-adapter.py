from sys import stdin 

trimmed_tab = {}

step=0

def reverse(seq):
	ret=''
	for b in seq:
		if b == 'N': ret='N'+ret
		if b == 'A': ret='T'+ret
		if b == 'C': ret='G'+ret
		if b == 'G': ret='C'+ret
		if b == 'T': ret='A'+ret
	return ret

missing_tseq = 0
missing_rname = 0
ambiguous = 0
for fl in stdin:
	if "PN:subread" in fl:
		step+=1 
		continue
	if fl[0]=='@': continue

	fli=fl.strip().split()
	rname=fli[0]
	flags=int(fli[1])
	seq=fli[9]
	cigar=fli[5]

	is_R2 = (flags & 0x80)>0

	# step1 = trimmed BAM
	# step2 = raw BAM 
	if step == 1:
		trimmed_tab[(rname, is_R2)] = (flags, cigar, seq)
	elif step == 2:
		if (rname, is_R2) in trimmed_tab:
			tflags, tcigar, tseq = trimmed_tab[(rname, is_R2)]
			if (tflags & 16) != (flags & 16): tseq = reverse(tseq)
			if not ( tseq in seq ): missing_tseq+=1
			trimseq_start = -1
			trimseq_hits = 0
			for ii in range( len(seq)-len(tseq)+1 ):
				assumed_trimseq = seq[ii: (ii+len(tseq)) ]
				if assumed_trimseq == tseq:
					trimseq_hits+=1
					trimseq_start=ii

			trimmed = []
			if trimseq_hits ==1:
				if trimseq_start > 0: trimmed.append([ 0 , trimseq_start ])
				if trimseq_start + len(tseq) < len(seq): trimmed.append([ trimseq_start + len(tseq), len(seq)])
			else:
				ambiguous +=1
				trimmed=[[99999999,99999999]]
		else:
			missing_rname +=1
			trimmed = [[0, len(seq)]]
		print( "%s/%s/%d\t"%(rname, is_R2, flags) + ("\t".join( [ "%d-%d"%(x[0],x[1]) for x in trimmed ] )))
	elif step > 2: raise("step > 2 ")
	

print("## Missing_Tseq=%d  Missing_Rname=%d  Ambiguous_Pos=%d"% (missing_tseq, missing_rname, ambiguous))
