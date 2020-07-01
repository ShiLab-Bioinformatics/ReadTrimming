from sys import argv

tool1=argv[1]
tool2=argv[2]

tool1_res = {}
tool2_res = {}
all_keys = {}

def loadsam(f):
	global all_keys
	ret={}
	for fl in open(f):
		if fl[0]=='@': continue
		fli=fl.split("\t")
		rname=fli[0]
		flag=int(fli[1])
		if flag & 2048 > 0 : continue
		assert(flag & 256 == 0)

		chro = fli[2]
		pos = int(fli[3])

		if (flag & 1)== 1:
			is_R2 = (flag & 128) >0
		else:
			# MUST be paired-end mapping
			assert(0)

		all_keys[(rname, is_R2)]=1

		if (flag & 4) == 0:
			genes=fl.strip().split("\tXT:Z:")
			assert(len(genes)<=2)
			if len(genes) ==2: genes=genes[1].split(",")
			else: genes=None
			ret[( rname, is_R2 )] = (genes, chro, pos)

	return ret

tool1_res = loadsam(argv[1])
tool2_res = loadsam(argv[2])

MM_match=0
MM_unmatch=0
MU=0
UM=0
UU=0

for k in all_keys.keys():
	if k in tool1_res and k in tool2_res:
		g1,c1,p1 = tool1_res[k]
		g2,c2,p2 = tool2_res[k]
		if False:
			print(g1)
			print(g2)
			if g1 and g2 :print( [xx for xx in g1 if xx in g2] )
			print( g1 and g2 and len( [xx for xx in g1 if xx in g2])>0 )
			print()

		if g1 and g2 and len( [xx for xx in g1 if xx in g2])>0:
			MM_match+=1
		elif c1 == c2 and abs(p1-p2)<=100: MM_match+=1
		else: MM_unmatch+=1
	elif k in tool1_res: MU+=1
	elif k in tool2_res: UM+=1
	else: UU+=1


print("%.4f"%( (UU+MM_match)*1./(MM_match+ MM_unmatch+ MU+ UM+ UU)))
