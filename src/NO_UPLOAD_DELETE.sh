for f in `ls |grep -Ev "(sh|py)$"|grep -v TruSeq[0-9]|grep -v ^hg38|grep -v ^SEQC`
do
	rm -f $f
done
