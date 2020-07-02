echo 
echo "This script generates results for Suppl Tab S3 (concordance of reads)"
echo "columns: SEQC-A and SEQC-B"
echo

fc=featureCounts
annot="hg38_RefSeq_exon.txt -F SAF"
for trimmer1 in RAW galore maticWindow maticInfo
do
	for test in SEQC-A SEQC-B
	do
		f1=Mapped-$test-Normal-$trimmer1.bam
		nohup $fc -M -O -a $annot -T6 -o del4-nomatter.FC -R SAM $f1 &>/dev/null
	done
done

for trimmer1 in RAW galore maticWindow maticInfo
do
	for trimmer2 in RAW galore maticWindow maticInfo
	do
		if [[ $trimmer1 > $trimmer2 ]]
		then 
			printf "$trimmer1\t$trimmer2\t"
			for test in SEQC-A SEQC-B
			do
				f1=Mapped-$test-Normal-$trimmer1.bam
				f2=Mapped-$test-Normal-$trimmer2.bam
				xx=$( python SupplTab3-Match-Read-Pos.py $f1.featureCounts.sam $f2.featureCounts.sam  )
				printf "%s\t" $xx
			done 
			echo
		fi
	done
done

