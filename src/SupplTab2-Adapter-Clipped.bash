echo
echo This script calculates the fractions of clipped bases and adapter bases
echo "Every two columns are for a data set: SEQC-A SEQC-B Simu:0.1% Simu:0.5% Simu:1%"
echo "Within the two columns for each dataset, the values are Adapter-bases% amd Clipped-bases-in-Adapter%"
echo

for trimmer in maticWindow maticInfo galore 
do
  printf $trimmer
  printf "\t"
  for t in  SEQC-A SEQC-B  Simu0010 Simu0050 Simu0100 
  do
	ada_trimmer=galore
    if [[ $trimmer =~ matic ]]
    then
       ada_trimmer=matic
    fi

	rawbam=Mapped-$t-Normal-RAW.bam
    trimmed_bam=Mapped-$t-Adapter-$trimmer.bam
    if ! test -f trimmed_bam
    then
		trimmed_bam=Mapped-$t-Normal-$trimmer.bam
    fi
	trr=$( python SupplTab2-adapter-clip.py <( samtools view $rawbam ) <( samtools view  $trimmed_bam ) SupplTab1-$t-$ada_trimmer-adapter.txt )
	printf "%s" $trr
  done
  echo
done
