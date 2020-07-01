
echo
echo This script calculates the fractions of trimmed bases and adapter bases
echo "Every three columns are for a data set: SEQC-A SEQC-B Simu:0.1% Simu:0.5% Simu:1%"
echo "Within the three columns for each dataset, the values are Trimmed-bases%, Clipped-bases-in-Trimed%, Adapter-bases-in-Trim-and-Clip%"
echo

for test in SEQC-A SEQC-B Simu0010 Simu0050 Simu0100 
do
	cat <( samtools view -h Mapped-$test-Normal-galoreAda.bam ) <(samtools view -h Mapped-$test-Normal-RAW.bam ) | python SupplTab1-pure-adapter.py > SupplTab1-$test-galore-adapter.txt
	cat <( samtools view -h Mapped-$test-Normal-maticAda.bam ) <(samtools view -h Mapped-$test-Normal-RAW.bam ) | python SupplTab1-pure-adapter.py > SupplTab1-$test-matic-adapter.txt
done

for trimmer in maticWindow maticInfo galore
do
  printf $trimmer
  printf "\t"
  for t in SEQC-A SEQC-B Simu0010 Simu0050 Simu0100 
  do
	ada_trimmer=galore
    if [[ $trimmer =~ matic ]]
    then
       ada_trimmer=matic
    fi

	rawbam=Mapped-$t-Normal-RAW.bam
    trimmed_bam=Mapped-$t-Adapter-$trimmer.bam
	if ! test -f $trimmed_bam
	then
    	trimmed_bam=Mapped-$t-Normal-$trimmer.bam
    fi
	trr=$( python SupplTab1-trim-clip.py <( samtools view $rawbam ) <( samtools view  $trimmed_bam ) SupplTab1-$t-$ada_trimmer-adapter.txt )
	printf "%s" $trr
  done
  echo
done
