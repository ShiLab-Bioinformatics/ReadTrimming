if ! test -f Idx-hg38-full.00.b.tab
then
	nohup Rscript Pre2-Subread-Index.R &> /dev/null
fi

for test in Simu0010 Simu0050 Simu0100 SEQC-A SEQC-B
do
  for trimmer in maticWindow maticInfo galore RAW galoreAda maticAda
  do
	 for purp in Normal Adapter
	 do
		r1=Trimmed-MaticWindow-$purp-${test}_1P.fq.gz
		r2=Trimmed-MaticWindow-$purp-${test}_2P.fq.gz
		if [[ maticAda ==  $trimmer ]]
		then
		  r1=Trimmed-MaticAda-${test}_1P.fq.gz
		  r2=Trimmed-MaticAda-${test}_2P.fq.gz
		elif [[ maticInfo ==  $trimmer ]]
		then
		  r1=Trimmed-MaticInfo-$purp-${test}_1P.fq.gz
		  r2=Trimmed-MaticInfo-$purp-${test}_2P.fq.gz
		elif [[ galoreAda == $trimmer ]]
		then
		  r1=Trimmed-GaloreAda-${test}_R1_val_1.fq.gz
		  r2=Trimmed-GaloreAda-${test}_R2_val_2.fq.gz
		elif [[ galore == $trimmer ]]
		then
		  r1=Trimmed-Galore-${test}_R1_val_1.fq.gz
		  r2=Trimmed-Galore-${test}_R2_val_2.fq.gz
		elif [[  RAW == $trimmer ]]
		then
		  r1=Simulation_R1.fastq.gz
		  r2=Simulation_R2.fastq.gz
		  if [[ $test =~ SEQC- ]]
		  then
		    sc=`echo $test |cut -f2 -d'-' `
		    r1=SEQC_ILM_AGR_${sc}_1-15M_R1.fastq.gz
		    r2=SEQC_ILM_AGR_${sc}_1-15M_R2.fastq.gz
		  elif [[ $test =~ SEQC50- ]]
		  then
		    sc=`echo $test |cut -f2 -d'-' `
		    r1=SEQC50-${sc}_R1.fastq.gz
		    r2=SEQC50-${sc}_R2.fastq.gz
		  elif [[ $test =~ Simu0 ]]
		  then
				bv=$( echo $test |sed 's/Simu//g' )
				r1=hg38-Adaptor-${bv}_R1.fastq.gz
				r2=hg38-Adaptor-${bv}_R2.fastq.gz
		  fi
		fi
		if test -f Mapped-$test-$purp-$trimmer.bam.indel.vcf
		then
			continue
		fi
		nohup Rscript Pre2-Subread-Run.R $r1 $r2 Mapped-$test-$purp-$trimmer.bam  &>/dev/null
	done
  done
done

purp=Normal
for test in Chopped-SEQC-A Chopped-SEQC-B
do
  for trimmer in maticWindow maticInfo galore RAW
  do
	R=Trimmed-Galore-$test-Normal_trimmed.fq.gz
	if [[ $trimmer == maticWindow ]]
	then
		R=Trimmed-MaticWindow-$test-Normal.fq.gz
	elif [[ $trimmer == maticInfo ]]
	then
		R=Trimmed-MaticInfo-$test-Normal.fq.gz
	elif [[ $trimmer == RAW ]]
	then
		R=$test.fastq.gz
	fi

    if test -f Mapped-$test-$purp-$trimmer.bam.indel.vcf
    then
  	  continue
    fi
    nohup Rscript Pre2-Subread-Run.R $R NA Mapped-$test-$purp-$trimmer.bam  &>/dev/null
#    Rscript Pre2-Subread-Run.R $R NA Mapped-$test-$purp-$trimmer.bam  
  done
done
