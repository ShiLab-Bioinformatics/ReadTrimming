galore=trim_galore
matic="java -jar trimmomatic-0.39.jar "
GALORE_OPT=" --illumina  -j 8 --paired "
GALORE_ADAPTOR="   --illumina  -j 8 --paired -q 0 "
MARIC_OPT=" PE -threads 8 -phred33"


for test in  SEQC-A SEQC-B
do
	
	R1=SEQC_ILM_AGR_A_1-15M_R1.fastq.gz
	if [[ $test == SEQC-B ]]
	then
		R1=SEQC_ILM_AGR_B_1-15M_R1.fastq.gz
	fi

	zcat $R1 | awk ' NR%4==1 || NR%4==3 {print;next} {print substr($0,51)}  ' |gzip -c > Chopped-$test.fastq.gz
done


for test in SEQC-A SEQC-B Simu0010 Simu0050 Simu0100 
do
	ADA_MA=10
	ADA_FA=TruSeq3-PE.fa
	if [[ $test =~ SEQC50 ]]
	then
		ADA_FA=TruSeq2-PE.fa
	fi

	for usage in Normal Adapter
	do
		R1=SEQC_ILM_AGR_A_1-15M_R1.fastq.gz
		R2=SEQC_ILM_AGR_A_1-15M_R2.fastq.gz
		if [[ $test == SEQC-B ]]
		then
			R1=SEQC_ILM_AGR_B_1-15M_R1.fastq.gz
			R2=SEQC_ILM_AGR_B_1-15M_R2.fastq.gz
		elif [[ $test =~ "Simu0" ]]
		then
			bv=` echo $test | sed "s/Simu//g"`
			R1=hg38-Adaptor-${bv}_R1.fastq.gz
			R2=hg38-Adaptor-${bv}_R2.fastq.gz
		fi

		if [[ $usage == "Normal" ]] 
		then
			MATIC_TRIMMERS_WINDOW=" ILLUMINACLIP:$ADA_FA:2:30:$ADA_MA  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15  MINLEN:36 "
			MATIC_TRIMMERS_INFO=" ILLUMINACLIP:$ADA_FA:2:30:$ADA_MA  LEADING:3 TRAILING:3  MAXINFO:50:0.5  MINLEN:36 "
			MATIC_TRIMMERS_ADAPTOR=" ILLUMINACLIP:$ADA_FA:2:30:$ADA_MA:8:true "
		else
			MATIC_TRIMMERS_WINDOW=" ILLUMINACLIP:$ADA_FA:2:30:$ADA_MA:8:true  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15  MINLEN:36 "
			MATIC_TRIMMERS_INFO=" ILLUMINACLIP:$ADA_FA:2:30:$ADA_MA:8:true  LEADING:3 TRAILING:3  MAXINFO:50:0.5  MINLEN:36 "
			MATIC_TRIMMERS_ADAPTOR=" ILLUMINACLIP:$ADA_FA:2:30:$ADA_MA:8:true "
		fi

		$matic $MARIC_OPT $R1 $R2 -baseout Trimmed-MaticWindow-$usage-$test.fq.gz  $MATIC_TRIMMERS_WINDOW
		$matic $MARIC_OPT $R1 $R2 -baseout Trimmed-MaticInfo-$usage-$test.fq.gz  $MATIC_TRIMMERS_INFO
	done

	$galore $GALORE_OPT $R1 $R2 --basename Trimmed-Galore-$test
	$galore $GALORE_ADAPTOR $R1 $R2 --basename Trimmed-GaloreAda-$test &> Trimmed-GaloreAda-$test.log
	$matic $MARIC_OPT $R1 $R2 -baseout Trimmed-MaticAda-$test.fq.gz  $MATIC_TRIMMERS_ADAPTOR &> Trimmed-MaticAda.fq.gz.log
done


galore=trim_galore
matic="java -jar trimmomatic-0.39.jar "
ADAFA=TruSeq3-PE-A2R.fa
ADAMA=10
MATIC_TRIMMERS_WINDOW=" ILLUMINACLIP:$ADAFA:2:30:$ADAMA  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15  MINLEN:25 "
MATIC_TRIMMERS_INFO=" ILLUMINACLIP:$ADAFA:2:30:$ADAMA  LEADING:3 TRAILING:3  MAXINFO:50:0.5  MINLEN:25 "
GALORE_OPT=" --illumina  -j 8   "
MARIC_OPT=" SE -threads 8 -phred33"
GALORE_ADAPTOR="   --illumina  -j 8 -q 0 "
purp=Normal

for test in SEQC-A SEQC-B
do
	R=Chopped-$test.fastq.gz
	$galore $GALORE_OPT $R --basename Trimmed-Galore-Chopped-$test-$purp
	$matic $MARIC_OPT $R Trimmed-MaticWindow-Chopped-$test-$purp.fq.gz  $MATIC_TRIMMERS_WINDOW
	$matic $MARIC_OPT $R Trimmed-MaticInfo-Chopped-$test-$purp.fq.gz  $MATIC_TRIMMERS_INFO
done
