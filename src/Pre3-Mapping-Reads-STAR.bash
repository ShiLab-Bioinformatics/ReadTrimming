st=STAR
index=STAR273a-hg38-hasAnno
purp=Normal

for test in SEQC-A SEQC-B 
do
	trimmer=RAW 
	sc=`echo $test |cut -f2 -d'-' `
	r1=SEQC_ILM_AGR_${sc}_1-15M_R1.fastq.gz
	r2=SEQC_ILM_AGR_${sc}_1-15M_R2.fastq.gz

	# Normal STAR
	$st --runThreadN 8 --outFilterMultimapNmax 1 --genomeDir $index --readFilesIn \
	$r1 $r2 --readFilesCommand zcat --outFileNamePrefix \
	del4-STAR-$test-$trimmer.bam --outSAMtype BAM Unsorted;\
	$st --runThreadN 8 --outFilterMultimapNmax 1 --genomeDir $index --readFilesIn \
	$r1 $r2 --readFilesCommand zcat --outFileNamePrefix \
	STAR-$test-$trimmer.bam  --outSAMtype BAM Unsorted --sjdbFileChrStartEnd  \
	del4-STAR-$test-$trimmer.bamSJ.out.tab

	# "EndToEnd"
	$st --runThreadN 8 --outFilterMultimapNmax 1 --genomeDir $index --readFilesIn \
	$r1 $r2 --readFilesCommand zcat --outFileNamePrefix \
	del4-STAR-$test-$trimmer.bam --outSAMtype BAM Unsorted;\
	$st --runThreadN 8 --outFilterMultimapNmax 1 --genomeDir $index --readFilesIn \
	$r1 $r2 --readFilesCommand zcat --outFileNamePrefix \
	STARe2e-$test-$trimmer.bam  --outSAMtype BAM Unsorted --sjdbFileChrStartEnd  \
	del4-STAR-$test-$trimmer.bamSJ.out.tab --alignEndsType EndToEnd 
done
