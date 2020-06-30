fc=featureCounts
anno=hg38_RefSeq_exon.txt
purp=Normal

echo 
echo  ==== Subread ALIGNER TESTS ====
echo 

for trimmer in RAW maticWindow maticInfo galore
do
  printf $trimmer
  printf "\t"
  for t in  SEQC-A SEQC-B Simu0010 Simu0050 Simu0100
  do
	fn=Mapped-$t-$purp-$trimmer.bam
    nohup $fc -a $anno -o del4-genes.txt0 -T6 -p -M -F SAF $fn &>/dev/null
	cat del4-genes.txt0|grep -v ^# |grep -v Leng |cut -f 1,7>del4-genes.txt

    rrv=$(if [[ $t =~ Simu ]]
	then
		bv=`echo $t|sed 's/Simu//g'`
    	python eval*py genes-Simulation del4-genes.txt hg38-Adaptor-$bv.trueCountsGene |cut -f12 -d' '
	elif [[ $t =~ "-A" ]]
	then
    	python eval*py genes-UHRR del4-genes.txt TaqMan-RTPCR-data.txt  |cut -f12 -d' '
	else
    	python eval*py genes-HBRR del4-genes.txt TaqMan-RTPCR-data.txt  |cut -f12 -d' '
	fi)
    printf "%.5f\t" $rrv
  done
  echo
done

echo 
echo  ==== OTHER ALIGNER TESTS ====
echo 

for mode in  RAW-STAR RAW-STAR-E2E
do
  printf $mode
  printf "\t"
  for t in  SEQC-A SEQC-B Simu0010 Simu0050 Simu0100
  do
	if [[ $mode =~ "STAR-E2E" ]]
	then
		fn=STARe2e-$t-RAW.bamAligned.out.bam
	elif [[ $mode =~ "STAR" ]]
	then
		fn=STAR-$t-RAW.bamAligned.out.bam
	fi

    nohup $fc -a $anno -o del4-genes.txt0 -T6 -p -M -F SAF $fn &>/dev/null
	cat del4-genes.txt0|grep -v ^# |grep -v Leng |cut -f 1,7>del4-genes.txt


    rrv=$(if [[ $t =~ Simu ]]
	then
		bv=`echo $t|sed 's/Simu//g'`
    	python eval*py genes-Simulation del4-genes.txt hg38-Adaptor-$bv.trueCountsGene |cut -f12 -d' '
	elif [[ $t =~ "-A" ]]
	then
    	python eval*py genes-UHRR del4-genes.txt TaqMan-RTPCR-data.txt  |cut -f12 -d' '
	else
    	python eval*py genes-HBRR del4-genes.txt TaqMan-RTPCR-data.txt  |cut -f12 -d' '
	fi)

    printf "%.5f\t" $rrv

  done
  echo
done
