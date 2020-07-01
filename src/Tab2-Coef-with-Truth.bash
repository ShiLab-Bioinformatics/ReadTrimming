fc=featureCounts
anno=hg38_RefSeq_exon.txt
purp=Normal

echo 
echo  ==== Subread ALIGNER TESTS ====
echo "The columns are SEQC-A SEQC-B Chopped-SEQC-A (50bp) Chopped-SEQC-B (50bp) Simu:0.1% Simu:0.5% Simu:1%"
echo 

for trimmer in RAW maticWindow maticInfo galore
do
  printf $trimmer
  printf "\t"
  for t in  SEQC-A SEQC-B Chopped-SEQC-A Chopped-SEQC-B Simu0010 Simu0050 Simu0100
  do
    if [[ $t =~ Chop ]]
	then
	  fn=Mapped-$t-$trimmer.bam
	else
	  fn=Mapped-$t-$purp-$trimmer.bam
	fi
    #$fc -a $anno -o del4-genes.txt0 -T6 -p -M -F SAF $fn 
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
echo "The columns are SEQC-A and SEQC-B "
echo 

for mode in  RAW-STAR RAW-STAR-E2E
do
  printf $mode
  printf "\t"
  for t in  SEQC-A SEQC-B 
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
