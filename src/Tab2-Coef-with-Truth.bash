fc="Rscript Tab2-featureCounts.R"
anno=hg38_RefSeq_exon.txt
purp=Normal

echo 
echo "This script calculates the correlations between read-mapping and counting results with known truths."
echo "The results are presented in Tab 2 (real data, Subread) and Suppl Tab S4 (simulation) and S6 (reaal data, STAR)."
echo 
echo  ==== Subread TESTS ====
echo "The columns are SEQC-UHRR, SEQC-HBRR, Chopped-SEQC-A (50bp), Chopped-SEQC-B (50bp), Simu:0.1%, Simu:0.5% and Simu:1%"
echo 

for trimmer in RAW maticWindow maticInfo galore
do
  printf $trimmer
  printf "\t"
  for t in  SEQC-A SEQC-B Chopped-SEQC-A Chopped-SEQC-B Simu0010 Simu0050 Simu0100
  do
	fn=Mapped-$t-$purp-$trimmer.bam
	rm -f del4-genes.txt
    nohup $fc del4-genes.txt $fn &>/dev/null
#    $fc del4-genes.txt $fn 

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
echo  ==== STAR TESTS ====
echo "The columns are SEQC-UHRR and SEQC-HBRR "
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

	rm -f del4-genes.txt
    nohup $fc del4-genes.txt $fn &>/dev/null

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
echo "The results for Subread-UHRR and Subread-HBRR are in the Subread result table."
