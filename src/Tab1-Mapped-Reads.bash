echo
echo "This script gives fractions of mapped reads and mapped read-bases."
echo "The results are in Tab 1 (Percentages of mapped read bases with or without read trimming)."
echo "Column order: SEQC-A SEQC-B Simu:0.1% Simu:0.5% Simu:1%"
echo
for fii in READ-BASES
do
		echo
		echo "====== MAPPED $fii ======"
		purp=Normal
		for trimmer in RAW maticWindow maticInfo galore 
		do
			for aligner in Subread STAR
			do
				printf $trimmer
				printf "\t"
				printf $aligner
				printf "\t"
				for dset in SEQC-A SEQC-B Simu0010 Simu0050 Simu0100 
				do
					if [[ $aligner == Subread ]]
					then
						fn=Mapped-$dset-$purp-$trimmer.bam 
					elif [[ $aligner == STAR ]]
					then
						fn=STAR-$dset-$trimmer.bamAligned.out.bam
					fi

					if test -f $fn 
					then

						if [[ $fii == READS ]]
						then
							dset_size=` samtools view  Mapped-$dset-$purp-RAW.bam | wc -l`
						else
							dset_size=` samtools view  Mapped-$dset-$purp-RAW.bam | awk '{a+=length($10)}END{print a}' `
						fi

						if [[ $fii == READS ]]
						then
							map_read_pct=$( samtools view $fn | awk -v ds_size=$dset_size ' and($2,4)==0{mm++} END{ printf("%.3f", mm*100./ds_size) } ')
						else
							map_read_pct=$( samtools view $fn | awk -v ds_size=$dset_size ' {inpreads ++; if(and($2,4)==0) mapped_reads++; inpbases=length($10); rawcigar=$6; gsub( /[*]|[0-9]+[IMDHN]/,"",$6 ) ; split( $6, aa,"S" ); clipbases=0; for(si in aa){clipbases+=aa[si]};    INPB += inpbases; if(and($2,4)==0 ){MAPB += inpbases - clipbases ; CLIPB += clipbases  } ; if(0)if(clipbases >0 || rawcigar~/S/){ print rawcigar, clipbases } } END{print MAPB * 100.0/ds_size}' )
						fi
					else
						map_read_pct="NA"
					fi
					printf "%s" $map_read_pct
					printf "\t"
				done
				echo
			done
		done
done
