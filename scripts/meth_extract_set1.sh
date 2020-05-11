for sample_set in AA_F00 AA_F25 AH_F25

do

    echo starting $sample_set

    for lane_input in lane1 lane2

    do

        echo starting ${sample_set} ${lane_input}

        ~/bin/Bismark-0.22.3/bismark_methylation_extractor --scaffolds --gzip --buffer_size 30% \
            --comprehensive \
            --no_header \
            --parallel 1 \
            --ignore 2 --ignore_r2 2 \
            --output ~/tonsa_epigenetics/analysis/methylation_extract/${lane_input} \
            --genome_folder /data/copepods/tonsa_genome/ \
            ~/tonsa_epigenetics/analysis/aligned/${lane_input}/${sample_set}*.bam

            #~/tonsa_epigenetics/analysis/aligned/${lane_input}/${sample_set}*pe.bam

        echo Finished ${sample_set} ${lane_input}

    done

        echo Finished ${sample_set} ALL SAMPLES DONE

done
