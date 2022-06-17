cd /home/arcasHLA

./arcasHLA genotype /mnt/possorted_genome_bam.extracted.fq.gz \
	-g A,B,C,DPB1,DQB1,DQA1,DRB1 \
	--log /out/genotype_log.txt \
	--single \
	-o /out/ \
	-t 4
