# Read in data
import pandas as pd
samples_table = pd.read_csv("samples.csv").set_index("sample",drop=False)

# Definitions for rule all
gex_table = samples_table[samples_table.sample_gex == True]
SAMPLES_GEX = gex_table.loc[:,"sample"].values.tolist()

citeseq_table = samples_table[samples_table.sample_citeseq == True]
SAMPLES_CITESEQ = citeseq_table.loc[:,"sample"].values.tolist()

vdj_table = samples_table[samples_table.sample_vdj == True]
vdj_library_table = vdj_table[~vdj_table['fastq_vdj'].astype(str).str.contains("False")]
vdj_no_library_table = vdj_table[vdj_table['fastq_vdj'].astype(str).str.contains("False")]
SAMPLES_VDJ_library = vdj_library_table.loc[:,"sample"].values.tolist()
SAMPLES_VDJ_library = [s + "_TCR" for s in SAMPLES_VDJ_library]
SAMPLES_VDJ_no_library = vdj_no_library_table.loc[:,"sample"].values.tolist()
# SAMPLES_VDJ = SAMPLES_VDJ_library + SAMPLES_VDJ_no_library

velocyto_table = samples_table[samples_table.sample_velocity == True]
SAMPLES_VELO = velocyto_table.loc[:,"sample"].values.tolist()

arcas_table = samples_table[samples_table.sample_arcas == True]
SAMPLES_ARCAS = arcas_table.loc[:,"sample"].values.tolist()

# Define local rule
localrules:
	all

# Rule to define output directories
rule all:
	input:
		expand("cellranger/{final_gex}",final_gex=SAMPLES_GEX),
		expand("citeseq/{final_citeseq}",final_citeseq=SAMPLES_CITESEQ),
		expand("vdj_cellranger/{final_vdj}",final_vdj=SAMPLES_VDJ_library),
		expand("vdj_trust4/{final_vdj}",final_vdj=SAMPLES_VDJ_no_library),
		expand("velocyto/{final_velocity}",final_velocity=SAMPLES_VELO),
		expand("arcas/{final_arcas}/genotype/",final_arcas=SAMPLES_ARCAS)

# Function to read in gex fastq paths from samples_table
def get_gex_input(wildcards):
	return gex_table.loc[wildcards.sample,"fastq_gex"]

# Cellranger count rule
rule count:
	input:
		get_gex_input
	output:
		directory("cellranger/{sample}/")
	threads:
		8
	resources:
		runtime="03-00:00:00",
		mem_mb=128000
	shell:
		"""
		module load cellranger/6.1.2
		"""
		"""
		mkdir -p cellranger
		cd cellranger
		cellranger count --id={wildcards.sample} \
			--transcriptome=/bgfs/genomics/refs/CellRanger/refdata-cellranger-GRCh38-3.0.0 \
			--fastqs={input} \
			--sample={wildcards.sample} \
			--localcores=8 \
			--localmem=128
		"""

# Function to read in citeseq fastq paths from samples_table
def get_citeseq_input(wildcards):
	return citeseq_table.loc[wildcards.sample,"fastq_citeseq"]

# CITEseq count rule
rule citeseq:
	input:
		fastq=get_citeseq_input,
		barcodes="cellranger/{sample}/"
	output:
		directory("citeseq/{sample}")
	threads:
		4
	resources:
		runtime="12:00:00",
		mem_mb=40000
	shell:
		"""
		module load gcc/8.2.0 r/3.6.0
		"""
		"""
		find {input.fastq} -name "*_R1_*" | xargs cat > {input.fastq}/r1_merged_fastq.gz
		find {input.fastq} -name "*_R2_*" | xargs cat > {input.fastq}/r2_merged_fastq.gz
		"""
		"""
		Rscript cell_barcode_id.R {input.barcodes}/outs/filtered_feature_bc_matrix {wildcards.sample}
		EXPECTED_CELLS=$(wc -l {wildcards.sample}_whitelist.csv | cut -d' ' -f 1)
		"""
		"""
		module purge
		module load cite-seq-count/1.4.3
		"""
		"""
		CITE-seq-Count -R1 {input.fastq}/r1_merged_fastq.gz \
			-R2 {input.fastq}/r2_merged_fastq.gz \
			-t citeseq_totalseqC_10tags_human_murine.csv \
  			-cbf 1 -cbl 16 -umif 17 -umil 28 \
			-trim 10 \
			-o {output} \
			-cells $EXPECTED_CELLS \
			-wl {wildcards.sample}_whitelist.csv \
			-T 4
		"""

# Function to define vdj samples with libraries
def get_vdj_library_input(wildcards):
	return vdj_library_table.loc[wildcards.sample,"fastq_vdj"]

# VDJ with cellranger with VDJ library
rule vdj_lib:
	input:
		get_vdj_library_input
	output:
		directory("vdj_cellranger/{sample}_TCR")
	threads:
		4
	resources:
		runtime="01-00:00:00",
		mem_mb=60000
	shell:
		"""
		module load cellranger/6.1.2
		"""
		"""
		mkdir -p vdj_cellranger
		cd vdj_cellranger
		cellranger vdj --id={wildcards.sample}_TCR \
			--fastqs={input} \
			--reference=/bgfs/genomics/refs/CellRanger/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0 \
			--sample={wildcards.sample}_TCR \
			--localcores=4 \
			--localmem=59
		"""

# VDJ with TRUST without VDJ library
rule vdj_no_lib:
	input:
		"cellranger/{sample}"
	output:
		directory("vdj_trust4/{sample}")
	threads:
		1
	resources:
		runtime="12:00:00",
		mem_mb=16000
	shell:
		"""
		/zfs1/tbruno/arc85/TRUST4/run-trust4 \
			-f /zfs1/tbruno/arc85/TRUST4/hg38_bcrtcr.fa \
			--ref /zfs1/tbruno/arc85/TRUST4/human_IMGT+C.fa \
			-b {input}/outs/possorted_genome_bam.bam \
			--barcode CB \
			--od {output} \
			-t 1
		"""

# Velocyto rule
rule velocity:
	input:
		"cellranger/{sample}"
	output:
		directory("velocyto/{sample}")
	threads:
		4
	resources:
		runtime="12:00:00",
		mem_mb=64000
	shell:
		"""
		module load gcc/8.2.0
		module load samtools/1.9
		module load velocyto/0.17
		"""
		"""
		velocyto run10x -m /zfs1/tbruno/arc85/00_INBOX/18071_for_velocity/mat/grch38_repeat_mask.gtf \
			{input} \
			/bgfs/genomics/refs/CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
		mkdir -p velocyto/{wildcards.sample}
		mv cellranger/{wildcards.sample}/velocyto velocyto/{wildcards.sample}
		"""

# ARCAS extract rule
rule arcas_extract:
	input:
		"cellranger/{sample}"
	output:
		directory("arcas/{sample}/fastq/")
	threads:
		4
	resources:
		runtime="02:00:00",
		mem_mb=64000
	shell:
		"""
		mkdir -p {output}
		module load singularity/3.9.6
		singularity exec --bind {input}:/mnt \
			--bind /ihome/tbruno/arc85/arcasHLA:/home/arcasHLA \
			--bind {output}:/out \
			--bind .:/script \
			/ihome/tbruno/arc85/arcashla.sif /bin/bash /script/arcas_extract.sh
		"""

# ARCAS align rule
rule arcas_align:
	input:
		"arcas/{sample}/fastq"
	output:
		directory("arcas/{sample}/genotype/")
	threads:
		4
	resources:
		runtime="02:00:00",
		mem_mb=64000
	shell:
		"""
		mkdir -p {output}
		module load singularity/3.9.6
		singularity exec --bind {input}:/mnt \
			--bind /ihome/tbruno/arc85/arcasHLA:/home/arcasHLA \
			--bind {output}:/out \
			--bind .:/script \
			/ihome/tbruno/arc85/arcashla.sif /bin/bash /script/arcas_align.sh
		"""
