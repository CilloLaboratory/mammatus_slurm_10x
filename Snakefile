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
SAMPLES_VDJ_no_library = vdj_no_library_table.loc[:,"sample"].values.tolist()
SAMPLES_VDJ_no_library = [s + "_TCR" for s in SAMPLES_VDJ_no_library]

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
		expand("citeseq_cellranger/{final_citeseq}",final_citeseq=SAMPLES_CITESEQ),
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
		runtime="3d",
		mem_mb=128000
	shell:
		"""
		module load cellranger/7.0.1
		"""
		"""
		mkdir -p cellranger
		cd cellranger
		cellranger count --id={wildcards.sample} \
			--transcriptome=/ix1/acillo/arc85/references/cellranger_ref_230418/GRCh38 \
			--fastqs={input} \
			--sample={wildcards.sample} \
			--localcores=8 \
			--localmem=128 \
			--include-introns=false
		"""

# Function to read in citeseq library CSV path from samples_table
def get_citeseq_library(wildcards):
	return citeseq_table.loc[wildcards.sample,"citeseq_library"]

# Cellranger FB rule
rule citeseq_cellranger:
	input:
		get_citeseq_library
	output:
		directory("citeseq_cellranger/{sample}")
	threads:
		8
	resources:
		runtime="12h",
		mem_mb=62000
	shell:
		"""
		module purge
		module load cellranger/7.0.1
		"""
		"""
		mkdir -p citeseq_cellranger
		cd citeseq_cellranger
		cellranger count --id={wildcards.sample} \
   		--libraries={input} \
		--transcriptome=/ix1/acillo/arc85/references/cellranger_ref_230418/GRCh38 \
		--feature-ref=/ix1/acillo/arc85/references/citeseq_references/citeseq_reference_list_cellranger.csv \
		--localcores=8 \
		--localmem=62 
		"""

# Function to define vdj samples with libraries
def get_vdj_library_input(wildcards):
	return vdj_library_table.loc[wildcards.sample,"fastq_vdj"]

# VDJ with cellranger with VDJ library
rule vdj_lib:
	input:
		get_vdj_library_input
	output:
		directory("vdj_cellranger/{sample}")
	threads:
		4
	resources:
		runtime="1d",
		mem_mb=60000
	shell:
		"""
		module load cellranger/7.0.1
		"""
		"""
		mkdir -p vdj_cellranger
		cd vdj_cellranger
		cellranger vdj --id={wildcards.sample} \
			--fastqs={input} \
			--reference=/ix/tbruno/arc85/00_INBOX/cellranger_ref_230418/vdj_reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
			--sample={wildcards.sample} \
			--localcores=4 \
			--localmem=59
		"""

# VDJ with TRUST without VDJ library
rule vdj_no_lib:
	input:
		"cellranger/{sample}"
	output:
		directory("vdj_trust4/{sample}_TCR")
	threads:
		1
	resources:
		runtime="12h",
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
		runtime="12h",
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
		runtime="2h",
		mem_mb=64000
	shell:
		"""
		mkdir -p {output}
		module load singularity/3.9.6
		singularity exec --bind {input}:/mnt \
			--bind /ihome/acillo/arc85/arcasHLA:/home/arcasHLA \
			--bind {output}:/out \
			--bind .:/script \
			/ihome/acillo/arc85/arcashla.sif /bin/bash /script/arcas_extract.sh
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
		runtime="2h",
		mem_mb=64000
	shell:
		"""
		mkdir -p {output}
		module load singularity/3.9.6
		singularity exec --bind {input}:/mnt \
			--bind /ihome/acillo/arc85/arcasHLA:/home/arcasHLA \
			--bind {output}:/out \
			--bind .:/script \
			/ihome/acillo/arc85/arcashla.sif /bin/bash /script/arcas_align.sh
		"""
