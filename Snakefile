#---------------------------
# @author: Kedi Cao
# @email: kedi.cao@ukr.de
# @date: 28.08.2024
#---------  CONFIG set up  --------------- 
# This allows snakemake load yaml for configuration
import glob
import os
import re

configfile: os.path.join(workflow.basedir,"config","snake_config.yaml")
configfile: os.path.join(workflow.basedir,"config","cluster_config.yaml")
configfile: os.path.join(workflow.basedir, "config", "container_config.yaml")

# Input files and directories  
INPUT_DIR = config["INPUT_DIR"]
REFERENCE = config["REFERENCE"]
ALTREF = config["ALTREF"]
kraken_db_path = config["kraken_db_path"]
OUTPUT_PATH = config["OUTPUT_PATH"]
WORKFLOW_PATH = config["WORKFLOW_PATH"]
ASSEMBLER = config["ASSEMBLER"]
resources = config["resources"]
threads = config["threads"]
containers = config["containers"]
POPINS4SNAKE = config["POPINS4SNAKE"]
SICKLE = config["SICKLE_BIN"]
GATB = config["MINIA_BIN"]
VELVET = config["VELVET_BIN"]


relative_paths = {
    WORK_DIR := os.path.join(OUTPUT_PATH, config["WORK_DIR"]),
    RESULTS_DIR := os.path.join(OUTPUT_PATH, config["RESULTS_DIR"])    
    # POPINS4SNAKE := os.path.join(WORKFLOW_PATH,config["POPINS4SNAKE"]),
    # SICKLE := os.path.join(WORKFLOW_PATH,config["SICKLE_BIN"]),
    # GATB := os.path.join(WORKFLOW_PATH,config["MINIA_BIN"]),
    # VELVET := os.path.join(WORKFLOW_PATH, config["VELVET_BIN"])
    }

# Natural sort helper
def natural_key(s):
    # Splits on any digits (\d+) and tries to convert them to integers,
    # ensuring numeric chunks are sorted as numbers, not strings.
    return [int(chunk) if chunk.isdigit() else chunk for chunk in re.split(r'(\d+)', s)]

################################################################################
# Automatically grabbing sample names from INPUT_DIR
selected_samples = config["SAMPLES"]
# Recursively find all BAM files in INPUT_DIR (including subdirectories)
bam_files = glob.glob(os.path.join(INPUT_DIR, "**", "*.bam"), recursive=True)
# Extract sample names (assuming the sample name is the filename without the .bam extension)
auto_samples = [os.path.splitext(os.path.basename(bam))[0] for bam in bam_files]


# If the user specified a list of samples, filter to those; otherwise, take all
if selected_samples:
    auto_samples = set(auto_samples).intersection(selected_samples)

# Assuming the sample name is the filename without the .bam extension.
SAMPLES = sorted(auto_samples, key=natural_key)

sample_num = len(SAMPLES)
################################################################################

input_files = [
    RESULTS_DIR + "/supercontigs.fa",
    RESULTS_DIR + "/locations.txt",
    RESULTS_DIR + "/insertions_genotypes.vcf.gz"
    ]


if config["remove_contamination"]=="no":
    input_files.append(expand("{p}/{s}/{f}", p=WORK_DIR, s=SAMPLES, f=["POPINS_SAMPLE_INFO"]))
elif config["remove_contamination"]=="yes":
    input_files.append(expand("{p}/{s}/contaminate_info/paired_report.txt", p=WORK_DIR, s=SAMPLES))
    input_files.append(expand("{p}/{s}/contaminate_info/single_report.txt", p=WORK_DIR, s=SAMPLES))

if config["ANALYSIS"]=="yes":
    input_files.append(expand("{p}/{f}", p=RESULTS_DIR, f=["analysis_table.pdf","gc_content_dist.pdf","contig_length_dist.pdf","heatmap_coverage.pdf"]))
    if config["SICKLE"]=="no":
        input_files.append(expand("{p}/{f}", p=RESULTS_DIR, f=["read_numbers.pdf"]))
        if config["remove_contamination"]=="yes":
            input_files.append(expand("{p}/{f}", p=RESULTS_DIR, f=["read_numbers_no_sickle_after_clean.pdf"]))
    elif config["SICKLE"]=="yes":
        if config["remove_contamination"]=="no":
            input_files.append(expand("{p}/{f}", p=RESULTS_DIR, f=["read_numbers_before_filter.pdf","read_numbers_after_filter.pdf"]))
        elif config["remove_contamination"]=="yes":
            input_files.append(expand("{p}/{f}", p=RESULTS_DIR, f=["read_numbers_before_clean_filter.pdf","read_numbers_after_clean.pdf","read_numbers_after_clean_filter.pdf"]))



# these files are all temporary files now and will be removed after workflow is finished
# expand("{p}/{s}/{f}", p=WORK_DIR, s=SAMPLES, f=["non_ref_new.bam", "non_ref_new.bam.bai"]),
# expand("{p}/{s}/{assemblr}.contigs.fa", p=WORK_DIR, s=SAMPLES, assemblr=ASSEMBLER),
# expand("{p}/{s}/locations.txt", p=WORK_DIR, s=SAMPLES),
# "single.fastq","paired.1.fastq","paired.2.fastq","non_ref.bam",
# elif config["remove_contamination"]=="yes":
#   input_files.append(expand("{p}/{s}/{f}",p=WORK_DIR, s=SAMPLES, f=["single_clean.fastq","paired_1_clean.fastq","paired_2_clean.fastq"]))
#   input_files.append(expand("{p}/{s}/remapped_non_ref.bam", p=WORK_DIR, s=SAMPLES))

################################################################################
if config["ANALYSIS"]=="yes":
    include: os.path.join(WORKFLOW_PATH, "snakemodules","analysis.smk")

if config["remove_contamination"]=="no":
    include: os.path.join(WORKFLOW_PATH, "snakemodules", "crop_unmapped.smk")
elif config["remove_contamination"]=="yes":
    include: os.path.join(WORKFLOW_PATH, "snakemodules", "crop_remapped.smk")
    include: os.path.join(WORKFLOW_PATH, "snakemodules", "kraken.smk")

include: os.path.join(WORKFLOW_PATH, "snakemodules", "contigmap.smk")
include: os.path.join(WORKFLOW_PATH, "snakemodules", "merge.smk")
include: os.path.join(WORKFLOW_PATH, "snakemodules", "position.smk")
include: os.path.join(WORKFLOW_PATH, "snakemodules", "genotype.smk")
 
if config['ASSEMBLER'] == 'minia':
    include: os.path.join(WORKFLOW_PATH, "snakemodules", "minia.smk")
elif config["ASSEMBLER"] == 'velvet':
    include: os.path.join(WORKFLOW_PATH, "snakemodules", "kraken.smk")



rule all:
    input:
        input_files
