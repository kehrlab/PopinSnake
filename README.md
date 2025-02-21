# PopinSnake workflow

This Git repository contains a generic [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for running the program [*PopIns4Snake*](https://gitlab.informatik.hu-berlin.de/fonda_a6/popins4snake.git).
It is supposed to facilitate your PopIns run and provides examples for a streamlined workflow. 

## Contents

1. [Environment Setup](#Environment-Setup)
1. [Program Installation](#program-installation)
1. [Workflow Configuration](#setting-up-your-workflow)
1. [Workflow execution](#workflow-execution)
1. [Example data](#example-data)
1. [References](#references)


## Environment Setup

As a minimum, an installation of the [Conda package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (Miniconda) is required for running the workflow.
In addition, we recommend to install [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) in the Conda base environment:

```
$ conda install mamba -n base -c conda-forge
```

If you do *not* have [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed already, create and activate a Conda environment for it:

```
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
$ conda activate snakemake
```

Most dependencies of the PopIns workflow will be automatically installed via Conda/Mamba when running the workflow for the first time.

Dependencies that are not available in Conda are included as additional programs in this repository and need to be installed before you can run the workflow.
These installations require a C++14 capable compiler, e.g. [GCC](https://gcc.gnu.org) version &ge;4.9.2, [CMake](https://cmake.org) version &ge;2.8.12 and [Zlib](http://www.zlib.net/) to be available on your system.
If *not* already available on your system, they can be installed via Conda/Mamba:

```
$ conda activate base
$ mamba create -c conda-forge -n gcc cxx-compiler zlib cmake
$ conda activate gcc
```

## Program Installation

Download the PopinSnake workflow from this repository via

```
$ git clone --recursive https://gitlab.informatik.hu-berlin.de/fonda_a6/popinSnake.git
$ cd popinSnake
```

The flag `--recursive` is required to download the required programs.

Before you can run the workflow, you need to install the programs `popins4snake` and `sickle`.
The `gatb-minia-pipeline` program requires no separate installation.

#### Installing the *PopIns4Snake* program

[PopIns4Snake](https://gitlab.informatik.hu-berlin.de/fonda_a6/popins4snake) is the main program executed by the workflow.
Before we can compile *PopIns4Snake*, we need to compile and install its dependency [Bifrost](https://github.com/pmelsted/bifrost).
For this, navigate to the `bifrost` folder and compile *Bifrost* using the flag `MAX_KMER_SIZE=64`:

```
$ cd submodules/popins4snake/external/bifrost
$ mkdir local
$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=../local -DMAX_KMER_SIZE=64
$ make
$ make install
```

This installs *Bifrost* in the folder `submodules/popins4snake/external/bifrost/local`, which is where the *PopIns4Snake* Makefile will look for it.

Now we can go back to the main `popins4snake` folder and compile *PopIns4Snake*:

```
$ cd ../../../
$ mkdir build
$ make
```

A binary `popins4snake` should now be available in `submodules/popins4snake/`.
The workflow is pre-configured to search for it in this folder.


#### Installing the *Sickle* program

[Sickle](https://github.com/najoshi/sickle) is a small program for filtering and trimming sequencing reads.
All you need to do to install *Sickle* is to navigate to its folder and run make:

```
$ cd submodules/sickle
$ make
```

This should create the binary `sickle` in `submodules/sickle/`.
It is not required to copy the binary to a location available in your path as described in the *Sickle* installation instructions.
The workflow is pre-configured to search for the `sickle` binary in the `submodules/sickle/` folder.

_If running into environment issue while compiling `sickle`, try_
```
$ make CFLAGS+="-I${CONDA_PREFIX}/include"
```

</br>

In case you created the GCC Conda environment, remember to now activate the Snakemake environment again:

```
$ conda activate snakemake
```

## Workflow Configuration

You can find all configurable parameters of the workflow in the file `snake_config.yaml`.
You will need to adjust them according to your input data, preferred output folders, etc. as described in the following.

The pre-configured values allow running the workflow on the [example data](#example-data) included in this repository.

#### Input File Paths

Input to the workflow is a set of BAM files along with BAM index (BAI) files, and the corresponding reference genome in FASTA format.

```
# Define input directories here
REFERENCE:
   ../path/to/hg38.fa
INPUT_DIR:
   ../path/to/my_bam_files/
```

The example contains separate subfolders for each sample (input BAM file and index file) in the `example_data`:

```
INPUT_DIR
├── genome1
│   ├── genome1.bam
│   └── genome1.bam.bai
├── genome2
│   ├── genome2.bam
│   └── genome2.bam.bai
└── genome3
    ├── genome3.bam
    └── genome3.bam.bai
```
However, the workflow is also able to recognize a list of bam files and index files directly under the `INPUT_DIR`, there's no need to follow the exact file structure as the example.
<!--
Further, the workflow requires you to specify all sample names:

```
# Provide sample names here
SAMPLES:
   - genome1
   - genome2
   - genome3
```
-->

#### Workflow and Output File Paths
```
WORKFLOW_PATH:
   ..prefix/of/workflow/path
OUTPUT_PATH:
   ..prefix/of/output/path
```
`OUTPUT_PATH` should indicate the full path of where the workflow is located, this path serves as a prefix of path for all the sub-directories within the workflow.
e.g.: `/home/user/popinSnake`\
The workflow automatically creates two output folders `WORK_DIR` and `RESULTS_DIR` within the workflow directory.\
`WORK_DIR` is intended for intermediate data per sample, i.e., the workflow will create a subfolder per sample in this directory.
Depending on the number of BAM files input to the workflow, the space requirements of this folder may grow quite large.
You may specify it in a temporary location or choose to delete it after successful completion of the workflow.\
`RESULTS_DIR` will hold the final VCF and FASTA output files describing the called insertions as well as the figures created during intermediate data analyses.
You probably want to keep most files written to this folder.


#### Other Workflow Behaviours
You can influence the behaviour of the workflow by changing the parameters in this section:

Use sickle filtering or the default popins4snake filtering:
```
SICKLE:
   "yes" or "no"
```
The `SICKLE` parameter choose which method to use for trimming the reads by defining the quality and length thresholds.\
By selecting `"yes"` it will use [`sickle`](https://github.com/najoshi/sickle) for quality trimming and filtering. The default quality is 20 and length is 60.\
By selecting `"no"` it will use the default `popins4snake` parameters for quality trimming and filtering. The default quality is 20 and length is 60. \
The trimming quality and length for either methods can be adjusted in the paramter section in the configuration file using `min-qual` and `min-read-len`.

Choose assembly tool:
```
ASSEMBLER:
   "minia" or "velvet"
```
The `ASSEMBLER` parameter specifies the program to use for assembling reads without alignment to the reference genome into contiguous sequences (contigs).
The default workflow installation supports two different assembly programs, [`"velvet"`](https://doi.org/10.1101/gr.074492.107) and [`"minia"`](https://doi.org/10.1186/1748-7188-8-22). 

Choose to view the data analysis with jupyter notebooks or not:
```
ANALYSIS:
   "yes" or "no"
```
The `ANALYSIS` parameter takes in either `"yes"` or `"no"` as inputs to activate a sub-module for data analysis with jupyter notebooks. \
The *`unmapped_analysis`* returns figures of unmapped reads before assembly.
The *`assembly_analysis`* returns figures showing the quality of the assembly.
The *`converage_analysis`* returns figures showing contig coverage before positioning and genotyping.
To use the jupyter notebook analysis, it requires an installation of snakemake version 5.32.0 and the snakemake flag `--edit-notebook`.
Please refer to the snakemake documentation on how to use [jupyter notebooks](https://snakemake.readthedocs.io/en/v6.12.0/snakefiles/rules.html?highlight=jupyter%20notebook#jupyter-notebook-integration) for further information.

#### Contamination Removal Module
The *popinSnake* workflow can identify and remove contaminations using [`kraken2`](https://github.com/DerrickWood/kraken2).\
Configure the contamination removal module here:
```
remove_contamination:
   "yes" or "no"
ALTREF:
   ../path/to/altref.fa
kraken_db_path:
   /example_data/kraken_db
```
- The `remove_contamination` lets user choose to identify and remove sample contamination with or without using `kraken2`. The parameter takes in either `"yes"` or `"no"` as inputs.
- An "alternative reference" in FASTA format for contamination removal can be specified.
- The `kraken_db_path` parameter is used to define a path of which the kraken2 database is stored. 
- In order to use this submodule, user must choose their preferred [kraken database](https://benlangmead.github.io/aws-indexes/k2), and follow the instruction below to download the database and store it under specified `kraken_db_path` e.g. /example_data/kraken_db:
```
cd /example_data/kraken_db/
mkdir -p {kraken_db_path}
wget {link_to_your_selected_kraken_lib} -O {kraken_db_path}/kraken_db.tar.gz
tar -xzvf {kraken_db_path}/kraken_db.tar.gz -C {kraken_db_path}/
rm {kraken_db_path}/kraken_db.tar.gz
```

<!--
- The `install_kraken_db` can automatically install the selected kraken database by link defined in the following parameter `kraken_lib` if user don't wish to download the database by themselves. The parameter takes in either `"yes"` or `"no"` as inputs.
- The `kraken_lib` takes a link of the chosen [kraken database](https://benlangmead.github.io/aws-indexes/k2) and combine with `install_kraken_db` set to `"yes"`, the workflow will automatically download the database from the link and store it under /example_data/kraken_db. The default database in the config file is Standard-8GB.
- The `python_script` contains the python code for parsing the contaminated genome information to create files for remapping.

install_kraken_db:
   "yes" or "no"
kraken_lib:
   "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz"
-->


<!-- 
The `REMAP` parameter can be used to activate a sub-workflow for remapping unaligned reads to the alternative reference (specified in the [input section](#configuring-locations-of-the-input-data) as `ALTREF`) by setting it to `REMAP="yes"`.
The purpose of the sub-workflow is to remove potential contamination in your read data.
More specifically, after cropping unaligned reads from the input BAM files, the remapping sub-workflow will attemt to align the unaligned reads to the alternative reference that consists of e.g. microbial genomes.
Any reads with an alignment to a sequence in the alternative reference will not be passed on to the following assembly steps. -->


#### Configuring workflow parameters

The main parameters of the programs called by the workflow can be configured in this section.
In most cases you will not need to change the default values.

Available parameters include:
<!-- - `memory`: The maximum amount of memory `samtools sort` is allowed to use per thread. Default: `2G`.
- `threads`: The number of threads allowed for parallelizing rules calling `samtools sort`, `bwa mem`, `minia`, and `popins4snake merge-contigs`. Default: `16`. -->
- `receive_email_notifications`: Choose between "yes" or "no" to receive email notifications before each jupyter notebook execution
- `email`: Specify the email address you wish to receive notifications from.
- `kmerlength`: The k-mer length for Velvet assembly. Default: `29`.
- `k_merge`: The k-mer length for merging contigs across samples using `popins4snake merge-contigs`. Default: `63`.
- `readlen`: The read length required by `popins4snake place-refalign` and `popins4snake place-splitalign`. Default: `150`.
- `min-qual`: For rules `popins2_crop_unmapped`, `popins2_crop_remapped` and `popins2_sickle`, minimum average quality value in windows for read quality trimming. Default: `30`.
- `min-read-len`: For rules `popins2_crop_unmapped`, `popins2_crop_remapped` and `popins2_sickle`, minimum read length to keep the read after quality trimming. Default: `60`.


#### [Optional] Configuring alternative installation directories

Define paths to program binaries. If you strictly followed the installation instructions above, you **will not need** to change the preconfigured values in this section. The default path to `submodules` is under `WORKFLOW_PATH`, which was appended to the binary paths in the Snakefile.
If you install any of the submodules in another location, it is possible to change the default paths to submodules binaries, however it requires you to remove the appended paths from the main Snakefile:
```
WORK_DIR:
   workdir
RESULTS_DIR:
   results
POPINS2_BIN:
   path/to/submodules/popins4snake/popins4snake
MINIA_BIN:
   path/to/submodules/gatb-minia-pipeline/gatb
SICKLE_BIN:
   path/to/submodules/sickle/sickle
VELVET_BIN:
   path/to/submodules/velvet/velvet
python_script:
   snakemodules/scripts/remap_classified_human.py
```


## Workflow Execution

After the initial setup, you can execute the workflow from the main directory of `popinSnake/` which contains the main Snakefile:

```
$ snakemake --use-conda --cores 1
```

The tag `--use-conda` is required for successful workflow execution unless you have all workflow dependencies (and their correct versions) specified in Conda environments (see `workflow/envs/` folder) installed in a location available on your path.

__Note: When you run the workflow for the first time using `--use-conda`, it will create the conda environments.
This process may take some time. The next time you run the workflow, this process is not necessary again.__

If you chose not to install `mamba` (see [Prerequisites](#prerequisites)), then you now have to add `--conda-frontend conda` to the snakemake command.

Similarly, other Snakemake options are available.
For example, you can specify `--cores all` to use all available cores instead of just a single core or specify any other number of cores.

The Snakemake command includes execution of all Jupyter notebooks we implemented for intermediate data analysis.
The output of these analyses can be found in the `results` folder.

If from the config file, user chose to use the data analysis module with jupyter notebooks, the command for execute the workflow is as below:
```
$ snakemake --use-conda --cores 1 --edit-notebook [OUTPUT_PATH]/popinSnake/results/insertions_genotypes.vcf.gz
```
@NOTICE: If your workflow doesn't wait for the jupyter notebooks to be opened halfway through the analysis, we suggest switching to snakemake version 5.32.0 for a better experience. 
<!-- `[WORKFLOW_PATH]` is as defined in the config file where the main popinSnake workflow is.  -->

## Example data

We provide an `example_data` folder containing examplary data generated for testing purposes.
It also shows you how to prepare your data for the workflow.
The `snake_config.yaml` file is preconfigured to run the workflow on the example data.
You can simply execute Snakemake to run the workflow on this data (see also [Workflow execution](#workflow-execution)):
```
$ snakemake --use-conda --cores 1
```

__Note: When you run the workflow for the first time using `--use-conda`, it will create the conda environments.
This process may take some time. The next time you run the workflow, this process is not necessary again.__

In this example, the FASTA file `chr21_ins.fa` is used as a reference genome.
A subfolder for each sample is present and its name matches the BAM file name.
You can find three BAM files (`.bam`) and the corresponding index files (`.bam.bai`) in the three sample folders `S0001`, `S0002`, `S0003`.

<!-- In addition, a BWA-indexed alternative reference for contamination removal is provided in the subfolder `example_data/virus_ref/`.
However, the default configuration will not use it as the `REMAP` parameter is set to `"no"`.
Change the configuration to `REMAP="yes"` in order to test the remapping sub-workflow on the example data. -->

## Cluster Setup
The popinSnake workflow also supports cluster execution and SLURM scheduling.
The memory and time for each rules can be configured in `cluster_config.yaml`

This configuration file defines resource profiles and threading parameters for various tasks in the workflow. Under the `resources` section, `memory` (in megabytes) and `time` limits (in HH:MM:SS format) are specified for standard setups, as well as dynamic and sample-based allocations that scale with the number of samples being processed. Additional specialized resource configurations are given for rules like email, analysis, kraken, bwa, samtools, and a high-memory variant of samtools, ensuring that each tool or step has appropriate computational limits. Under the threads section, default single-thread usage is defined, as well as parallelization settings for multi-threaded tools, enabling efficient utilization of available CPU cores during resource-intensive steps. This configuration allows the workflow manager to schedule and run tasks efficiently on a cluster or high-performance computing environment.

For the rules use `sample_based` scheduling, user need to define the allowed memory for each sample, the total number of samples will be obtained automatically from the input folder and the final memory will be calculated and defined.

For the rules use `dynamic_schedule` to handle memory-intensive tasks efficiently. Users define the initial memory allocation with the `mem_init` parameter in each rule. If a rule fails due to insufficient memory, Snakemake retries the rule up to the number of times specified by `--restart-times`, e.g. 
```
$ snakemake --restart-times 3 --cores 10
```
doubling the memory allocation on each retry. This ensures tasks with unpredictable memory requirements can adjust dynamically without manual intervention. Make sure your cluster supports increasing memory requests, and verify logs to confirm failures are memory-related.


## References

Cao K., Elfaramawy N., Weidlich M., Kehr B. (2022)
__From Program Chains to Exploratory Workflows: PopinSnake for Insertion Detection in Genomics.__
In: 2023 IEEE 19th International Conference on e-Science (e-Science).

Krannich T., White W. T. J., Niehus S., Holley G., Halldórsson B. V., Kehr B. (2022)
__Population-scale detection of non-reference sequence variants using colored de Bruijn graphs.__
[Bioinformatics, 38(3):604–611](https://academic.oup.com/bioinformatics/article/38/3/604/6415820).

Kehr B., Helgadóttir A., Melsted P., Jónsson H., Helgason H., Jónasdóttir Að., Jónasdóttir As.,	Sigurðsson Á., Gylfason A., Halldórsson G. H., Kristmundsdóttir S., Þorgeirsson G., Ólafsson Í., Holm H., Þorsteinsdóttir U., Sulem P., Helgason A., Guðbjartsson D. F., Halldórsson B. V., Stefánsson K. (2017).
__Diversity in non-repetitive human sequences not found in the reference genome.__
[Nature Genetics,](http://rdcu.be/pDbJ) [49(4):588–593](https://www.nature.com/articles/ng.3801).

Kehr B., Melsted P., Halldórsson B. V. (2016).
__PopIns: population-scale detection of novel sequence insertions.__
[Bioinformatics, 32(7):961-967](https://academic.oup.com/bioinformatics/article/32/7/961/2240308).
