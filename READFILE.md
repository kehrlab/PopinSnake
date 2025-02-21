# popinSnake file documentation

**Popins2 Snakemake Workflow** can discover and characterize non-reference-insertions of 100bp or longer on a population scale.


The workflow takes in a reference file of a DNA sequence `(.fasta/.fa)` and it's index, and many sample files on a population scale. 
The sample files`(.bam)` are assumed to have been aligned to the same reference genome with bwa.

By following the preset rules, the workflow takes as input a reads-to-reference alignment, 
assembles unaligned reads using a standard assembly tool,
merges the contigs of different individuals into high-confidence sequences, 
anchors the merged sequences into the reference genome, 
and finally genotypes all individuals for the discovered insertions.

It is possible to assign alternative reference files for REMAPPING, if you would like to remove the possible contamination (virus/bacteria) in the samples.

**FASTA file** 

`.fasta/.fa` : contains DNA sequences

`.fa.fai` : index file for fasta file

`.fa.amb/.fa.ann/.fa.bwt/.fa.pac/.fa.sa` : more index files for fasta file, required for bwa alignment

**SAM file** 

`.sam` : a tab-delimited text file that contains sequence alignment data.

**BAM file** 

`.bam` : the binary version of a SAM file
            
`.bam.bai` : the index file for BAM file

**YAML file**

`.yaml/.yml` : configuration file or environment file


---

## Initial file structure:

```
$ tree .
.
├── README.md
├── ref
│   ├── chr21_ins.fa
│   ├── chr21_ins.fa.fai
│   ├── S0001
│   │   ├── S0001.bam
│   │   └── S0001.bam.bai
│   ├── S0002
│   │   ├── S0002.bam
│   │   └── S0002.bam.bai
│   ├── S0003
│   │   ├── S0003.bam
│   │   └── S0003.bam.bai
│   └── virus_ref
│       ├── altref_virus_clean.fa
│       ├── altref_virus_clean.fa.amb
│       ├── altref_virus_clean.fa.ann
│       ├── altref_virus_clean.fa.bwt
│       ├── altref_virus_clean.fa.pac
│       └── altref_virus_clean.fa.sa
├── results
├── unmapped
└── workflow
    ├── envs
    │   ├── environment.yaml
    │   └── env_py2.yaml
    ├── popins2_snake_structure.py
    ├── snake_config.yaml
    ├── Snakefile
    └── subworkflow_remapping
        ├── snake_config.yaml
        └── Snakefile
     
```
---


## Snakefile - Main Workflow

Along with the snakefile always comes with a configfile where we create the configuration to make the workflow more flexible. In our configfile we define a dictionary of configuration parameters and their values.

```
configfile: "snake_config.yaml"
```
in `configfile` we define the following parameters:

`POPINS2_BIN`: the directory where your popins2 program is located in
        
    /home/user/popins4snake/popins2

`REFERENCE`: the directory where your reference fo choice is located in
   
    /home/user/popins2_snake/ref/chr21_ins.fa

`WORK_DIR`: working directory, where you wish to run your workflow and store all the output files

    /home/user/popins2_snake

`ASSEMBLER`: use the assembler of your choice, possible options are "velvet" or "minia"

`REMAP`: whether you wish to include the **REMAPPING** subworkflow, "yes" or "no"

`SAMPLES`: Your choice of sample names according to sample files

<br/>

----
<br/>
Snakemake’s language is similar to that of standard Python syntax and works backwards by requesting output files and defining each step required to produce them

from the input_files we define all the files we expect to see when the workflow finishes execution:

Files to be expected in `{WORK_DIR}/unmapped/{SAMPLES}/`:
- contig_mapped_unsorted.sam, * *temporary file*
- contig_mapped_unsorted.bam,
- contig_mapped.bam,
- non_ref_new.bam,
- non_ref_new.bam.bai,
- {assemblr}.contigs.fa
- locations.txt,
- locations_unplaced.txt,
- locations_placed.txt

If no REMAPPING, the additional files are:
- single.fastq,
- paired.1.fastq,
- paired.2.fastq,
- non_ref.bam,
- POPINS_SAMPLE_INFO

Files to be expected in `{WORK_DIR}/results/`:    
- supercontigs.fa,
- supercontigs.fa.amb,
- supercontigs.fa.ann,
- supercontigs.fa.bwt,
- supercontigs.fa.pac,
- supercontigs.fa.sa,
- locations.txt,
- insertions.vcf,
- groups.txt,
- place-finish.done,
- genotype-finish.done

</br>

---

</br>

## Introducing the REMAPPING sub-workflow

</br>

How to include the REMAPPING subworkflow into the snakefile:

If in the `configfile` we define REMAP parameter as yes, the main workflow will try to include the REMAPPING subworkflow in the pipeline. This main workflow will run the subworkflow before its own execution. the purpose of the subworkflow is to map samples to an alternative reference to eliminate potential contaminations from the sample files. This implement can reduce the error of NRV calling by reducing the noise and avoid inaccurate mapping.

The `Snakefile` of the subworkflow for remapping is located on filepath `/home/user/popins2_snake/workflow/subworkflow_remapping/` along with its `configfile`

</br>

----

</br>

## Snakemake Rules:

<br/>

### \# **Rule \<popins2_crop_unmapped>** (no REMAP) :

**Input file** located on filepath `{WORK_DIR}}/ref/{sample}/`: 
- `{sample}.bam`

    BAM files of samples, containing alignment data. assume they are pre-aligned to the reference file
    bam files contains 
    - Header
    - Alignment Information: QName Flag Chr Pos MQ Cigar MChr MPos ISize Seq Qual Tags
 
**Output files** located on filepath`{WORK_DIR}/unmapped/{sample}/` :
 - `single.fastq`
    - One read from a read-pair that's not mapped to the reference genome but has its paired read mapped to the reference stored in `mates.bam`
 - `paired.1.fastq`             
    - one sequence from the the read-pairs that's not mapped to the reference but has its paired read  in `paired.2.fastq`
 - `paired.2.fastq`
    - one sequence from the the read-pair that's not mapped to the reference but has its paired read in `paired.1.fastq`
 - `mastes.bam`
    - the reference-aligned end of a read-pair, with its unaligned end in `single.fastq` 
 - `POPINS_SAMPLE_INFO`          

</br>

### \# **Rule \<popins2_sort>** :
**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `mates.bam`

**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `non_ref.bam`

    - This output comes from `mates.bam` after being sorted by Qname using **samtools**
        
</br>     

### \# Rule **\<popins2_sickle>**:

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:

- `single.fastq`, 
- `paired.1.fastq`, 
- `paired.2.fastq`

**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `filtered.single.fastq`, 
- `filtered.paired.1.fastq`, 
- `filtered.paired.2.fastq`,
- `filtered.single2.fastq` \* *temporary file*

    - Most modern sequencing technologies produce reads that have deteriorating quality towards the 3'-end and some towards the 5'-end as well. We use **sickle** to filter and trim reads to improve assemble quality

</br>

### \# Rule **\<popins2_sickle>** (after REMAP) :

**Input files** located on `{WORK_DIR}}/unmapped/{sample}/`: 

- `remapped_single.fastq`, 
- `remapped_paired.1.fastq`, 
- `remapped_paired.2.fastq`
    - These new input files come from the REMAPPING subworkflow

**Output files** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `filtered.single.fastq`, 
- `filtered.paired.1.fastq`, 
- `filtered.paired.2.fastq`,
- `filtered.single2.fastq`, \* *temporary file*
    - Quality filtering/trimming with **sickle**
 
</br>      

### \# Rule **\<popins2_minia>** or **\<popins2_velvet>** :
**Input files** located on `{WORK_DIR}}/unmapped/{sample}/`:     
- `filtered.single.fastq`, 
- `filtered.paired.1.fastq`, 
- `filtered.paired.2.fastq`

**Output files** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `{assemblr}.contigs.fa`
    - The `{assemblr}` could be either minia or velvet, which is defined previously in the `configfile`
    - Assemble all unmapped reads and read-pairs into contigs using `minia` or `velvet`. The assembly tool of selection computes a set of contigs when given an individual’s set of unaligned reads

</br>

### \# Rule **\<popins2_merge_contigs>** : 
**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `{assemblr}.contigs.fa`

**Output files** located on `{WORK_DIR}/results/`:
- `supercontigs.fa`,
- `supercontigs.gfa`,
- `supercontigs.bfg_colors`,
- `popins2.merge.log`
    - This rule merges sets of contigs into single, long contigs, which we call supercontigs. A supercontig represents a contig if the full length of the contig aligns to the supercontig. Ideally supercontigs satisfy two points: 
         1. Every input contig should be represented by a supercontig in the set
         2. Any two supercontigs in the set cannot be further merged 

</br>

### \# Rule **\<index_supercontigs>**:
**Input file** located on `{WORK_DIR}/results/`:
- `supercontigs.fa`

**Output files** located on `{WORK_DIR}/results/`:
- `supercontigs.fa.amb`,
- `supercontigs.fa.ann`,
- `supercontigs.fa.bwt`,
- `supercontigs.fa.pac`,
- `supercontigs.fa.sa`
    - This rule creates index files for supercontigs. The index files are for bwa alignment laThe input files located in ter in the workflow.

</br>

### \# Rule **\<map_supercontigs>** (after REMAP):
**Input file** located on `{WORK_DIR}/results/`:
- `supercontigs.fa.amb`,
- `supercontigs.fa.ann`,
- `supercontigs.fa.bwt`,
- `supercontigs.fa.pac`,
- `supercontigs.fa.sa`,
- `supercontigs.fa`,

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `remapped_single.fastq`,
- `remapped_paired.1.fastq`,
- `remapped_paired.2.fastq`
    * These files come from the REMAPPING subworkflow
        
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `contig_mapped_unsorted.sam`, \**temporary file*
    - Align supercontigs with the unmapped (also remapped) read-pairs and store alignment information   
       
</br>

### \# Rule **\<map_supercontigs>** (no REMAP):
**Input file** located on `{WORK_DIR}/results/`:
- `supercontigs.fa.amb`,
- `supercontigs.fa.ann`,
- `supercontigs.fa.bwt`,
- `supercontigs.fa.pac`,
- `supercontigs.fa.sa`,
- `supercontigs.fa`,

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `single.fastq`,
- `paired.1.fastq`,
- `paired.2.fastq`   
        
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `contig_mapped_unsorted.sam`
    - Align unmapped read-pairs to the supercontigs and store alignment information. The output SAM file is only temporarily stored so it won't take up too much space            

</br>

### \# Rule **\<name_sort_unsorted**:

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:        
- `contig_mapped_unsorted.sam`
      
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `contig_mapped_unsorted.bam`,
    - convert SAM file into BAM file
- `contig_mapped.bam`
    - sort the aligned reads by read names



</br>

### \# Rule **\<remap_merge_set_mate>** (after REMAP):

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`: 
- `contig_mapped.bam`
- `remapped_non_ref.bam`
    - This file comes from REMAP subworkflow
   
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `merged.bam`


</br>

### \# Rule **\<popins2_merge_set_mate>** (no REMAP):

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`: 
- `contig_mapped.bam`,
- `non_ref.bam`
  
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `merged.bam`
    - As one end of the read-pairs that are aligned to reference genome and the other end of the read-paris are aligned to the supercontigs. It is time to anchor the two end of the read-pairs with each other and put the pairs together. 

</br>

### \# Rule **\<coordinate_sort_unsorted>**:
   
**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `merged.bam`
   
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `non_ref_new.bam`
    - sort the paired reads by their relative location on the reference or supercontigs
    
</br>

### \# Rule **\<index_sorted>**:
   
**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `non_ref_new.bam`

**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `non_ref_new.bam.bai`
    - generate index files for the bam file of location sorted reads


</br>

### \# Rule **\<remap_find_locations>** (after REMAP):

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `non_ref_new.bam`,
- `remapped_non_ref.bam`
    - This file comes from REMAP subworkflow
        
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `locations.txt`
    - find contig locations relative to the reference

</br>

### \# Rule **\<popins2_find_locations>**(no REMAP):
        
**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `non_ref_new.bam`
- `non_ref.bam`
   
**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `locations.txt`
    - find contig locations relative to the reference
    - the files contain information about: location, contig information, forward or reverse strand ...

</br>

### \# Rule **\<popins2_place_refalign>**:
 
**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `non_ref_new.bam`,
- `non_ref_new.bam.bai`,
- `locations.txt`

**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `locations_unplaced.txt`

**Output file** located on `{WORK_DIR}/results/`:
- `locations.txt`,
- `insertions.vcf`,
- `groups.txt`,
    - This rule use the alignment information and relative locations of the contigs to place contigs along the reference genome. If a contig is able to align one end with the reference and not the other end, such contigs are stored as insertions



</br>

### \# Rule **\<popins2_place_splitalign>**:

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `locations_unplaced.txt`

**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `locations_placed.txt`
    - Use split reads information to anchor unplaced contigs to their precise location on reference genome


</br>

### \# Rule **\<popins2_place_finish>**:

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:        - `locations_placed.txt`

**Output file** located on `{WORK_DIR}/results/`:
- `place-finish.done`
    - No actual file created by this rule. Snakemake provides a ***touch*** function to check whether the previous rules are completed.
    

</br>

### \# Rule **\<popins2_genotype>**:

**Input file** located on `{WORK_DIR}/results/`:        
- `place-finish.done`

**Output file** located on `{WORK_DIR}}/unmapped/{sample}/`:
- `insertions.vcf`
    - based on the alignment score, the genotype of each insertion can be calculated. The file contains each insertion variation, their position, and genotype information

</br>

### \# Rule **\<popins2_check_genotype_status>**:

**Input file** located on `{WORK_DIR}}/unmapped/{sample}/`:        
- `insertions.vcf`

**Output file** located on `{WORK_DIR}/results/`:
- `genotype-finish.done`