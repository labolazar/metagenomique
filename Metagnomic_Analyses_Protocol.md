---
title: "Protocol for metagenomic analyses"
author: "Karine Villeneuve"
date: "`r Sys.Date()`"
output: rmdformats::robobook
---

# Pre-Assembly 

## Unzip 
Open a terminal window and move to the folder with your raw `fastq.gz`. 
```{bash, eval=FALSE}
gunzip *.gz
```

## Interleaving 
I used the script <font color='green'>interleave.py</font> from this [GitHub gist](https://gist.github.com/ngcrawford/2232505). The script is currently available in the folder `Script Métagénomique` of our Teams (see `General`). Transfer the script to the folder with the fastq files on the server, make the script executable and run the script using python. 

Make the script executable and chose option **a** or **b** depending on the amount of fastq you have. 
```{bash, eval=FALSE}
chmod +x interleave.py
```

a. If you have only a few fastq use this simple command line. 
```{bash, eval=FALSE}
python3 interleave.py file_R1.fastq file_R1_fastq > interleaved.fastq
```

b. If you have multiple fastq use this loop with nohup and a bash script.
```{bash, eval=FALSE}
for R1 in *_R1_001.fastq ; do python3 interleave_fastq.py $R1 "${R1/R1/R2}"  > $R1.interleave.fastq ; done 
```

View the top of the new interleaved files to make sure your reads alternate between R1 and R2. 
```{bash, eval=FALSE}
grep @M interleaved.fastq | head
```


I recommend moving all the interleaved fastq to a new folder (I called the new folder `interleaved`)

<font color='red'>*I am looking for a better way to rename the interleave file from the loop. But this will do for now.*</font>

## Sickle 

Sickle is a tool that uses sliding windows along with quality and length thresholds to determine when quality is sufficiently low to trim the 3’-end of reads and also determines when the quality is sufficiently high enough to trim the 5’-end of reads. It generates two files for each fastq you input : `.trim.fastq` and `.singles.fastq`. We only need the `.trim.fastq` and therefor you can move all the others to another folder. 

Running sickle with a bash script. 

```{bash, eval=FALSE}
#!/bin/bash
for i in *.fastq
  do sickle pe -c $i -t sanger -m $i.trim.fastq -s $i.singles.fastq
done
```
*sickle pe (paired end) -c (inputfile) -t sanger (from illumina) -m (outputfilename) -s (exclutedreadsfilename)*


## Quality check with Fastqc 

Run fastqc on the output (trimmed) file and the non-trimmed file. The process can take a lot of time so a recommend using `nohup`. 

```{bash, eval=FALSE}
#!/bin/bash 
fastqc *.fastq --outdir=/home/kvilleneuve/Metagenomic_analyses/interleaved/fastqc
``` 
Transfer the HTML files from the fastqc directory to your computer in order to view them in your a browser. 

## Fastq to Fasta 

Gzip the combined_trimmed.fastq file using and then convert them to fasta. I recommend using nohup and a bash script for both command since it can be very long. 
```{bash, eval=FALSE}
#!/bin/bash
gzip *.fastq
```

```{bash, eval=FALSE}
#!/bin/bash
for i in *.gz ; 
  do seqtk seq -a $i > $i.fa ; 
done 
```

***

# Assembly 

For the moment, there is no way of knowing which assembly program is best suited for your sample. Therefor, I try the different assemblers and look at (1) the output number of contigs, (2) the number of contigs > 2000 bp and (3) the number of contigs > 1000 bp. Because the assembly process is quite long I always use nohup. 

## IDBA 

```{bash, eval=FALSE}
#!/bin/bash
for i in *.fa
  do idba_ud -l $i -o $i.assembly --pre_correction --mink 65 --maxk 115 --step 10 --seed_kmer 55 i--num_threads 20
done
```

The output is a new folder ending in `assembly` for each of our file. In this folder we are interested in the file called `contig.fa`. 

## Megahit 
[Github](https://github.com/voutcn/megahit)
```{bash, eval=FALSE}
#!/bin/bash
for i in *.fa
  do megahit --12 $i --k-list 21,33,55,77,99,121 --min-count 2 --verbose -t 25 -o ~/Metagenomic_analyses/fasta/meg_$i --out-prefix megahit_$i  
done
``` 

## SPAdes
[Github](https://github.com/voutcn/megahit)

The lastest Python version accepted is 3.5, therefore you have to create a virtual environment if you Python version is newer to use SPAdes v 3.9.0. The process takes a few minutes, nohup can be used but is optional. 

```{bash, eval=FALSE}
git clone https://github.com/ablab/spades.git
git checkout remotes/origin/spades_3.9.0

conda create -n envPython python=3.5
path/to/python3.5 -m venv path/to/virtual/env
source path/to/virtual/env/bin/activate  

spades.py --12 contig.fa --meta -o /home/.../.../SPAdes

deactivate

``` 
*For SPAdes v3.9.0, there is an error in ./spades_compile.sh when compiling the source code, but it is easily fixed by adding the following line in spades/assembler/src/modules/dev_support/segfault_handler.hpp

```{bash, eval=FALSE}
#include <functional> 
``` 


***


# Post-Assembly stats 

## Number of contigs 
```{bash, eval=FALSE}
grep -c ">" contig.fa 
```
## Lenght of contigs and histogram
Use seqkit to extract the length of every contigs. Remove the first row of the document lengths (word length) using vi and use pipeline to create histogram with the lengths file. 
```{bash, eval=FALSE}
seqkit fx2tab --length --name --header-line contig.fa > length.tab 
# Open length.tab  in vi or nano and delete the first row 
cut -f 2 length.tab > lengths
```

<font color='red'>*Eventually remove this*</font> Replace with a script giving : **1**. The total amount of sequence **2**. Sequence > 2kbp **3**. Sequences > 1kbp **4**. Threshold for GC content ? 
```{bash, eval=FALSE}
less lengths | Rscript -e 'data=abs(scan(file="stdin")); png("seq.png"); hist(data,xlab="sequences", breaks=250, xlim=c(0, 5000))'
```
The output is a png file called “seq.png”. If x axis of the histogram is not right change the xlim=c(x,x) values. Copy the file to your local computer to view it (in your local computer terminal navigate to the local directory where you want the file to be copied)
## Lenght and GC 
High GC organisms tend not to assemble well and may have an uneven read coverage distribution. I used this modified script `length+GC.pl`
```{bash, eval=FALSE}
perl length+GC.pl contig.fa > contig_GC.txt
```
## Keeping sequences above 2000 base pairs (2kbp -kilobase pairs)  
```{bash, eval=FALSE}
perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if(length($fa{$s})>2000) }}' contig.fa > 2000bp.fa
``` 

# Binning 
Coverage-based binning approaches will require you to map reads to assembled contigs

## Mapping 
[Github](https://github.com/imrambo/genome_mapping)
This wrapper maps **FASTQ** reads against an assembly (e.g. genome) in **FASTA** format using BWA-MEM
For each sample, create a folder and copy into this folder these two files : the **long contig.fa** (1000 or 2000 bp) and the **trim and interleaved fastq** (fastq after sickle)
1. Go to your home directory and clone the git repository 
```{bash, highlight=TRUE, eval=FALSE}
git clone https://github.com/imrambo/genome_mapping.git
```
2. Go to the new directory created `genome_mapping° and create a conda environment 
```{bash, highlight=TRUE, eval=FALSE}
conda env create -f environment_linux64.yml
```
3. Activate the conda environement
```{bash, highlight=TRUE, eval=FALSE}
conda activate scons_map
```
3. Run a dry run to ensure everything runs smoothly  
```{bash, highlight=TRUE, eval=FALSE}
scons  --dry-run --fastq_dir=/home/karine/VP/megahit/concat --assembly=/home/karine/VP/megahit/concat/concatenated_1kb.fa --outdir=/home/karine/VP/megahit/concat/map --sampleids=fastq_concat.fastq --align_thread=5 --samsort_thread=5 --samsort_mem=768M --nheader=8 --tmpdir=/home/karine/tmp --logfile=concat.log
```
4. Run the script 
```{bash, highlight=TRUE, eval=FALSE}
scons  --fastq_dir=/home/karine/VP/SV10 --assembly=/home/karine/VP/SV10/SV10_1kb.fa --outdir=/home/karine/VP/SV10_map --sampleids=SV10_combined --align_thread=5 --samsort_thread=5 --samsort_mem=768M --nheader=8 --tmpdir=/home/karine/tmp --logfile=SV10log
```
In the outpir directory created  you will find a file that ends with `.sorted.bam` which we will need to create the depth file 
## Depth file {.tabset}
The depth allows you to know how many sequence you can align with certain sections of your contigs. Section with very little depth (few sequences) are not reputable to use
## MetaBAT 
[MetaBAT](https://peerj.com/articles/1165/)
Efficient tool for accurately reconstructing single genomes from complex microbial communities
a. Create `vi` file named `run_metabat.sh` with the following command
```{bash, highlight=TRUE, eval=FALSE}
#!/bin/bash
metabat2 -i 2000kb.fa -a depth.txt -o bins_dir/bin -t 20 --minCVSum 0 --saveCls -d -v --minCV 0.1 -m 2000
```
`minCVsum` : assigning number of tetranucleotide frequency graphs, don’t grab negative numbers 
`-m` : min size of contig to be considered for binning
b. Run Metabat using `nohup ./`and `&`
The output is a folder called `bins_dir` containing all the bins created 


# Bin quality 

# Bin cleaning 

# Taxonomy 

# Phylogenetic tree 

# Metabolic pathway 

# Other information 

## Workflow (alternative)

[iMetAMOS](https://metamos.readthedocs.io/en/v1.5rc3/content/workflows.html)

## Articles 

New approaches for metagenome assembly with short reads. [Article](https://academic.oup.com/bib/article/21/2/584/5363831) by Ayling *et al.*, 2019. 




