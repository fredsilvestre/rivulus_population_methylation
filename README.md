# 1° Background

This repository contains the scripts, the workflow and the explanations to reproduce the data obtained from a research project about mangrove rivulus (Kryptolebias marmoratus) wild populations. This project is a part of the PhD thesis of Valentine Chapelle which aims to study the epigenetic diversity of DNA methylation in wild populations of this self-fertilizing fish species. Four populations have been sampled in Emerson Point Preserve (EPP - 42 fish) and Long Key (LK - 30 fish) in Florida ; in Long Caye (LC - 28 fish) and Twin Caye (TC - 40 fish) in Belize, in May and July 2019. These four populations are characterized by a gradient of hermaphrodites/males ratio (EPP > LK > LC > TC) and thus a gradient of genetic diversity. The main scientific question is to know whether wild populations with very low genetic diversity can compensate with a higher epigenetic diversity to maintain a variation of phenotypes. The phenotypes have been analyzed at the behavioral level by measuring boldness and exploratory personality traits during the days following sampling. Brain, liver and gonads have then been removed for DNA methylation analysis using Reduced Representation Bisulfite Sequencing (RRBS). 

Directly after dissection, organs were stored in RNA Later Stabilization Solution (Invitrogen), brought back in Belgium and kept at -80°C until DNA extraction. DNA of the organs from the 140 individuals were extracted using Nucleospin Tissue XS kit (Macherey-Nagel). Libraries were prepared from 100 ng of purified DNA using Premimum RRBS kit (Diagenode) and purified twice using Agencourt AMPure XP beds (Beckman Coulter). Libraries were sequenced on two NextSeq550 High Output flow cells (Illumina) as single-end 100 bp reads with a minimum sequencing depth of 10 millions reads per sample. Quality was assessed using FastQC (Babraham Bioinformatics) ; trimming was performed using TrimGalore (Babraham Bioinformatics) ; genome alignment and methylation calling was performed using Bismark (Babraham Bioinformatics). We used the reference genome of Kryptolebias marmoratus assemble version 2 (RefSeq assembly accession: GCF_001649575.2). The Bioconductor R package methylKit was used to identify differentially methylated regions (DMRs) and CpGs (DMCs). Visualization and annotation was processed into SeqMonk (Babraham Bioinformatics). Functional enrichment was performed using the R package gprofiler2. 

## References

Reference book: Computational genomics with R: https://compgenomr.github.io/book/

The workflow is based on Laine et al. (2022) An ecologist's guide for studying DNA methylation variation in wild vertebrates. Mol. Ecol. Resour.

Main references for wild populations:
- van Oers et al (2020) Epigenetics of Animal Personality: DNA Methylation Cannot Explain the Heritability of Exploratory Behavior in a Songbird. Int. Comp. Biol.
- Liu et al. (2020). A comprehensive evaluation of computational tools to identify differential methylation regions using RRBS data. Genomics.
- Wreczycka et al/ (2017). Strategies for analyzing bisulfite sequencing data. J. Biotechnol.
- Akalin et al (2012) methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. Genome Biol.
- Liu et al (2020) A comprehensive evaluation of computational tools to identify di!erential methylation regions using RRBS data. Genomics.
- Blanc et al. (2021) The insecticide permethrin induces transgenerational behavioral changes linked to transcriptomic and epigenetic alterations in zebrafish (Danio rerio). Sci. Tot. Env.
- Watson et al (2021) Urbanization is associated with modifications in DNA methylation in a small passerine bird. Evol Applications.
- Johnson and Kelly (2019) Population epigenetic divergence exceeds genetic divergence in the Eastern oyster Crassostrea virginica in the Northern Gulf of Mexico. Evol Applications.
- Gavery et al (2018) Characterization of Genetic and Epigenetic Variation in Sperm and Red Blood Cells from Adult Hatchery and Natural-Origin Steelhead, Oncorhynchus mykiss. G3.
- Aluru et al (2018) Role of DNA methylation in altered gene expression patterns in adult zebrafish (Danio rerio) exposed to 3, 3’, 4, 4’, 5-pentachlorobiphenyl (PCB 126). Env Epigenetics.
- Wogan et al (2019). Genome-wide epigenetic isolation by environment in a widespread Anolis lizard. Mol Ecol.
- Chapelle (2023). Adaptation and evolution in the mangrove with low genetic diversity: a combined field and laboratory study on DNA methylation variation in the mangrove rivulus Kryptolebias marmoratus. PhD thesis, UNamur.

Additional references:
- Voisin et al. (2022). Genome-wide DNA methylation of the liver reveals delayed effects of early-life exposure to 17-α-ethinylestradiol in the self-fertilizing mangrove rivulus. Epigenetics.
- Chapelle and Silvestre (2022). Population Epigenetics: The Extent of DNA Methylation Variation in Wild Animal Populations. Epigenomes.
- Vogt (2022). Environmental Adaptation of Genetically Uniform Organisms with the Help of Epigenetic Mechanisms—An Insightful Perspective on Ecoepigenetics. Epigenomes.
- Vernaz et al. (2021)	Mapping epigenetic divergence in the massive radiation of Lake Malawi cichlid fishes. Nat Commun. 


# 2° Experimental design

## First questions

- Create an experimental file with all the metadata.
- Think to make the project FAIR (findable, accessible, interoperable and reusable) > provide the data and the codes (Github or Zenodo).
- Know the DNA methylation pattern in the species of interest.
- Check if there is a reference genome.
- Choose the tissues.
- Know the genetic structure of the studied population.
- Choose to pool or not (only when there is low DNA or as exploratory study)

## Choosing the method

- WGBS or RRBS (or other technique ?).
- Choose the sequencing depth (= read Length * number of sequencing reads / genome length) (to have a sufficant coverage per CpG) and the number of replicates > find the good balance. Minimum 15x depth is recommended.
- Single-end or paired-end (paired-end is better to align but might not be necessary if the alignment is good enough with single-end)
- Remove the overlap (when the fragments are bigger than the reads length)
- Directional (only original top and bottom strands in single-end mode) or non-directional libraries (also the complementary strands). TruSeq adaptors are automatically directional. To find out:  https://www.seqanswers.com/forum/applications-forums/sample-prep-library-generation/15699-did-i-make-directional-or-non-directional-bs-seq-libraries?t=18908  

# 3° Preparation of the data

## 1° Cloning Github repository to your local folder

Be sure to note the exact pathway to the folder where you downloaded the clone.Do not forget to return to the parent directory where you cloned the repository before you run each step.

```{r}

git clone https://github.com/fredsilvestre/rivulus_population_methylation

```

Then create the following folders within the cloned repository in order to organize the data.

```{r}

mkdir ./rivulus_population_methylation/data/
        
mkdir ./rivulus_population_methylation/data/fastq

mkdir ./rivulus_population_methylation/data/trimgalore_output

mkdir ./rivulus_population_methylation/data/bismark_output

mkdir ./rivulus_population_methylation/data/csv

mkdir ./rivulus_population_methylation/data/fastqc

mkdir ./rivulus_population_methylation/data/genome

```


## 3° Download the fastq files from bisulfite sequencing


## 4° Quality control of raw sequencing data

Run FASTQC and MULTIQC for multiple comparison on the .fasta.gz files. BSEQC can also be used for more details.
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


## 5° Trimming

Run TRIM GALORE (which is using CUTADAPT) in order to:
	- remove adapters
	- remove low quality nucleotides (< 20)
	- remove the end-repaired C at the extremity (--rrbs)

Possible options are: 
options:  --rrbs ; 
possible add: 
--non_directional
--paired
--clip_r2 2 (to remove 2 bases at 5' from R2)
--output_directory #to redirect to another directory
--quality <INT> #Phread scrore of 20 by deffault or change to 30
--fastqc #to run fastqc after trimming
--phred33 # by deffault
--illumina # if not, will detect adaptors automatically
--stringency <INT> #deffault is 1 > will decrease the A at the end because it is recognized as adaptor. We can use 2 instead
--hardtrim5 <int> # to trim a number of bp at 5'
--length <INT> #by deffault it discard reads smaller than 20bp ; we can change to 40 > it will increase the alignment

```{r}

~/TrimGalore-0.6.6/trim_galore --rrbs --quality 28 --illumina --stringency 2 --length 40 --output_dir ./rivulus_population_methylation/TrimGalore_output ./rivulus_population_methylation/fastq/*.gz

```


## 6° Download the reference genome from NCBI

The last genome assembly of the mangrove rivulus (RefSeq assembly accession: GCF_001649575.2 ) can be found here @ https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_001649575.2/
You can either download it from this website and unzip it, or use the 3 command lines together.

```{r}

cd ./rivulus_population_methylation/data/genome/

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Kryptolebias_marmoratus/latest_assembly_versions/GCF_001649575.2_ASM164957v2/GCF_001649575.2_ASM164957v2_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Kryptolebias_marmoratus/latest_assembly_versions/GCF_001649575.2_ASM164957v2/GCF_001649575.2_ASM164957v2_genomic.fna.gz

gunzip GCF*
```


## 7° Prepare the genome in Bismark

Installation information for Bismark can be found @ https://github.com/FelixKrueger/Bismark 
"Bismark is a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step. The output can be easily imported into a genome viewer, such as SeqMonk, and enables a researcher to analyse the methylation levels of their samples straight away. It's main features are:
Bisulfite mapping and methylation calling in one single step
Supports single-end and paired-end read alignments
Supports ungapped, gapped or spliced alignments
Alignment seed length, number of mismatches etc. are adjustable
Output discriminates between cytosine methylation in CpG, CHG and CHH context"

Change the extension of the genome file into .fa to be recognized by Bismark.
Run the lines below to prepare the genome in Bismark.
```{r}

/Applications/Bismark-0.22.3/bismark_genome_preparation --path_to_aligner /Applications/bowtie/bin/ --verbose ~/rivulus_population_methylation/data/genome/
        
```


## 8° Align sequences to the reference genome using Bismark


```{r}

/Applications/Bismark-0.22.3/bismark --genome ~/rivulus_population_methylation/data/genome/ --score_min L,0,-0.6 --output_dir ~/rivulus_population_methylation/data/genome/bismark_output ~/rivulus_population_methylation/data/genome/trimGalore_output/*.fq.gz

```


## 9° Samtools

The output are sam.gz files that must first be uncompressed.
Bismark writes a report in .txt where we can find many infos such as the average % methylation.

We must first reorganize the  files after Bismark with the samtools package from the terminal. When you align FASTQ files with all current sequence aligners, the alignments produced are in random order with respect to their position in the reference genome. In other words, the BAM file is in the order that the sequences occurred in the input FASTQ files.

Now we read a file from Bismark with the format .SAM
SAM files must be sorted by chromosome and read position columns, using ‘sort’ command in unix-like machines will accomplish such a sort easily. BAM files should be sorted and indexed. This could be achieved with samtools (http://www.htslib.org/doc/samtools.html).
Follow samtools tutorial: http://quinlanlab.org/tutorials/samtools/samtools.html#:~:text=samtools%20%E2%80%9Csort%E2%80%9D&text=In%20other%20words%2C%20the%20BAM,in%20the%20input%20FASTQ%20files.&text=Doing%20anything%20meaningful%20such%20as,occur%20in%20%E2%80%9Cgenome%20order%E2%80%9D.
List of terminal commands: https://www.makeuseof.com/tag/mac-terminal-commands-cheat-sheet/


```{r}

cd ~
mkdir ~/samtools
cd ~/samtools
git clone https://github.com/samtools/htslib
git clone https://github.com/samtools/samtools
cd samtools
make

./samtools/samtools sort .bam -o .sorted.bam

```

