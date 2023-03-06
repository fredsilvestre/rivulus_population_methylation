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
