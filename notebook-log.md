# Notebook for Botany 563
(Written by Ari Maurer in Spring 2026) 

## Dataset for class final project
Molecular Phylogenetics and Evolutionary Diversification of Labyrinth Fishes (Perciformes: Anabantoidei) by Lukas Rüber , Ralf Britz , Rafael Zardoya (2006)

Full text PDF: https://www.researchgate.net/publication/6925298_Molecular_Phylogenetics_and_Evolutionary_Diversification_of_Labyrinth_Fishes_Perciformes_Anabantoidei

Dataset available from: 
https://datadryad.org/dataset/doi:10.5061/dryad.59zw3r2c0 

Dataset summary:
Mitochondrial and nuclear nucleotide sequence data was used to construct a phylogeny of anabantoids. The mitochondrial data set included the complete cytochrome b, partial 12S rRNA, complete tRNA Val, and partial 16S rRNA genes (3332 bp) of 57 species representing all 19 anabantoid genera. The nuclear data set included the partial RAG1 gene (1494 bp) of 21 representative species

## All Software Used 
| Software | Description | Strengths | Weaknesses | Assumptions | User Choices |
|---|---|---|---|---|---|
| Muscle | Program for creating multiple alignments of protein sequences. Includes fast distance estimation using k mer counting, progressive alignment, and refinement using tree-dependent restricted partitioning. | Fast and accurate for large datasets | MUSCLE's accuracy drops when aligning highly divergent sequences | MUSCLE assumes using k-mer counting (a fast distance measurement) and then iterative refinment can produce results similiar to computationally expensive methods but faster | User choices TBD |

## Data alignment with MUSCLE
### Installation 
1. Install Miniconda (https://www.anaconda.com/docs/getting-started/miniconda/main)
- Miniconda is a minimal installer for Anaconda that includes only Conda (environment manager), Python, and a few essentail packages
- Miniconda will be used to install MUSCLE
2. Install MUSCLE using Miniconda

```
 conda install muscle
```

3. Confirm MUSCLE installation

```
muscle -version
```
### Run MUSCLE on the anabantoid data

```
muscle -align anabantoid.fasta -output aligned_anabantoid.fasta
```

** Data not given in unaligned FASTA format, data given in already aligned Nexus files, so this step not run. Data given as UCE_loci_75p_complete.gzip   > gzip directory of nexus-formatted alignments for individual UCE loci. All alignments in this directory contain at least 75% of taxa in our total dataset. **

## Distance Tree Calculations
### Installing neccessary packages
```
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
```

### Loading the packages
```
library(ape)
library(adegenet)
library(phangorn)
```
### Loading the data

```
uce-2 <- read.nexus.data("uce-2.nex")
uce-3 <- read.nexus.data("uce-3.nex")
uce-4 <- read.nexus.data("uce-4.nex")
uce-5 <- read.nexus.data("uce-5.nex")
uce-6 <- read.nexus.data("uce-6.nex")
uce-7 <- read.nexus.data("uce-7.nex")
uce-8 <- read.nexus.data("uce-8.nex")
uce-9 <- read.nexus.data("uce-9.nex")
uce-10 <- read.nexus.data("uce-10.nex")
uce-11 <- read.nexus.data("uce-11.nex")
uce-12 <- read.nexus.data("uce-12.nex")
```
### Convert list to DNAbin
```
Uce2bin <- nexus2DNAbin(Uce2)
```

### Computing the genetic distances
Tamura and Nei 1993 model chosen, which allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate (more on Models of Evolution).

```
Uce2D <- dist.dna(Uce2bin, model="TN93")
A <- dist.dna(uce-2, model="TN93")
```

### Get the NJ tree
```
tre2 <- nj(Uce2D)
```

### Ladderize
Before plotting, we can use the ladderize function which reorganizes the internal structure of the tree to get the ladderized effect when plotted
```
tre2L <- ladderize(tre2)
```

### Plot the tree
```
plot(tre2L, cex=.6)
title("UCE-2")
```
