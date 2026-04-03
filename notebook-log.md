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
|RAxML|RAxML is an open-source software that estimates maximum likelihood phylogenies|RAxML typically finds higher-scoring trees than IQTree|IQTree typically finds trees with better stability|RAxML uses randomized stepwise addition parsimony trees. RAxML does model support for you, and supports 22 classical GTR models. RAxML can support DNA and protein data, as well as binary, multi-state morphological and RNA secondary structure data. RAxML can correct for ascertainment bias.|All model parameters can either be optimized or fixed to user-specified values. User needs to check what model RAxML chooses. |
|IQtree|Software for phylogenetic inference using maximum likelihood|Utilizes fast model selection via ModelFinder. IQtree faster than RAxML in optimization of model parameters|---|Two starting trees computed using maximum parsimony and neighbor joining which are then optimized by hill climbing nearest neighbor interchange moves. Uses stochastic NNI to escape local optima.|User choices TBD|

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
--> I chose the first 10 UCE files (2, 3, 6, 7, 8, 9, 11, 13, 14, 15) from the data set to make trees for 

### Set working directory in terminal then open R 

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
Uce2 <- read.nexus.data("uce-2.nexus")
Uce3 <- read.nexus.data("uce-3.nexus")
Uce6 <- read.nexus.data("uce-6.nexus")
Uce7 <- read.nexus.data("uce-7.nexus")
Uce8 <- read.nexus.data("uce-8.nexus")
Uce9 <- read.nexus.data("uce-9.nexus")
Uce11 <- read.nexus.data("uce-11.nexus")
Uce13 <- read.nexus.data("uce-13.nexus")
Uce14 <- read.nexus.data("uce-14.nexus")
Uce15 <- read.nexus.data("uce-15.nexus")
```
### Convert list to DNAbin
```
Uce2bin <- nexus2DNAbin(Uce2)
Uce3bin <- nexus2DNAbin(Uce3)
Uce6bin <- nexus2DNAbin(Uce6)
Uce7bin <- nexus2DNAbin(Uce7)
Uce8bin <- nexus2DNAbin(Uce8)
Uce9bin <- nexus2DNAbin(Uce9)
Uce11bin <- nexus2DNAbin(Uce11)
Uce13bin <- nexus2DNAbin(Uce13)
Uce14bin <- nexus2DNAbin(Uce14)
Uce15bin <- nexus2DNAbin(Uce15)
```

### Computing the genetic distances
Tamura and Nei 1993 model chosen, which allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate (more on Models of Evolution).

```
Uce2D <- dist.dna(Uce2bin, model="TN93")
Uce3D <- dist.dna(Uce3bin, model="TN93")
Uce6D <- dist.dna(Uce6bin, model="TN93")
Uce7D <- dist.dna(Uce7bin, model="TN93")
Uce8D <- dist.dna(Uce8bin, model="TN93")
Uce9D <- dist.dna(Uce9bin, model="TN93")
Uce11D <- dist.dna(Uce11bin, model="TN93")
Uce13D <- dist.dna(Uce13bin, model="TN93")
Uce14D <- dist.dna(Uce14bin, model="TN93")
Uce15D <- dist.dna(Uce15bin, model="TN93")
```

### Get the NJ tree
```
tre2 <- nj(Uce2D)
tre3 <- nj(Uce3D)
tre6 <- nj(Uce6D)
tre7 <- nj(Uce7D)
tre8 <- nj(Uce8D)
tre9 <- nj(Uce9D)
tre11 <- nj(Uce11D)
tre13 <- nj(Uce13D)
tre14 <- nj(Uce14D)
tre15 <- nj(Uce15D)
```

### Ladderize
Before plotting, we can use the ladderize function which reorganizes the internal structure of the tree to get the ladderized effect when plotted
```
tre2L <- ladderize(tre2)
tre3L <- ladderize(tre3)
tre6L <- ladderize(tre6)
tre7L <- ladderize(tre7)
tre8L <- ladderize(tre8)
tre9L <- ladderize(tre9)
tre11L <- ladderize(tre11)
tre13L <- ladderize(tre13)
tre14L <- ladderize(tre14)
tre15L <- ladderize(tre15)
```

### Plot the tree
```
plot(tre2L, cex=.6)
title("UCE-2")

plot(tre3L, cex=.6)
title("UCE-3")

plot(tre6L, cex=.6)
title("UCE-6")

plot(tre7L, cex=.6)
title("UCE-7")

plot(tre8L, cex=.6)
title("UCE-8")

plot(tre9L, cex=.6)
title("UCE-9")

plot(tre11L, cex=.6)
title("UCE-11")

plot(tre13L, cex=.6)
title("UCE-13")

plot(tre14L, cex=.6)
title("UCE-14")

plot(tre15L, cex=.6)
title("UCE-15")
```

## Parsimony Tree Calculations
--> I chose the first 10 UCE files (2, 3, 6, 7, 8, 9, 11, 13, 14, 15) from the data set to make trees for 

### Set working directory in terminal then open R 

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
Uce2 <- read.nexus.data("uce-2.nexus")
Uce3 <- read.nexus.data("uce-3.nexus")
Uce6 <- read.nexus.data("uce-6.nexus")
Uce7 <- read.nexus.data("uce-7.nexus")
Uce8 <- read.nexus.data("uce-8.nexus")
Uce9 <- read.nexus.data("uce-9.nexus")
Uce11 <- read.nexus.data("uce-11.nexus")
Uce13 <- read.nexus.data("uce-13.nexus")
Uce14 <- read.nexus.data("uce-14.nexus")
Uce15 <- read.nexus.data("uce-15.nexus")
```

### Convert list to DNAbin
```
Uce2bin <- nexus2DNAbin(Uce2)
Uce3bin <- nexus2DNAbin(Uce3)
Uce6bin <- nexus2DNAbin(Uce6)
Uce7bin <- nexus2DNAbin(Uce7)
Uce8bin <- nexus2DNAbin(Uce8)
Uce9bin <- nexus2DNAbin(Uce9)
Uce11bin <- nexus2DNAbin(Uce11)
Uce13bin <- nexus2DNAbin(Uce13)
Uce14bin <- nexus2DNAbin(Uce14)
Uce15bin <- nexus2DNAbin(Uce15)
```

### Convert to phanghorn object
```
Uce2P <- as.phyDat(Uce2bin)
Uce3P <- as.phyDat(Uce3bin)
Uce6P <- as.phyDat(Uce6bin)
Uce7P <- as.phyDat(Uce7bin)
Uce8P <- as.phyDat(Uce8bin)
Uce9P <- as.phyDat(Uce9bin)
Uce11P <- as.phyDat(Uce11bin)
Uce13P <- as.phyDat(Uce13bin)
Uce14P <- as.phyDat(Uce14bin)
Uce15P <- as.phyDat(Uce15bin)
```
### Create a starting tree for the search on tree space and compute the parsimony score of this tree 
```
tre2.ini <- nj(dist.dna(Uce2bin,model="raw"))
parsimony(tre2.ini, Uce2P)

tre3.ini <- nj(dist.dna(Uce3bin,model="raw"))
parsimony(tre3.ini, Uce3P)

tre6.ini <- nj(dist.dna(Uce6bin,model="raw"))
parsimony(tre6.ini, Uce6P)

tre7.ini <- nj(dist.dna(Uce2bin,model="raw"))
parsimony(tre7.ini, Uce7P)

tre8.ini <- nj(dist.dna(Uce8bin,model="raw"))
parsimony(tre8.ini, Uce8P)

tre9.ini <- nj(dist.dna(Uce9bin,model="raw"))
parsimony(tre9.ini, Uce9P)

tre11.ini <- nj(dist.dna(Uce11bin,model="raw"))
parsimony(tre11.ini, Uce11P)

tre13.ini <- nj(dist.dna(Uce13bin,model="raw"))
parsimony(tre13.ini, Uce13P)

tre14.ini <- nj(dist.dna(Uce14bin,model="raw"))
parsimony(tre14.ini, Uce14P)

tre15.ini <- nj(dist.dna(Uce15bin,model="raw"))
parsimony(tre15.ini, Uce15P)
```

### Search for the tree with maximum parsimony 
```
tre2.pars <- optim.parsimony(tre2.ini, Uce2P)
tre3.pars <- optim.parsimony(tre3.ini, Uce3P)
tre6.pars <- optim.parsimony(tre6.ini, Uce6P)
tre7.pars <- optim.parsimony(tre7.ini, Uce7P)
tre8.pars <- optim.parsimony(tre8.ini, Uce8P)
tre9.pars <- optim.parsimony(tre9.ini, Uce9P)
tre11.pars <- optim.parsimony(tre11.ini, Uce11P)
tre13.pars <- optim.parsimony(tre13.ini, Uce13P)
tre14.pars <- optim.parsimony(tre14.ini, Uce14P)
tre15.pars <- optim.parsimony(tre15.ini, Uce15P)
```

### Plot tree
```
plot(tre2.pars, cex=0.6)
plot(tre3.pars, cex=0.6)
plot(tre6.pars, cex=0.6)
plot(tre7.pars, cex=0.6)
plot(tre8.pars, cex=0.6)
plot(tre9.pars, cex=0.6)
plot(tre11.pars, cex=0.6)
plot(tre13.pars, cex=0.6)
plot(tre14.pars, cex=0.6)
plot(tre15.pars, cex=0.6)
```

## Maximum Likelihood Tree Calculations
--> I chose the first 10 UCE files (2, 3, 6, 7, 8, 9, 11, 13, 14, 15) from the data set to make trees for 

## RAxML 
###  Download RAxML
https://github.com/amkozlov/raxml-ng 

### Convert Data
RAxML cannot read Nexus format. Nexus files can be converted to PHYLIP format, which RAxML can read using R:

```
library (ape)
alignment2 <- read.nexus.data("uce-2.nexus")
write.dna(alignment, file = "uce-2.phy", format = "sequential")

alignment3 <- read.nexus.data("uce-3.nexus")
write.dna(alignment, file = "uce-3.phy", format = "sequential")

alignment6 <- read.nexus.data("uce-6.nexus")
write.dna(alignment, file = "uce-6.phy", format = "sequential")

alignment7 <- read.nexus.data("uce-7.nexus")
write.dna(alignment, file = "uce-7.phy", format = "sequential")

alignment8 <- read.nexus.data("uce-8.nexus")
write.dna(alignment, file = "uce-8.phy", format = "sequential")

alignment9 <- read.nexus.data("uce-9.nexus")
write.dna(alignment, file = "uce-9.phy", format = "sequential")

alignment11 <- read.nexus.data("uce-11.nexus")
write.dna(alignment, file = "uce-11.phy", format = "sequential")

alignment13 <- read.nexus.data("uce-13.nexus")
write.dna(alignment, file = "uce-13.phy", format = "sequential")

alignment14 <- read.nexus.data("uce-14.nexus")
write.dna(alignment, file = "uce-14.phy", format = "sequential")

alignment15 <- read.nexus.data("uce-15.nexus")
write.dna(alignment, file = "uce-15.phy", format = "sequential")
```

### Run RAxML to calculate the ML tree

```
raxml-ng --msa uce-2.phy --model LG+G8+F

raxml-ng --msa uce-3.phy --model LG+G8+F

raxml-ng --msa uce-6.phy --model LG+G8+F

raxml-ng --msa uce-7.phy --model LG+G8+F

raxml-ng --msa uce-8.phy --model LG+G8+F

raxml-ng --msa uce-9.phy --model LG+G8+F

raxml-ng --msa uce-11.phy --model LG+G8+F

raxml-ng --msa uce-13.phy --model LG+G8+F

raxml-ng --msa uce-14.phy --model LG+G8+F

raxml-ng --msa uce-15.phy --model LG+G8+F
```

### Run both ML and Bootstrap

```
raxml-ng --all --msa uce-2.fasta --model LG+G8+F --bs-trees 10 --uce-2-raxml-boostrap
```

## IQtree
### Download IQtree 
https://iqtree.github.io/ 

### Input data 
```
iqtree3 -s uce-2.nexus
iqtree3 -s uce-3.nexus
iqtree3 -s uce-6.nexus
iqtree3 -s uce-7.nexus
iqtree3 -s uce-8.nexus
iqtree3 -s uce-9.nexus
iqtree3 -s uce-11.nexus
iqtree3 -s uce-13.nexus
iqtree3 -s uce-14.nexus
iqtree3 -s uce-15.nexus
```

### Check plot in R
```
library(ape)
tre2 = read.tree(file="uce-2.nexus.treefile")
plot(tre2)

library(ape)
tre3 = read.tree(file="uce-3.nexus.treefile")
plot(tre3)

library(ape)
tre6 = read.tree(file="uce-6.nexus.treefile")
plot(tre6)

library(ape)
tre7 = read.tree(file="uce-7.nexus.treefile")
plot(tre7)

library(ape)
tre8 = read.tree(file="uce-8.nexus.treefile")
plot(tre8)

library(ape)
tre9 = read.tree(file="uce-9.nexus.treefile")
plot(tre9)

library(ape)
tre11 = read.tree(file="uce-11.nexus.treefile")
plot(tre11)

library(ape)
tre13 = read.tree(file="uce-13.nexus.treefile")
plot(tre13)

library(ape)
tre14 = read.tree(file="uce-14.nexus.treefile")
plot(tre14)

library(ape)
tre15 = read.tree(file="uce-15.nexus.treefile")
plot(tre15)
```

### Root tree
```
plot(tre2)
nodelabels()
rtre = root(tre2, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre3)
nodelabels()
rtre = root(tre3, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre6)
nodelabels()
rtre = root(tre6, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre7)
nodelabels()
rtre = root(tre7, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre8)
nodelabels()
rtre = root(tre8, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre9)
nodelabels()
rtre = root(tre9, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre11)
nodelabels()
rtre = root(tre11, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre13)
nodelabels()
rtre = root(tre13, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre14)
nodelabels()
rtre = root(tre14, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)

plot(tre15)
nodelabels()
rtre = root(tre15, node=151, resolve.root=TRUE)
plot(rtre, type = "phylogram", cex = 0.3, no.margin = TRUE)
```

### Close R
### Quantify support for the estimated tree- check what best fit model is for each on previous step
```
iqtree -s uce-2.nexus -m TPM3u+I+R4 -b 10 -pre uce-2.nexus-iqtree-bootstrap

iqtree -s uce-3.nexus -m XXXXXXXX T -b 10 -pre uce-3.nexus-iqtree-bootstrap

iqtree -s uce-6.nexus -m XXXXXXXX -b 10 -pre uce-6.nexus-iqtree-bootstrap

iqtree -s uce-7.nexus -m XXXXXXXX -b 10 -pre uce-7.nexus-iqtree-bootstrap

iqtree -s uce-8.nexus -m XXXXXXXX -b 10 -pre uce-8.nexus-iqtree-bootstrap

iqtree -s uce-9.nexus -m XXXXXXXX -b 10 -pre uce-9.nexus-iqtree-bootstrap

iqtree -s uce-11.nexus -m XXXXXXXX -b 10 -pre uce-11.nexus-iqtree-bootstrap

iqtree -s uce-13.nexus -m XXXXXXXX -b 10 -pre uce-13.nexus-iqtree-bootstrap

iqtree -s uce-14.nexus -m XXXXXXXX -b 10 -pre uce-14.nexus-iqtree-bootstrap

iqtree -s uce-15.nexus -m XXXXXXXX -b 10 -pre uce-15.nexus-iqtree-bootstrap
```

### Open R
### Plot tree again now with bootstrap support
```
library(ape)

tre2b = read.tree(file="uce-2.nexus-iqtree-bootstrap.treefile")
plot(tre2b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre2b = root(tre2b, node=151, resolve.root=TRUE)
plot(rtre2b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre2b$node.label)

tre3b = read.tree(file="uce-3.nexus-iqtree-bootstrap.treefile")
plot(tre3b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre3b = root(tre3b, node=151, resolve.root=TRUE)
plot(rtre3b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre3b$node.label)

tre6b = read.tree(file="uce-6.nexus-iqtree-bootstrap.treefile")
plot(tre6b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre6b = root(tre6b, node=151, resolve.root=TRUE)
plot(rtre6b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre6b$node.label)

tre7b = read.tree(file="uce-7.nexus-iqtree-bootstrap.treefile")
plot(tre7b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre7b = root(tre7b, node=151, resolve.root=TRUE)
plot(rtre7b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre7b$node.label)

tre8b = read.tree(file="uce-8.nexus-iqtree-bootstrap.treefile")
plot(tre8b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre8b = root(tre8b, node=151, resolve.root=TRUE)
plot(rtre8b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre8b$node.label)

tre9b = read.tree(file="uce-9.nexus-iqtree-bootstrap.treefile")
plot(tre9b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre9b = root(tre9b, node=151, resolve.root=TRUE)
plot(rtre9b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre9b$node.label)

tre11b = read.tree(file="uce-11.nexus-iqtree-bootstrap.treefile")
plot(tre11b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre11b = root(tre11b, node=151, resolve.root=TRUE)
plot(rtre11b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre11b$node.label)

tre13b = read.tree(file="uce-13.nexus-iqtree-bootstrap.treefile")
plot(tre13b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre13b = root(tre13b, node=151, resolve.root=TRUE)
plot(rtre13b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre13b$node.label)

tre14b = read.tree(file="uce-14.nexus-iqtree-bootstrap.treefile")
plot(tre14b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre14b = root(tre14b, node=151, resolve.root=TRUE)
plot(rtre14b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre14b$node.label)

tre15b = read.tree(file="uce-15.nexus-iqtree-bootstrap.treefile")
plot(tre15b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels()
rtre15b = root(tre15b, node=151, resolve.root=TRUE)
plot(rtre15b, type = "phylogram", cex = 0.3, no.margin = TRUE)
nodelabels(rtre15b$node.label)
```
