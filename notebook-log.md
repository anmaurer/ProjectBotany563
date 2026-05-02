# Notebook for Botany 563
(Written by Ari Maurer in Spring 2026) 

## Dataset for class final project
Dispersal sweepstakes: Biotic interchange propelled air-breathing fishes across the globe by Harrington et al., 2024

Full text PDF: https://deepblue.lib.umich.edu/items/e98808ae-5575-453f-99dc-f92c9ec2bc87 

Dataset available from: 
https://datadryad.org/dataset/doi:10.5061/dryad.59zw3r2c0 

Dataset summary:
The authors obtained UCE DNA sequence data through hybrid enrichment of genomic libraries. DNA samples were prepared using Kapa HyperPrep kits and UCE probes were designed to target ~1300 ultraconserved elements from the clade Acanthomorpha. Illumina HiSeq4000 was used to perform sequencing. Data was processed and aligned by the authors using the software Phyluce before download. 

## All Software Used 
| Software | Description | Strengths | Weaknesses | Assumptions | User Choices |
|---|---|---|---|---|---|
| Muscle | Program for creating multiple alignments of protein sequences. Includes fast distance estimation using k mer counting, progressive alignment, and refinement using tree-dependent restricted partitioning. | Fast and accurate for large datasets | MUSCLE's accuracy drops when aligning highly divergent sequences | MUSCLE assumes using k-mer counting (a fast distance measurement) and then iterative refinment can produce results similiar to computationally expensive methods but faster | User specifies what kind of sequencing data they are inputting |
|MAFFT|Multiple sequence alignment program designed for amino acid or nucleotide sequences|Balances speed and accuracy, flexible in handling diverse dataset sizes|For very long DNA sequences, can be time-consuming. Can over-align. |Assumes order of residues is conserved.|Users can select alignment strategy and adjust parameters for gap penalties.|
|RAxML|RAxML is an open-source software that estimates maximum likelihood phylogenies|RAxML typically finds higher-scoring trees than IQTree|IQTree typically finds trees with better stability|RAxML uses randomized stepwise addition parsimony trees. RAxML does model support for you, and supports 22 classical GTR models. RAxML can support DNA and protein data, as well as binary, multi-state morphological and RNA secondary structure data. RAxML can correct for ascertainment bias.|All model parameters can either be optimized or fixed to user-specified values. User needs to check what model RAxML chooses. |
|IQtree|Software for phylogenetic inference using maximum likelihood|Utilizes fast model selection via ModelFinder. IQtree faster than RAxML in optimization of model parameters|Sensitive to model mis-specification, can be slow|Two starting trees computed using maximum parsimony and neighbor joining which are then optimized by hill climbing nearest neighbor interchange moves. Uses stochastic NNI to escape local optima.|User choices TBD|
|MrBayes|The program performs Bayesian inference of phylogeny using a variant of Markov chain monte carlo|Bayesian inference allows for the incorporation of prior knowledge, no longer need homogenous data, and the user can account for complexity. Bayesian phylogenetic inference promises users direct probability statements about trees/clades, integration over uncertainty, and use of prior information. MCMC is used to approximate the posterior probabilities of trees. MrBayes allows for mixed models and higher computation efficiency. Can unlink topology and branch lengths. Allows for analysis of more than 350 taxa.|MCMC has a high computational cost, and exact computation of the posterior is impossible. Slow.|Assumes all parameters are random variables. |User must choose substitution model, the prior, and the details of the MCMC analysis (number of generations, how convergence was assessed, effective sample size etc.)|
|Astral|Coalescent method for reconstructing species tree after inferring a set of gene trees.|Statistically consistent under the multi-species coalescent model, is scalable, and has shown high accuracy in simulated and empirical studies. |Astral can be sensitive to gene tree estimation error. Astral is statistically inconsistent under models of gene evolution that include gene flow.|Assumes Multi-species coalescent model, maxmizes the number of shared quartet trees with the input gene trees which must be unrooted.|User must provide Newick format gene trees|

## Data alignment with MUSCLE
*I did not do this step since the data was already aligned before download.*
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
--> I chose the first 3 UCE files (2, 3, 6) from the data set to make trees for 

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
```

### Get the NJ tree
```
tre2 <- nj(Uce2D)
tre3 <- nj(Uce3D)
tre6 <- nj(Uce6D)
```

### Ladderize
Before plotting, we can use the ladderize function which reorganizes the internal structure of the tree to get the ladderized effect when plotted
```
tre2L <- ladderize(tre2)
tre3L <- ladderize(tre3)
tre6L <- ladderize(tre6)
```
### Calculate bootstrap support values
```
Uce2D <- as.matrix(Uce2D)
bs2 <- boot.phylo(tre2L, Uce2D, function(x) nj(dist.dna(x)), B = 100)

Uce3D <- as.matrix(Uce3D)
bs <- boot.phylo(tre3L, Uce3D, function(x) nj(dist.dna(x)), B = 100)

Uce6D <- as.matrix(Uce6D)
bs <- boot.phylo(tre6L, Uce6D, function(x) nj(dist.dna(x)), B = 100)
```
### Plot the tree (simple plot to check)
```
plot(tre2L, cex=.6)
title("UCE-2")

plot(tre3L, cex=.6)
title("UCE-3")

plot(tre6L, cex=.6)
title("UCE-6")
```

### Plot the tree more legibly (vertical plot)
```
pdf("UCE2_tree_vertical.pdf", width = 10, height = 20)
plot(tre2L, cex = 0.5)
nodelabels(bs2, cex=0.6)
title("UCE-2 Distance Tree")
dev.off() 

pdf("UCE3_tree_vertical.pdf", width = 10, height = 20)
plot(tre3L, cex = 0.5)
nodelabels(bs3, cex=0.6)
title("UCE-3 Distance Tree")
dev.off() 

pdf("UCE6_tree_vertical.pdf", width = 10, height = 20)
plot(tre6L, cex = 0.5)
nodelabels(bs6, cex=0.6)
title("UCE-6 Distance Tree")
dev.off() 
```

## Parsimony Tree Calculations
--> I chose the first 3 UCE files (2, 3, 6) from the data set to make trees for 

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
```

### Search for the tree with maximum parsimony 
```
tre2.pars <- optim.parsimony(tre2.ini, Uce2P)
tre3.pars <- optim.parsimony(tre3.ini, Uce3P)
tre6.pars <- optim.parsimony(tre6.ini, Uce6P)
```
### Bootstrap support calculations
```
bs2 <- bootstrap.phyDat(Uce2P, FUN = function(x) pratchet(x),bs = 100)
bs3 <- bootstrap.phyDat(Uce3P, FUN = function(x) pratchet(x),bs = 100)
bs6 <- bootstrap.phyDat(Uce6P, FUN = function(x) pratchet(x),bs = 100)
```

### Plot tree (basic)
```
plot(tre2.pars, cex=0.6)
plot(tre3.pars, cex=0.6)
plot(tre6.pars, cex=0.6)
```

### Plot tree more legibly (Vertical)
```
pdf("UCE2_Ptree_verticalBS.pdf", width = 10, height = 20)
plot(tre2.pars, cex = 0.5)
nodelabels(prop.clades(tre2.pars, bs2), cex=0.6)
title("UCE-2 Parsimony Tree")
dev.off()

pdf("UCE3_Ptree_vertical.pdf", width = 10, height = 20)
plot(tre3.pars, cex = 0.5)
nodelabels(prop.clades(tre2.pars, bs3), cex=0.6)
title("UCE-3 Parsimony Tree")
dev.off()

pdf("UCE6_Ptree_vertical.pdf", width = 10, height = 20)
plot(tre6.pars, cex = 0.5)
nodelabels(prop.clades(tre2.pars, bs6), cex=0.6)
title("UCE-6 Parsimony Tree")
dev.off()
```

## Maximum Likelihood Tree Calculations
--> I chose the first 10 UCE files (2, 3, 6) from the data set to make trees for 

## RAxML 
###  Download RAxML
https://github.com/amkozlov/raxml-ng 

### Convert Data
RAxML cannot read Nexus format. Nexus files can be converted to PHYLIP format, which RAxML can read using R:

```
library (ape)
alignment2 <- read.nexus.data("uce-2.nexus")
write.dna(alignment2, file = "uce-2.phy", format = "sequential")

alignment3 <- read.nexus.data("uce-3.nexus")
write.dna(alignment3, file = "uce-3.phy", format = "sequential")

alignment6 <- read.nexus.data("uce-6.nexus")
write.dna(alignment6, file = "uce-6.phy", format = "sequential")
```

### Run RAxML to calculate the ML tree

```
raxml-ng --msa uce-2.phy --model LG+G8+F

raxml-ng --msa uce-3.phy --model LG+G8+F

raxml-ng --msa uce-6.phy --model LG+G8+F
```

### Run both ML and Bootstrap (!! Takes forever !!) 

```
raxml-ng --all --msa uce-2.phy --model LG+G8+F --bs-trees 100 --prefix uce-2-raxml-boostrap

raxml-ng --all --msa uce-3.phy --model LG+G8+F --bs-trees 100 --prefix uce-3-raxml-boostrap

raxml-ng --all --msa uce-6.phy --model LG+G8+F --bs-trees 100 --prefix uce-6-raxml-boostrap
```

### Root and plot tree (Open R)

```
library(ape)

tre2 = read.tree(file="uce-2-raxml-boostrap.raxml.support")
rtre2 = root(tre2, node=151, resolve.root=TRUE)
plot(rtre2)
nodelabels(rtre2$node.label)

tre3 = read.tree(file="uce-3-raxml-boostrap.raxml.support")
rtre3 = root(tre3, node=151, resolve.root=TRUE)
plot(rtre3)
nodelabels(rtre3$node.label)

tre6 = read.tree(file="uce-6-raxml-boostrap.raxml.support")
rtre6 = root(tre6, node=151, resolve.root=TRUE)
plot(rtre6)
nodelabels(rtre6$node.label)
```

### Plot tree more legibly (vertical)
```
tre2 = read.tree(file="uce-2-raxml-boostrap.raxml.support")
rtre2 = root(tre2, node=151, resolve.root=TRUE)
pdf("UCE2_MLtree_verticalX.pdf", width = 10, height = 20)
plot(rtre2, cex = 0.5)
nodelabels(rtre2$node.label)
title("UCE-2 Maximum Likelihood Tree")
dev.off()

tre3 = read.tree(file="uce-3-raxml-boostrap.raxml.support")
rtre3 = root(tre3, node=151, resolve.root=TRUE)
pdf("UCE3_MLtree_vertical.pdf", width = 10, height = 20)
plot(tre3, cex = 0.5)
title("UCE-3 Maximum Likelihood Tree")
dev.off()

tre6 = read.tree(file="uce-6-raxml-boostrap.raxml.support")
rtre6 = root(tre6, node=151, resolve.root=TRUE)
pdf("UCE6_MLtree_vertical.pdf", width = 10, height = 20)
plot(rtre6, cex = 0.5)
title("UCE-6 Maximum Likelihood Tree")
dev.off()
```

## IQtree
*Another way to conduct a maximum likelihood analysis. This is not the analysis I used in my paper, however, I did include the models that IQtree determined in downstream Bayesian analysis* 
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

tre3 = read.tree(file="uce-3.nexus.treefile")
plot(tre3)

tre6 = read.tree(file="uce-6.nexus.treefile")
plot(tre6)
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
```

### Close R
### Quantify support for the estimated tree- check what best fit model is for each on previous step
```
iqtree -s uce-2.nexus -m TPM3u+I+R4 -b 10 -pre uce-2.nexus-iqtree-bootstrap

iqtree -s uce-3.nexus -m XXXXXXXX T -b 10 -pre uce-3.nexus-iqtree-bootstrap

iqtree -s uce-6.nexus -m XXXXXXXX -b 10 -pre uce-6.nexus-iqtree-bootstrap
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
```
## Bayesian Inference Phylogenetic Method (MrBayes) 
### Download MrBayes

```
conda install -c bioconda mrbayes
```
### Create a mrbayes block in a separate text file called mbblock.txt 

```
begin mrbayes;
    set autoclose=yes;

    prset brlenspr=unconstrained:exp(10.0);
    prset shapepr=exp(1.0);
    prset tratiopr=beta(1.0,1.0);
    prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);

    lset nst=2 rates=gamma ngammacat=4;

    mcmcp ngen=1000000 samplefreq=100 printfreq=1000 nruns=2 nchains=4 savebrlens=yes;

    outgroup Channa_bleheri_48;

    mcmc;

    sumt burnin=2500 contype=allcompat;

end;
```

### Append the MrBayes block to the end of each nexus file 

```
cat uce-2.nexus mmmblock.txt > uce-2-mb.nex

cat uce-3.nexus mbblock.txt > uce-3-mb.nex

cat uce-6.nexus mbblock.txt > uce-6-mb.nex
```

### Run MrBayes

```
mb uce-2-mb.nex

mb uce-3-mb.nex

mb uce-6-mb.nex
```

### Uce Tracer to access if the chain converged and had good mixing 
Go to File --> Import Trace file and select uce-X-mb.nex.p

Look at the trace plot:

Good behavior:
The trace looks like a fuzzy horizontal band
No obvious trends up or down after burn-in
Bad behavior:
Strong trends
Long flat stretches
Sudden jumps followed by no mixing
ESS Rules of thumb:

ESS > 200 → acceptable
ESS > 500 → good
ESS < 100 → problematic

### Basic Plots in R (Open R)

```
library(ape)

tre2 = read.nexus(file="uce-2-mb.nex.con.tre")
plot(tre2)

tre3 = read.nexus(file="uce-3-mb.nex.con.tre")
plot(tre3)

tre6 = read.nexus(file="uce-6-mb.nex.con.tre")
plot(tre6)
```
### Plot trees more legibly (vertical)
```
tre2 = read.nexus(file="uce-2-mb.nex.con.tre")
pdf("UCE2_Btree_verticalBS.pdf", width = 10, height = 20)
plot(tre2, cex = 0.5)
title("UCE-2 Bayesian Tree")
dev.off()

tre3 = read.nexus(file="uce-3-mb.nex.con.tre")
pdf("UCE3_Btree_vertical.pdf", width = 10, height = 20)
plot(tre3, cex = 0.5)
title("UCE-3 Bayesian Tree")
dev.off()

tre6 = read.nexus(file="uce-6-mb.nex.con.tre")
pdf("UCE6_Btree_vertical.pdf", width = 10, height = 20)
plot(tre6, cex = 0.5)
title("UCE-6 Bayesian Tree")
dev.off()
```

## Coalescent Methods (Astral)
### Astral only takes Newick files, so concatenate previous treefiles from IQTree in Anabantoid WD
*I used the UCE2, UCE3, and UCE6 files*
```
cat *.treefile > gene_trees.tre
```

### Move file to ASTRAL WD
### Run ASTRAL
```
java -jar astral.5.7.8.jar -i gene_trees.tre -o out.tre
```

### plot tree in R
```
library(ape)
tre = read.tree(file="out.tre")
plot(tre)
```

### Plot tree more legibly (vertical)
```
tre = read.tree(file="out.tre")
pdf("UCE_Atree_vertical.pdf", width = 10, height = 20)
plot(tre, cex = 0.5)
nodelabels(tre$node.label, cex=0.6)
title("UCE-2, UCE-3, UCE-6 Astral Tree")
dev.off()
```


