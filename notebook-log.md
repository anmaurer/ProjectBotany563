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
muscle -anabantoid.fasta -output aligned_anabantoid.fasta
```

** Data not given in unaligned FASTA format, data given in already aligned Nexus files, so this step not run. Data given as UCE_loci_75p_complete.gzip   > gzip directory of nexus-formatted alignments for individual UCE loci. All alignments in this directory contain at least 75% of taxa in our total dataset. **
