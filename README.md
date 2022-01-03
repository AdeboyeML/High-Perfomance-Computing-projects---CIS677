## High-Perfomance-Computing-projects-CIS677
This repo includes projects from high performance computing course - **CIS 677 in  GVSU**

### MPI for Genomics: Genome-wide Association Study (GWAS)
- Problem Description:
The aim of this project is to develop an **MPI-based message-passing application/program** that analyses raw microarray data and identifies top differentially expressed genes between the ***Renal and Control groups*** of patient samples.
----------------------------------------------------------------------------------------------------------------
#### - Random T-Test distribution and D scores of some of the top discriminant genes i.e., GENE377X and GENE2107X:

<p float="left">
  <img src="./genome_wide_t-test_project/ttest_imgs/ttest_dist.png" width="700" height="600">
  <img src="./genome_wide_t-test_project/ttest_imgs/ttest_dist2.png" width="700" height="600">
</p>

#### - Top 20 Discriminant Genes from both Sequential and MPI Parallel Programs:

<p float="left">
  <img src="./genome_wide_t-test_project/ttest_imgs/topgenes_seq.png" width="700" height="600">
  <img src="./genome_wide_t-test_project/ttest_imgs/topgenes_mpi.png" width="700" height="600">
</p>

#### - Open MPI Speed Up (~100X):

<p float="left">
  <img src="./genome_wide_t-test_project/ttest_imgs/ttest_speedup.png" width="700" height="600">
</p>
