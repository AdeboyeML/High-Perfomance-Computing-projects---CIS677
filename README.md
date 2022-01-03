## High-Perfomance-Computing-projects-CIS677
This repo includes projects from high performance computing course - **CIS 677 in  GVSU**

### MPI for Genomics: Genome-wide Association Study (GWAS)
- Problem Description:
The aim of this project is to develop an **MPI-based message-passing application/program** that analyses raw microarray data and identifies top differentially expressed genes between the ***Renal and Control groups*** of patient samples.
----------------------------------------------------------------------------------------------------------------
#### - Random T-Test distribution and D scores of some of the top discriminant genes i.e., GENE377X and GENE2107X:

<p float="left">
  <img src="./genome_wide_t-test_project/ttest_imgs/ttest_dist.png" width="500" height="400">
  <img src="./genome_wide_t-test_project/ttest_imgs/ttest_dist2.png" width="500" height="400">
</p>

#### - Top 20 Discriminant Genes from both Sequential and MPI Parallel Programs:

<p float="left">
  <img src="./genome_wide_t-test_project/ttest_imgs/topgenes_seq.png" width="500" height="400">
  <img src="./genome_wide_t-test_project/ttest_imgs/topgenes_mpi.png" width="500" height="400">
</p>

#### - Open MPI Speed Up (~100X):

<p float="left">
  <img src="./genome_wide_t-test_project/ttest_imgs/ttest_speedup.png" width="500" height="400">
  <img src="./genome_wide_t-test_project/ttest_imgs/perf_analysis.png" width="500" height="400">
</p>


-----------------------------------------------------------------------------------------------------------------------------------


### GPU-Accelerated computation of Voronoi Diagram

- The main idea is to write a CUDA-based massively multi-threaded application that when given a set of seeds and a 2-dimensional canvas, computes a discrete approximation of the corresponding Voronoi diagram.
----------------------------------------------------------------------------------------------------------------

#### - Voronoi images for different seed sizes (1024 & 10,000)

<p float="left">
  <img src="./voronoi_project/imgs_vor/image_1024_1024.png" width="400" height="400">
  <img src="./voronoi_project/imgs_vor/image_10000_2048.png" width="400" height="400">
</p>

<p float="left">
  <img src="./voronoi_project/imgs_vor/img_1024.gif" width="400" height="400">
</p>

#### - Performance Analysis
<p float="left">
  <img src="./voronoi_project/imgs_vor/perf_analysis.png" width="550" height="400">
  <img src="./voronoi_project/imgs_vor/perf_analysis2.png" width="550" height="400">
</p>
