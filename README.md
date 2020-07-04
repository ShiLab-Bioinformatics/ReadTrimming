# Read-Trimming
Yang Liao and Wei Shi. "Read trimming is not required for mapping and quantification of RNA-seq reads at the gene level". Under Review.

## Data
All the data included in this study can be downloaded from the link:

https://latrobeuni-my.sharepoint.com/:f:/g/personal/yliao_ltu_edu_au/EgBbclblJpFPgcGZO3W-cKkB8XFgFtNbr63B8GfJitaAPQ?e=PvHBq0

## Reproduce the results
All the code can be found in the /src directory. To reproduce the results reported in the paper, you will need to download all the data and code and save them to the same directory. You also need to install the following programs:

1. [The R environment](https://www.r-project.org/)
2. [Rsubread (v1.32.3)](https://bioconductor.org/packages/3.8/bioc/src/contrib/Archive/Rsubread/)
3. [Trimmomatic (v0.39)](http://www.usadellab.org/cms/?page=trimmomatic)
4. [Trim Galore (v0.6.2)](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
5. [Python (v3.4.0 or newer)](https://www.python.org/downloads/)
6. [NumPy](https://numpy.org/) and [SciPy](https://www.scipy.org/) for your Python installation
7. [Bash (any modern version should be OK)](https://www.gnu.org/software/bash/)
8. [Samtools (v0.1.15+)](http://www.htslib.org/)

You need to make sure that the path to executables in these installed programs should be included in your PATH environment variable. The JAR file of Trimmomatic should be included in your current working directory. More details can be found in the shell scripts.

After setting up the environment and downloading all the data, you can simply run **Rerun-Analysis.bash** to reproduce all the results included in the paper. This shell script invokes other shell and R scripts, and Python programs.

It is recommended to run the analyses on an x86-64 Linux computer with at least 1TB of free disk space, 8 CPU cores and 64GB of memory.
