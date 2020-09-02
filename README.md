# Read-Trimming
Yang Liao and Wei Shi. "Read trimming is not required for mapping and quantification of RNA-seq reads at the gene level". Accepted for publication by [NAR Genomics and Bioinformatics](https://academic.oup.com/nargab).

## Data
All the data included in this study can be downloaded from the link:

https://latrobeuni-my.sharepoint.com/:f:/g/personal/yliao_ltu_edu_au/EgBbclblJpFPgcGZO3W-cKkB8XFgFtNbr63B8GfJitaAPQ?e=PvHBq0

## Software
The following programs should be installed to reproduce the results reported in the paper:

1. [R (v3.6.1)](https://www.r-project.org/)
2. [Rsubread (v2.0.0)](https://bioconductor.org/packages/3.10/bioc/html/Rsubread.html)
3. [Trimmomatic (v0.39)](http://www.usadellab.org/cms/?page=trimmomatic)
4. [Trim Galore (v0.6.2)](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
5. [Python (v3.4.0 or newer)](https://www.python.org/downloads/)
6. [NumPy](https://numpy.org/) and [SciPy](https://www.scipy.org/) for your Python installation
7. [Bash (any modern version should be OK)](https://www.gnu.org/software/bash/)
8. [Samtools (v0.1.15 or newer)](http://www.htslib.org/)
9. [STAR (v2.7.3a)](https://github.com/alexdobin/STAR/releases/tag/2.7.3a)

You need to make sure that the path to executables in these installed programs should be included in your PATH environment variable. The JAR file of Trimmomatic should be included in your current working directory.

## Analysis code
All the code can be found in the /src directory. To reproduce the results reported in the paper, you will need to download all the data and code and save them to the same directory. 

After setting up the environment and downloading all the data, you can simply run **Rerun-Analysis.bash** to reproduce all the results included in the paper. This shell script will invoke other shell scripts.

It is recommended to run the analyses on an x86-64 Linux computer with at least 1TB of free disk space, 8 CPU cores and 64GB of memory.
