# oisst_warming_trends
Accompanying code repository with scripts and programs for our paper:

> Reduced nearshore warming associated with Eastern Boundary Upwelling Systems.
> Rui Seabra, Rubén Varela, António M. Santos, Moncho Gomez-Gesteira, Claudia Meneghesso, David S. Wethey, Fernando P. Lima.
> Frontiers in Marine Science, 2019.

The repository is structured as follows:
- `code`: C++ code to be used with our genesis library. The directory contains all prototypes and helper programs that were used to evaluate the data and make the plots reported in the papers.
- `data`: Contains information and some processed files about the empirical datasets that were used as a basis for the analyses and evaluations.
- `results`: A large collection of scripts, data, results and figures that were used in the papers. The directory contains all analyses that we used for the papers.
- `software`: Information about the extrnal software that was used.

The code provided allows the reproduction of the analysis from scratch, and up to the production of all figures and tables with the necessary stats.

**IMPORTANT REMARK:** The code was designed to run on a 160 thread local server with 1.5 TB of RAM. It makes frequent use of parallel computations and deals with several "dificult" parts of the computation by  using brute-force computation. Please bear in mind that if this code is to be run in a regular machine more efficient memory management and a lot of patience may be required...