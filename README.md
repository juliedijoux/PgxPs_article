# Requirements
All statistical analyses were performed using R version 4.2.0. The **`biostat_analysis_PgxPs.rmd`** file must be located in the directory (e.g., **`biostat_analysis`**), along with the **`sample_table.csv`** file, which contains the metadata, and the **`0_DADA2_output`** folder, which includes the output files from the DADA2 pipeline (see the *Bioinformatics-processing* branch of the repository) for the 16S and ITS libraries: **`seqtab.nochim_16S.rdata`** and **`seqtab.nochim_ITS.rdata`** for the abundance tables, and **`taxo_tab_16S.rdata`** and **`taxo_tab_ITS.rdata`** for the taxonomic affiliations. These files are necessary to run the entire R script.

# References
Csárdi, G., Nepusz, T., Traag, V., Horvát, S., Zanini, F., Noom, D., et al. (2024). igraph: Network analysis and visualization in R. https://cran.r-project.org/package=igraph

Dusa, A. (2024). venn: Draw Venn diagrams. https://cran.r-project.org/package=venn

Liu, C., Cui, Y., Li, X., & Yao, M. (2021). microeco: An R package for data mining in microbial community ecology. FEMS Microbiology Ecology, 97, fiaa255. https://doi.org/10.1093/femsec/fiaa255

Mangiafico, S. S. (2024). rcompanion: Functions to support extension education program evaluation. New Brunswick, New Jersey. https://cran.r-project.org/package=rcompanion/

Oksanen, J., Simpson, G. L., Blanchet, F. G., Kindt, R., Legendre, P., Minchin, P. R., et al. (2022). vegan: Community ecology package. https://cran.r-project.org/package=vegan

R Core Team. (2022). R: A language and environment for statistical computing. Vienna, Austria. https://www.r-project.org/

RStudio Team. (2020). RStudio: Integrated development environment for R. Boston, MA. http://www.rstudio.com/

Revelle, W. (2024). psych: Procedures for psychological, psychometric, and personality research. Evanston, Illinois. https://cran.r-project.org/package=psych

Sievert, C. (2020). Interactive web-based data visualization with R, plotly, and shiny. Chapman and Hall/CRC. https://plotly-r.com
