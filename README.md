# MG_Growth_Room

R code used for data manipulation, statistical analyses and figure generation for the manuscript "Intermittent water stress favors microbial generalists that better help wheat under drought"

If you intend to re-run the analyses, please download the entire folder structure (click on "Code", then "Download ZIP" )

The raw data needed to run the scripts can be found on Zenodo: [https://doi.org/10.5281/zenodo.10140593](https://can01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fdoi.org%2F10.5281%2Fzenodo.10140593&data=05%7C01%7CEtienne.Yergeau%40inrs.ca%7C7af38c8c2dfc41ae71d108dbe652027e%7C80c9267a177f49439958b599ca0e54c2%7C0%7C0%7C638357012125865438%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&sdata=IxkiOpbsrPO1lgEZZYczhdIIwmQaMqzRk7YO0fzC9Q8%3D&reserved=0). It should be copied to the /data/raw folder.

The scripts should be run in order (01 to 08).

**01-LoadPackages.R**: Download, install and load necessary packages

**02-LoadRawDataNormalise.R**: Load raw data, normalize (if needed), order and save as intermediate files (.RDS) in the /data/intermediate folder.

**03-PlantBiomass.R**: ANOVA and figure generation for plant root and shoot biomass and leaf water content

**04-PCoAPermanova.R**: Generate PCoA ordination based on Bray-Curtis dissimilarity of gene relative abundance and test effect of treatments on the gene relative abundance table.

**05-Traits.R**: Look for genes related to resistance to drought, run ANOVA tests for their summed relative abundance and permanova for their gene table, and generate figures.

**06-Taxonomy.R**: Generate stack bar chart figure for the phylum relative abundance.

**07-MAGs.R**: ANOVA test for the effect of the treatments on the relative abundance of the MAGs. For the significant MAGs, look for occurence, prevalence and diversity of the functional genes identified in 05-Traits.R. Generate figures.

**08-Figures.R**: Create multi-pannel figures for publication.
