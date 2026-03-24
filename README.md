# BINF-6110-Assignment-3

## Introduction 

The human gut microbiome is a complex community of microorganisms that plays a central role in host metabolism, immune regulation, and disease susceptibility. Diet is one of the most influential modulators of gut microbial composition, with long-term dietary patterns capable of driving substantial shifts in community structure and function [1]. Plant-based diets are characterised by high fibre intake, which selectively enriches taxa capable of complex carbohydrate fermentation and short-chain fatty acid production, metabolites with well-documented anti-inflammatory properties. In contrast, omnivore diets, which include red meat and animal products, have been associated with enrichment of taxa involved in protein fermentation and bile acid metabolism. Notably, the genus *Prevotella* has been repeatedly linked to fibre-rich, plant-based diets, while Bacteroides predominates in individuals consuming high-fat, protein-rich diets [2]. However, the health implications of *Prevotella* abundance remain contested, with associations reported for both beneficial metabolic outcomes and inflammatory conditions, a discrepancy that may be explained by strain-level genomic diversity within the species [2].
Shotgun metagenomics provides a culture-independent, whole-genome sequencing approach that enables taxonomic and functional characterisation of microbial communities at greater resolution than 16S rRNA amplicon sequencing. Unlike amplicon-based methods, shotgun sequencing captures the full genetic complement of a sample, allowing species-level and strain-level classification without the primer bias inherent to targeted approaches. This study re-analyses publicly available shotgun metagenomic data from Italian omnivore and vegan donors originally deposited by De Filippis et al. [2], retrieved from the NCBI. Raw reads were processed using Kraken2 for taxonomic classification and Bracken for abundance re-estimation at the species level. Kraken2 uses an exact k-mer matching algorithm against a reference database to assign reads to taxonomic nodes with high speed and accuracy [3], while Bracken applies a Bayesian re-estimation step to redistribute higher-level classifications into robust species-level abundance estimates [4].
Downstream community analysis was performed in R using the phyloseq framework for data management and visualisation, vegan for ecological distance metrics and multivariate statistics, and ANCOMBC2 for differential abundance testing [5]. Beta diversity was assessed using Bray-Curtis dissimilarity with both principal coordinates analysis (PCoA) and non-metric multidimensional scaling (NMDS) ordination, and group-level differences were tested using PERMANOVA. Differential abundance between omnivore and vegan groups was evaluated using ANCOMBC2, which accounts for unequal sampling fractions and structural zeros, a common challenge in compositional microbiome data [5]. Together, this pipeline provides a systematic characterisation of how diet shapes gut microbial community structure and identifies taxa specifically enriched or depleted in vegan relative to omnivore samples.

## Methods

### Data Acquisition

Shotgun metagenomic sequencing data were obtained from De Filippis et al. [2], originally deposited in the NCBI Sequence Read Archive under BioProject PRJNA421881. Six samples were retrieved (SRR8146935, SRR8146936, SRR8146938, SRR8146944, SRR8146951, SRR8146952), corresponding to three omnivore and three vegan donors. Raw SRA files were downloaded using prefetch from the SRA Toolkit on the Narval high-performance computing cluster (Digital Research Alliance of Canada). FASTQ files were extracted from SRA format using fasterq-dump with paired-end splitting (--split-files, 32 threads), executed as a SLURM array job across all six samples. Output files were compressed with gzip and integrity-verified prior to downstream processing.

### Taxonomic Classification

Taxonomic classification of paired-end reads was performed using Kraken2 [3] against the Standard 8GB database kraken database, with a confidence threshold of 0.15 (--confidence 0.15) and 32 threads per sample. Classification was run as a SLURM array job. Kraken2 outputs per-read classification files and per-sample taxonomic report files. Species-level abundance re-estimation was subsequently performed using Bracken [4], which applies a Bayesian approach to redistribute reads assigned to higher taxonomic nodes back to species level, producing corrected abundance reports per sample. The six Bracken species-level report files were combined into a single BIOM-format abundance table using kraken-biom (v1.0.1), executed via Apptainer on Narval. The resulting kraken_table.biom file served as input for all downstream R analyses.

### Downstream Analysis in R

All downstream analyses were performed in R using the phyloseq framework (v1.46) for data import, management, and visualisation [6]. The BIOM file was imported using the biomformat package and converted to a phyloseq object. Sample metadata, including BioSample accession, organism, and diet group assignment, were retrieved programmatically from the NCBI SRA and BioSample databases using the rentrez package and merged with the phyloseq object.
Sequencing adequacy was assessed by rarefaction curve analysis using the rarecurve function from the vegan package [7]. Taxonomic composition was visualised at the genus, family, and species levels by aggregating taxa with tax_glom(), transforming to relative abundance with transform_sample_counts(), and plotting stacked bar charts using ggplot2 [8].
Alpha diversity was estimated using observed species richness, Chao1, Shannon, and Simpson indices via the estimate_richness() function in phyloseq. Differences in Shannon diversity between diet groups were assessed using a Welch two-sample t-test. Beta diversity was assessed using Bray-Curtis dissimilarity. Principal coordinates analysis (PCoA) and non-metric multidimensional scaling (NMDS) ordinations were computed using the ordinate() function in phyloseq. The statistical significance of community composition differences between diet groups was tested using PERMANOVA from the vegan package [7]. Differential abundance analysis between omnivore and vegan groups was performed using ANCOMBC2 from the ANCOMBC package [5], with Diet as the fixed effect, Holm multiple testing correction, structural zero detection enabled (struc_zero = TRUE), and a minimum library size of 1,000 reads (lib_cut = 1000). Results were visualised as a ranked differential abundance plot showing log fold change with standard error bars for all 66 tested taxa.

## Results

### Sequencing Depth and Quality Control

Rarefaction curves were generated to assess whether sequencing depth was sufficient to capture the microbial diversity present in each sample (Figure 1). Sequencing depth ranged from 348,708 reads to 676,187 reads. All six rarefaction curves reached a clear plateau, indicating that sequencing depth was adequate to capture the observed species richness in each sample and that additional sequencing would yield minimal new taxa.

<img width="924" height="490" alt="rare" src="https://github.com/user-attachments/assets/06a1541e-3649-46c7-8d32-e169085d4657" />

Figure 1: Rarefaction curves showing the accumulation of observed species as a function of sequencing depth for each sample. All samples approach a clear asymptote, indicating sufficient sequencing depth to capture the majority of microbial diversity present. This suggests that downstream diversity comparisons are unlikely to be confounded by undersampling.

### Taxonomic Composition

Taxonomic composition was characterised at the genus, family, and species levels (Figures 2, 3, 4). At the genus level, *Segatella* was the dominant taxon across vegan samples, while omnivore samples displayed greater variability, with *Alistipes*, *Bacteroides*, and *Faecalibacterium* more prominently represented. At the family level, *Prevotellaceae* was visually dominant in vegan donors, while omnivore samples showed higher relative abundance of *Bacteroidaceae* and *Rikenellaceae*. At the species level, the top 15 most abundant species included *Segatella copri*, *Faecalibacterium prausnitzii*, and *Alistipes onderdonkii*, with notable compositional differences between diet groups.

<img width="924" height="559" alt="genus" src="https://github.com/user-attachments/assets/b6c1ab60-649e-436b-8ac3-9540e94a377a" />

Figure 2: Relative abundance of dominant microbial genera across samples, aggregated at the genus level and stratified by diet. Patterns of genus-level composition vary substantially between individuals, with no clear separation between omnivore and vegan groups, further supporting the absence of strong diet-driven differences in community structure.

<img width="924" height="490" alt="family comp" src="https://github.com/user-attachments/assets/6895cd0f-5688-4745-a6b7-18ca440ceac9" />

Figure 3: Family-Level Taxonomic Composition by Diet Group. Relative abundance of dominant microbial families across samples, aggregated at the family level and stratified by diet. While differences in family-level composition are observed between individual samples, no consistent diet-associated clustering is evident.

<img width="924" height="490" alt="rel_taxa" src="https://github.com/user-attachments/assets/fb347005-8cc6-4bdc-b918-4cf6ceac27df" />

Figure 4: Relative abundance of the 15 most abundant microbial species across samples, stratified by diet group. Taxa were aggregated at the species level following Bracken re-estimation and normalised to relative abundance. Considerable inter-individual variability is evident within both diet groups, with no consistent pattern of species dominance distinguishing omnivore and vegan donors.

### Alpha Diversity

Alpha diversity was estimated using Observed species richness, Chao1, Shannon, and Simpson indices (Figure 5). Omnivore samples showed higher mean Shannon diversity (mean = 2.81) compared to vegan samples (mean = 2.32), as well as higher mean Simpson index (mean = 0.845 vs 0.717). However, vegan samples showed marginally higher mean observed species richness (mean = 112.7 vs 100.7). A Welch two-sample t-test comparing Shannon diversity between diet groups was not statistically significant (t = 0.849, df = 3.61, p = 0.449), which is expected given the small sample size of three donors per group.

<img width="924" height="490" alt="Alpha" src="https://github.com/user-attachments/assets/e6cf39dd-9e29-4917-95e9-8d233019ec54" />

Figure 5: Alpha diversity metrics comparing gut microbial communities between omnivore and vegan donors. While vegan samples exhibited marginally higher observed richness and Chao1 estimates, omnivore samples showed consistently higher Shannon and Simpson diversity, indicating greater community evenness. Substantial overlap between groups across all metrics suggests no statistically significant differences in alpha diversity.

### Beta Diversity

Beta diversity was assessed using Bray-Curtis dissimilarity visualised through PCoA and NMDS ordination (Figures 6, 7). The NMDS solution converged with a near-zero stress value (stress = 1.78 x 10^-5), indicating low-dimensional representation of the dissimilarity structure. Visual inspection of both ordinations suggested some separation between omnivore and vegan samples, though with considerable spread within groups. PERMANOVA did not detect a statistically significant difference in community composition between diet groups (R^2 = 0.181, F = 0.887, p = 0.500), indicating that diet explained approximately 18% of the variance in Bray-Curtis dissimilarity but that this effect was not significant given the sample size.

<img width="924" height="362" alt="PcoA" src="https://github.com/user-attachments/assets/966e7f19-ec87-42ce-a4b0-f9fb02ef692c" />

Figure 6: Principal coordinates analysis (PCoA) of Bray–Curtis dissimilarity illustrating differences in microbial community composition between omnivore and vegan samples. Each point represents an individual sample, coloured by diet group, with dashed lines indicating convex hulls encompassing each group. While some separation along the second axis is observed, substantial overlap between groups is evident, consistent with the non-significant PERMANOVA result and indicating limited diet-associated structuring of community composition.

<img width="924" height="362" alt="NMDS" src="https://github.com/user-attachments/assets/59cc7cc5-85a2-47c5-aa80-b8bbbf4595e2" />

Figure 7. Non-metric multidimensional scaling (NMDS) ordination based on Bray–Curtis dissimilarity showing the relative similarity of microbial communities across samples. Points represent individual samples coloured by diet group. Despite minor clustering tendencies, considerable overlap between omnivore and vegan samples is observed, supporting the conclusion that diet does not drive a statistically significant shift in overall community composition in this dataset.

### Differential Abundance

Differential abundance analysis was performed using ANCOMBC2 with Holm multiple testing correction across 66 taxa (Figure 8). No taxa reached statistical significance after correction (q < 0.05). The species with the largest positive log fold change in vegan donors were *Blautia wexlerae* (LFC = 3.64), *Anaerostipes hadrus* (LFC = 3.17), and *Ruminococcus bicirculans* (LFC = 3.09), all with q-values of 1 following correction. The species most enriched in omnivore donors included *Alistipes onderdonkii* (LFC = −3.96), *Butyricimonas virosa* (LFC = −2.69), and *Alistipes ihumii* (LFC = −2.52). The absence of statistically significant findings is consistent with the limited statistical power provided by a sample size of three per group.

<img width="924" height="490" alt="diff_abd" src="https://github.com/user-attachments/assets/422d0b59-c2e3-4ee2-ab46-22f6b22637e7" />

Figure 8: Differential abundance analysis of microbial taxa between vegan and omnivore groups using ANCOMBC2. Log fold change values represent enrichment in vegan relative to omnivore samples, with error bars indicating standard error. The vertical dashed line denotes no difference. Although several taxa display large effect sizes, no taxa were significantly different following Holm-adjusted multiple testing correction (q < 0.05), reflecting high variability and limited statistical power.

## Discussion

This study re-analysed shotgun metagenomic data from De Filippis et al. [2] to characterise gut microbial community differences between omnivore and vegan donors. While no statistically significant differences were detected in alpha diversity, beta diversity, or differential abundance after multiple testing correction, the directional trends observed across analyses are consistent with the broader literature on diet-microbiome interactions.
The dominance of *Segatella copri* in vegan samples is consistent with the results of De Filippis et al. [2], who demonstrated that fibre-rich plant-based diets select for distinct *P. copri* strains with enhanced capacity for complex carbohydrate degradation. The high relative abundance of *Prevotellaceae* in vegan donors observed at the family level further supports this pattern. *Segatella copri* is canonically associated with agrarian and plant-based diets, where high dietary fibre availability promotes its proliferation through carbohydrate fermentation pathways [2].
The enrichment of *Blautia wexlerae* and *Anaerostipes hadrus* in vegan donors, while not reaching significance, is also biologically plausible. Both are butyrate and acetate-producing members of the Lachnospiraceae family, and their abundance has been positively associated with plant-based dietary patterns and anti-inflammatory metabolic profiles [1]. Conversely, the enrichment of *Alistipes onderdonkii* and *Alistipes ihumii* in omnivore donors is consistent with the known association of *Alistipes* species with protein fermentation and animal-based diets, as members of this genus have been linked to the metabolism of amino acids derived from dietary protein [2].
The omnivore group showed higher Shannon and Simpson diversity compared to vegans, contrary to the common assumption that plant-based diets universally increase microbial diversity. This may reflect the compositional dominance of *Segatella copri* in vegan samples, which can suppress evenness even when richness is maintained. However, given the sample size of three donors per group, none of these findings can be interpreted as statistically conclusive. The PERMANOVA result (R^2 = 0.181, p = 0.500) indicates a potentially meaningful effect size that was underpowered to detect significance. Future studies with larger cohorts and matched covariates such as age, BMI, and antibiotic history would be necessary to confirm these directional trends.

## References

[1] Fackelmann, G., Manghi, P., Carlino, N., Heidrich, V., Piccinno, G., Ricci, L., Piperni, E., Arrè, A., Bakker, E., Creedon, A. C., Francis, L., Capdevila Pujol, J., Davies, R., Wolf, J., Bermingham, K. M., Berry, S. E., Spector, T. D., Asnicar, F., & Segata, N. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. Nature Microbiology, 10(1), 41–52. https://doi.org/10.1038/s41564-024-01870-z

[2] De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct genetic and functional traits of human intestinal Prevotella copri strains are associated with different habitual diets. Cell Host & Microbe, 25(3). https://doi.org/10.1016/j.chom.2019.01.004

[3] Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1). https://doi.org/10.1186/s13059-019-1891-0

[4] Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: Estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. https://doi.org/10.7717/peerj-cs.104

[5] Lin, H., & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. Nature Communications, 11, 3514. https://doi.org/10.1038/s41467-020-17041-7

[6] McMurdie, P. J., & Holmes, S. (2013). Phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217

[7] Oksanen, J., Blanchet, F. G., Friendly, M., Kindt, R., Legendre, P., McGlinn, D., Minchin, P. R., O’Hara, R. B., Simpson, G. L., Solymos, P., Stevens, M. H. H., Szoecs, E., & Wagner, H. (2022). vegan: Community Ecology Package (R package version 2.6-4). https://CRAN.R-project.org/package=vegan

[8] Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer.

