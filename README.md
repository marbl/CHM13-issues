# CHM13-issues
CHM13 human reference genome issue tracking

For any downstream analysis, please use the following files:
* Possible consensus or mis-assembly issue: <ver.>_issues.bed
* Het sites: <ver.>/chm13.draft_<ver.>.curated_sv.20210612.vcf, <ver.>/chm13.draft_<ver.>.hets_combined.20210615.bed

## Releases
* 2022-12-02 Issues added for v2.0. X and Y were simultaneously used in T2T-HG002XYv2.7, and issues found on the Y are appended to v1.1_issues.bed. Note the sequencing data used is from HG002
* 2021-10-13 Het regions lifted over from v1.0 to v1.1
* 2021-06-23 Updating 3 additional issues and adding error k-mers in v1.0 and v1.1
* 2021-06-15 Validated het SVs and clusters of heterozygous sites in v1.0 assembly
* 2021-04-28 Issues track for HiFi and ONT read alignments from Winnowmap 2.01
* 2021-03-08 Combined low coverage and clipped regions
* 2021-02-23 Low coverage regions for HiFi, CLR, and ONT read alignments

## Issues.bed file format
| Label | Description | R,G,B | Color|
| :--- | :--- | :---: | :---: |
| Low | Low coverage | 204,0,0 | red |
| Low_Qual | Low coverage from lower consensus quality | 204,0,0 | red |
| Error_Kmer | K-mers identified as errors from the Illumina-HiFi hybrid 21-mers | 0,0,0 | black |
| Collapse | Approximate region conatining sequence collapse | 204,0,0 | red |
| Chimeric_Hap | Chimeric consensus of two haplotypes | 204,0,0 | red

## Methods
Brief descriptions are provided for
* [Coverage Analysis](https://github.com/marbl/CHM13-issues/blob/main/coverage.md)
* [Error Detection](https://github.com/marbl/CHM13-issues/blob/main/error_detection.md)
* [Heterozygous variants in CHM13](https://github.com/marbl/CHM13-issues/blob/main/het_variants.md)

More details for the polishing and evaluation methods applied on CHM13 is available in [T2T-Polish](https://github.com/arangrhie/T2T-Polish). For the methods used for polishing a nd evaluating the Y, see this [preprint](https://doi.org/10.1101/2022.12.01.518724) for more details.

## Citation

Please cite the papers below if any of the materials posted on this github are used:
- Issues found in HG002XYv2.7, especially the Y: Rhie A, Nurk S, Cechova M, Hoyt S, Taylor DJ et al., The complete sequence of a human Y chromosome. bioRxiv (2022) https://doi.org/10.1101/2022.12.01.518724
- General methods for finding issues and hets: Mc Cartney AM, Shafin K, Alonge M et al., Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. Nat Methods (2022) https://doi.org/10.1038/s41592-022-01440-3 (bioRxiv version: https://doi.org/10.1101/2021.07.02.450803)
- Issues for v0.9-v1.1: Nurk S, Koren S, Rhie A, and Rautiainen M et al., The complete sequence of a human genome. Science (2022) https://doi.org/10.1126/science.abj6987