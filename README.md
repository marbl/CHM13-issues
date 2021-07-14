# CHM13-issues
CHM13 human reference genome issue tracking

For any downstream analysis, please use the following files:
* Possible consensus or mis-assembly issue: <ver.>_issues.bed
* Het sites: <ver.>/chm13.draft_<ver.>.curated_sv.20210612.vcf, <ver.>/chm13.draft_<ver.>.hets_combined.20210615.bed

## Releases
* 2021-06-23 Updating 3 additional issues and adding error k-mers in v1.0 and v1.1
* 2021-06-15 Validated het SVs and clusters of heterozygous sites in v1.0 assembly
* 2021-04-28 Issues track for HiFi and ONT read alignments from Winnowmap 2.01
* 2021-03-08 Combined low coverage and clipped regions
* 2021-02-23 Low coverage regions for HiFi, CLR, and ONT read alignments

## Alignments
HiFi and ONT reads were aligned to the assembly using [Winnowmap 2.01](https://github.com/marbl/Winnowmap/releases). From the reads, top 0.02% of the repetitive 15-mers were collected using Meryl and downweighted in the process of minimizer sampling. Platform specific mapping parameters were provided; map-pb for hifi and map-ont for ONT reads, respectively. Reads were sorted and indexed using samtools, and filtered for primary reads with `samtools view -F0x100`.

*Note*: We did not perform this analysis on CLR reads as the information added was marginal or noisier compared to HiFi and ONT.


```shell
# Collect repetitive 15-mers
meryl count k=15 $ref output merylDB
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

# Align
winnowmap --MD -W repetitive_k15.txt -ax map-$platform -t$cpus $ref $reads > $output.sam

# Sort and filter
samtools sort -@$cpus -m2G -O bam -o $output.bam $output.sam
samtools view -F0x104 -hb $output.srt.bam > $output.primary.bam
```

These alignments were converted to pairwise alignment format (paf) using [paftools from minimap2 v2.17](https://github.com/lh3/minimap2/tree/master/misc) and processed to collect low coverage and excessive clipped regions.
Spurious alignment blocks shorter than 1kb or lower than 85% identity were removed from the ONT primary alignments.

```shell
# Convert to paf
samtools view -h -@$cpus $output.primary.bam | k8 paftools.js sam2paf - | cut -f 1-16 - >  $output.primary.paf

# Filter ONT for alignment blocks < 1kb and identity < 85%
awk '$11>1000 && $10/$11>0.85' ont.primary.paf > ont_pri.len1k_idy85.paf
```

## Low coverage region
Supportive regions were collected for HiFi (>7x) and ONT (>10x) reads with [asset commit v102, 0133f268eebf308a1c3eb356b564550526465157](https://github.com/dfguan/asset).

```shell
# Collect supportive regions from hifi
ast_pb -m 7 -M 10000000 -l 0 $paf > $pre.bed

# Collect supportive regions from ONT
ast_pb -m 10 -M 10000000 $paf > $pre.bed
```

Alignment blocks were trimmed 300bp on both sides to avoid spurious mapping while collecting coverage in ONT. Low coverage regions were obtained with `bedtools subtract`, and merged when regions were 5kb apart with `bedtools merge -d 5000`. Regions overlapping collapsed rDNA region (only applied in v1.0 assembly) and 3kb around the end of chromosomes were excluded.

To distinguish the cause of low coverage, we collected % GA/TC, GC, and AT micro satellite repeats in every 128 bp window in the assembly that appears as dimers when compressing homopolymers. Any low coverage region overlapping with >80% of any micro satellite repeat within 10kb distance was marked accordingly.

In addition, regions overlapping with known consensus issues obtained with [Merqury](https://github.com/marbl/merqury) Illumina-HiFi hybrid 21-mers were marked as `Low_Qual`.

## Clipped region
Total number of reads with more than 100 bp soft-clipped or hard-clipped bases were collected in every 1024 bp window for both HiFi and ONT reads. Clipped regions were identified when 1) >10x HiFi or ONT reads were clipped; or 2) >10% of the HiFi or >15% of the ONT total aligned reads were clipped.

### Associated files
* <ver.>_issues.bed: Low coverage not associated with sequencing biases and other issues

| Label | Description | R,G,B | Color|
| :--- | :--- | :---: | :---: |
| Low | Low coverage | 204,0,0 | red |
| Low_Qual | Low coverage from lower consensus quality | 204,0,0 | red |
| Error_Kmer | K-mers identified as errors from the Illumina-HiFi hybrid 21-mers | 0,0,0 | black |
| Collapse | Approximate region conatining sequence collapse | 204,0,0 | red |
| Chimeric_Hap | Chimeric consensus of two haplotypes | 204,0,0 | red


* issues_raw/\<platform\>.issues.bed : Low coverage + Flanking clipped regions

| Label | Description | R,G,B | Color|
| :--- | :--- | :---: | :---: |
| Low | Low coverage | 204,0,0 | red |
| Low_Qual | Low coverage from lower consensus quality | 204,0,0 | red |
| Low_GA/TC, Low_AT | Low coverage from sequencing bias | 153,153,255 | light purple |
| Clipped | Region with excessive read clipping | 153,153,153 | gray |


#### Other files under issues_raw:
* clipped.bed
* Absolute num. reads with clipping: hifi_pri.w1k.clip_abs.wig, ont_pri.len1k_idy85.w1k.clip_abs.wig
* Relative fraction of clipped reads compared to all reads: hifi_pri.w1k.clip_norm.wig, ont_pri.len1k_idy85.w1k.clip_norm.wig
* Value of -1 represents windows with no read alignments

## Het sites

Heterozygous sturcutral variants and clusters of heterozygous regions collected on v1.0 release.

A manuscript describing the details of each variant calling methods are in preparation and will be posted soon.

### Clusters of heterozygous regions: chm13.draft_<ver.>.hets_combined.YYYYMMDD.bed

| Abbreviation | Description |
| :---: | :--- |
| ID | Numbered ID |
| SC | Source of detection method |
| AF | Allele frequency
| PL | Platform support |
| LEN | SV size |
| (TYPE) | INS, DEL, DUP, INV for Insertion, Deletion, Duplication, Inversion |

| SC | Source |
| :---: | :--- |
| Curated | Hets from SV calling methods, manually curated for the position, LEN, AF, and SV type (manusript in prep.) |
| Clipped | Hets from read clippings in HiFi and ONT, curated for position, LEN, AF, and SV type (manusript in prep.) |
| NucFreq | Clusters of heterozygous small variants called with [NucFreq](https://github.com/mrvollger/NucFreq)
| tQ | Heterozygous SVs and AF detected from [TandemMapper and TandemQUAST](https://doi.org/10.1093/bioinformatics/btaa440) |


### Curated SVs: chm13.draft_<ver.>.curated_sv.YYYYMMDD.vcf

SVs manually validated for breakpoints, CCS AF, and SV type. CLR and ONT AF may be slightly lower than the actual AF due to mapping biases.

## Citation

Please cite the paper below if any of the materials posted on this github are used:

Mc Cartney AM, Shafin K, Alonge M et al., Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. bioRxiv (2021) https://doi.org/10.1101/2021.07.02.450803
Nurk, Koren, Rhie, and Rautiainen et al., The complete sequence of a human genome. bioRxiv (2021) https://doi.org/10.1101/2021.05.26.445798
