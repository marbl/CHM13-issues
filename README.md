# CHM13-issues
CHM13 human reference genome issue tracking

For any downstream analysis, please use the following files:
* Possible consensus or mis-assembly issue: issues.bed
* Validated het sites (TBA)

## Releases
* 2021-03-08 Combined low coverage and clipped regions
* 2021-02-23 Low coverage regions for HiFi, CLR, and ONT read alignments

## Alignments
HiFi, ONT, CLR reads were aligned to the assembly using [Winnowmap 1.11](https://github.com/marbl/Winnowmap/releases). From the reads, top 0.02% of the repetitive 15-mers were collected using Meryl and provided as “bad mers” to avoid seeding. Platform specific mapping parameters were provided; map-pb for hifi, map-pb-clr for CLR, and map-ont for ONT reads, respectively. Reads were sorted and indexed using samtools, and filtered for primary reads with `samtools view -F256`. Spurious alignment blocks shorter than 1kbp or lower than 85% identity were removed from the ONT primary alignments.

These alignments were converted to pairwise alignment format (paf) and processed to collect low coverage and excessive clipped regions.

## Low coverage region
Low coverage (<10x) regions were collected for HiFi, CLR, ONT reads with [asset commit ver. 0133f268eebf308a1c3eb356b564550526465157](https://github.com/dfguan/asset) `ast_pb -M 10000000`, which collects coverage from alignment blocks trimmed 500bp on both sides of the block ends to avoid spurious mapping. Regions over 10x coverage are reported as ‘supportive’. Low coverage regions were obtained with `bedtools subtract`, and merged when regions were 1kbp apart with `bedtools merge -d 1000`. Regions overlapping rDNA gaps +/- 10bp were excluded.

### Associated files
* low_coverage.filt.bed : bed file for the region described above
* low_coverage.bed      : Using ONT primary reads with no alignment block or identity threshold
* low_coverage.30k.bed  : Same as low_coverage.bed, except using ONT reads > 30 kbp
* low_coverage.100k.bed : Same as low_coverage.bed, except using ONT reads > 100 kbp

## Clipped region
Total number of reads with more than 100 bp soft-clipped or hard-clipped bases were collected in every 1024 bp window for both HiFi and ONT reads. Clipped regions were identified when 1) >10x HiFi or ONT reads were clipped; or 2) >10% of the HiFi or >15% of the ONT total aligned reads were clipped.

### Associated files
* clipped.bed
* Absolute num. reads with clipping: hifi_pri.w1k.clip_abs.wig, ont_pri.len1k_idy85.w1k.clip_abs.wig
* Relative fraction of clipped reads compared to all reads: hifi_pri.w1k.clip_norm.wig, ont_pri.len1k_idy85.w1k.clip_norm.wig
* Value of -1 represents windows with no read alignments

### Color codes
| Found in | R,G,B | Color|
| :---: | :---: | :--- |
| HIFI only | 204,0,0 | red |
| CLR only | 255,153,0 | dark yellow |
| ONT only | 0,102,255 | dark blue |
| HIFI & ONT | 102,0,204 | dark purple |
| CLR & HIFI | 153,0,51 | dark red |
| CLR & ONT | 204,0,204 | dark pink |
| CLR & HiFi & ONT | 0,0,0 | black |
| Clipped | 153,153,153 | gray |

## het sites

TBA

