# Heterozygous variants in CHM13

More details are described in _Mc Cartney AM, Shafin K, Alonge M et al., Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies._ bioRxiv (2021) https://doi.org/10.1101/2021.07.02.450803.

## Curated SVs
SV calling procedure was performed on the T2T-CHM13v0.9 and T2T-CHM13v1.0 assembly using HiFi and ONT read alignments as described [here](). During manual inspection of SV calls on T2T-CHM13v0.9 assembly with Winnowmap v1.11, we found excessive read clippings in highly repetitive regions, which was fixed in v2.01. Manual inspection was performed on the newly called SVs compared to SVs called from T2T-CHM13v0.9, excluding SVs called only by ONT. The final curated set of SVs were lifted over to T2T-CHM13v1.1 using [picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard-) `LiftoverVCF`.

* Resulting curated heterozygous variants are listed as <ver.>/hets/chm13.draft_<ver.>.curated_sv.YYYYMMDD.vcf.

## Clusters of heterozygous region

Heterozygous SVs and clusters of heterozygous regions were collected on T2T-CHM13v1.0 release and lifted over to T2T-CHM13v1.1 using the command line tool of [UCSC liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).

* Resulting regions are listed as <ver.>/hets/chm13.draft_<ver.>.hets_combined.YYYYMMDD.bed

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
| Curated | Hets from SV calling methods, manually curated for the position, LEN, AF, and SV type |
| Clipped | Hets from read clippings in HiFi and ONT, curated for position, LEN, AF, and SV type |
| NucFreq | Clusters of heterozygous small variants called with [NucFreq](https://github.com/mrvollger/NucFreq)
| tQ | Heterozygous SVs and AF detected from [TandemMapper](https://github.com/ablab/TandemMapper2) |

### Curated SV region
Hets from SV calling methods, manually curated for the position, LEN, AF, and SV type.

### Clipped
Heterozygous regions from read clippings in HiFi and ONT, curated for position, LEN, AF, and SV type. `Clipped` regions from issues_raw/\<platform\>.issues.bed were manually examined.
See [Coverage analysis]() for more details.

### NucFreq
Clusters of heterozygous small variants called with [NucFreq](https://github.com/mrvollger/NucFreq) (Vollger et al., Nat Methods, 2018). 
The frequency of the first and second most common bases was calculated from the HiFi read [alignments](). using NucFreq with the following command:

```
NucPlot.py --obed {output.bed} --bed {region.bed} --minobed 2 {input.bam} {output.png}.
```

The fourth column in the resulting bed file lists the frequency of the most common base, and the fifth column lists the frequency of the second most common base. The bed files were converted to wiggle files using the following commands: 

```
echo "variableStep chrom={chromosome}" > {chromosome.wig} \
    cat {output.bed} | awk -v OFS="\t" '{print $2, $4}' >> {chromosome.wig} \
    sed '/start/d' {chromosome.wig} >> {allChromosomes_mostCommonBase.wig}
```
for the first most common base and
```
echo "variableStep chrom={chromosome}" > {chromosome.wig} \
    cat {output.bed} | awk -v OFS="\t" '{print $2, $5}' >> {chromosome.wig} \
    sed '/start/d' {chromosome.wig} >> {allChromosomes_secondMostCommonBase.wig}
```
for the second most common base.

Wiggle files are available to download or directly load into IGV or UCSC genome browser.
* v1.0
  * [most common bases](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13v1.0_hifi_20k_wm_2.01/chm13v1.0.mostCommonBase.bw)
  * [2nd most common bases](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13v1.0_hifi_20k_wm_2.01/chm13v1.0.2ndmostCommonBase.bw)

* v1.1
  * [most common bases](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13v1.1_hifi_20k_wm_2.01/chm13v1.1.mostCommonBase.bw)
  * [2nd most common bases](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13v1.1_hifi_20k_wm_2.01/chm13v1.1.2ndmostCommonBase.bw)

Heterozygous variants were also detected in the NucFreq bed file using the [hetDetection R script](https://github.com/mrvollger/NucFreq/blob/master/HetDetection.R). The resulting table indicates regions where the second most common base is present in at least 10% of reads in at least 5 positions within a 500 bp region.

### TandemMapper
Heterozygous SVs and AF were detected with [TandemMapper (developmental version, now released as v2.0.0-alpha](https://github.com/ablab/TandemMapper2).
Quality assessment with both accurate PacBio HiFi and error-prone ONT reads revealed positions in the assembly with consistent discrepancies
between mapped reads and the assembly in distances between consecutive rare k-mers (default k = 301 for HiFi reads and k = 51 for ONT reads),
a k-mer was considered rare if its frequency in the assembly did not exceed 10.

First, cenSat regions were extracted from the assembly using coordinates from 
[Nurk et al., 2021, Table S4](https://www.biorxiv.org/content/10.1101/2021.05.26.445798v1).
```
seqtk subseq $ref $bed > censat.fasta
```

Reads corresponding to these regions were extracted from Winnowmap alignments
```
samtools view -b -L $bed $output.primary.bam > $output.censat.bam
bedtools bamtofastq -i $output.censat.bam -fq $censat.reads.fastq
seqtk seq -A $censat.reads.fastq > $censat.reads.fasta
```

Then, TandemMapper was run on cenSat regions. 
```
tandemmapper2 $censat.fasta --reads $censat.reads.fasta -t $cpus -o $output_dir
cat $output_dir/$censat_kmers_dist_diff.bed > tandemmapper_issues.bed
```
Only issues with minimal AF 20 were reported. All found issues were manually validated. 
