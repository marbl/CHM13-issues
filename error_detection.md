# Error detection

Error detection and polishing was performed to remove any remaining errors detectable from read alignment and variant calling based methods. Below are provided for a practical guidance with the caveat that versions might be older than what is available today. Using the most latest and stable tool is always recommended. 

More details are described in _Mc Cartney AM, Shafin K, Alonge M et al., Chasing perfection: validation and polishing strategies for
telomere-to-telomere genome assemblies. bioRxiv (2021)_ https://doi.org/10.1101/2021.07.02.450803.

## T2T-CHM13v0.9 to T2T-CHM13v1.0
Short and long variants were called on Illumina, HiFi, CLR, and ONT alignments.
* [v0.9_patch.vcf.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v0.9_to_v1.0/v0.9_patch.vcf.gz): Edits on CHM13v0.9 coordinates, with the REF bases representing error sequences and the ALT bases being the replacement incorporated into CHM13v1.0 assembly.
* [v0.9 -> v1.0 chain](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v0.9_to_v1.0/v0.9_to_v1.0.chain)

## Alignments
In addition to the long-read alignments used in the [Coverage analysis](https://github.com/marbl/CHM13-issues/blob/main/coverage.md), PCR-free Illumina reads were aligned with [bwa mem v0.7.15](https://github.com/lh3/bwa) and removed PCR duplicate-like redundancies using [biobambam2 v2.0.87](https://github.com/gt1/biobambam2) `bamsormadup`.

## Variant calling
**S**hort **N**ucleotide **V**ariants (SNV) and *S*tructural *V*ariants (SV) were called using 
[DeepVariant v0.10~v1.2](https://github.com/google/deepvariant), [PEPPER-DeepVariant v1.0](https://github.com/kishwarshafin/pepper),
[Parliament2 v0.1.11](https://github.com/fritzsedlazeck/parliament2) and [Sniffles v1.0.12](https://github.com/fritzsedlazeck/Sniffles).

### SNV-like error detection
Both Illumina and HiFi read alignments were used to call SNPs and indels with the `hybrid` model of DeepVariant. ONT alignments were used to call SNPs using PEPPER-DeepVariant.
In addition, DeepVariant results using `wgs` or `pacbio (ccs)` model were produced but not used in the final analysis.

* All variant calling results are available to download [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/CHM13/).


DeepVariant Illumina wgs and hybrid mode was run using the following:
```
# Merge and index bams
samtools merge -@$CPUS hifi_illumina.bam hifi.bam illumina.bam
samtools index hifi_illumina.bam

# DeepVariant
run_deepvariant \
    --model_type=$MODEL \
    --ref=${INPUT_DIR}/$REF \
    --reads=${INPUT_DIR}/${BAM} \
    --output_vcf=${OUTPUT_DIR}/$OUT.vcf.gz \
    --num_shards=${THREADS}
```
Using models downloaded from [here](https://console.cloud.google.com/storage/browser/deepvariant/models/DeepVariant).

PEPPER-DeepVariant was run using the following:
```
run_pepper_margin_deepvariant call_variant \
    -b "${INPUT_DIR}/${BAM}" \
    -f "${INPUT_DIR}/${REF}" \
    -o "${OUTPUT_DIR}" \
    -p "${OUTPUT_PREFIX}" \
    -t ${THREADS} \
    --ont
```

Latest models and code versions are always recommended to use. For more information, see [DeepVariant](https://github.com/google/deepvariant/blob/r1.2/docs/deepvariant-quick-start.md).

To exclude potentially spurious variant calls, variants with low allele fraction support or low genotype quality (VAF<=0.5, GQ<=30 for Illumina/HiFi, and GQ<=25 for ONT) were removed using [bcftools v1.9~v1.10.2](https://github.com/samtools/bcftools).

```
bcftools view -f "PASS" -e 'FORMAT/VAF<=0.5 | FORMAT/GQ<=30' -Oz $HYBRID.vcf.gz > $HYBRID.PASS.vaf0.5.gq30.vcf.gz
bcftools view -f "PASS" -e 'FORMAT/VAF<=0.5 | FORMAT/GQ<=25' -Oz $ONT.vcf.gz > $ONT.PASS.vaf0.5.gq25.vcf.gz
```

The Illumina/HiFi hybrid and ONT variant calls were combined using [this custom script](https://github.com/kishwarshafin/T2T_polishing_scripts/blob/master/polishing_merge_script/vcf_merge_t2t.py).

Finally, we filtered small polishing edits using a development version of [Merfin](https://github.com/arangrhie/merfin) to ensure all retained edits did not introduce any false 21-mers that were absent from the Illumina or HiFi reads. 

```
# Merfin was run twice using 21-mers obtained from Illumina and HiFi 20kb library
merfin -vmer -memory1 30 -memory2 200 -sequence ${INPUT_DIR}/${REF} -seqmers ${INPUT_DIR}/${REF}.meryl \
  -readmers ${INPUT_DIR}/${READ}.meryl -vcf ${INPUT_DIR}/${VCF} -peak 106.8 -output illumina

merfin -vmer -memory1 30 -memory2 200 -sequence ${INPUT_DIR}/${REF} -seqmers ${INPUT_DIR}/${REF}.meryl \
  -readmers ${INPUT_DIR}/${READ}.meryl -vcf ${INPUT_DIR}/${VCF} -peak 31.8 -output hifi20k

# Only the intersection of the output ${OUTPUT_PREFIX}.polish.vcf was used.
bcftools isec -p ${OUT_DIR} hifi20k.polish.vcf.gz illumina.polish.vcf.gz
```

In Merfin v1.0 release version, we recommend to perform Merfin using Illumina k-mers to maximize error correction in HiFi based consensus
as we found prevalent SNV-like error types in HiFi (homopolymer and 2-mer microsatellite length biases) were better corrected with Illumina
reads. For more details, see this [wiki](https://github.com/arangrhie/merfin/wiki/Best-practices-for-Merfin) page and this 
[paper by McCartney et al., 2021](https://www.biorxiv.org/content/10.1101/2021.07.02.450803v1) and this 
[paper by Formenti et al., 2021](https://www.biorxiv.org/content/10.1101/2021.07.16.452324v1).

### SV-like error detection
SV callings were performed separately on short-read alignments and long-read alignments.

(1) For short-read SV calling, we used Illumina alignments as input to [Parliament2 v0.1.11](https://github.com/fritzsedlazeck/parliament2)
using default settings.
```
parliament2 --bam $BAM -r $REF_FAGZ \
    --prefix $VCF_PREF --bai $BAI --fai REF_FAI --filter_short_contigs \
    --breakdancer --breakseq --manta --cnvnator --lumpy --delly_deletion \
    --genotype --svviz_only_validated_candidates
```
SV calls longer than 30 bp and supported by at least two SV calling methods were selected for further validation.

(2) For long-read SV calling, we relied on HiFi, CLR, and ONT alignments to call SVs with Sniffles.
```
sniffles -m $BAM -v $VCF_PREF.raw.vcf -d 500 -n -1 -s 3
```

All SVs with less than 30% of reads supporting the ALT allele were removed with a [custom script](https://github.com/malonge/CallSV/blob/master/filter.py).
Remaining SV calls were refined for its ALT insertion and deletion sequences with [Iris v1.0.3](https://github.com/mkirsche/Iris) using 
[Minimap2](https://github.com/lh3/minimap2) and [Racon](https://github.com/lbcb-sci/racon) for aligning and polishing, respectively.
```
python3 filter.py $VCF_PREF.raw.vcf > $VCF_PREF.f.vcf
iris --pacbio --keep_long_variants --keep_files minimap_path=$MM2 \
    samtools_path=$SAMTOOLS \
    racon_path=$RACON genome_in=$REF \
    vcf_in=$VCF_PREF.f.vcf reads_in=$BAM \
    vcf_out=$VCF_PREF.i.vcf \
    out_dir=$VCF_PREF.iris_out
```

Our approach yielded three independent technology-specific call sets that we merged using Jasmine [v1.0.2](https://github.com/mkirsche/Jasmine).
```
jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 --output_genotypes
```
SV calls longer than 30 bp and supported by at least two technologies were selected for further validation.

SV calls from (1) and (2) were validated through manual inspection on IGV. SVs were categorized as `true heterozygous variant in CHM13`, `error with possible correction`, and `not enough evidence to determine`.

(3) Missing telomere patch
Canonical telomeric 6-mers (TTAGGG) were identified as in the [Vertebrate Genomes Project](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere)
with the commands:
```
find asm_v0.9.fasta > telomere
java -cp telomere.jar FindTelomereWindows telomere 99.9 > windows
cat windows |awk '{if ($NF > 0.5) print $2"\t"$4"\t"$5"\t"$3"\t"$NF}'|sed s/\>//g|bedtools merge -d 500 -i - -c 4 -o distinct > telomere
```

This analysis revealed one missing telomere on the p-arm of Chromosome 18. We used telomeres identified in ONT reads >50kb from 
[Miga et al, Nature, 2020](https://doi.org/10.1038/s41586-020-2547-7). Five ONT reads with telomere sequence mapping to this locus based on 
[GraphAligner](https://github.com/maickrau/GraphAligner) alignments crossing loop nodes 
(m12_2123, m12_2124, m15_26, m16_215, m20_400, utg012711l, utg032985l, utg053820l, utg053823l, utg054799l)
were manually identified. The longest was selected as template (7916a377-12c1-4251-b4fc-36b9520adb4f), all others aligned to it and polished with 
[Medaka v1.0.3](https://github.com/nanoporetech/medaka):
```
medaka -v -i reads.fasta -d template.fasta -o medaka.fasta
```

Telomere signal in all HiFi reads was identified with the commands:
```
find hifi_20kb.fasta > telomere
java -cp telomere.jar FindTelomereWindows telomere 99.9 > windows
cat windows |awk '{if ($NF > 0.5) print $2"\t"$4"\t"$5"\t"$3"\t"$NF}'|sed s/\>//g|bedtools merge -d 500 -i - -c 4 -o distinct > telomere
```

Seven additional HiFi reads were recruited from a manual analysis of the simplified version of the string graph. 
We looked for trimmed tips that could extend the end of our traversal (m15_26) and found two nodes (utg012710l and utg053821l)
consting of seven reads. All seven reads had telomere signal and were aligned to the medaka consensus and polished with Racon with the commands:
```
minimap2 -t16 -ax map-pb medaka.fasta hifi_tel.fastq > medaka.sam
racon hifi_tel.fastq medaka.sam medaka.fa > racon.fasta
```

Finally, the polished result was patched into the assembly with [PatchPolish](https://github.com/malonge/PatchPolish). 

## Consensus v1.0

The SNV-like edits, `error with possible correction` SV edits, and the Chr. 18 telomere patch from the CHM13v0.9 assembly were incorporated into one vcf file. The vcf was filtered one more time to exclude any edits in the known collapsed rDNA locus to prevent over-polishing.

Collapsed rDNA locus in CHM13v0.9:
```
chr13	5749744	10823810
chr14	2100016	2875404
chr15	2506440	5292239
chr21	3108296	6349726
chr22	4793756	5749371
```

The CHM13v1.0 assembly was generated using [bcftools v1.10](https://github.com/samtools/bcftools).
```
bcftools consensus -H 1 --chain v0.9_to_v1.0.chain -f chm13.draft_v0.9.fasta v0.9_patch.vcf.gz > chm13.draft_v1.0.fasta
```

In addition, the mitochondrial sequence (chrM) was manually rotated to start with the conventional Phenylalanine tRNA sequence
with [MitoVGP](https://github.com/VGP/vgp-assembly/tree/master/mitoVGP).


## T2T-CHM13v1.0 to T2T-CHM13v1.1
Although most of the errors were corrected, remaining errors in the telomeres were found undetected. 
Low HiFi coverage at chromosomal ends due to natural sampling dropouts resulted in lower consensus quality,
and the repetitiveness in telomeric repeats prevented reliable mapping and error detection with Illumina reads.
Therefore, a special model was trained to handle strand biases in ONT and a custom approach was taken to recover
canonical k-mer repeats. Additionally, while annotating genes, we found one true heterozygous SNP causing premature stop in _PTGS1_, on chr9:134589884 (v1.0).

* [v1.0_patch.vcf.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.0_patch.vcf.gz)
* [v1.0  -> v1.1 chain](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.0_to_v1.1_rdna_merged.chain)
* [v1.0 <-  v1.0 chain](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.1_to_v1.0_rdna_merged.chain)

### Telomere polishing
TBD

### Racon polishing
See [this page](https://github.com/arangrhie/T2T-Polish/tree/master/automated_polishing) for more details.
