# CHM13-issues
CHM13 human reference genome issue tracking

## low_coverage.bed

This file contains regions with low coverage support in HiFi, CLR, and ONT platforms.

Reads were aligned with [Winnowmap 1.11](https://github.com/marbl/Winnowmap/releases), with the read alignments downloadable [here](https://github.com/nanopore-wgs-consortium/CHM13#downloads).

Primary read alignments have been processed with [asset commit ver. 0133f268eebf308a1c3eb356b564550526465157](https://github.com/dfguan/asset), and regions meeting the following criteria are reported in each platform:
* Less than 10 reads aligned
* Merged when regions are less than 1kbp apart
* Exclude low coverage at both 1kbp ends of each chromosome due to natural coverage drop
* Exclude low coverage overlapping rDNA gaps + 10bp around

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

Two additional files are provided using the same criteria above, except the ONT reads were prefiltered from the primary alignments.
* low_coverage.30k.bed  : Using ONT reads > 30 kbp
* low_coverage.100k.bed : Using ONT reads > 100 kbp

## het sites

TBA

