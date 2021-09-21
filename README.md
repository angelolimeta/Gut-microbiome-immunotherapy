# Meta analysis of gut microbiome composition in patients undergoing immunotherapy

This repository contains scripts for meta-analysis of four studies which observed that the gut microbiome composition plays a large role in predicting the treatment outcome of patients undergoing cancer immunotherapy.

* Last update: 2021-09-21

This repository is administered by Angelo Limeta ([@angelolimeta](https://github.com/angelolimeta)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Study information
This project concerns five published papers:
* [Gopalakrishnan et al. 2018, _Science_](http://science.sciencemag.org/content/359/6371/97)
* [Matson et al. 2018, _Science_](http://science.sciencemag.org/content/359/6371/104)
* [Routy et al. 2018, _Science_](http://science.sciencemag.org/content/359/6371/91)
* [Frankel et al. 2018, _Neoplasia_](https://www.sciencedirect.com/science/article/pii/S1476558617302385)
* [Peters et al. 2019, _Genome Medicine_](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0672-4)

All studies utilize metagenomic sequencing through either **Whole Genome Shotgun (WGS)** or **16S rRNA**.

### WGS Data

| Study                 | Responders | Non-responders | Total | Sequencing Method       | Database | ID         | Usage in meta analysis |
|-----------------------|-----------:|---------------:|------:|:-----------------------:|----------|:----------:|:------------------:|
| Gopalakrishnan et al. |         14 |             11 |    25 |   Illumina HiSeq 2000   | ENA      | PRJEB22893 | Discovery |
| Matson et al.         |         16 |             26 | 42(39)|  Illumina NextSeq 500   | SRA      |PRJNA399742 | Discovery |
| (Routy et al.)        |          ? |              ? |   127 | Thermofisher Ion Proton | ENA      | PRJEB22863 | Not used |
| Frankel et al.        |         24 |              9 |    33 |   Illumina HiSeq 2000   | SRA      |PRJNA397906 | Discovery |
| Peters et al.        |         18 |              9 |    27 |   Illumina HiSeq 2500   | SRA      |PRJNA397906 | Validation |

In the Peters et al. data set, response is defined as PFS >= 6 mo. for R and PFS < 6 mo. for NR. 
Only WGS samples taken prior to immune checkpoint inhibition are used. For example, the Routy et al. data set included several WGS samples from the same patient at different timepoints during treatment (*T0*,*T1*,*T2*). However, only the *T0* samples are used.

### Response Criteria

Treatment response is evaluated through the Response Evaluation Criteria in Solid Tumors 1.1 (RECIST 1.1). Patients with **complete** or **partial** response, according to RECIST 1.1 criteria, are classified as responders; whereas patients with **stable** or **progressive** disease are classified as non-responders.
