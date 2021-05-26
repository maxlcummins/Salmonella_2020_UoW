# Salmonella genomics analysis

## Outline of collection & primary QC
211 putative Salmonella genomes were recovered by University of Wollongong collaborators from gulls and humans. Kraken analysis later revealed many of these to be Genera other than Salmonella.

| Genus           | counts |
| --------------- | ------ |
| Salmonella      | 145    |
| Salmonella\*     | 7      |
| Proteus         | 28     |
| Citrobacter     | 14     |
| Hafnia          | 7      |
| Lelliottia      | 4      |
| Enterobacter    | 2      |
| Escherichia     | 2      |
| Obesumbacterium | 1      |
| Serratia        | 1      |

*Note: * indicates low confidence in Salmonella ID - were ID'd as Salmonella but came with a FAIL or WARNING flag from Kraken.*

All strains which were Salmonella or Salmonella\* were uploaded to Enterobase for additional QC.

Six of the Salmonella* were accepted by Enterobase:

* ISLHD_2017_1_S76
* SIML_2017_12
* SIML_2017_17
* W_1_C11
* W_1_E11
* ISLHD_2019_26

One of these Salmonella\* were rejected by Enterobase:

* SIML_2017_20

Of those Salmonella which passed our Kraken analysis, some were rejected by enterobase:

* Seagull_18_89_S8
* Seagull_18_154_S34
* SIML_2017_11

All other strains which passed Kraken also passed Enterobase, leaving us with 148.

A list of our strains in the final subset can be foud at:
 `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/delims/all_Salm_assemblies.txt`

*Note that the strain name for SG19_116.05_E.coli in enterobase is "SG19_116.05_E_coli"*

This highlights the need to use multiple approaches for QC.

## Genomic analysis

Firstly, pipelord2.0 was used to assemble and analyse our samples. This was firstly run to generate pMLST, genotype, AMR-snps, MLST, kraken and genomic annotation data. The config file used for this run is available at `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/snakemake/Salmonella_UoW_config.yaml`, as can the snakemake workflow at `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/snakemake/workflow`

**Note this was performed for ALL isolates, Salmonella or otherwise.**

```
nohup snakemake -j all --use-conda -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_out.log 2> /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_err.log &

JOBID = 26407
```

# Pangenomic analysis and tree building

Pangenomic analysis was performed on multiple subsets of strains, including the entire cohort and a subset of *Salmonella enterica* serovar Typhimurium.

## Typhimurium strains

Enterobase's SISTR serovar determinations were used to select a subset of *S.* Typhimurium strains. This list of sample names was then exported to a text file: `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/delims/Typhimurium_Salm_assemblies.txt`


Snakemake was ran again using a different snakefile which subsets the collection based on a given list of sample names. `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/snakemake/workflow`

```
nohup snakemake -j all -s workflow/Treebuild.smk -p --use-conda > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_Typhimurium_Roary_out.log 2> /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_Typhimurium_Roary_err.log &

JOBID = 208155
```

Upon producing a phylogenetic tree it became apparent that one of the strains was particularly divergent from the Typhimurium cohort, SIML_2019_11. Analysis of the MLST data indicated it was of ST 5060. There was also a strain of ST 2089 but this sat within the phylogeny neatly, indicating it is perhaps of the same CC as ST19 and the rest of the Typhimurium cohort.

|   ST |counts|
|------|------|
|   19 |    45|
| 2089 |     1|
| 5060 |     1|

___

![Figure 1](../misc/Virotype_first_typhi.png "Figure 1")
Figure 1 - Typhimurium phylogenetic tree adjacent to VAG data. Shown to the left of the tree is ST data, while tip labels are shown in blue and red representing human and gull sourced isolates, respectively.

The branching in this tree is obscured by the presence of SIML_2019_11, therefore a new tree needed to be generated. An additional fofn was generated wherein this strain was removed and the pangenomic analysis and tree building was performed again. `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/delims/Typhi_minus_outlier.txt`


```
nohup snakemake -j all -s workflow/Treebuild.smk -p --use-conda > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_Typhimurium_new_Roary_out.log 2> /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_Typhimurium_new_Roary_err.log &

JOBID = 150892
```

## All Salmonella strains

The tree building snakemake workflow was run again after modifying the config file. Within the `Salmonella_UoW_config.yaml` file, Lines containing the path to the file of file names \(fofn\) and the output prefix were hashed out and new lines containing pertaining to a fofn for the Salmonella strains was added along with a new prefix

```
nohup snakemake -j all -s workflow/Treebuild.smk -p --use-conda > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_roary_out.log 2> /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/Salm_UoW_roary_err.log &

#JOBID = 267521
```


## Rerun in May

I had to remove SIML_2019_03 (because it was a subspecies diarizonae and a massive outlier), SIML_2017_17 and SIML_2017_12 because they were low quality assemblies \(as determined by failed SISTR serovaring and crappy genotype data. Note these strains were still accepted by enterobase!\)

Therefore I edited the previous FOFNs to remove these lines and wrote the Salmonella full cohort fofn to `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/delims/May_all_Salm_assemblies.txt` and the Salmonella Typhimurium FOFN to `/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/delims/May_Typhimurium.txt`. I decided to run the pangenome on the full cohort first as it is the slowest to run.

```
nohup snakemake -j all -s workflow/Treebuild.smk -p --use-conda > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/May_Salm_UoW_roary_out.log 2> /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/May_Salm_UoW_roary_err.log &
[1] 21219
```