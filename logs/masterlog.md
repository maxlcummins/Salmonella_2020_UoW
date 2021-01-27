#  University of Wollongong Salmonella collection analysis

## Collection description
211 strains of putative Salmonella sourced from gulls and clinical settings near Wollongong.

These were identified as *Salmonella* using a combination of microbiological tests and MALDI-TOF analysis. Here I have documented the analytical workflow ranging from quality control and assembly of WGS data to its bioinformatic analysis.

---

## Configuration
Here we setup some variables which will be inherited by future steps.

```
# Define the pipelord directory where our pipelines are kept
PIPELORD_DIR=~/Data/pipelord/

# Define the project directory where we will send our outputs
PROJ_DIR=/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all

# Make the appropriate subdirectories for our analysis
cd $PROJ_DIR
mkdir output
mkdir logs
mkdir snakemake
mkdir scripts

#
dt=$(date '+%d/%m/%Y %H:%M:%S')

#
dt=$(date '+%d/%m/%Y %H:%M:%S')
MASTER_CONF=/projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/snakemake/masterconfig.yaml
```

You must also modify the masterconfig.yaml configuration file before proceeding, however.

* PLACEHOLDER


___
# Job Running
## FastP - Quality control of reads
fastp version: 

Snakefile:

```
# Set the task variable for our output names
TASK=fastp

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/fastp.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

## Shovill - Genome assembly
Shovill version: 1.0.4

Shovill yaml: /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/Snakemake/shovill/shovill.yaml

```
# Set the task variable for our output names
TASK=shovill

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/shovill.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

## Kraken2 - Species identification
Prior to any further analysis we will screen our samples to try and confirm their genus and species identification. This prevents future downstream issues, particularly with determining phylogenetic distance.

kraken2 version: 2.1.1

kraken2 yaml: /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/snakemake/kraken2/kraken.yaml

```
# Set the task variable for our output names
TASK=kraken2

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/kraken2.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

## Sistr - Serovar/ Species identification
We will couple sistr analysis with kraken data to determine species ID and qc our assemblies.

sistr version: 

sistr yaml: /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/snakemake/sistr/sistr.yaml

```
# Set the task variable for our output names
TASK=sistr

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/*.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

## Genotyping with abricate
We can start out genotyping as we can filter out non-Salmonella later anyway quite readily and computationally it doesnt hurt to do the extra strains.

abricate version: 

abricate yaml: /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/snakemake/abricate/abricate.yaml

```
# Set the task variable for our output names
TASK=abricate

# Change to the pipelord directory
cd ${PIPELORD_DIR}/${TASK}lord

# Create a directory for our task output log
mkdir ${PROJ_DIR}/snakemake/${TASK}

#Change Snakefile config variable to our master config
perl -p -i -e "s@^configfile.*@configfile: \"${MASTER_CONF}\"@g" Snakefile

# Copy the Snakefile and Environment yaml/s to our project directory
cp Snakefile config/*.yaml ${PROJ_DIR}/snakemake/${TASK}

# Create a directory for our task output log
mkdir ${PROJ_DIR}/logs/${TASK}

# Run task
nohup snakemake -j --use-conda  -p > ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.err 2> ${PROJ_DIR}/logs/${TASK}/nohup_${TASK}.out &

# Save our job ID
PROCESS_ID=$!
echo "$dt" "$PWD" "JOB_ID =" "$PROCESS_ID" >> ${PROJ_DIR}/logs/${TASK}/JOB_IDs
```

# Downstream analysis and species ID




Strains removed from analysis

#removed due to assembly size (Cut off at 8.2 mbp)
-rw-r--r-- 1 malcummi AusGEM 9578573 Jan  8 13:42 Seagull_18_124_2_S27.fasta
-rw-r--r-- 1 malcummi AusGEM 9669035 Jan  8 13:42 Seagull_18_116_S20.fasta
-rw-r--r-- 1 malcummi AusGEM 9750944 Jan  8 16:09 SIML_2017_12.fasta
-rw-r--r-- 1 malcummi AusGEM 9769862 Jan  8 13:42 Seagull_18_127_S28.fasta

#Many Salmonella were found to be Proteus!
#These include:
Seagull_18_119_S25
Seagull_18_229_S54
Seagull_18_242_S60
Seagull_18_256_S69
Seagull_18_213_S50
Seagull_18_234_S55
Seagull_18_262_S113
Seagull_18_254_S112
Seagull_18_49_S5
Seagull_18_160_S36
Seagull_18_261_S73
Seagull_18_265_S115
Seagull_18_196_S44
Seagull_18_210_S49
Seagull_18_263_S114
Seagull_18_8_S1
Seagull_18_205_S108
Seagull_18_92_S11
Seagull_18_244_S61
Seagull_18_245_S62
Seagull_18_249_S65
Seagull_18_246_S63
Seagull_18_202_1_S45
Seagull_18_40_S3
Seagull_18_202_2_S46
Seagull_18_94_S13
Seagull_18_217_S52
Seagull_18_216_S51
Seagull_18_87_S7
Seagull_18_118_S24
Seagull_18_258_S70
Seagull_18_238_S57
Seagull_18_260_S72
Seagull_18_34_S2
Seagull_18_266_S74
Seagull_18_70_S6
Seagull_18_183_S41
Seagull_18_156_S35
Seagull_18_241_S111
Seagull_18_166_S39
Seagull_18_235_S56
Seagull_18_208_S48
Seagull_18_248_S64
Seagull_18_252_S67
Seagull_18_267_S75
Seagull_18_136_S31

/projects/AusGEM/i3_genomes/UoW_Salmonella_all/rm_proteus.sh


################################################################################################################################################
########################       Shuvlord - Assembly                                             #################################################
################################################################################################################################################

#Shovill
nohup snakemake -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/shovill/nohup_shuvlord.out &
#JOBID = 23171

################################################################################################################################################
########################       Pangenomelord - Pangenomic analysis                               ###############################################
################################################################################################################################################

#Some of the samples had to be rerun as the wrong reads were being used - these were created with read index files rather than sequence read data.
nohup snakemake -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/pangenomelord/nohup_pangenomelord.out &
#JOBID = 1724


cd /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/output/pangenomelord/Roary.out
mkdir ../initial_pangenome_211
cp core_gene_alignment.aln ../initial_pangenome_211





## WENT TO RUN THE step below but then found this in the logs when I couldnt locate the aln file...
#Number of clusters (72887) exceeds limit (50000). Multifastas not created. Please check the spreadsheet for contamination from different species or increase the --group_limit parameter.
#2021/01/12 14:37:04 Exiting early because number of clusters is too high
## Run SNP-sites
#source deactivate
#source activate /home/malcummi/Data/pipelord/snplord/.snakemake/conda/88bf0609 # snp_sites from snplord
#snp-sites -c core_gene_alignment.aln > ../initial_pangenome_211/core_gene_alignment_snp_sites.aln
#source deactivate

#Deleted genomes above 8.2MB in size, deleted ROARY output and reran pangenomelord
cd /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/output/shovill/final_assemblies
mkdir ../junk_assemblies
mv Seagull_18_124_2_S27.fasta Seagull_18_116_S20.fasta SIML_2017_12.fasta Seagull_18_127_S28.fasta ../junk_assemblies/

cd ../../shovill
rm -rf Seagull_18_124_2_S27.out Seagull_18_116_S20.out SIML_2017_12.out Seagull_18_127_S28.out Roary.out

cd ~/Data/pipelord/pangenomelord
source deactivate
source activate snakemake

nohup snakemake -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/pangenomelord/nohup_pangenomelord2.out &
#JOBID = 28433

#Same error:
##Use of uninitialized value in require at /data/malcummi/pipelord/pangenomelord/.snakemake/conda/2a112bc4/lib/site_perl/5.26.2/x86_64-linux-thread-multi/Encode.pm line 61.
##Number of clusters (65766) exceeds limit (50000). Multifastas not created. Please check the spreadsheet for contamination from different species or increase the --group_limit parameter.
##2021/01/12 16:06:59 Exiting early because number of clusters is too high

#Reran the command with a new snakefile which has a different --group_limit parameter (70000) (Snakefile_Salmonella_group_lim_up)

nohup snakemake -s Snakefile_Salmonella_group_lim_up -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/pangenomelord/nohup_pangenomelord3.out &
#JobID = 29500


## Ran SNP-sites
snp-sites -c core_gene_alignment.aln > core_gene_alignment_snp_sites.aln
source deactivate

source activate /home/malcummi/Data/pipelord/snplord/.snakemake/conda/4ead7ed0 # snp_dists from snplord
snp-dists -c core_gene_alignment_snp_sites.aln > core_gene_alignment_snp_sites_snp_dists.csv
source deactivate

#source activate /home/malcummi/Data/pipelord/snplord/.snakemake/conda/89835887 # fasttree from snplord
#fasttree -gtr -nt ../output/core_gene_alignment_snp_sites.aln > ../output/core_gene_alignment_snp_sites.tree

### Cam has suggested using IQtree on the core_gene_alignment.aln instead
conda create -n iqtree -c bioconda iqtree
source activate iqtree
cd ../output

#This was run again on the snp_sites
nohup iqtree -s core_gene_alignment_snp_sites.aln -m MFP -bb 1000 -nt AUTO >iqtree_core_genome_snp_sites_aln.out 2>iqtree_core_genome_snp_sites_aln.err &
#JOBID=13708

#This was run again on the full alignment (rather than the snp_sites)
nohup iqtree -s core_gene_alignment.aln -m MFP -bb 1000 -nt AUTO >iqtree_core_genome_aln.out 2>iqtree_core_genome_aln.err &
#JOBID=32647

################################################################################################################################################
########################       sistrlord - serovar analysis                                      ###############################################
################################################################################################################################################


nohup snakemake -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/sistrlord/nohup_sistrlord.out &
#JOBID = 15401

nohup snakemake -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/sistrlord/nohup_sistrlord2.out &
#JOBID =12906

nohup snakemake -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/sistrlord/nohup_exploration_wild_animals_enterobase_sistrlord.out &
#JOBID =9142



################################################################################################################################################
########################       abricatelord - abricate analysis                                      ###############################################
################################################################################################################################################


nohup snakemake -j --use-conda  -p > /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/logs/abricate/nohup_abricatelord.out &
#JOBID = 23304

cd /projects/AusGEM/Users/Max/Manuscripts/UoW_Salmonella_all/output/abricate/Salmonella

for f in *; do cat ${f}/*.tab > ${f}.txt; done

cat *.txt > Salmonella_abricate.txt