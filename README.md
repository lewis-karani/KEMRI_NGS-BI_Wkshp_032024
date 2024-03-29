# **Assembling the SARS-CoV-2 genome: A WALK THROUGH**


Download the contents of the shared google drive folder. It should contain all your required data.

Get raw fastq.gz files here https://drive.google.com/drive/folders/1RMkGdnlRxW5rVCHIfCjuPHXJU2a1lvBU?usp=drive_link or use your own data
The SARS-CoV-2 genome reference is contained in the “genome_refs” folder
The bed file containing the ARTIC V3/4 primer coordinates are contained in the “genome_ref/ARTIC_V3”
You will work within a Linux environment, either Mac OSX, Ubuntu or the Windows Subsystem Linux (WSL) 
Download all the necessary tools you will need for this walkthrough

## Setting up the necessary tools
It is highly recommended to set up the necessary tools in an environment. This guarantees stability of the pipeline and its dependencies. 

create an environment named covmap-bsl by running the following commands using the yaml file provided i.e. _covmap-bsl.yml_

```
conda env create --file covmap-bsl.yml
```
Tools required for this walkthrough:
- fastqc - https://github.com/s-andrews/FastQC
- trimmomatic - http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf 
- fastp  - https://github.com/OpenGene/fastp
- bwa - http://bio-bwa.sourceforge.net/bwa.shtml 
- samtools - http://www.htslib.org/doc/samtools.html 
- ivar - https://github.com/andersen-lab/ivar#readme 


**Viewing the data:**
You can view the quality of the raw fastq files using fastqc by running the following commands

```
fastqc *.fastq.gz
```
This creates .html files that you can open using your web browser.

**Trimming:**
Once you determine the quality of your sequence data, proceed to trimming adapters and bad quality sequences using trimmomatic PE. Ensure your adapter file is in the path of the software in the ILLUMINACLIP option. You can view the trimmed files with fastqc.

_collapse all the adapter fasta files into one called adapters.fa_ edit <user> to your terminal username
```
cat /home/<user>/mambaforge/envs/covmap-bsl/share/trimmomatic-0.39-2/adapters/*.fa > /home/<user>/mambaforge/envs/covmap-bsl/share/trimmomatic-0.39-2/adapters/adapters.fa
``` 

Run trimmomatic on the data
```
trimmomatic PE -threads 4 WHO_1_R1_001.fastq.gz WHO_1_R2_001.fastq.gz WHO_1_R1_pair
.fastq.gz WHO_1_R1_unpair.fastq.gz WHO_1_R2_pair.fastq.gz WHO_1_R2_unpair.fastq.gz ILLUMINACLIP:/home/<user>/mambaforge/envs/covmap-bsl/share/trimmomatic-0.39-2/adapters/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70
```

**_trim multiple samples at once:_**
use the  _cov-trim.sh_ bash script under software and tools folder…

Edit the path to you adapter file with nano/vim editor  _/home/<user>/miniconda3/envs/covmap-bsl/share/trimmomatic/adapters/adapters.fa_
Save the file and exit
Make it executable by running the following command chmod +x cov-trim.sh
To run it, just type: ./cov-trim.sh *.fastq.gz assuming your script is in the same folder as the samples


Create a directory called genome_ref and download reference genome into it
```
mkdir ./genome_ref
cd ./genome_ref
wget -O sarcov2-Wu1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text"
cd ..
```
index the reference fasta
```
bwa index ./genome_ref/sarcov2-Wu1.fasta
```
**Map reads against reference:**
In this step, bwa mem is used to map the trimmed fastq files against the reference SARS-CoV-2 genome.
The output of bwa is piped “|” into samtools as input to sort the alignments, which are then output as bam files.
```
bwa mem -t 4 ./genome_ref/sarcov2-Wu1.fasta WHO_3_R1_pair.fastq.gz WHO_3_R2_pair.fastq.gz | samtools sort | samtools view -F 4 -o WHO_3_sorted.bam
```
*Remove primers:*
In this step, we use ivar to remove primers from the alignment map.
Make a directory artic_v3 inside genome_ref and download the primer bed file into the directory

```
cd genome_ref
mkdir artic_v3
cd artic_v3
wget https://github.com/lewis-karani/KEMRI_NGS-BI_Wkshp_032024/blob/main/genome_ref/artic_v3/ARTIC-V3.bed
cd ..
 
```
ivar trim -e -i WHO_3_sorted.bam -b ./genome_ref/artic_v3/ARTIC-V3.bed -p WHO_3_ptrm.bam
```
Sorting bams:
```
samtools sort WHO_3_ptrm.bam -o WHO_3_ptrm_sorted.bam
```
Calling the consensus sequence:
samtools creates the read pileups which ivar uses to call consensus at a minimum depth quality of 10.
```
samtools mpileup -A -d 1000 -B -Q 0 --reference ./genome_ref/sarscov2-Wu1.fasta WHO_3_ptrm_sorted.bam.| ivar consensus -p WHO_3.consensus -n N -m 10
```
Visualizing the read mapping with igv
convert WHO_1.sorted.bam to .sam
```
samtools view -h WHO_3_sorted.bam > WHO_3_sorted.sam
```
**GenomeQC and classification with nextclade**
Nextclade Web is a useful tool to curate your consensus genome. It also classifies your SARS-CoV-2 genome into clades and pango lineages, as well as catalogues the nucleotide and amino acid mutations in your consensus genome as compared to the Wu-H1 reference. Drag and drop your genome into nexclade for analysis.

https://clades.nextstrain.org/  

**Lineage calling with pangolin**
Pangolin is a lineage classification tool that you can use to classify your consensus sequence to help track the circulating SARS-CoV-2 lineages. Note: Pangolin and Pangolearn versions used to classify the sequence must always be documented as new data keeps refining the classification system.
 
https://pangolin.cog-uk.io/ 



## Looping through assembly steps...
step 1: mapping reada againt the reference

```
for i in `ls *.fastq.gz | cut -f1 -d'_'`; do bwa mem -t 22 ./genome_ref/sarscov2-Wu1.fasta ${i}*_R1_pair.fastq ${i}*_R2_pair.fastq | samtools sort | samtools view -F 4 -o ${i}.sorted.bam; done &
```
step 2: remove primers

```
for i in *.bam; do ivar trim -e -i ${i} -b ./genome_ref/artic_v3/ARTIC-V3.bed -p ${i}.primertrim; done &
```

step 3: sort bams

```
for i in *.primertrim.bam; do samtools sort ${i} -o ${i}.sorted;done &
```
step 4:

```
for i in *.bam.sorted; do samtools mpileup -A -d 1000 -B -Q 0 --reference ~/_SARSCOV2/genome_ref/sarscov2-Wu1.fasta ${i} | ivar consensus -p ${i}.consensus -n N; done &
```

**execute covmap-bsl.py to run multiple samples together**
```
python3 covmap-bsl.py
```




## Bonus Section: 
# phylogenetics

**Methods**
- Distance based (calculate genetic distances between sequences and construct a tree based on these distances) i.e.  Neighbour Joining (NJ)
- Character based (analyze the characters (e.g., nucleotide or amino acid sequences) directly to infer evolutionary relationships) i.e. Maximum Parsimony (MP)
- Maximum Likelihood (L): estimate the likelihood of different evolutionary models and selects the tree that maximizes this likelihood.
- Bayesian Methods: use a probabilistic model to estimate the posterior distribution of trees, taking into account prior knowledge and the likelihood of the data given the model
- Phylogenetic Network Methods: In cases where the evolutionary history of a set of taxa cannot be fully represented by a tree (e.g., due to hybridization or recombination), phylogenetic network methods are used to represent these complex evolutionary relationships.

**Core steps**

*_make alignment_** 
- mafft <- good for very large alignments
- muscle <-fast good for midsized alignments
- clustal o
- geneious
- clc genomics
- seaview
- Aliview - http://www.ormbunkar.se/aliview/downloads/windows/windows-version-1.28/without_installer_version/
- NextAlign - https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextalign-cli.html

*_estimate a substitution model_*
especially if you want to use Maximum Likelihood/ Bayesian methods 
- jmodeltest: https://github.com/ddarriba/jmodeltest2
- megaX: https://www.megasoftware.net/
- modelfinderplus 
  
*_make your tree_*
- phyml: https://github.com/stephaneguindon/phyml
- iqtree: http://www.iqtree.org/
- megaX: https://www.megasoftware.net/

*_vizualize tree_*
- figtree https://github.com/rambaut/figtree/releases/tag/v1.4.4
- itol 


Note
issues with java? download here: https://www.java.com/en/download/

Install Nextstrain CLI: https://docs.nextstrain.org/en/latest/install.html
Run Augur for Sequence Alignment:
```
augur align --sequences input.fasta --output aligned_sequences.fasta
```
Run Augur for Phylogenetic Analysis:
```
augur tree --alignment aligned_sequences.fasta --output tree.json
```
Visualize the Phylogenetic Tree: 
```
auspice view --flags

```


