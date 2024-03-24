import os
import subprocess
import glob

# Step 1: Run bwa mem and samtools sort
fastq_files = glob.glob('*_R1_pair.fastq.gz')
for fastq_file in fastq_files:
    sample_name = fastq_file.split('_')[0]
    r1_file = fastq_file
    r2_file = fastq_file.replace('_R1_', '_R2_')
    cmd1 = f"bwa mem -t 22 /home/gathii.kimita/_SARSCOV2/genome_ref/sarscov2-Wu1.fasta {r1_file} {r2_file}"
    cmd2 = "samtools sort"
    cmd3 = f"samtools view -F 4 -o {sample_name}.sorted.bam"
    cmd = f"{cmd1} | {cmd2} | {cmd3}"
    subprocess.run(cmd, shell=True, check=True)

# Step 2: Run ivar trim
bam_files = [f for f in os.listdir('.') if f.endswith('.bam')]
for bam_file in bam_files:
    cmd1 = f"ivar trim -e -i {bam_file}"
    cmd2 = f"-b /home/gathii.kimita/_SARSCOV2/genome_ref/artic_v3/ARTIC-V3.bed -p {bam_file}.primertrim"
    cmd = f"{cmd1} {cmd2}"
    subprocess.run(cmd, shell=True, check=True)

# Step 3: Sort bams
primertrim_bam_files = [f for f in os.listdir('.') if f.endswith('.primertrim.bam')]
for primertrim_bam_file in primertrim_bam_files:
    cmd = f"samtools sort {primertrim_bam_file} -o {primertrim_bam_file}.sorted"
    subprocess.run(cmd, shell=True, check=True)

# Step 4: Run samtools mpileup and ivar consensus
sorted_bam_files = [f for f in os.listdir('.') if f.endswith('.bam.sorted')]
for sorted_bam_file in sorted_bam_files:
    cmd1 = f"samtools mpileup -A -d 1000 -B -Q 0"
    cmd2 = f"--reference /home/gathii.kimita/_SARSCOV2/genome_ref/sarscov2-Wu1.fasta {sorted_bam_file}"
    cmd3 = f"| ivar consensus -p {sorted_bam_file}.consensus -n N"
    cmd = f"{cmd1} {cmd2} {cmd3}"
    subprocess.run(cmd, shell=True, check=True)

