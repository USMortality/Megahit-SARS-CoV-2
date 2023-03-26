# Install deps
sudo yum update -y
sudo amazon-linux-extras install -y epel

# Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
chmod +x Miniconda3-py39_4.12.0-Linux-x86_64.sh
./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p $HOME/miniconda
sudo ln -s /home/ec2-user/miniconda/condabin/conda /usr/local/bin/conda
rm -rf Miniconda3-py39_4.12.0-Linux-x86_64.sh

# Install Megahit
# https://github.com/voutcn/megahit
conda install -y -c bioconda megahit

# Install SRA tools
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/setup-yum.sh
chmod +x setup-yum.sh
./setup-yum.sh
source /etc/profile.d/sra-tools.sh
vdb-config --simplified-quality-scores yes
printf '/LIBS/GUID = "%s"\n' $(uuidgen) >>${HOME}/.ncbi/user-settings.mkfg
rm setup-yum.sh

# Download SARS-CoV-2 contigs
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7094943/
# https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR10971381&display=data-access
if ! [ -f "SRR10971381" ]; then
  wget https://sra-pub-sars-cov2.s3.amazonaws.com/run/SRR10971381/SRR10971381
fi

# Convert SRA format data to fastq format
if ! [ -f "SRR10971381_1.fastq" ]; then
  fasterq-dump --progress SRR10971381
fi

# Install Trimmomatric for quality trimming
sudo yum install -y java
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE\
SRR10971381_{1,2}.fastq {1,2}{,un}paired.fastq.gz AVGQUAL:20 HEADCROP:12\
LEADING:3 TRAILING:3 MINLEN:75 -threads 4

# Run megahit assembly (no trimming)
# /home/ec2-user/miniconda/bin/megahit -1 SRR10971381_1.fastq -2 SRR10971381_2.fastq -o out |& tee output.txt

# Run megahit assembly (trimmed)
/home/ec2-user/miniconda/bin/megahit -1 1paired.fastq.gz -2 2paired.fastq.gz -o out |& tee output.txt

# Zip log & output.
zip -r9 out.zip out/final.contigs.fa output.txt
