# Install deps
sudo yum install -y https://repo.ius.io/ius-release-el7.rpm
sudo yum update -y
sudo yum install -y python36u python36u-libs python36u-devel python36u-pip
sudo amazon-linux-extras install -y epel

# Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
chmod +x Miniconda3-py39_4.12.0-Linux-x86_64.sh
./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p $HOME/miniconda
sudo ln -s /home/ec2-user/miniconda/condabin/conda /usr/local/bin/conda
rm -rf Miniconda3-py39_4.12.0-Linux-x86_64.sh

# Install Megahit
# https://github.com/voutcn/megahit
conda create --name megahit python=3.6 pip -y
conda init bash
source /home/ec2-user/.bashrc
conda activate megahit
conda install -y -c bioconda megahit=1.1.3

# Install SRA tools
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/setup-yum.sh
chmod +x setup-yum.sh
./setup-yum.sh
source /etc/profile.d/sra-tools.sh
vdb-config --simplified-quality-scores yes
printf '/LIBS/GUID = "%s"\n' `uuidgen` >> ${HOME}/.ncbi/user-settings.mkfg
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

# Run megahit assembly
megahit -1 SRR10971381_1.fastq -2 SRR10971381_2.fastq -o out |& tee output.txt

# Zip log & output.
zip -r9 out.zip out/final.contigs.fa output.txt
