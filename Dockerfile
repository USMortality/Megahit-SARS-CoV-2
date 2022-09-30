FROM amazonlinux

RUN mkdir /home/ec2-user
WORKDIR /home/ec2-user

# Install deps
RUN yum update -y
RUN amazon-linux-extras install -y epel

# Install conda
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh --silent --output miniconda.sh
RUN sh miniconda.sh -b -p miniconda
RUN ln -s /home/ec2-user/miniconda/condabin/conda /usr/local/bin/conda
RUN rm -rf Miniconda3-py39_4.12.0-Linux-x86_64.sh

# Install Megahit
# https://github.com/voutcn/megahit
RUN conda install -y -c bioconda megahit

# Install SRA tools
# https://github.com/ncbi/sra-tools
RUN yum install -y tar gzip util-linux zip
RUN curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/setup-yum.sh --silent --output setup-yum.sh
RUN chmod +x setup-yum.sh
RUN ./setup-yum.sh
RUN source /etc/profile.d/sra-tools.sh && vdb-config --report-cloud-identity no
RUN printf '/LIBS/GUID = "%s"\n' `uuidgen` >> ${HOME}/.ncbi/user-settings.mkfg
RUN rm setup-yum.sh

# Download SARS-CoV-2 contigs
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7094943/
# https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR10971381&display=data-access
RUN curl https://sra-pub-sars-cov2.s3.amazonaws.com/run/SRR10971381/SRR10971381 --output SRR10971381

# Convert SRA format data to fastq format
RUN source /etc/profile.d/sra-tools.sh && fasterq-dump --progress SRR10971381

# Run megahit assembly
RUN source /etc/profile.d/sra-tools.sh && megahit -1 SRR10971381_1.fastq -2 SRR10971381_2.fastq -o out > megahit.log

# Zip log & output.
RUN zip -r9 out.zip out/final.contigs.fa output.txt