#!/bin/bash

sudo apt-get update
sudo apt-get install -y default-jre mafft samtools seqtk bedtools ncbi-blast+ sra-toolkit

# Dependencies
if [ ! -d "sratoolkit" ]; then
  wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz
  mkdir sratoolkit
  tar zvxf sratoolkit.3.1.1-ubuntu64.tar.gz --strip-components=1 -C sratoolkit
fi

if [ ! -d "megahit" ]; then
  wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
  mkdir megahit
  tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz --strip-components=1 -C megahit
fi

if [ ! -d "minimap2" ]; then
  wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2
  mkdir minimap2
  tar jxvf minimap2-2.28_x64-linux.tar.bz2 --strip-components=1 -C minimap2
fi

if [ ! -d "trimmomatic" ]; then
  wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
  unzip Trimmomatic-0.39.zip -d trimmomatic
fi

export PATH=$PWD/sratoolkit/bin:$PWD/megahit/bin:$PWD/minimap2:$PATH

if ! [ -f "/usr/bin/python" ]; then
  sudo ln -s /usr/bin/python3 /usr/bin/python
fi
