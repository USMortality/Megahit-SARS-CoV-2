# Megahit SARS-CoV-2

This repo contains a 1-click script to configure an Amazon (Amazon Linux) Server with all necessary tools, to run the de-novo assembly of the claimed SARS-CoV-2 sequence with Megahit.

# Usage
Either run the `./megahit.sh` file on a Amazon EC2 Linux instance, or run the `Dockerfile` via `docker build . -t megahit -f ./Dockerfile --platform linux/amd64`.

# Output
`out/` folder contains the output of a run, that was carried out on a 32 core Amazon c5 instance, which took about ~15 minutes or so to assemble.

# Result
Megahit result:
`29463 contigs, total 14438186 bp, min 200 bp, max 29802 bp, avg 490 bp, N50 458 bp`

The longest sequence generated by Megahit is 29,802bp long.

Here we are not able to generate the exact claimed SARS-CoV-2 sequences, that Wu et al. (2020) had published (30,473bp/29,875bp/29,903bp).

# Run

```
aws ec2 run-instances --image-id ami-026b57f3c383c2eec --count 1 --instance-type c5a.8xlarge --key-name ben --ebs-optimized --block-device-mapping "[ { \"DeviceName\": \"/dev/xvda\", \"Ebs\": { \"VolumeSize\": 100 } } ]"
```

`vim ~/.ssh/config`

```
Host megahit
Hostname 18.232.52.37
Port 22
User ec2-user
IdentityFile ~/.ssh/id_rsa
IdentitiesOnly yes
```

```
scp megahit.sh megahit:/home/ec2-user/. 
ssh megahit  
```

```
chmod +x megahit.sh
screen
./megahit.sh
```


scp megahit:/home/ec2-user/out.zip .

# Download original genomes

```
curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=MN908947.1'>out/MN908947.1.fasta
curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=MN908947.2'>out/MN908947.2.fasta
curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=MN908947.3'>out/MN908947.3.fasta
```

# Show differences
```
cat MN908947.1.fasta MN908947.3.fasta|mafft -|awk '{printf(/^>/?"\n%s\n":"%s"),$0}'|grep .>temp.aln;paste <(sed -n 2p temp.aln|grep -o .) <(sed -n 4p temp.aln|grep -o .)|awk '$1!=$2{print NR,$1,$2}'
```

# Compare the produced sequence to the original MN908947.3 Wuhan-Hu-1 Isolate
```
# Make Blast DB of target sequence
makeblastdb -in MN908947.3.fasta -dbtype nucl  

# Invert the sequence (end to front)
seqtk seq -r k141_13590.fasta > k141_13590r.fasta

# Compare thes genomes
blastn -query k141_13590r.fasta -db MN908947.3.fasta -evalue 1 -task megablast -outfmt 6 > k141_13590r_MN908947.3.crunch
```