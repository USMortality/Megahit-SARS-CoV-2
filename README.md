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