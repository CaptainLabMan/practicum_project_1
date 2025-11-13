mkdir -p refs
mkdir -p reads/fastqc
mkdir -p reads/trimmed
mkdir -p reads/trimmed/fastqc
mkdir -p alignments
mkdir -p mpileup
mkdir -p vcf
mkdir -p data/k12
mkdir -p scripts
wget -P reads https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/23769689/amp_res_1.fastq.gz
wget -P reads https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/23769692/amp_res_2.fastq.gz
wget -P refs https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget -P refs https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
wget -P refs https://github.com/timflutre/trimmomatic/blob/master/adapters/NexteraPE-PE.fa
wget -P data/k12 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz
wget -P scripts https://raw.githubusercontent.com/pcingola/SnpEff/refs/heads/master/scripts/_OLD/vcfEffOnePerLine.pl
chmod +x scripts/vcfEffOnePerLine.pl