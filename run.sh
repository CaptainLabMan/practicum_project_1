gunzip -c reads/amp_res_1.fastq.gz | wc -l > reads/reads_wc_stats.txt
gunzip -c reads/amp_res_2.fastq.gz | wc -l >> reads/reads_wc_stats.txt
seqkit stats reads/amp_res_1.fastq.gz > reads/seqsit_stats_output.txt  
seqkit stats reads/amp_res_2.fastq.gz >> reads/seqsit_stats_output.txt  
fastqc -o ./reads/fastqc reads/amp_res_1.fastq.gz reads/amp_res_2.fastq.gz
echo "Trimmomatic running..."
trimmomatic PE -phred33 reads/amp_res_1.fastq.gz reads/amp_res_2.fastq.gz reads/trimmed/amp_res_1.fastq_1P.gz reads/trimmed/amp_res_1.fastq_1U.gz reads/trimmed/amp_res_2.fastq_1P.gz reads/trimmed/amp_res_2.fastq_1U.gz ILLUMINACLIP:refs/NexteraPE-PE.fa:2:30:10:2:True LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20 2> reads/trimmed/trimmomatic.log
fastqc -o ./reads/trimmed/fastqc reads/trimmed/amp_res_1.fastq_1P.gz reads/trimmed/amp_res_2.fastq_1P.gz
echo "Trimmomatic running..."
trimmomatic PE -phred33 reads/amp_res_1.fastq.gz reads/amp_res_2.fastq.gz reads/trimmed/amp_res_1.fastq_1.2P.gz reads/trimmed/amp_res_1.fastq_1.2U.gz reads/trimmed/amp_res_2.fastq_1.2P.gz reads/trimmed/amp_res_2.fastq_1.2U.gz ILLUMINACLIP:refs/NexteraPE-PE.fa:2:30:10:2:True LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:20 2> reads/trimmed/trimmomatic_2.log  
fastqc -o ./reads/trimmed/fastqc reads/trimmed/amp_res_1.fastq_1.2P.gz reads/trimmed/amp_res_2.fastq_1.2P.gz
bwa index refs/GCF_000005845.2_ASM584v2_genomic.fna.gz
bwa mem refs/GCF_000005845.2_ASM584v2_genomic.fna.gz reads/trimmed/amp_res_1.fastq_1P.gz reads/trimmed/amp_res_2.fastq_1P.gz > alignments/alignment.sam 2> alignments/bwa_mem.log
samtools view -Sb alignments/alignment.sam > alignments/alignment.bam 2> alignments/samtools_sam_to_bam.log  
samtools flagstat alignments/alignment.bam > alignments/samtools_flagstat.txt 2> alignments/samtools_flagstat.log  
samtools sort alignments/alignment.bam -o alignments/alignment_sorted.bam 2> alignments/samtools_sort.log  
samtools index alignments/alignment_sorted.bam 2> alignments/samtools_index.log
gunzip -c refs/GCF_000005845.2_ASM584v2_genomic.fna.gz > refs/GCF_000005845.2_ASM584v2_genomic.fna
samtools mpileup -f refs/GCF_000005845.2_ASM584v2_genomic.fna alignments/alignment_sorted.bam > mpileup/my.mpileup 2> mpileup/mpileup.log
varscan mpileup2snp mpileup/my.mpileup --min-var-freq 0.8 --variants --output-vcf 1 > vcf/VarScan_results.vcf 2> vcf/VarScan_results.log
touch snpEff.config
echo "k12.genome : ecoli_K12" > snpEff.config
gunzip -c data/k12/GCF_000005845.2_ASM584v2_genomic.gbff.gz > data/k12/genes.gbk
snpeff build -genbank -v k12 > data/k12/snpeff_build.txt 2> data/k12/snpeff_build.log
snpeff ann k12 vcf/VarScan_results.vcf > vcf/VarScan_results_annotated.vcf
cat vcf/VarScan_results_annotated.vcf | ./scripts/vcfEffOnePerLine.pl | snpsift extractFields - CHROM POS REF ALT ID FILTER ADP WT HET HOM NC \
"ANN[*].ALLELE" \
"ANN[*].EFFECT" \
"ANN[*].IMPACT" \
"ANN[*].GENE" \
"ANN[*].GENEID" \
"ANN[*].FEATURE" \
"ANN[*].FEATUREID" \
"ANN[*].BIOTYPE" \
"ANN[*].RANK" \
"ANN[*].HGVS_C" \
"ANN[*].HGVS_P" \
"ANN[*].CDNA_POS" \
"ANN[*].CDNA_LEN" \
"ANN[*].CDS_POS" \
"ANN[*].CDS_LEN" \
"ANN[*].AA_POS" \
"ANN[*].AA_LEN" \
"ANN[*].DISTANCE" \
"ANN[*].ERRORS" \
> VarScan_results_annotated.tsv
awk -F'\t' 'NR==1 || $29 == "0"' VarScan_results_annotated.tsv > VarScan_results_annotated_main.tsv
python3 scripts/tsv2md.py VarScan_results_annotated_main.tsv > VarScan_results_annotated_main.md
python3 scripts/tsv2md.py reads/seqsit_stats_output.txt > reads/seqsit_stats_output.md
tree -h > project_tree.txt