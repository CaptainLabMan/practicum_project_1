# Prerequisites
### To install all dependencies, you must have [Mamba](https://github.com/conda-forge/miniforge) installed on your system.  

🟢 **Create the environment with the following command:**
```bash
mamba env create -f environment.yml -n practicum_project_1
```  

# 1. Where to get the data.
🟢 **Automatic/semi-automatic installation of all components:**  
 - Run **setup.sh** file.  
 ```bash
    sh setup.sh  
```  
 - Donwload [**reads**](https://figshare.com/articles/dataset/amp_res_2_fastq_zip/10006541/3) manually and move them to the **/reads** directory, navigate to the folder and run the **unzip 10006541.zip** command. ⚠️ **After that, return to the main project directory.**  

# 🟢 Now you can run the **run.sh** file for automatic commands execution.  
```bash
sh run.sh  
```  

# 🟡 All subsequent commands will be executed individually.  
⚠️ **(if you ran run.sh, you don't need to execute them).**  

# 2. Inspect raw sequencing data manually.  
<details> 
<summary>Show code</summary>
</details>

<details> 
<summary>Show code</summary>

```bash
gunzip -c reads/amp_res_1.fastq.gz | wc -l  
gunzip -c reads/amp_res_2.fastq.gz | wc -l   

seqkit stats reads/amp_res_1.fastq.gz > seqsit_stats_output.txt  
seqkit stats reads/amp_res_2.fastq.gz >> seqsit_stats_output.txt  
```
</details>

🤔 **Task:** *From the line count, use what you know about the fastq format to calculate the number of reads in each file, and record in your lab notebook.*  
✅ **Answer.** Reads count. (To get the number of reads, you should to dovode these numbers by 4.):  
 - amp_res_1.fastq.gz - 1823504  
 - amp_res_2.fastq.gz - 1823504  

# 3. Inspect raw sequencing data with FastQC. Filtering the reads.  
<details> 
<summary>Show code</summary>

``` bash
fastqc -o ./reads/fastqc reads/amp_res_1.fastq.gz reads/amp_res_2.fastq.gz  
```
</details>

🤔 **Task:** *Do the basic statistics match what you calculated for the number of reads last time?*  
✅ **Answer:** Yes  

🤔 **Task:** *On the left, you’ll see a navigation window with green (normal), yellow (slightly abnormal), and red (very unusual) circles for several kinds of data analysis. If you have any red circles, record them in your notebook:*  
✅ **Answer:**  
> amp_res_1.fastq.gz: Per base sequence quality, Per tile sequence quality  
> amp_res_2.fastq.gz: Per base sequence quality  

🤔 **Task:** *Mention the QC results in your lab report.*  
✅ **Answer. FastQC results:**  
> **Per base sequence quality:**  
>> The FastQC analysis revealed a decrease in read quality towards the end of both forward and reverse reads.

> **Per tile sequence quality:**  
>> A significant drop in quality was detected for forward reads in specific tiles of the flow cell (red and yellow areas). Reverse reads also showed regions of reduced quality, but the result is acceptable. Overall, the general result remains acceptable due to the limited and localized nature of these deviations. These issues may have been caused by the presence of a bubble or edge effects in the flow cell.

🤔 **Task:** *What do you think we should do about anything FastQC identified as unusual?*  
✅ **Answer:**
> It depends on the results. In this case, we can remove reads or their parts that do not meet the quality criteria.

# 4. (Optional, 1 bonus point) Filtering the reads. 
<details> 
<summary>Show code</summary>

```bash 
trimmomatic PE -phred33 reads/amp_res_1.fastq.gz reads/amp_res_2.fastq.gz reads/trimmed/amp_res_1.fastq_1P.gz reads/trimmed/amp_res_1.fastq_1U.gz reads/trimmed/amp_res_2.fastq_1P.gz reads/trimmed/amp_res_2.fastq_1U.gz ILLUMINACLIP:refs/NexteraPE-PE.fa:2:30:10:2:True LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20 2> reads/trimmed/trimmomatic.log

fastqc -o ./reads/trimmed/fastqc reads/trimmed/amp_res_1.fastq_1P.gz reads/trimmed/amp_res_2.fastq_1P.gz
```  
</details>

📈 **Trimmomatic output:** Input Read Pairs: 455876 Both Surviving: 430758 (94,49%) Forward Only Surviving: 9340 (2,05%) Reverse Only Surviving: 527 (0,12%) Dropped: 15251 (3,35%)

🤔 *What happens if we increase the quality score at all steps to 30? Try to modify the previous command (be sure to name them something distinct, so as not to overwrite your data).*  

<details> 
<summary>Show code</summary>

```bash
trimmomatic PE -phred33 reads/amp_res_1.fastq.gz reads/amp_res_2.fastq.gz reads/trimmed/amp_res_1.fastq_1.2P.gz reads/trimmed/amp_res_1.fastq_1.2U.gz reads/trimmed/amp_res_2.fastq_1.2P.gz reads/trimmed/amp_res_2.fastq_1.2U.gz ILLUMINACLIP:refs/NexteraPE-PE.fa:2:30:10:2:True LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:20 2> reads/trimmed/trimmomatic_2.log 

fastqc -o ./reads/trimmed/fastqc reads/trimmed/amp_res_1.fastq_1.2P.gz reads/trimmed/amp_res_2.fastq_1.2P.gz
```
</details>

✅ **Answer:**
> **Total Sequences**
> **Before trimming:**
>> amp_res_1.fastq.gz - 455876  
>> amp_res_2.fastq.gz - 455876  

> **After trimming (qual=20):**
>> amp_res_1.fastq_1P.gz - 430758  
>> amp_res_2.fastq_1P.gz - 430758

> **After trimming (qual=30):**
>> amp_res_1.fastq_1.2P.gz - 363413  
>> amp_res_2.fastq_1.2P.gz - 363413

# 5 . Aligning sequences to reference
## 5.1 Index the reference file 
```bash
mamba install bioconda::bwa  
bwa index refs/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

*Record the command you used in your lab notebook.*  
Смотреть выше

## 5.2 Align your reads
```bash
bwa mem refs/GCF_000005845.2_ASM584v2_genomic.fna.gz reads/trimmed/amp_res_1.fastq_1P.gz reads/trimmed/amp_res_2.fastq_1P.gz > alignments/alignment.sam 2> alignments/bwa_mem.log
```

## 5.3. Compress SAM file
```bash
 mamba install bioconda::samtools  
 samtools view -Sb alignments/alignment.sam > alignments/alignment.bam 2> alignments/samtools_sam_to_bam.log  
 samtools flagstat alignments/alignment.bam > alignments/samtools_flagstat.txt 2> alignments/samtools_flagstat.log
```  

*What percentage of reads are mapped?*  
 - 860682 + 0 mapped (99.87% : N/A)  
 - 860426 + 0 primary mapped (99.87% : N/A)  

 ## 5.4 Sort and index BAM file
 ```bash
 samtools sort alignments/alignment.bam -o alignments/alignment_sorted.bam 2> alignments/samtools_sort.log  
 samtools index alignments/alignment_sorted.bam  2> alignments/samtools_index.log
 ```

 ```bash
 mamba install bioconda::igv  
 gunzip -c refs/GCF_000005845.2_ASM584v2_genomic.fna.gz > refs/GCF_000005845.2_ASM584v2_genomic.fna
 ```

 #  6. Variant calling
 ```bash
  samtools mpileup -f refs/GCF_000005845.2_ASM584v2_genomic.fna alignments/alignment_sorted.bam > mpileup/my.mpileup 2> mpileup/mpileup.log  
  mamba install bioconda::varscan  
  varscan mpileup2cns mpileup/my.mpileup --min-var-freq 0.8 --variants --output-vcf 1 > vcf/VarScan_results.vcf 2> vcf/VarScan_results.log
 ```  

 # 7. Variant effect prediction
PASS

# 8. Automatic SNP annotation
```bash
mamba install bioconda::snpeff
mamba install bioconda::snpsift

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
```