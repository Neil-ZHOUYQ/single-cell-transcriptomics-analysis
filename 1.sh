
prefetch -h
#Download SRA data
prefetch -X 200GB --option-file SRR_Acc_List.txt


parallel-fastq-dump -h

#Convert sra file to fastq file
parallel-fastq-dump --sra-id SRR11955372 --threads 100 --outdir out/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11955373 --threads 100 --outdir out/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11955374 --threads 100 --outdir out/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11955375 --threads 100 --outdir out/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11955376 --threads 100 --outdir out/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11955377 --threads 100 --outdir out/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11955378 --threads 100 --outdir out/ --split-files --gzip
parallel-fastq-dump --sra-id SRR11955379 --threads 100 --outdir out/ --split-files --gzip

#Or write a for loop, shell script
for i in SRR1*
do  
echo $i  
parallel-fastq-dump  --sra-id $i --threads 100 --outdir out/ --split-files --gzip
done


#Merge sequencing files for the same sample, pay attention to the output filename format
cat SRR11955372_1.fastq.gz SRR11955373_1.fastq.gz SRR11955374_1.fastq.gz SRR11955375_1.fastq.gz SRR11955376_1.fastq.gz SRR11955377_1.fastq.gz SRR11955378_1.fastq.gz SRR11955379_1.fastq.gz > BC22_S1_L001_R1_001.fastq.gz
cat SRR11955372_2.fastq.gz SRR11955373_2.fastq.gz SRR11955374_2.fastq.gz SRR11955375_2.fastq.gz SRR11955376_2.fastq.gz SRR11955377_2.fastq.gz SRR11955378_2.fastq.gz SRR11955379_2.fastq.gz > BC22_S1_L001_R2_001.fastq.gz

cat SRR11955375_1.fastq.gz SRR11955377_1.fastq.gz > BC22_S1_L001_R1_001.fastq.gz
cat SRR11955375_2.fastq.gz SRR11955377_2.fastq.gz > BC22_S1_L001_R2_001.fastq.gz


cellranger -h

#Human reference file
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" --no-check-certificate
#Extract
tar -zxvf refdata-gex-GRCh3D-2024-A.tar.gz

#Mouse reference file
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz" --no-check-certificate
#Extract
tar -zxvf refdata-gex-GRCm39-2024-A.tar.gz

#Run cellranger
cellranger count --id=BC22 --sample=BC22 --transcriptome=/home/20221028/1YJ/11cellranger/refdata-gex-GRCh38-2024-A --fastqs=/home/20221028/1YJ/11cellranger/data/out --nosecondary --create-bam=true

cellranger count --id=BC22 --sample=BC22 --transcriptome=/home/20221028/1YJ/11cellranger/refdata-gex-GRCh38-2024-A --fastqs=/home/20221028/1YJ/11cellranger/data2/out --nosecondary --create-bam=true