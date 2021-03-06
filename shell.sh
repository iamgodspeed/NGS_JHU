########## section 1 ##########

data : apple.genome, apple.genes, apple.A, apple.B, apple.C
#1.how many chromosomes
grep ">" apple.genome | wc -l
grep -c ">" apple.genome

#2.how many genes and transcript variants
cut -f1 apple.genes | sort -u | wc -l # genes
cut -f2 apple.genes | sort -u | wc -l # transcript variants

#3.single splice variant
cut -f1 apple.genes | sort | uniq -c | grep " 1 " | wc -l

#4.multiple splice variants
cut -f1 apple.genes | sort | uniq -c | grep -v " 1 " | wc -l

#5.how many genes on "+" or "-" strand
cut -f1,4 apple.genes | sort -u | cut -f2 | sort | uniq -c

#6.how many genes on each chr
cut -f1,3 apple.genes | sort -u | cut -f2 | sort | uniq -c

#7.how many transcript variants on each chr
cut -f2,3 apple.genes | sort -u | cut -f2 | sort | uniq -c

#8.common between A and B
cut -f1 apple.A | sort -u > apple.A.genes
cut -f1 apple.B | sort -u > apple.B.genes
comm -1 -2 apple.A.genes apple.B.genes | wc -l

#9.specific to A, specific to B
comm -2 -3 apple.A.genes apple.B.genes | wc -l
comm -1 -3 apple.A.genes apple.B.genes | wc -l

#10.common among A B C
cat apple.genes.{A,B,C} | sort | uniq -c | grep " 3 " | wc -l

########## section 2 ##########

#1.how many alignment in bam
samtools flagstat x.bam
samtools view x.bam | wc -l
samtools view x.bam | cut -f3 | grep -v "*" | wc -l

#2.how many alignment shows mate unmapped (mapped)
samtools view x.bam | cut -f7 | grep "*" | wc -l ("=")

#3.how many alignments contain Deletion(spliced)
samtools view x.bam | cut -f6 | grep D | wc -l (N)

extract chromosome position:

samtools sort x.bam x.sorted (generate x.sorted.bam)
samtools index x.sorted.bam (generate x.sorted.bam.bai)
samtools view –b x.sorted.bam “Chr3:1000-2000” > x.subset.bam  (* -b is very important *)

#11-15
#11.how many sequences
samtools view -H x.bam | grep "SN:" | wc -l

12.length of the first sequence
samtools view -H x.bam | grep "SN:" | more

13.aligner was used
samtools view –H x.bam | grep “^@PG” (^start of the line)

14.id of the first sequence
samtools view x.bam | head -1 | cut -f1

15.start position of mate
samtools view x.bam | head -1 | cut -f7,8

16-20
16.how many overlaps
bedtools intersect -abam x.bam -b x_annot.gtf -bed -wo > overlaps.bed
-abam
-bed

17.how many these overlaps are longer than 10 bases
cut –f22 overlaps.bed | sort –nrk1 > lengths

18.how many alignments overlap annotation
cut -f1-12 overlaps.bed | sort -u | wc -l
The minimum information needed to define the alignments is contained in columns 1-5
which include the read ID and the flag
specifying whether this is read 1 or read 2 in a pair with the same read ID).

19.how many exons overlap alignments
cut -f13-21 overlaps.bed | sort -u | wc -l

20.how many transcripts in annotation.gtf (in BED file, each line = 1 transcript)
cut -f9 x_annot.gtf | cut -d ' ' -f1,2 |  sort -u

########## section 3 ##########

reference : x.fas
sequence : x.fastq

# 1.how many sequences in genome (Ref)
grep "^>" x.fas | wc -l

# 2.name of the third sequence
grep “^>” x.fas | head -3 | tail -1

# 4.how many index bowtie2 creates
bowtie2-build x.fas ./index/index

# 5.extension of index files
bt2

# 6.how many reads in fastq file
wc -l x.fastq

# 7.full match vs local match
bowtie2 –x ./index –U x.fastq –S out.full.sam
bowtie2 –x ./index –U x.fastq –S out.local.sam --local

# filter mapped only
cat out.full.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l
cat out.local.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l

# 8.reads mapped
a read may have 0 (unmapped), 1 (unique), or multiple alignments in a SAM file
default bowtie2 reports only one match per read
= alignment

# 10.how many alignments contain in dels
cut –f6 out.full.sam | grep –c “[I,D]”
cut –f6 out.local.sam | grep –c “[I,D]”

# 11.how many "chr3" in vcf
samtools view –bT x.fas out.full.sam > out.full.bam
samtools sort out.full.bam out.full.sorted
samtools mpileup –f x.fas –uv out.full.sorted.bam > out.full.mpileup.vcf
cat out.full.mpileup.vcf | grep –v “^#” | cut –f1 | grep –c “^Chr3”

# 12.
cat out.full.mpileup.vcf | grep –v “^#” | cut –f4 | grep –p “^A$”

# 13. supporting reads
cat out.full.mpileup.vcf | grep –v “^#” | grep –c “DP=20;”

# 14.
cat out.full.mpileup.vcf | grep –v “^#” | grep –c INDEL

# 15.
cat out.full.mpileup.vcf | grep –v “^#” | cut –f1,2 | grep Chr1 | grep 175672

# 16-20 BCF
# 16.
samtools mpileup –f x.fas –g out.full.sorted.bam > out.full.mpileup.bcf
bcftools call –m –v –O v out.full.mpileup.bcf > out.final.vcf
cat out.final.vcf | grep –v “^#” |  cut –f1 | sort | uniq –c | grep “Chr3”

# 17. A->T
cat out.final.vcf | grep –v “^#” | cut –f4,5 | grep –p “^A\tT$” | wc -l

# 18.
cat out.final.vcf | grep –v “^ #” | grep –c INDEL

# 19.
cat out.final.vcf | grep –v “^#” | grep –c “DP=20;”

# 20.
cat out.final.vcf | grep –v “^#” | cut -f1-5 | grep Chr3 | grep 11937923

########## section 4 ##########

data: xy.fa(reference) xy.gtf(annotation) x.fastq y.fastq
# 1-5 mapping with tophat
mkdir Tophat
mkdir Tophat/X
mkdir Tophat/Y

mkdir index
bowtie2-build xy.fa index/index
cp xy.fa index/

tophat -o Tophat/X/ index/index x.fastq # splice mapping x
tophat -o Tophat/Y/ index/index y.fastq # splice mapping y

check align_summary.txt # check alignment in results

# 6-10 transcripts assemble with cufflinks
mkdir Cufflinks
mkdir Cufflinks/X
mkdir Cufflinks/Y

#(specified labels as prefix for naming the assembled transcripts)
cufflinks -o Cufflinks/X -L x Tophat/X/accepted_hits.bam
cufflinks -o Cufflinks/Y -L y Tophat/Y/accepted_hits.bam

#how many genes
cut -f9 Cufflinks/X/transcripts.gtf | cut -d ' ' -f2 | sort -u | wc -l
cut -f9 Cufflinks/Y/transcripts.gtf | cut -d ' ' -f2 | sort -u | wc -l

#how many transcripts
cut -f9 Cufflinks/X/transcripts.gtf | cut -d ' ' -f4 | sort -u | wc -l
cut -f9 Cufflinks/Y/transcripts.gtf | cut -d ' ' -f4 | sort -u | wc -l

#single transcript genes
cut -f9 Cufflinks/X/transcripts.gtf | cut -d ' ' -f2,4 | sort -u | cut -d ' ' -f1 | sort | uniq -c | grep -c "1"
cut -f9 Cufflinks/Y/transcripts.gtf | cut -d ' ' -f2,4 | sort -u | cut -d ' ' -f1 | sort | uniq -c | grep -c "1"

#multi exon transcripts
transcripts - single exon transcripts

# 11-15 compare with cuffcompare
cd Cufflinks/X
cuffcompare -r ../../xy.gtf -R transcripts.gtf

cd Cufflinks/Y
cuffcompare -r ../../xy.gtf -R transcripts.gtf

cut -f3 cuffcmp.transcripts.gtf.tmap | sort | uniq -c

# = :fully reconstructed
# c :partial
# j :novel
# i :intron

grep gene_name cuffcmp.transcripts.gtf.tmap

# 16-20 differential expression with cuffdiff
#create GTFs.txt
#Cufflinks/*/transcripts.gtf

cuffmerge -g xy.gtf GTF.txt

#genes
cut -f9 merged_asm/merged.gtf | cut -d ' ' -f2 | sort -u | wc -l

#transcripts
cut -f9 merged_asm/merged.gtf | cut -d ' ' -f4 | sort -u | wc -l

cuffdiff -o Cuffdiff merged_asm/merged.gtf Tophat/X/accepted_hits.bam Tophat/Y/accepted_hits.bam

wc -l gene_exp.diff

#genes differentlly expressed
grep -c yes gene_exp.diff

#transcripts differentlly expressed
grep -c yes isoform_exp.diff
