srun -N 1 --cpus-per-task=10 --time=5:00:00 --partition=compute --pty bash

mkdir GATK_tutorial
cd GATK_tutorial

# download data
wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
  tar xvzf samples.tar.gz
  rm samples.tar.gz


mkdir -p dbSNP
mkdir -p SAM
mkdir -p BAM
mkdir -p sortedBAM
mkdir -p SNPs
mkdir -p SNPs/{MERGED,persample}

GATK="${PWD}/"
dbSNP="${GATK}/dbSNP/"
SAM="${GATK}/SAM/"
BAM="${GATK}/BAM/"
sortedBAM="${GATK}/sortedBAM/"
SNPs="${GATK}/SNPs/"
MERGED="${SNPs}MERGED/"
outdir="${SNPs}persample/"

## module use /gpfs/shared/modulefiles_local/bio  to run these modules
module load samtools
module load bwa
module load R

# Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
#
# Rscript -e "install.packages('gplots', contriburl=contrib.url('http://cran.r-project.org/'))"
#
# Rscript -e "install.packages('reshape', contriburl=contrib.url('http://cran.r-project.org/'))"
#
# Rscript -e "install.packages('gsalib', contriburl=contrib.url('http://cran.r-project.org/'))"
#
# Rscript -e "install.packages('Biobase', contriburl=contrib.url('http://bioconductor.org/packages/release/bioc/'))"
#
# install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))
# install.packages('gplots', contriburl=contrib.url('http://cran.r-project.org/'))
# install.packages('reshape', contriburl=contrib.url('http://cran.r-project.org/'))
# install.packages('gsalib', contriburl=contrib.url('http://cran.r-project.org/'))
# install.packages('Biobase', contriburl=contrib.url('http://bioconductor.org/packages/release/bioc/'))

# mkdir -p $HOME/.rlibs
# export R_LIBS_USER=$HOME/.rlibs

# download reference genome
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-47/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz 
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
bwa index Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

###indexing the yeast Genome  skip for now
bwa mem -M -R '@RG\tID:ABC123.LANE3\tLB:LIB-12JaimeYeast\tPL: ILLUMINA\tSM:12'  ../genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz  12JaimeYeast-05272020-12_S5_L001_R1_001.fastq.gz 12JaimeYeast-05272020-12_S5_L001_R2_001.fastq.gz > 12JaimeYeast_aligned.sam


module load fastqc
mkdir fastqc

### Trim adaptor :  this is to remove the adaptor added by the library #######   cut adapt need to know the adaptors used
module load bio/cutadapt/2.0
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o 1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq -p 1JaimeYeast-05272020-1-sq_S4_L001_R2_001_trimmed.fastq 1JaimeYeast-05272020-1-sq_S4_L001_R1_001.fastq 1JaimeYeast-05272020-1-sq_S4_L001_R2_001.fastq
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o 1JaimeYeast-05272020-1-sq_S4_L002_R1_001_trimmed.fastq -p 1JaimeYeast-05272020-1-sq_S4_L002_R2_001_trimmed.fastq 1JaimeYeast-05272020-1-sq_S4_L002_R1_001.fastq 1JaimeYeast-05272020-1-sq_S4_L002_R2_001.fastq
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o 1JaimeYeast-05272020-1-sq_S4_L003_R1_001_trimmed.fastq -p 1JaimeYeast-05272020-1-sq_S4_L003_R2_001_trimmed.fastq 1JaimeYeast-05272020-1-sq_S4_L003_R1_001.fastq 1JaimeYeast-05272020-1-sq_S4_L003_R2_001.fastq
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o 1JaimeYeast-05272020-1-sq_S4_L004_R1_001_trimmed.fastq -p 1JaimeYeast-05272020-1-sq_S4_L004_R2_001_trimmed.fastq 1JaimeYeast-05272020-1-sq_S4_L004_R1_001.fastq 1JaimeYeast-05272020-1-sq_S4_L004_R2_001.fastq

## bbduk to automatically trim adaptors
bbduk in1=1JaimeYeast-05272020-1-sq_S4_L001_R1_001.fastq.gz in2=1JaimeYeast-05272020-1-sq_S4_L001_R2_001.fastq.gz out1=1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq.gz out2=1JaimeYeast-05272020-1-sq_S4_L001_R2_001_trimmed.fastq.gz ref=/gpfs/shared/apps_local/bbtools/38.79/resources/nextera.fa.gz k=23 mink=11 hdist=1 tbo tpe ktrim=l
#########RUN THIS: ADD read Group!!!!#############

for R1 in *_R1_001_trimmed.fastq;do
  SM=$(echo $R1 | cut -d"_" -f1)                                          ##sample ID
  LB=$(echo $R1 | cut -d"_" -f1,2)                                        ##library ID
  PL="Illumina"                                                           ##platform (e.g. illumina, solid)
  RGID=$(cat $R1 | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)       ##read group identifier
  PU=$RGID.$LB                                                            ##Platform Unit
  echo -e "SM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

  R2=$(echo $R1 | sed 's/_R1_/_R2_/')
  echo $R1 $R2
  bwa mem -t 40 -M -R '@RG\tID:'${RGID}'\tSM:'${SM}'\tPL:'${PL}'\tLB:'${LB}'\tPU:'${PU} Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa $R1 $R2 > ${SAM}${R1%_R1_001.fastq}.sam
done



# Generate BAM file
for samfile in ${SAM}*.sam;do
  sample=$(basename "${samfile%.*}")
  echo "Doing: " $sample
  samtools view -bS -o ${BAM}${sample}.bam $samfile
  echo "Created: " ${BAM}${sample}.bam
  samtools sort ${BAM}${sample}.bam -o ${sortedBAM}${sample}.sorted.bam
  echo "Created: " ${sortedBAM}${sample}.sorted.bam
done







wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar
chmod u+x picard.jar




wget https://de.cyverse.org/dl/d/6177B1E0-718A-4F95-A83B-C3B88E23C093/GenomeAnalysisTK-3.7-0.tar.bz2
tar xjf GenomeAnalysisTK-3.7-0.tar.bz2



java -Xmx10g -jar ${GATK}picard.jar CreateSequenceDictionary R=Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa O=Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.dict
 samtools faidx Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa


# Merge BAM replicates
#samtools view -H ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq.sorted.bam | grep -v "^@RG" | samtools reheader - ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq.sorted.bam > ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq.sorted.bam
#samtools view -H ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L002_R1_001_trimmed.fastq.sorted.bam | grep -v "^@RG" | samtools reheader - ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L002_R1_001_trimmed.fastq.sorted.bam > ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L002_R1_001_trimmed.fastq.sorted.bam
#samtools view -H ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L003_R1_001_trimmed.fastq.sorted.bam | grep -v "^@RG" | samtools reheader - ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L003_R1_001_trimmed.fastq.sorted.bam > ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L003_R1_001_trimmed.fastq.sorted.bam
#samtools view -H ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L004_R1_001_trimmed.fastq.sorted.bam | grep -v "^@RG" | samtools reheader - ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L004_R1_001_trimmed.fastq.sorted.bam > ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L004_R1_001_trimmed.fastq.sorted.bam

java  -jar ${GATK}picard.jar MergeSamFiles I="${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq.sorted.bam" I="${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L002_R1_001_trimmed.fastq.sorted.bam" I="${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L003_R1_001_trimmed.fastq.sorted.bam" I="${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L004_R1_001_trimmed.fastq.sorted.bam" OUTPUT="${sortedBAM}1_trimmed_merged.sorted.bam"




#Mark duplicates


for sample in ${sortedBAM}*.sorted.bam;do
  #name=${sample%.sorted.bam}
  #name=$(basename "${sample%.*}")
  name=$(basename "${sample%.sorted.bam}")
  echo "Doing: " $name
  java  -Xmx10g -jar ${GATK}picard.jar MarkDuplicates INPUT=$sample OUTPUT=${sortedBAM}${name}.dedup.bam METRICS_FILE=$name.metrics.txt;
done

# cd sortedBAM
#samtools view -H ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq.sorted.bam
#samtools view -H ${sortedBAM}BD143_TGACCA_L006.sorted.bam
#samtools view -H ${sortedBAM}BD143_TGACCA_merged.sorted.bam
#samtools view -H ${sortedBAM}BD143_TGACCA_merged.dedup.bam

mv ${sortedBAM}1JaimeYeast-05272020-1-sq_S4_L00*_R1_001_trimmed.fastq.sorted* .



http://m.ensembl.org/Saccharomyces_cerevisiae/Info/Annotation#assembly

# wget 'ftp://ftp.ensemblgenomes.org/pub/fungi/release-47/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz'
wget 'ftp://ftp.ensemblgenomes.org/pub/fungi/release-47/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz' -O ${dbSNP}saccharomyces_cerevisiae.vcf.gz

# mv Saccharomyces_cerevisiae.vcf.gz Saccharomyces_cerevisiae.vcf.gz

gunzip -c ${dbSNP}saccharomyces_cerevisiae.vcf.gz > ${dbSNP}saccharomyces_cerevisiae.vcf
## to extract chromosome 15 for Top1
grep "^#" ${dbSNP}saccharomyces_cerevisiae.vcf > ${dbSNP}saccharomyces_cerevisiae_chr15.vcf
grep "^15" ${dbSNP}saccharomyces_cerevisiae.vcf | sed 's/^15/chr15/' >> ${dbSNP}saccharomyces_cerevisiae_chr15.vcf
# Run Recalibration
# BQSR stands for Base Quality Score Recalibration.
samtools view -H ${sortedBAM}1_trimmed_merged.sorted.bam | grep -v "^@RG" | samtools reheader - ${sortedBAM}1_trimmed_merged.sorted.bam > ${sortedBAM}1_trimmed_merged.sorted.bam


for sample in ${sortedBAM}*.dedup.bam;do
  #name=${sample%.dedup.bam}
  name=$(basename "${sample%.dedup.bam}")
  echo "Doing: " $name
  samtools index $sample
  java -Xmx10g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -I $sample -knownSites ${dbSNP}saccharomyces_cerevisiae.vcf -o ${SNPs}${name}.1st.table
  java -Xmx10g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -I $sample -knownSites ${dbSNP}saccharomyces_cerevisiae.vcf -BQSR ${SNPs}${name}.1st.table -o ${SNPs}${name}.2nd.table
  java -Xmx10g -jar GenomeAnalysisTK.jar -T PrintReads -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -I $sample -BQSR ${SNPs}${name}.2nd.table -o ${SNPs}${name}.recal.bam
  java -Xmx10g -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -before ${SNPs}${name}.1st.table -after ${SNPs}${name}.2nd.table -plots ${SNPs}${name}.BQSR.pdf
done



# HERE
for sample in ${SNPs}*.recal.bam;do
  name=$(basename "${sample%.recal.bam}")
  java -Xmx10g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --dbsnp ${dbSNP}saccharomyces_cerevisiae.vcf -I $sample --emitRefConfidence GVCF -nct 3 -o ${outdir}${name}.g.vcf
done


#Combine call
java -Xmx10g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --dbsnp ${dbSNP}saccharomyces_cerevisiae.vcf  --variant ${outdir}1JaimeYeast-05272020-1-sq_S4_L001_R1_001_trimmed.fastq.g.vcf --variant ${outdir}1JaimeYeast-05272020-1-sq_S4_L002_R1_001_trimmed.fastq.g.vcf -o ${MERGED}raw_variants.vcf


 # split variants into SNPs and indels
mkdir ${MERGED}/SNPs
mkdir ${MERGED}/INDELs
java -Xmx10g -jar ${GATK}GenomeAnalysisTK.jar -T SelectVariants -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V ${MERGED}raw_variants.vcf -selectType SNP -o "${MERGED}/SNPs/"raw_SNP.vcf
java -Xmx10g -jar ${GATK}GenomeAnalysisTK.jar -T SelectVariants -R ${GATK}Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V ${MERGED}raw_variants.vcf -selectType INDEL -o "${MERGED}/INDELs/"raw_INDEL.vcf


# Distribution of variants
cd $MERGED
mkdir both
cp INDELs/* ./both/
cp SNPs/* ./both/
cd "$MERGED/both"
wget https://raw.githubusercontent.com/drtamermansour/angus/2017/densityCurves.R
for var in "SNP" "INDEL";do
 for ann in "QD" "MQRankSum" "FS" "SOR" "ReadPosRankSum";do
  annFile=$var.$ann; echo $annFile;
  awk -v k="$ann=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' raw_$var.vcf > $annFile
  grep -v "^\." $annFile > known.$annFile
  grep "^\." $annFile > novel.$annFile
  Rscript densityCurves.R "$annFile"
  rm $annFile known.$annFile novel.$annFile
done; done

# Apply filters
cd ../../.. # Move back to genome
java -Xmx10g -jar GenomeAnalysisTK.jar -T VariantFiltration -R Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V $MERGED/both/raw_SNP.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filterName "snp_filter" -o $MERGED/both/filtered_SNP.vcf

java -Xmx10g -jar GenomeAnalysisTK.jar -T VariantFiltration -R Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V $MERGED/both/raw_INDEL.vcf --filterExpression "QD < 2.0 || FS > 200.0" --filterName "indel_filter" -o $MERGED/both/filtered_INDEL.vcf


#NEXT>>>>>> R Programming
GWAS
Family-Based/Gene-based Study
Ethnicity Mapping
Cancer gene discovery/Clonal evolution
Personalized medicine
