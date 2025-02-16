##Analysis
###1.Run Skeleton script from (Analysis script) /mnt/imperial_data/Imperial_bioinformatics_pipeline/data/results/NextSeq_Cvg(Platform)/Reads/181015_NB551088_0013_AHY7LFAFXX_Cvg(RunName)
##QC
###1.Run CoverageSummaryScript (QC script) from /mnt/imperial_data/Imperial_bioinformatics_pipeline/data/results/NextSeq_Cvg(Platform)/Reads
###2.Add Sample Sheet for CoverageSummaryScript (QC script) in  /mnt/imperial_data/Imperial_bioinformatics_pipeline/data/results/NextSeq_Cvg(Platform)/Reads/181015_NB551088_0013_AHY7LFAFXX_Cvg(RunName)/Run_51-Sample_Sheet_151018.csv(SampleSheet)


##These parameters will modify by Master script
RunDate=181015_NB551088_0013_AHY7LFAFXX_PandaQC_20AE07731_3
PoolName=ICCNexteraV4_169
Platform=NextSeq
mode=Live
nt=15
dcovg=1000
memory="-Xmx25g"
Skip_size_Check=1
TopDir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/results/$Platform
StoreDir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/Store/Scripts/Research

# # # Testing
# echo -e "TopDir: $TopDir"
# echo -e "StoreDir: $StoreDir"
###############################################################################
#Here I added two brackets instead of 1 of the original script
#Don't forget to change the paths
if [[ $mode == "Test" ]]
then
	TopDir=/data/Research/PipelineDevelopments/Test_data/$Platform
	StoreDir=/data/Research/PipelineDevelopments/Research
	echo -e "\nTest mode selected. \n"
fi
##############################################################################
###Target Files 
###Target file from Ensembl70 & refSeq -- Overhangs exon(+/-) 40bp into introns (169 genes)
TargetFile=$StoreDir/TargetFiles/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.mergeBed.bed
###ProteinCodingregion
CDSFile=$StoreDir/TargetFiles/ICC_169Genes_Nextera_V4_ProteinCodingExons.mergeBed.bed
###Canonical Transcript 
CanonTranFile=$StoreDir/TargetFiles/ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.mergeBed.bed


##############################################################################
mkdir -p $TopDir/$RunDate
mkdir -p $TopDir/$RunDate/$PoolName

SampleDir=$TopDir/$RunDate
runDir=$TopDir/Reads/$RunDate

# # # #Testing
# echo -e "SampleDir: $SampleDir"
# echo -e "runDir: $runDir"

## create Output dir
outdir="gatk_snp_indel"
mkdir -p $SampleDir/$PoolName/$outdir
MyOutDir=$TopDir/$RunDate/$PoolName/$outdir

# #Testing
# echo -e "MyOutDir: $MyOutDir"
###############################################################################

#The dict and fai files are also available in the same folder.
#Therefore, I'll keep the code that creates them commented.
Ref=/mnt/imperial_data/Imperial_bioinformatics_pipeline/Ref_files_to_Aswan/UCSChg19/allchrom.Chr1ToChrM.validated.fa
dbSNP=/mnt/imperial_data/Imperial_bioinformatics_pipeline/Ref_files_to_Aswan/dbSNP138/hg19_dbSNP138.vcf
omni=/mnt/imperial_data/Imperial_bioinformatics_pipeline/Ref_files_to_Aswan/GATK_bundle_2.8/1000G_omni2.5.hg19.vcf
TenK=/mnt/imperial_data/Imperial_bioinformatics_pipeline/Ref_files_to_Aswan/GATK_bundle_2.8/1000G_phase1.snps.high_confidence.hg19.vcf
hapmap=/mnt/imperial_data/Imperial_bioinformatics_pipeline/Ref_files_to_Aswan/GATK_bundle_2.8/hapmap_3.3.hg19.vcf
mills=/mnt/imperial_data/Imperial_bioinformatics_pipeline/Ref_files_to_Aswan/GATK_bundle_2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf
TenKIndel=/mnt/imperial_data/Imperial_bioinformatics_pipeline/Ref_files_to_Aswan/GATK_bundle_2.8/1000G_phase1.indels.hg19.vcf
snpeff=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/Install/snpEff/
group=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/Install/filo2/bin/
VEP=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/Install/ensembl-tools-release-83/scripts/variant_effect_predictor/
LOFTEE=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/Install/loftee-master/src
			  
PATH=$PATH ; export PATH

###########################################################################
###########################################################################
## FastQC 
###########################################################################
###########################################################################
	echo -e "FastQC started"
	$fastqc/fastqc $runDir/$Lib/$Read1
	$fastqc/fastqc $runDir/$Lib/$Read2
	echo -e "FastQC ended"	
#############################################################################
## Get Read name
	echo -e "Get Read name started"
	Read1aTR=`ls $runDir/$Lib | grep _"$Lane"_ | grep R1_001 | grep -v prinseq | grep -v fastqc | awk 'NR==1'`
	Read1bTR=`ls $runDir/$Lib | grep _"$Lane"_ | grep R2_001 | grep -v prinseq | grep -v fastqc | awk 'NR==1'`
	#Test
	ReadLenTR=`less $runDir/$Lib/$Read1aTR | head | awk 'NR==2' | wc -m`
	# #Testing	
	# echo -e "Read1aTR: $Read1aTR"
	# echo -e "Read1bTR: $Read1bTR"
	# echo -e "ReadLenTR: $ReadLenTR"
	echo -e "Get Read name ended"
#############################################################################	
# ## Unzip Reads
	echo -e "Unzip Reads started"
	gzip -d $runDir/$Lib/$Read1aTR
	gzip -d $runDir/$Lib/$Read1bTR
	echo -e "Unzip Reads ended"
##############################################################################
## Get Unzipped read name
	echo -e "Get Unzipped read name started"
	Read1aTRa=`ls $runDir/$Lib | grep _"$Lane"_ | grep R1_001 | grep -v prinseq | grep -v fastqc | awk 'NR==1'`
	Read1bTRa=`ls $runDir/$Lib | grep _"$Lane"_ | grep R2_001 | grep -v prinseq | grep -v fastqc | awk 'NR==1'`
	# echo -e "Read1aTRa: $Read1aTRa"
	# echo -e "Read1bTRa: $Read1bTRa"
	echo -e "Get Unzipped read ended"
##############################################################################
### Get Sample ID 
	echo -e "Get Sample ID started"
	SampleID=`ls $runDir/$Lib | grep _"$Lane"_ | grep R1_001 | grep -v prinseq | grep -v fastqc | awk 'NR==1'| awk 'BEGIN {FS="_";} {print $1}'`
	Barcode=`ls $runDir/$Lib | grep _"$Lane"_ | grep R1_001 | grep -v prinseq | grep -v fastqc | awk 'NR==1' | awk 'BEGIN {FS="_";} {print $2}'`
	#LaneID=`ls $runDir/$Lib | grep $Read1aTRa | grep -v prinseq | grep -v fastqc | awk 'BEGIN {FS="_";} {print $3}'`
	# echo -e "SampleID: $SampleID"
	# echo -e "Barcode: $Barcode"
	# echo -e "LaneID: $LaneID"
	echo -e "Get Sample ID ended"
##############################################################################	






###########################################################################
###########################################################################
# Trim low qual bases < 20 from both ends using Prinseq
###########################################################################
###########################################################################

## Trim low qual bases < 20 from both ends using #http://prinseq.sourceforge.net/manual.html
	echo -e "Trimming started"
	perl $prinseq/prinseq-lite.pl \
		-fastq $runDir/$Lib/$Read1aTRa \
		-fastq2 $runDir/$Lib/$Read1bTRa \
		-trim_qual_right 20 \
		-trim_qual_left 20 \
		-trim_qual_window 5 \
		-min_len 35 \
		-out_good $runDir/$Lib/"$SampleID"_"$Lane"_trimqual_prinseq_good_Read
	
	Step_Check $? $runDir/$Lib/"$SampleID"_"$Lane"_trimqual_prinseq_good_Read_1.fastq 100 	#min 100K [750 reads]
	echo -e "Trimming Step Check ended"

# ## Get read name after trimmed
	# echo -e "Get read name after trimmed started"
	Read1a=`ls $runDir/$Lib | grep _"$Lane"_ | grep Read_1 | grep -v fastqc | awk 'NR==1'`
	Read1b=`ls $runDir/$Lib | grep _"$Lane"_ | grep Read_2 | grep -v fastqc | awk 'NR==1'`
	# echo -e "Read1a: $Read1a"
	# echo -e "Read1b: $Read1b"
	echo -e "Get read name after trimmed ended"




###########################################################################
###########################################################################
# FastQC after trimmed
###########################################################################
###########################################################################

	$fastqc/fastqc $runDir/$Lib/$Read1a
	$fastqc/fastqc $runDir/$Lib/$Read1b


## Zip trimmed read
	gzip $runDir/$Lib/$Read1a
	gzip $runDir/$Lib/$Read1b


# Zip orginal read

	gzip $runDir/$Lib/$Read1aTRa
	gzip $runDir/$Lib/$Read1bTRa


# ##Get Read name after trim&zip
	# echo -e "Get Read name after trim&zip started"
	Read1=`ls $runDir/$Lib | grep _"$Lane"_ | grep Read_1 | grep -v fastqc | awk 'NR==1'`
	Read2=`ls $runDir/$Lib | grep _"$Lane"_ | grep Read_2 | grep -v fastqc | awk 'NR==1'`

## Create directory for sample
	mkdir $MyOutDir/"$SampleID"





###########################################################################
###########################################################################
## BWA alignment
###########################################################################
###########################################################################
## Get BWA version
	$bwa 2> /mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/bwa.txt
	bwaver=`cat /mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/bwa.txt | grep Version | awk '{print $2}'`
	# echo -e "bwaver : $bwaver"

	## -M	 Mark shorter split hits as secondary (for Picard compatibility). -t number of threads
	echo -e "BWA mapping started"
	
	echo -e "bwa mem started"
	$bwa mem -M -t $nt $Ref $runDir/$Lib/$Read1 $runDir/$Lib/$Read2 > $MyOutDir/"$SampleID"/$SampleID."$Lane".sam
	Step_Check $? $MyOutDir/$SampleID/$SampleID."$Lane".sam 250 	#min 250K
	echo -e "Step check bwa mem ended"

	echo -e "samtools view started"
	$samtools view -bS $MyOutDir/"$SampleID"/"$SampleID"."$Lane".sam > $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam
	Step_Check $? $MyOutDir/$SampleID/$SampleID."$Lane".bam 50	#min 50K
	echo -e "Step check samtools view ended"

	echo -e "samtools sort started"
	$samtools sort $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam $MyOutDir/"$SampleID"/"$SampleID"."$Lane"1
	echo -e "samtools sort ended"
	
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/"$SampleID"."$Lane"1.bam $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam
	# echo -e "mv ended"

	echo -e "samtools index started"
	$samtools index $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam
	Step_Check $? $MyOutDir/$SampleID/$SampleID."$Lane".bam.bai 100	#min 100K
	echo -e "Step check samtools index ended"
	
	date

	echo -e "BWA run completed"






###########################################################################
###########################################################################	
## Add Read Groups
###########################################################################
###########################################################################

	echo -e "Add RG tag started"
	$Java7 $memory -jar -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp $picard/AddOrReplaceReadGroups.jar \
	INPUT=$MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam \
	OUTPUT=$MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam.bam \
	SORT_ORDER=coordinate \
	RGID=$RunDate \
	RGLB="$ReadLen"PEBC \
	RGPL=ILLUMINA \
	RGPU=$LaneID \
	RGSM=$SampleID \
	RGCN=BRU-RBH \
	RGDS=BWA-MEM."$bwaver" \
	VALIDATION_STRINGENCY=SILENT
	
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam.bam 50
	echo -e "Step check Add RG tag ended"
	
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam.bam $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam
	# echo -e "mv ended"
	
	echo -e "samtools index started"
	$samtools index $MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam
	Step_Check $? $MyOutDir/$SampleID/$SampleID."$Lane".bam.bai 100
	echo -e "Step check samtools index started"
###########################################################################







###########################################################################
###########################################################################
## Mark duplicates
###########################################################################
###########################################################################
	echo -e "mark duplicates started"
	$Java7 -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/ \
	$memory -jar \
	$picard/MarkDuplicates.jar \
	INPUT=$MyOutDir/"$SampleID"/"$SampleID"."$Lane".bam \
	OUTPUT=$MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam \
	METRICS_FILE=$MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam.metrics.txt \
	VALIDATION_STRINGENCY=LENIENT
	
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam 50
	echo -e "Step check mark duplicates ended"
	
	echo -e "samtools index started"
	$samtools index $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam.bai 100
	echo -e "Step check samtools index ended"
	
	# echo -e "markDup_bam_file assignment started"
	markDup_bam_file=$MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam
	# echo -e "markDup_bam_file: $markDup_bam_file"
	# echo -e "markDup_bam_file assignment ended"






	
#############################################################################	
# ## Get GATK version
	# echo -e "GATK version started"
	GATK=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/Install/GenomeAnalysisTK/GenomeAnalysisTK.jar
	GATK_version=`$Java7 $memory -jar -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp $GATK -version`
	# echo -e "GATK_version : $GATK_version"
	# echo -e "GATK version ended"
################################################################################
	# # #DONT FORGET TO FIX THE EXTRAXTION OF THE GATK VERSION!!!!!!!!!!!!!!!!
	echo -e "\nStarting Analysis for GATK_version:$GATK_version : $Lib"
###############################################################################
#get default values of variables to reset later
	SampleDir_default=$SampleDir
	MyOutDir_default=$MyOutDir
	GATKs_default=$GATKs
	# # #Testing
	# echo -e "SampleDir_default: $SampleDir_default"
	# echo -e "MyOutDir_default: $MyOutDir_default"
	# echo -e "GATKs_default: $GATKs_default"
##############################################################################








###########################################################################
###########################################################################
## Target stuff
###########################################################################
###########################################################################

# Create directory by SampleName/Target/ and SampleName/Exon
	
	Target=Target
	Exon=ProteinCodingTarget
	CanonTran=CanonTranCodingTarget

	 mkdir $MyOutDir/"$SampleID"/"$Target"
	 mkdir $MyOutDir/"$SampleID"/"$Exon"
	 mkdir $MyOutDir/"$SampleID"/"$CanonTran"
##############################################################################	
# Create directory  for EnrichmentReport and make symbolic link of files
	
	Enrich_Report="Enrichment_Report"
	
	 mkdir -p $SampleDir/"$Enrich_Report"
	 mkdir $SampleDir/"$Enrich_Report"/"$SampleID"
	 mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$Target"
	 mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"
	 mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"

	Exon20x=ProteinCodingTarget_20x

	 mkdir $MyOutDir/"$SampleID"/"$Exon20x"
	 mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$Exon20x"

	Exon30x=ProteinCodingTarget_30x

	 mkdir $MyOutDir/"$SampleID"/"$Exon30x"
	 mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$Exon30x"

	CanonTran20x=CanonTranCodingTarget_20x
	
	 mkdir $MyOutDir/"$SampleID"/"$CanonTran20x"
	 mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran20x"

	CanonTran30x=CanonTranCodingTarget_30x

	 mkdir $MyOutDir/"$SampleID"/"$CanonTran30x"
	 mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran30x"






###########################################################################
###########################################################################
## Creating Intervals (Depricated)
###########################################################################
###########################################################################
		 echo -e "Creating Intervals started"
		$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-nt $nt \
		-I $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam \
		-R $Ref \
		-L $TargetFile \
		--out $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam.intervals \
		-known $mills \
		-known $TenKIndel
			
		Step_Check $? $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam.intervals 10
		echo -e "Step_Check Creating Intervals ended"






###########################################################################
###########################################################################
## Indel Realigning (Depricated)
###########################################################################
###########################################################################

	# ##Indel Realigning -model USE_READS	
	## If GATK_2_3, use -rf NotPrimaryAlignment [Ignore all non-primary alignments]
		
		echo -e "Indel Realigning started"
		$Java7 -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/ $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-I $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam \
		-R $Ref \
		-T IndelRealigner \
		-targetIntervals $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.bam.intervals \
		--out $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam \
		-known $mills \
		-known $TenKIndel \
		-model USE_READS
		
		Step_Check $? $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam 50
		echo -e "Step_Check Indel Realigning ended"

	echo -e "\nDone Indel Realigning ($GATK_version) : $Lib"

	echo -e "samtools index started"
	$samtools index $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam
	echo -e "samtools index ended"









###########################################################################
###########################################################################
## BaseRecalibrator
###########################################################################
###########################################################################

##CountCovariates changed to BaseRecalibrator in v2  (For 454  include  -cov HomopolymerCovariate; For SOLiD -cov PrimerRoundCovariate )
	
	echo -e "BaseRecalibrator started"
	$Java7 -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/ $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-I $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam \
		--knownSites $dbSNP \
		--knownSites $mills \
		--knownSites $TenKIndel \
		-L $TargetFile \
		-T BaseRecalibrator \
		-nct $nt \
		-cov ReadGroupCovariate \
		-cov QualityScoreCovariate \
		-cov CycleCovariate \
		-cov ContextCovariate \
		--out $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam.recalibTable.grp
	
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam.recalibTable.grp 100
	echo -e "Step_Check BaseRecalibrator ended"

	echo -e "\nDone : CountCovariates ($GATK_version) : $Lib"



###########################################################################
###########################################################################
## PrintReads (changed to ApplyBQSR???)
	
	echo -e "TableRecalibration started"
	$Java7 -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/ $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-I $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam \
		-T PrintReads \
		-nct $nt \
		--out $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.recalibrated.bam \
		-BQSR $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.bam.recalibTable.grp \
		-S LENIENT
	
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.recalibrated.bam 50
	echo -e "Step_Check TableRecalibration ended"
	
	# echo -e "Cp started"
	cp $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.recalibrated.bai $MyOutDir/"$SampleID"/"$SampleID"."$Lane".markDup.Realigned.recalibrated.bam.bai
	# echo -e "Cp ended"
# # ###############################################################################	
	done

###########################################################################
###########################################################################
## Merging BAM files by Lanes
###########################################################################
###########################################################################

	## Merging BAM files by Lanes
	NumLanes=`ls $runDir/$Lib | grep fastq | grep -v prinseq | grep -v fastqc  | awk 'BEGIN {FS="_";} {print $3}' |sort -u |wc -l | awk 'NR==1'| awk '{print $1}'`
	echo -e "\n Number of Lanes = $NumLanes"
	
	if [ $NumLanes -eq 4 ]
	then
		echo -e "Number of lanes = $NumLanes"
		Lane1=`ls $MyOutDir/"$SampleID"/ | grep 'markDup.Realigned.recalibrated.bam$' | awk 'NR==1'`
		Lane2=`ls $MyOutDir/"$SampleID"/ | grep 'markDup.Realigned.recalibrated.bam$' | awk 'NR==2'`
		Lane3=`ls $MyOutDir/"$SampleID"/ | grep 'markDup.Realigned.recalibrated.bam$' | awk 'NR==3'`
		Lane4=`ls $MyOutDir/"$SampleID"/ | grep 'markDup.Realigned.recalibrated.bam$' | awk 'NR==4'`
		
		### Merge BAM files
		$Java7 $memory -jar -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/ $picard/MergeSamFiles.jar \
		INPUT=$MyOutDir/"$SampleID"/$Lane1 \
		INPUT=$MyOutDir/"$SampleID"/$Lane2 \
		INPUT=$MyOutDir/"$SampleID"/$Lane3 \
		INPUT=$MyOutDir/"$SampleID"/$Lane4 \
		OUTPUT=$MyOutDir/"$SampleID"/"$SampleID".bam \
		CREATE_INDEX=TRUE \
		VALIDATION_STRINGENCY=SILENT
		
		mv $MyOutDir/"$SampleID"/"$SampleID".bai $MyOutDir/"$SampleID"/"$SampleID".bam.bai







###########################################################################
###########################################################################
## Mark duplicates for merged BAM
###########################################################################
###########################################################################
	echo -e "Mark duplicates on merged started"
	$Java7 -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/ \
	$memory -jar \
	$picard/MarkDuplicates.jar \
	INPUT=$MyOutDir/"$SampleID"/"$SampleID".bam \
	OUTPUT=$MyOutDir/"$SampleID"/"$SampleID".markDup.bam \
	METRICS_FILE=$MyOutDir/"$SampleID"/"$SampleID".markDup.bam.metrics.txt \
	VALIDATION_STRINGENCY=LENIENT
	
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID".markDup.bam 50
	echo -e "Step_Check Mark duplicates on merged ended"
	
	echo -e "samtools index on merged started"
	$samtools index $MyOutDir/"$SampleID"/"$SampleID".markDup.bam
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.bai 100
	echo -e "Step_Check samtools index on merged ended"
	
	# echo -e "markDup_bam_file assignment on merged started"
	markDup_bam_file=$MyOutDir/"$SampleID"/"$SampleID".markDup.bam
	# echo -e "markDup_bam_file : $markDup_bam_file"
	# echo -e "markDup_bam_file assignment on merged ended"
	
echo -e "\nStarting Analysis for $GATK_version : $Lib"

###########################################################################
###########################################################################
## Creating Intervals on merged BAM (Depricated)
###########################################################################
###########################################################################
	echo -e "Creating Intervals on merged started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-nt $nt \
		-I $MyOutDir/"$SampleID"/"$SampleID".markDup.bam \
		-R $Ref \
		-L $TargetFile \
		--out $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.intervals \
		-known $mills \
		-known $TenKIndel
	
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.intervals 10
	echo -e "Step_Check Creating Intervals on merged ended"
###########################################################################
###########################################################################
## Indel Realigning on merged BAM (Depricated)
###########################################################################
###########################################################################
##Indel Realigning  -model USE_READS
	echo -e "Indel Realigning on merged started"	
	$Java7 -Djava.io.tmpdir=/mnt/imperial_data/Imperial_bioinformatics_pipeline/data/tmp/ $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-I $MyOutDir/"$SampleID"/"$SampleID".markDup.bam \
		-R $Ref \
		-T IndelRealigner \
		-targetIntervals $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.intervals \
		--out $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam \
		-known $mills \
		-known $TenKIndel \
		-model USE_READS
	
	Step_Check $? $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam 50
	echo -e "Step_Check Indel Realigning on merged ended"
	echo -e "\nDone Indel Realigning ($GATK_version) : $Lib"
	
	# echo -e "samtools index on merged started"
	 $samtools index $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam
	# echo -e "samtools index on merged ended"






###########################################################################
###########################################################################
## Filter BAM file for MAP qual > 8 for OnTarget, OnProtCodingTarget, OnCanonTranFile
###########################################################################
###########################################################################

TargetFile=$StoreDir/TargetFiles/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.mergeBed.bed
CDSFile=$StoreDir/TargetFiles/ICC_169Genes_Nextera_V4_ProteinCodingExons.mergeBed.bed
CanonTranFile=$StoreDir/TargetFiles/ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.mergeBed.bed

	# Define orgFile, OnTargetFile, OnProtCodingTargetFile, OnCanonTranFile
	OrgFile="$SampleID".markDup.Realigned.recalibrated.bam
	OnTargetFile="$SampleID".markDup.Realigned.recalibrated.OnTarget.q8.bam
	OnProtCodingTargetFile="$SampleID".markDup.Realigned.recalibrated.ProteinCoding.OnTarget.q8.bam
	OnCanonTranFile="$SampleID".markDup.Realigned.recalibrated.CanonTranCoding.OnTarget.q8.bam
	##Testing
	# echo -e "OrgFile: $OrgFile"
	# echo -e "OnTargetFile: $OnTargetFile"
	# echo -e "OnProtCodingTargetFile: $OnProtCodingTargetFile"
	# echo -e "OnCanonTranFile: $OnCanonTranFile"
	
	## Generate read on $Target MAP qual > 8
	echo -e "samtools view OnTargetFile started"
	$samtools view -uq 8 $MyOutDir/"$SampleID"/"$OrgFile" | $bedtools/intersectBed -abam stdin -b $TargetFile -u > $MyOutDir/"$SampleID"/$Target/"$OnTargetFile"
	Step_Check $? $MyOutDir/"$SampleID"/$Target/"$OnTargetFile" 10
	echo -e "Step_Check samtools view OnTargetFile ended"
	echo -e "samtools index OnTargetFile started"
	$samtools index $MyOutDir/"$SampleID"/$Target/"$OnTargetFile"
	echo -e "samtools index OnTargetFile ended"
	
	# # Generate read on Protein coding target  MAP qual > 8
	echo -e "samtools view OnProtCodingTargetFile started"
	$samtools view -uq 8 $MyOutDir/"$SampleID"/"$OrgFile" | $bedtools/intersectBed -abam stdin -b $CDSFile -u > $MyOutDir/"$SampleID"/$Exon/"$OnProtCodingTargetFile"
	Step_Check $? $MyOutDir/"$SampleID"/$Exon/"$OnProtCodingTargetFile" 5
	echo -e "Step_Check samtools view OnProtCodingTargetFile ended"
	echo -e "samtools index OnProtCodingTargetFile started"	
	$samtools index $MyOutDir/"$SampleID"/$Exon/"$OnProtCodingTargetFile"
	echo -e "samtools index OnProtCodingTargetFile ended"
	
	# ## Generate read on Canonical Transcript coding target  MAP qual > 8
	echo -e "samtools view OnCanonTranFile started"
	$samtools view -uq 8 $MyOutDir/"$SampleID"/"$OrgFile" | $bedtools/intersectBed -abam stdin -b $CanonTranFile -u > $MyOutDir/"$SampleID"/$CanonTran/"$OnCanonTranFile"
	Step_Check $? $MyOutDir/"$SampleID"/$CanonTran/"$OnCanonTranFile" 5
	echo -e "Step_Check samtools view OnCanonTranFile ended"
	echo -e "samtools index OnCanonTranFile started"
	$samtools index $MyOutDir/"$SampleID"/$CanonTran/"$OnCanonTranFile"
	echo -e "samtools index OnCanonTranFile ended"







###########################################################################
###########################################################################
## Flagstat reports
###########################################################################
###########################################################################

##Generating flagstat reports 
##Take out chimeric reads to get exact number of reads in flagstat -F 0x100

	##Target File
	echo -e "samtools view -F OnTargetFile started"
	$samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$Target/$OnTargetFile | $samtools flagstat - > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.flagstat
	echo -e "samtools view -F OnTargetFile ended"
	##ProteinCodingTarget file
	echo -e "samtools view -F OnProtCodingTargetFile started"
	$samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile | $samtools flagstat - > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.flagstat
	echo -e "samtools view -F OnProtCodingTargetFile ended"
	##Canonical Transcript coding file
	echo -e "samtools view -F OnCanonTranFile started"
	$samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile | $samtools flagstat - > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.flagstat
	echo -e "samtools view -F OnCanonTranFile ended"

	echo -e "samtools view -F OrgFile started"
	$samtools view -F 0x100 -uq 8 $MyOutDir/"$SampleID"/$OrgFile | $samtools flagstat - > $MyOutDir/"$SampleID"/$OrgFile.q8.flagstat
	echo -e "samtools view -F OrgFile ended"
	echo -e "samtools view -F OrgFile started"
	$samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$OrgFile | $samtools flagstat - > $MyOutDir/"$SampleID"/$OrgFile.flagstat
	echo -e "samtools view -F OrgFile ended"
	






###########################################################################
###########################################################################
## Coverage Stats and sort for OnProtCodingTarget
###########################################################################
###########################################################################

	## On ProteinCoding Target
	## CoverageStas [ duplicateReads were removed to calculate stats]
	echo -e "CoverageStats OnProtCodingTargetFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile | $bedtools/coverageBed -abam stdin -b $CDSFile  > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats
	echo -e "CoverageStats OnProtCodingTargetFile ended"
	## CoveragePerBase [ duplicateReads were removed to calculate stats]
	echo -e "CoveragePerBase OnProtCodingTargetFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile | $bedtools/coverageBed -abam stdin -b $CDSFile -d > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoveragePerBase
	echo -e "CoveragePerBase OnProtCodingTargetFile ended"
	## CoverageHistogram [ duplicateReads were removed to calculate stats]
	echo -e "CoverageHistogram OnProtCodingTargetFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile | $bedtools/coverageBed -abam stdin -b $CDSFile -hist > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageHist
	echo -e "CoverageHistogram OnProtCodingTargetFile ended"
	
	# sort the CoverageBed, CoverageStats output by chr,position
	echo -e "sort the CoverageBed, CoverageStats output by chr,position started"
	sort -k 1,1 -k 2,2n $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats.sorted
	echo -e "sort the CoverageBed, CoverageStats output by chr,position ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats.sorted $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats
	# echo -e "mv ended"
	echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoverageStats completed ($GATK_version) : $Lib"
	
	## sort the CoverageBed,CoveragePerBase output by chr,position,coverage (Change colum 5 into 4 in Exome $Target)
	echo -e "sort the CoverageBed,CoveragePerBase output by chr,position,coverage started"
	sort -k 1,1 -k 2,2n -k 5,5n $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoveragePerBase  > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoveragePerBase.sorted
	echo -e " sort the CoverageBed,CoveragePerBase output by chr,position,coverage ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoveragePerBase.sorted $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoveragePerBase
	# echo -e "mv ended"
	echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoveragePerBase completed ($GATK_version) : $Lib"






###########################################################################
###########################################################################
## Callable Loci
###########################################################################
###########################################################################

	## No of callable Bases by Exon --minDepth 10
	echo -e "No of callable Bases by Exon --minDepth 10 started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T CallableLoci \
		-I $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile \
		-o $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable \
		--minMappingQuality 10 \
		--minBaseQuality 20 \
		--minDepth 10 \
		-l INFO -format BED \
		-L $CDSFile \
		-summary $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.summary
	echo -e "No of callable Bases by Exon --minDepth 10 ended"	
	echo -e " \n Done : CallableLoci  for ($GATK_version) : $Lib"


###########################################################################
###########################################################################
## Depth of coverage
###########################################################################
###########################################################################
## Depth of coverage in $Target
	echo -e "Depth of coverage in $Target started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T DepthOfCoverage \
		-I $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile \
		-o $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.DepthOfCvg.txt  \
		-L $CDSFile \
		--minMappingQuality 10 \
		--minBaseQuality 20 \
		--stop 2000 \
		--includeDeletions --summaryCoverageThreshold \
		--outputFormat table --summaryCoverageThreshold --omitDepthOutputAtEachBase
	echo -e "Depth of coverage in $Target ended"
	echo -e " \n Done : DepthOfCoverage for  : $Lib"



###########################################################################
###########################################################################
## Mean coverage per Exon
###########################################################################
###########################################################################

	# echo -e "cat started"
	cat $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.DepthOfCvg.txt.sample_interval_summary | sed 's/chr16:997401/chr16:997401-997401/' | awk 'BEGIN {FS=":"} {OFS="\t"} {print $1,$2}'| awk 'BEGIN {FS="-"} {OFS="\t"} {print $1,$2}' |  awk 'BEGIN {FS=" "} {OFS="\t"} {print $1,$2-1,$3,$4,$5,$9}' | sed '1d'| $bedtools/sortBed -i stdin > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.temp.depth1.bed
	# echo -e "cat ended"
	# echo -e "cut started"
	cut -f4 $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats | awk 'BEGIN {FS="|";} {print $1}' > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.temp.depth2.bed 
	# echo -e "cut ended"
	# echo -e "paste started"
	paste <(cut -f1-3 $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.temp.depth1.bed) <(cut -f4 $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.temp.depth2.bed ) <(cut -f4- $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.temp.depth1.bed) | cut -f1-4,6 > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed 
	# echo -e "paste ended"
	
	# echo -e "rm started"
	rm -rf $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.temp*
	# echo -e "rm ended"
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	# echo -e "ln ended"






###############################################################################	
	## CollectAlignmentSummaryMetrics in OnTarget.bam
	echo -e "CollectAlignmentSummaryMetrics in OnTarget.bam started"
	$Java7 $memory -jar $picard/CollectAlignmentSummaryMetrics.jar \
		INPUT=$MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile \
		OUTPUT=$MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.AlignSumMet.txt \
		REFERENCE_SEQUENCE=$Ref \
		ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT
	echo -e "CollectAlignmentSummaryMetrics in OnTarget.bam ended"
	echo -e " \n Done : CollectAlignmentSummaryMetrics for  ($GATK_version) : $Lib"





###############################################################################
	## Generate Bases callable by target file
	# echo -e "awk started"
	awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable | $bedtools/intersectBed -a $CDSFile -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | $group/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.byTarget
	# echo -e "awk ended"




# ##############################################################################
	### Perl script to generate this format 
	#Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
	#chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
	# echo -e "perl script started"
	perl $StoreDir/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.byTarget > $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.byTarget.temp
	# echo -e "perl script ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.byTarget
	# echo -e "mv ended"




# ###############################################################################
# # Generate  soft links in $TopDir/$RunDate/$Enrich_Report folder
	# echo -e "Generate  soft links in Enrich_Report folder started"
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$OrgFile.q8.flagstat $SampleDir/"$Enrich_Report"/$SampleID/
	ln -s $MyOutDir/"$SampleID"/$OrgFile.flagstat $SampleDir/"$Enrich_Report"/$SampleID/

	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.flagstat $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoveragePerBase $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageHist $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.bases.callable.summary $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.AlignSumMet.txt $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.DepthOfCvg.txt.sample_interval_summary $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
	# echo -e "Generate  soft links in Enrich_Report folder ended"
##############################################################################
##  No of callable Bases by Exon --minDepth 30
	echo -e "No of callable Bases by Exon --minDepth 30 started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
			-R $Ref \
			-T CallableLoci \
			-I $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile \
			-o $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable \
			--minMappingQuality 10 \
			--minBaseQuality 20 \
			--minDepth 30 \
			-l INFO -format BED \
			-L $CDSFile \
			-summary $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable.summary
	echo -e "No of callable Bases by Exon --minDepth 30 ended"
	echo -e " \n Done : CallableLoci  for ($GATK_version) : $Lib"
# ##############################################################################
	# echo -e "awk started"
	awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable | $bedtools/intersectBed -a $CDSFile -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | $group/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable.byTarget
	# echo -e "awk ended"
# ##############################################################################
	### Perl script to generate this format 
	#Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
	#chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
	# echo -e "perl script started"
	perl $StoreDir/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable.byTarget > $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable.byTarget.temp
	# echo -e "perl script ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable.byTarget
	# echo -e "mv ended"
# ###############################################################################	
	# # Generate  soft links in $TopDir/$RunDate/$Enrich_Report folder
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$Exon30x/$OnProtCodingTargetFile.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$Exon30x"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$Exon30x"/
	# echo -e "ln ended"
# ##############################################################################
# ##  No of callable Bases by Exon --minDepth 20
	echo -e "No of callable Bases by Exon --minDepth 20 started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
			-R $Ref \
			-T CallableLoci \
			-I $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile \
			-o $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable \
			--minMappingQuality 10 \
			--minBaseQuality 20 \
			--minDepth 20 \
			-l INFO -format BED \
			-L $CDSFile \
			-summary $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable.summary
	echo -e "No of callable Bases by Exon --minDepth 20 ended"
	echo -e " \n Done : CallableLoci  for ($GATK_version) : $Lib"
# ##############################################################################
	# echo -e "awk started"
	awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable | $bedtools/intersectBed -a $CDSFile -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | $group/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable.byTarget
	# echo -e "awk ended"
# ##############################################################################
	### Perl script to generate this format 
	#Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
	#chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
	# echo -e "perl started"
	perl $StoreDir/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable.byTarget > $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable.byTarget.temp
	# echo -e "perl ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable.byTarget
	# echo -e "mv ended"
# ###############################################################################
# # Generate  soft links in $TopDir/$RunDate/$Enrich_Report folder
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$Exon20x/$OnProtCodingTargetFile.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$Exon20x"/
	ln -s $MyOutDir/"$SampleID"/$Exon/$OnProtCodingTargetFile.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$Exon20x"/
	# echo -e "ln ended"
###############################################################################
# ##End of ProteinCoding Target
# ##############################################################################
# ##OnCanonical Transcript File
	### CoverageStas [ duplicateReads were removed to calculate stats]
	echo -e "CoverageStats OnCanonTranFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile | $bedtools/coverageBed -abam stdin -b $CanonTranFile  > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats
	echo -e "CoverageStats OnCanonTranFile ended"
	### CoveragePerBase [ duplicateReads were removed to calculate stats]
	echo -e "CoveragePerBase OnCanonTranFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile | $bedtools/coverageBed -abam stdin -b $CanonTranFile -d > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoveragePerBase
	echo -e "CoveragePerBase OnCanonTranFile ended"
	### CoverageHistogram [ duplicateReads were removed to calculate stats]
	echo -e "CoverageHistogram OnCanonTranFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile | $bedtools/coverageBed -abam stdin -b $CanonTranFile -hist > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageHist
	echo -e "CoverageHistogram OnCanonTranFile ended"
# # ##############################################################################
	# ## sort the CoverageBed, CoverageStats output by chr,position
	echo -e "sort the CoverageBed, CoverageStats output by chr,position started"
	sort -k 1,1 -k 2,2n $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats.sorted
	echo -e "sort the CoverageBed, CoverageStats output by chr,position ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats.sorted $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats
	# echo -e "mv ended"
	echo -e "\nDone : Sorting of Canonical Transcript coding OnTarget (Targets).CoverageStats completed ($GATK_version) : $Lib"
# # ##############################################################################
	# ### sort the CoverageBed,CoveragePerBase output by chr,position,coverage (Change colum 5 into 4 in Exome $Target)
	echo -e "sort the CoverageBed,CoveragePerBase output by chr,position,coverage (Change colum 5 into 4 in Exome $Target) started"
	sort -k 1,1 -k 2,2n -k 5,5n $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoveragePerBase  > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoveragePerBase.sorted
	echo -e "sort the CoverageBed,CoveragePerBase output by chr,position,coverage (Change colum 5 into 4 in Exome $Target) ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoveragePerBase.sorted $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoveragePerBase
	# echo -e "mv ended"
	echo -e "\nDone : Sorting of Canonical Transcript coding OnTarget (Targets).CoveragePerBase completed ($GATK_version) : $Lib"
# # ###############################################################################
# ##  No of callable Bases by Exon --minDepth 10
	echo -e "No of callable Bases by Exon --minDepth 10 started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T CallableLoci \
		-I $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile \
		-o $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable \
		--minMappingQuality 10 \
		--minBaseQuality 20 \
		--minDepth 10 \
		-l INFO -format BED \
		-L $CanonTranFile \
		-summary $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.summary
	echo -e "No of callable Bases by Exon --minDepth 10 ended"
	echo -e " \n Done : CallableLoci  for ($GATK_version) : $Lib"
# ##############################################################################
# ### Depth of coverage in $Target
	echo -e "Depth of coverage in $Target started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T DepthOfCoverage \
		-I $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile \
		-o $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.DepthOfCvg.txt  \
		-L $CanonTranFile \
		--minMappingQuality 10 \
		--minBaseQuality 20 \
		--stop 2000 \
		--includeDeletions --summaryCoverageThreshold \
		--outputFormat table --summaryCoverageThreshold --omitDepthOutputAtEachBase
	echo -e "Depth of coverage in $Target ended"
	echo -e " \n Done : DepthOfCoverage for  ($GATK_version) : $Lib"
# # ##############################################################################
	# ## Getting mean coverage per Exon
	# echo -e "cat started"
	cat $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.DepthOfCvg.txt.sample_interval_summary | sed 's/chr8:19822821/chr8:19822820-19822821/' | awk 'BEGIN {FS=":"} {OFS="\t"} {print $1,$2}'| awk 'BEGIN {FS="-"} {OFS="\t"} {print $1,$2}' |  awk 'BEGIN {FS=" "} {OFS="\t"} {print $1,$2-1,$3,$4,$5,$9}' | sed '1d'| $bedtools/sortBed -i stdin > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.temp.depth1.bed
	# echo -e "cat ended"
	
	# echo -e "cut started"
	cut -f4 $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats | awk 'BEGIN {FS="|";} {print $1}' > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.temp.depth2.bed 
	# echo -e "cut ended"
	
	# echo -e "paste started"
	paste <(cut -f1-3 $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.temp.depth1.bed) <(cut -f4 $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.temp.depth2.bed ) <(cut -f4- $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.temp.depth1.bed) | cut -f1-4,6 > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed 
	# echo -e "paste ended"
	
	# echo -e "rm started"
	rm -rf $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.temp*
	# echo -e "rm ended"
	
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	# echo -e "ln ended"

# # # ##############################################################################	
# ### CollectAlignmentSummaryMetrics in OnTarget.bam
	echo -e "CollectAlignmentSummaryMetrics in OnTarget.bam started"
	$Java7 $memory -jar $picard/CollectAlignmentSummaryMetrics.jar \
		INPUT=$MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile \
		OUTPUT=$MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.AlignSumMet.txt \
		REFERENCE_SEQUENCE=$Ref \
		ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT
	echo -e "CollectAlignmentSummaryMetrics in OnTarget.bam ended"
	echo -e " \n Done : CollectAlignmentSummaryMetrics for  ($GATK_version) : $Lib"
# # ##############################################################################
	# ### Generate Bases callable by target file
	# echo -e "Generate Bases callable by target file started"
	awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable | $bedtools/intersectBed -a $CanonTranFile -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | $group/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.byTarget
	# echo -e "Generate Bases callable by target file ended"
# # ##############################################################################
	# #### Perl script to generate this format 
	# ##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
	# ##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
	# echo -e "perl script started"
	perl $StoreDir/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.byTarget > $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.byTarget.temp
	# echo -e "perl script ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.byTarget
	# echo -e "mv ended"
# # # ##############################################################################
# ## Generate  soft links in $TopDir/$RunDate/$Enrich_Report folder
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$OrgFile.q8.flagstat $SampleDir/"$Enrich_Report"/$SampleID/
	ln -s $MyOutDir/"$SampleID"/$OrgFile.flagstat $SampleDir/"$Enrich_Report"/$SampleID/

	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.flagstat $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoveragePerBase $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageHist $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.bases.callable.summary $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.AlignSumMet.txt $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.DepthOfCvg.txt.sample_interval_summary $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran"/
	# echo -e "ln ended"
# # ##############################################################################
# ###  No of callable Bases by Exon --minDepth 30
# echo -e "No of callable Bases by Exon --minDepth 30 started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
			-R $Ref \
			-T CallableLoci \
			-I $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile \
			-o $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable \
			--minMappingQuality 10 \
			--minBaseQuality 20 \
			--minDepth 30 \
			-l INFO -format BED \
			-L $CanonTranFile \
			-summary $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable.summary
	echo -e "No of callable Bases by Exon --minDepth 30 ended"
	echo -e " \n Done : CallableLoci  for ($GATK_version) : $Lib"
# # ##############################################################################
	# echo -e "awk started"
	awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable | $bedtools/intersectBed -a $CanonTranFile -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | $group/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable.byTarget
	# echo -e "awk ended"
# # # ##############################################################################
	# ### Perl script to generate this format 
	# #Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
	# #chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
	# echo -e "perl started"
	perl $StoreDir/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable.byTarget > $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable.byTarget.temp
	# echo -e "perl ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable.byTarget
	# echo -e "mv ended"
# # # ##############################################################################
# ## Generate  soft links in $TopDir/$RunDate/$Enrich_Report folder
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$CanonTran30x/$OnCanonTranFile.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran30x"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran30x"/
	# echo -e "ln ended"
# # ##############################################################################
# ##  No of callable Bases by Exon --minDepth 20
	# echo -e "No of callable Bases by Exon --minDepth 20 started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
			-R $Ref \
			-T CallableLoci \
			-I $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile \
			-o $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable \
			--minMappingQuality 10 \
			--minBaseQuality 20 \
			--minDepth 20 \
			-l INFO -format BED \
			-L $CanonTranFile \
			-summary $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable.summary
	echo -e "No of callable Bases by Exon --minDepth 20 ended"
	echo -e " \n Done : CallableLoci  for ($GATK_version) : $Lib"
# # # ##############################################################################
	# echo -e "awk started"
	awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable | $bedtools/intersectBed -a $CanonTranFile -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | $group/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable.byTarget
	# echo -e "awk ended"
# # # ##############################################################################
	# ### Perl script to generate this format 
	# #Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
	# #chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
	# echo -e "perl started"
	perl $StoreDir/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable.byTarget > $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable.byTarget.temp
	# echo -e "perl ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable.byTarget
	# echo -e "mv ended"
# # # ##############################################################################
# ### Generate  soft links in $TopDir/$RunDate/$Enrich_Report folder
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$CanonTran20x/$OnCanonTranFile.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran20x"/
	ln -s $MyOutDir/"$SampleID"/$CanonTran/$OnCanonTranFile.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$CanonTran20x"/
	# echo -e "ln ended"
# # # ##############################################################################
# # #End of Canonical Transcript
# # ##############################################################################
# # ## On Target File
	# ### CoverageStas [ duplicateReads were removed to calculate stats]
	echo -e "CoverageStats OnTargetFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Target/$OnTargetFile | $bedtools/coverageBed -abam stdin -b $TargetFile  > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageStats
	echo -e "CoverageStats OnTargetFile ended"
	### CoveragePerBase [ duplicateReads were removed to calculate stats]
	echo -e "CoveragePerBase OnTargetFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Target/$OnTargetFile | $bedtools/coverageBed -abam stdin -b $TargetFile -d > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoveragePerBase
	echo -e "CoveragePerBase OnTargetFile ended"
	### CoverageHistogram [ duplicateReads were removed to calculate stats]
	echo -e "CoverageHistogram OnTargetFile started"
	$samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Target/$OnTargetFile | $bedtools/coverageBed -abam stdin -b $TargetFile -hist > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageHist
	echo -e "CoverageHistogram OnTargetFile ended"
# # ##############################################################################
	# ## sort the CoverageBed, CoverageStats output by chr,position
	echo -e "sort the CoverageBed, CoverageStats output by chr,position started"
	sort -k 1,1 -k 2,2n $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageStats > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageStats.sorted
	echo -e "sort the CoverageBed, CoverageStats output by chr,position ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageStats.sorted $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageStats
	# echo -e "mv ended"
	echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoverageStats completed ($GATK_version) : $Lib"
# # ##############################################################################
	# ### sort the CoverageBed,CoveragePerBase output by chr,position,coverage (Change colum 5 into 4 in Exome $Target)
	echo -e "sort the CoverageBed,CoveragePerBase output by chr,position,coverage (Change colum 5 into 4 in Exome $Target) started"
	sort -k 1,1 -k 2,2n -k 5,5n $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoveragePerBase  > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoveragePerBase.sorted
	echo -e "sort the CoverageBed,CoveragePerBase output by chr,position,coverage (Change colum 5 into 4 in Exome $Target) ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoveragePerBase.sorted $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoveragePerBase
	# echo -e "mv ended"
	echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoveragePerBase completed ($GATK_version) : $Lib"
# # ##############################################################################
# ##  No of callable Bases by Target
	echo -e "No of callable Bases by Target started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T CallableLoci \
		-I $MyOutDir/"$SampleID"/$Target/$OnTargetFile \
		-o $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable \
		--minMappingQuality 10 \
		--minBaseQuality 20 \
		--minDepth 10 \
		-l INFO -format BED \
		-L $TargetFile \
		-summary $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.summary
	echo -e "No of callable Bases by Target started"
	echo -e " \n Done : CallableLoci  for ($GATK_version) : $Lib"
# # ###############################################################################
# ## Depth of coverage by Target
	echo -e "Depth of coverage by Target started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T DepthOfCoverage \
		-I $MyOutDir/"$SampleID"/$Target/$OnTargetFile \
		-o $MyOutDir/"$SampleID"/$Target/$OnTargetFile.DepthOfCvg.txt \
		-L $TargetFile \
		--minMappingQuality 10 \
		--minBaseQuality 20 \
		--stop 2000 \
		--includeDeletions --summaryCoverageThreshold \
		--outputFormat table --summaryCoverageThreshold --omitDepthOutputAtEachBase
	echo -e "Depth of coverage by Target ended"
	echo -e " \n Done : DepthOfCoverage for ($GATK_version)  : $Lib"
# # ###############################################################################
	# # Getting mean coverage per Exon
	# echo -e "cat started"
	cat $MyOutDir/"$SampleID"/$Target/$OnTargetFile.DepthOfCvg.txt.sample_interval_summary | sed 's/chr16:997401/chr16:997401-997401/' | awk 'BEGIN {FS=":"} {OFS="\t"} {print $1,$2}'| awk 'BEGIN {FS="-"} {OFS="\t"} {print $1,$2}' |  awk 'BEGIN {FS=" "} {OFS="\t"} {print $1,$2-1,$3,$4,$5,$9}' | sed '1d'| $bedtools/sortBed -i stdin > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.temp.depth1.bed
	# echo -e "cat ended"
	
	# echo -e "cut started"
	cut -f4 $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageStats | awk 'BEGIN {FS="|";} {print $1}' > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.temp.depth2.bed
	# echo -e "cut ended"
	
	# echo -e "paste started"
	paste <(cut -f1-3 $MyOutDir/"$SampleID"/$Target/$OnTargetFile.temp.depth1.bed) <(cut -f4 $MyOutDir/"$SampleID"/$Target/$OnTargetFile.temp.depth2.bed ) <(cut -f4- $MyOutDir/"$SampleID"/$Target/$OnTargetFile.temp.depth1.bed) | cut -f1-4,6 > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed 
	# echo -e "paste ended"
	
	# echo -e "rm started"
	rm -rf $MyOutDir/"$SampleID"/$Target/$OnTargetFile.temp*
	# echo -e "rm ended"
	
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	# echo -e "ln ended"
# # ###############################################################################
# ## CollectAlignmentSummaryMetrics in OnTarget.bam
	echo -e "CollectAlignmentSummaryMetrics in OnTarget.bam started"
	$Java7 $memory -jar $picard/CollectAlignmentSummaryMetrics.jar \
		INPUT=$MyOutDir/"$SampleID"/$Target/$OnTargetFile \
		OUTPUT=$MyOutDir/"$SampleID"/$Target/$OnTargetFile.AlignSumMet.txt \
		REFERENCE_SEQUENCE=$Ref \
		ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT
	echo -e "CollectAlignmentSummaryMetrics in OnTarget.bam ended"
	echo -e " \n Done : CollectAlignmentSummaryMetrics for  ($GATK_version) : $Lib"
# # ###############################################################################
# # ## Generate Bases callable by target file
	# echo -e "awk started"
	awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable | $bedtools/intersectBed -a $TargetFile -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | $group/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.byTarget
	# echo -e "awk ended"
# # ###############################################################################
	# #### Perl script to generate this format
	# ##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
	# ##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
	# echo -e "perl started"
	perl $StoreDir/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.byTarget > $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.byTarget.temp
	# echo -e "perl ended"
	# echo -e "mv started"
	mv $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.byTarget
	# echo -e "mv ended"
# # ###############################################################################
# # ### Generate  soft links in $TopDir/$RunDate/$Enrich_Report folder

	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.flagstat $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.flagstat $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoveragePerBase $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.CoverageHist $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.bases.callable.summary $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.AlignSumMet.txt $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	ln -s $MyOutDir/"$SampleID"/$Target/$OnTargetFile.DepthOfCvg.txt.sample_interval_summary $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/


# # #######################################################################################################################################
# # # ************************************************GATK - Indels and SNPs calling******************************************************#
# # #######################################################################################################################################


mkdir $MyOutDir/$SampleID/$Target/"UnifiedGenotyper"


# # ###############################################################################
# ##Calling UnifiedGenotyper  
	# ## -stand_emit_conf 30.0 into 10; --min_base_quality_score 20 into 10; Recommended by Risha
	# ## For GATK_2_3, use -A DepthOfCoverage instead of -A Coverage
	echo -e "Calling UnifiedGenotyper started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar -l INFO \
			-R $Ref \
			-L $TargetFile \
			-I $MyOutDir/"$SampleID"/$Target/$OnTargetFile \
			-T UnifiedGenotyper \
			-baq CALCULATE_AS_NECESSARY \
			-minIndelCnt 4 \
			-glm BOTH \
			--dbsnp $dbSNP \
			-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.indel.vcf \
			-A Coverage \
			-A AlleleBalance \
			-G Standard \
			--min_base_quality_score 10 \
			--metrics_file $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.indel.metrics \
			--num_threads $nt \
			-stand_call_conf 30.0 \
			-stand_emit_conf 10 \
			--downsample_to_coverage $dcovg \
			--downsampling_type BY_SAMPLE \
			--computeSLOD
	
	Step_Check $? $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.indel.vcf 1
	echo -e "Step_Check Calling UnifiedGenotyper ended"
	
	echo -e "\nDone : UnifiedGenotyper - SNP & Indels calling ($GATK_version) : $Lib"
# # ###############################################################################
# ## print SNP variant
	echo -e "print snp variant started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-T SelectVariants \
		--downsample_to_coverage $dcovg \
		--downsampling_type BY_SAMPLE \
		-R $Ref \
		-nt $nt \
		--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.indel.vcf \
		--selectTypeToInclude INDEL \
		-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.indel.vcf
	
	Step_Check $? $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.indel.vcf 1
	echo -e "Step_Check print snp variant ended"

	echo -e "\nDone : UnifiedGenotyper - Indels calling ($GATK_version) : $Lib"
# # ###############################################################################
# ## Indel filter
	echo -e "Indel filter started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T VariantFiltration \
		--downsample_to_coverage $dcovg \
		--downsampling_type BY_SAMPLE \
		-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.indel.filtered.vcf \
		--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.indel.vcf \
		-baq CALCULATE_AS_NECESSARY \
		--filterExpression "QD < 2.0" \
		--filterExpression "ReadPosRankSum < -20.0" \
		--filterExpression "FS > 200.0" \
		--filterName QDFilter \
		--filterName ReadPosFilter \
		--filterName FSFilter 
		
	Step_Check $? $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.indel.filtered.vcf 1
	echo -e "Step_Check Indel filter ended"
# # ###############################################################################
## print Indel variant
	echo -e "print Indel variant started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $Ref \
		--downsample_to_coverage $dcovg \
		--downsampling_type BY_SAMPLE \
		-nt $nt \
		--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.indel.vcf \
		--selectTypeToInclude SNP \
		-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.vcf
	
	Step_Check $? $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.vcf 1
	echo -e "Step_Check print Indel variant ended"
	
	echo -e "\nDone : UnifiedGenotyper - Indel filtering ($GATK_version) : $Lib"
# # ###############################################################################
# ## SNPs Filter
	echo -e "SNPs Filter started"
	$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T VariantFiltration \
		--downsample_to_coverage $dcovg \
		--downsampling_type BY_SAMPLE \
		-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.filtered.vcf \
		--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.vcf \
		--mask $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.indel.vcf \
		--maskName InDel \
		--filterExpression "QD < 2.0" \
		--filterExpression "MQ < 40.0" \
		--filterExpression "FS > 60.0" \
		--filterExpression "MQRankSum < -12.5" \
		--filterExpression "ReadPosRankSum < -8.0" \
		--filterName QDFilter \
		--filterName MQFilter \
		--filterName FSFilter \
		--filterName MQRankSumFilter \
		--filterName ReadPosFilter 
	
	Step_Check $? $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.filtered.vcf 1
	echo -e "Step_Check SNPs Filter ended"

	echo -e "\nDone : SNP filtering ($GATK_version) : $Lib"
	echo -e "\nDone : GATK & SAMTools analysis were completed for ($GATK_version) : $Lib"
# # ###############################################################################
	# echo -e "ln started"
	##Generating final filtered SNP VCF file
	ln -s $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.filtered.vcf $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.final.UnifiedGenotyper.snp.vcf
	##Generating final filtered Indel VCF file
	ln -s $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.indel.filtered.vcf $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.final.UnifiedGenotyper.indel.vcf
	# echo -e "ln ended"
	
	echo -e "\nDone : GATK UnifiedGenotyper analysis were completed for ($GATK_version) : $Lib"
###############################################################################
##Remove markDup, SAM, realigned BAM files
	# echo -e "rm1 started"
	rm -rf $MyOutDir/"$SampleID"/*.sam
	rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.bam
	rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.bai
	rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.intervals
	rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.metrics.txt
	#No such files to remove
	# rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bai
	# rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam
	# rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam.bai
	# rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam.recalibTable.grp
	#
	rm -rf $MyOutDir/"$SampleID"/*.L00*

	rm -rf $MyOutDir/$SampleID/"$SampleID".bam
	rm -rf $MyOutDir/$SampleID/"$SampleID".bam.bai
	rm -rf $MyOutDir/$SampleID/"$SampleID".markDup.Realigned.recalibrated.bai
	# echo -e "rm1 ended"
	
	# echo -e "rm2 started"
	rm -rf $MyOutDir/$SampleID/ProteinCodingTarget/"$SampleID".markDup.Realigned.recalibrated.ProteinCoding.OnTarget.q8.bam
	rm -rf $MyOutDir/$SampleID/ProteinCodingTarget/"$SampleID".markDup.Realigned.recalibrated.ProteinCoding.OnTarget.q8.bam.bai
	# echo -e "rm2 ended"
	
	# echo -e "rm3 started"
	rm -rf $MyOutDir/$SampleID/CanonTranCodingTarget/"$SampleID".markDup.Realigned.recalibrated.CanonTranCoding.OnTarget.q8.bam
	rm -rf $MyOutDir/$SampleID/CanonTranCodingTarget/"$SampleID".markDup.Realigned.recalibrated.CanonTranCoding.OnTarget.q8.bam.bai
	# echo -e "rm3 ended"
	
	# echo -e "rm4 started"
	rm -rf $MyOutDir/$SampleID/Target/UnifiedGenotyper/"$SampleID".markDup.Realigned.recalibrated.OnTarget.q8.bam.UnifiedGenotyper.indel.vcf
	rm -rf $MyOutDir/$SampleID/Target/UnifiedGenotyper/"$SampleID".markDup.Realigned.recalibrated.OnTarget.q8.bam.UnifiedGenotyper.indel.vcf.idx
	rm -rf $MyOutDir/$SampleID/Target/UnifiedGenotyper/"$SampleID".markDup.Realigned.recalibrated.OnTarget.q8.bam.UnifiedGenotyper.snp.vcf
	rm -rf $MyOutDir/$SampleID/Target/UnifiedGenotyper/"$SampleID".markDup.Realigned.recalibrated.OnTarget.q8.bam.UnifiedGenotyper.snp.vcf.idx
	# echo -e "rm4 ended"
# # ###############################################################################
# #Ts/TV ratio for UnifiedGenotyer SNP
	echo -e "Ts/TV ratio for UnifiedGenotyer SNP started"
	$Java7 $memory -jar $snpeff/SnpSift.jar tstv any $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.filtered.vcf  > $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.filtered.vcf.all.snp.TsTv.txt
	echo -e "Ts/TV ratio for UnifiedGenotyer SNP ended"
# # ################################################################################	
	# ## Symbolic link
	# echo -e "ln started"
	ln -s $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$OnTargetFile.UnifiedGenotyper.snp.filtered.vcf.all.snp.TsTv.txt $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/
	# echo -e "ln ended"
# # ################################################################################
	# echo -e "chmod started"
	 chmod 775 -R $SampleDir
	# echo -e "chmod ended"
# # ################################################################################	
	#Exiting For Loop, reset variables to default 
	SampleDir=$SampleDir_default
	MyOutDir=$MyOutDir_default
	GATKs=$GATKs_default
# ################################################################################
# # ##HaplotypeCaller
# # ## GATK new variant Call method -T HaplotypeCaller; This will replace UnifiedGenotyper in future
# # ################################################################################


mkdir $MyOutDir/"$SampleID"/$Target/"HaplotypeCaller"


# # # ################################################################################
# # HaplotypeCaller ##-nct $nt  not support with --bamOutput
echo -e "HaplotypeCaller started"
$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-L $TargetFile \
	-I $MyOutDir/"$SampleID"/$Target/$OnTargetFile \
	-T HaplotypeCaller \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-A Coverage \
	-A AlleleBalance \
	-G Standard \
	--dbsnp $dbSNP \
	--min_base_quality_score 10 \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.vcf \
	-stand_call_conf 30 \
	-stand_emit_conf 10 \
	--bamOutput $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.bam

Step_Check $? $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.vcf 1
echo -e "Step_Check HaplotypeCaller ended"

# echo -e "mv started"
mv $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.bai $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.bam.bai
# echo -e "mv ended"
# # ################################################################################
# ## Variant Annotater for haplotypescore Allele Balance
echo -e "Variant Annotater for haplotypescore Allele Balance started"
$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantAnnotator \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-L $TargetFile \
	-I $MyOutDir/"$SampleID"/$Target/$OnTargetFile \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.varAnnotate.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.vcf \
	-A Coverage \
	-A AlleleBalance \
	-A HaplotypeScore \
	-A InbreedingCoeff \
	-A HomopolymerRun \
	-A HardyWeinberg \
	-A GCContent \
	--dbsnp $dbSNP \
	-nt $nt 

Step_Check $? $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.varAnnotate.vcf 1
echo -e "Step_Check Variant Annotater for haplotypescore Allele Balance ended"
# # ################################################################################
# ## print SNP variant
echo -e "print SNP variant started"
$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
	-T SelectVariants \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-R $Ref \
	-nt $nt \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.varAnnotate.vcf \
	--selectTypeToInclude INDEL \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.indel.vcf
	
Step_Check $? $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.indel.vcf 1
echo -e "Step_Check print SNP variant ended"	
echo -e "\nDone : HaplotypeCaller Indels calling : $Lib"
# # ################################################################################
# ##Indel filter
echo -e "Indel filter started"
$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantFiltration \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.indel.filtered.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.indel.vcf \
	-baq CALCULATE_AS_NECESSARY \
	--filterExpression "QD < 2.0" \
	--filterExpression "ReadPosRankSum < -20.0" \
	--filterExpression "FS > 200.0" \
	--filterName QDFilter \
	--filterName ReadPosFilter \
	--filterName FSFilter 
	
Step_Check $? $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.indel.filtered.vcf 1
echo -e "Step_Check Indel filter ended"
# # ################################################################################
# ## print Indel variant
echo -e "print Indel variant started"
$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $Ref \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-nt $nt \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.varAnnotate.vcf \
	--selectTypeToInclude SNP \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.vcf
	
Step_Check $? $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.vcf 1
echo -e "Step_Check print Indel variant ended"

echo -e "\nDone : HaplotypeCaller Indel filtering : $Lib"
# # ################################################################################
# # SNPs Filter
echo -e "SNPs Filter started"
$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantFiltration \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.filtered.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.vcf \
	--mask $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.indel.vcf \
	--maskName InDel \
	--filterExpression "QD < 2.0" \
	--filterExpression "MQ < 40.0" \
	--filterExpression "FS > 60.0" \
	--filterExpression "MQRankSum < -12.5" \
	--filterExpression "ReadPosRankSum < -8.0" \
	--filterName QDFilter \
	--filterName MQFilter \
	--filterName FSFilter \
	--filterName MQRankSumFilter \
	--filterName ReadPosFilter

Step_Check $? $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.filtered.vcf 1
echo -e "Step_Check SNPs Filter ended"
echo -e "\nDone :  HaplotypeCaller SNP filtering : $Lib"
# # ################################################################################
# ##Generating final filtered SNP VCF file
ln -s $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.filtered.vcf $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.final.HaplotypeCaller.snp.vcf
###Generating final filtered Indel VCF file
ln -s $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.indel.filtered.vcf $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.final.HaplotypeCaller.indel.vcf
# # ################################################################################
# ## HaplotypeCaller -  generate genomic VCFs (gVCFS) 
$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-L $TargetFile \
	-I $MyOutDir/"$SampleID"/$Target/$OnTargetFile \
	-T HaplotypeCaller \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-A Coverage \
	-A AlleleBalance \
	-G Standard \
	-nct $nt \
	--dbsnp $dbSNP \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.GVCF.vcf \
	-stand_call_conf 30 \
	-stand_emit_conf 10 \
	--emitRefConfidence GVCF \
	--variant_index_type LINEAR \
	--variant_index_parameter 128000

Step_Check $? $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.indel.GVCF.vcf 1

# # ################################################################################

# # ################################################################################
# ##Ts/Tv ratio for HaplotyeCaller SNP
$Java7 $memory -jar $snpeff/SnpSift.jar tstv any $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.filtered.vcf  > $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.filtered.vcf.all.snp.TsTv.txt
# # # ################################################################################
# ##Symbolic link
ln -s $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$OnTargetFile.HaplotypeCaller.snp.filtered.vcf.all.snp.TsTv.txt $SampleDir/"$Enrich_Report"/$SampleID/"$Target"/

# ##############################################################################

    mkdir $MyOutDir/"$SampleID"/$Target/"VEP"

# ###############################################################################
    for caller in HaplotypeCaller UnifiedGenotyper
     do
	   for variantType in snp indel
	    do
	     # echo -e "VEP file assignments started"
	     inVEP=$MyOutDir/"$SampleID"/$Target/$caller/$OnTargetFile.$caller.$variantType.filtered.vcf
	     outVEP=$MyOutDir/"$SampleID"/$Target/VEP/$OnTargetFile.$caller.$variantType.VEP.vcf
	     # echo -e "VEP file assignments ended"
	     ##Testing
	     # echo -e "inVEP: $inVEP"
	     # echo -e "outVEP: $outVEP"
	

	## Count lines that don't match '#' in $inVEP -> exclude metadata lines -> count variants
	  # echo -e "VariantCount started"
	  VariantCount=`grep -cv '#' $inVEP`
	  # echo -e "VariantCount ended"
	##Testing
	  # echo -e "VariantCount: $VariantCount"
	
	# echo -e "Check VariantCount started"
	 if [ $VariantCount == 0 ]
	  then
		echo -e "\nERROR : No Variants in this VCF : $inVEP\n"
		continue
	 fi
	
	 perl $VEP/variant_effect_predictor.pl \
	  --dir_cache /home/ahcauc/.vep \
	  --assembly GRCh37 \
	  --fasta /home/ahcauc/.vep/homo_sapiens/83_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
	  --cache --offline --no_progress --fork $nt -i $inVEP  -o $outVEP \
	  --allele_number --hgvs --check_existing --canonical --ccds --maf_exac --pubmed --protein --sift b --polyphen b \
	  --fields ALLELE_NUM,Consequence,SYMBOL,SYMBOL_SOURCE,Gene,Feature,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,STRAND,CANONICAL,CCDS,ENSP,SIFT,PolyPhen,ExAC_MAF,PUBMED  \
	  --vcf --force_overwrite
	
	
	 Step_Check $? $outVEP 1

	
	##Tablelise VEP Annotated VCFs
	 out_tableize=$MyOutDir/"$SampleID"/$Target/VEP/$OnTargetFile.$caller.$variantType.VEP.txt

	 python $LOFTEE/tableize_vcf.py \
	 --vcf $outVEP --out $out_tableize --split_by_transcript --all_csqs --do_not_minrep --include_id --info ABHet,ABHom,AC,AF,AN,DP,FS,HaplotypeScore,MLEAC,MLEAF,MQ,MQ0,QD \
	 --vep_info Consequence,SYMBOL,SYMBOL_SOURCE,Gene,Feature,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,STRAND,CANONICAL,CCDS,ENSP,SIFT,PolyPhen,ExAC_MAF,PUBMED \
	 --mysql

	  done
   done
# # ################################################################################
# # ################################################################################
done
# # ################################################################################