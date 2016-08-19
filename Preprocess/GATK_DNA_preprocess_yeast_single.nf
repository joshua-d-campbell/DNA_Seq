#!/usr/bin/env nextflow


//############################################################################################################################
//
// Josh Campbell
// 7/20/2016
// Peforms alignment and preprocessing of paired-end exome sequencing data.
// For all samples derived from the same individual, an indel co-cleaning step will be performed on all bams jointly
//
//############################################################################################################################


//############################################################################################################################
//
// GATK tutorials this pipeline is based on:
// Mapping: http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
// Marking Duplicates: http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
// Realignment around indels: http://gatkforums.broadinstitute.org/gatk/discussion/2800/howto-perform-local-realignment-around-indels
// 
//############################################################################################################################

// List of parameters that can be passed to this workflow
params.ref = "/restricted/projectnb/cbmhive/references/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_S288C_R64-2-1_20150113.fasta"
params.gatk_jar = "/share/pkg/gatk/3.5/install/GenomeAnalysisTK.jar"
params.picard_jar = "/share/pkg/picard/2.1.1/install/picard-tools-2.1.1/picard.jar"
params.output_dir = "./"

// Set up global variables for requried parameters:
PROJECT = params.project
inputFile = file(params.infile)

// Set up global variables for parameters with preset defaults:
REF = file(params.ref)
GATK = file(params.gatk_jar)
PICARD = file(params.picard_jar)
OUTDIR = file(params.output_dir)


// Necessary columns needed in input file:
// INDIVIDUAl_ID	SAMPLE_ID	LIBRARY_ID	RG_ID	PLATFORM_UNIT	PLATFORM	PLATFORM_MODEL	RUN_DATE	CENTER	R1	R2

//#############################################################################################################
//
// Read input file and save it into list of lists
//
//#############################################################################################################

file_params = []
is_header = 1
for( line in inputFile.readLines() ) {
  if(is_header == 1) {
    header = line.split("\t").flatten()
    is_header = 0
  } else {
    file_params << line.split("\t").flatten()
  }  
}


// Send FASTQ files to two processes: FastQC and FastqToSam
(readPairsFastQC, readPairsFastqToSam) = Channel.from(file_params).into(2)


//#############################################################################################################
//
// Preprocess reads
// 1) Convert to BAM
// 2) Mark Adapters
//
//#############################################################################################################

process runFastqToSam {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5G"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastqToSam/"
    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastqToSam
    
    output:
    set indivID, sampleID, libraryID, rgID, file(outfile) into runFastqToSamOutput

    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".unaligned.bam"
    
    """
    module load java/1.8.0_66
    
	java -Xmx5G -XX:ParallelGCThreads=1 -jar ${PICARD} FastqToSam \
		FASTQ=${fastqR1} \
		FASTQ2=${fastqR2} \
		OUTPUT=${outfile} \
		READ_GROUP_NAME=${rgID} \
		SAMPLE_NAME=${sampleID} \
		LIBRARY_NAME=${libraryID} \
		PLATFORM_UNIT=${platform_unit} \
		PLATFORM=${platform} \
		PLATFORM_MODEL=${platform_model} \
		SEQUENCING_CENTER=${center} \
		RUN_DATE=${run_date} \
		TMP_DIR=tmp
    """
}

process runMarkIlluminaAdapters {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=6G"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/MarkIlluminaAdapters/"
    
    input:
    set indivID, sampleID, libraryID, rgID, ubam from runFastqToSamOutput
    
    output:
    set indivID, sampleID, libraryID, rgID, ubam, file(outfile_bam), file(outfile_metrics) into runMarkIlluminaAdaptersOutput
	
    script:
    outfile_bam = sampleID + "_" + libraryID + "_" + rgID + ".adapters_marked.bam"
    outfile_metrics = sampleID + "_" + libraryID + "_" + rgID + "_adapters_metrics.txt"
            
    """
    module load java/1.8.0_66
    
	java -Xmx5G -XX:ParallelGCThreads=1 -jar ${PICARD} MarkIlluminaAdapters \
		I=${ubam} \
		O=${outfile_bam} \
		M=${outfile_metrics} \
		TMP_DIR=tmp

    """
}




//#############################################################################################################
//
// Run BWA to align reads to genome
//
//#############################################################################################################

process runBWA {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/"
	
    input:
    set indivID, sampleID, libraryID, rgID, ubam, ubamxt, metrics from runMarkIlluminaAdaptersOutput
    
    output:
    set indivID, sampleID, libraryID, rgID, ubam, ubamxt into deleteBWAInput
    set indivID, sampleID, file(outfile_bam) into runBWAOutput
    
    script:
    outfile_bam = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"
	
    """
    module load java/1.8.0_66
    module load bwa/0.7.12
        
	set -o pipefail
	java -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=1 -Xmx128m -jar ${PICARD} SamToFastq \
		I=${ubamxt} \
		FASTQ=/dev/stdout \
		CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
		TMP_DIR=tmp | \
	bwa mem -M -t 10 -p ${REF} /dev/stdin | \
	java -XX:ParallelGCThreads=1 -Xmx4G -jar ${PICARD} MergeBamAlignment \
		ALIGNED_BAM=/dev/stdin \
		UNMAPPED_BAM=${ubamxt} \
		OUTPUT=${outfile_bam} \
		R=${REF} CREATE_INDEX=true ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS \
		TMP_DIR=tmp
	"""	
}




//#############################################################################################################
//
// Combined libraries from the same Individual/Sample to send to MarkDuplicates
//
//#############################################################################################################

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])




//#############################################################################################################
//
// Run Picard MarkDuplicates
// This is used to merge different libraries from the same sample or the same library run on different lanes
// Requires a lot of memory
// Need to set "ParallelGCThreads" otherwise it will "grab" extra available threads without asking (and potentially be terminated by SGE)
//
//#############################################################################################################

process runMarkDuplicates {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=48:00:00 -l mem_total=94G -pe omp 5"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/MarkDuplicates"
	
    input:
    set indivID, sampleID, aligned_bam_list from runBWAOutput_grouped_by_sample
    
    output:
    set indivID, sampleID, aligned_bam_list into deleteMarkDuplicatesInput
    set indivID, sampleID, file(outfile_bam), file(outfile_metrics) into runMarkDuplicatesOutput
 
    script:
    outfile_bam = sampleID + ".dedup.bam"
    outfile_metrics = sampleID + "_duplicate_metrics.txt"	
	        
    """
    module load java/1.8.0_66
    
	java -Xmx25G -XX:ParallelGCThreads=5 -Djava.io.tmpdir=tmp/ -jar ${PICARD} MarkDuplicates \
		INPUT=${aligned_bam_list.join(" INPUT=")} \
		OUTPUT=${outfile_bam} \
		METRICS_FILE=${outfile_metrics} \
		CREATE_INDEX=true \
		TMP_DIR=tmp
	"""  
}





//#############################################################################################################
//
// Perform realignment around indels
// 1) Identify regions for realignement
// 2) Perform realignment
//
//#############################################################################################################

process runRealignerTargetCreator {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=96:00:00 -l mem_total=25G"
    publishDir "${OUTDIR}/${indivID}/Processing/RealignerTargetCreator/"
    
    input:
    set indivID, sampleID, dedup_bam, metrics from runMarkDuplicatesOutput
    
    output:
    set indivID, sampleID, dedup_bam, file(target_file) into runRealignerTargetCreatorOutput
 	
    script:
    target_file = sampleID + "_target_intervals.list"
	        
    """
    module load java/1.8.0_66

	java -Xmx15g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T RealignerTargetCreator \
		-R ${REF} \
		-I ${dedup_bam} \
		-o ${target_file}
	"""  
}

process runIndelRealigner {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=48:00:00 -l mem_total=96G"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/"
	    
    input:
    set indivID, sampleID, dedup_bam, target_file from runRealignerTargetCreatorOutput
 	    
    output:
    set indivID, file(outfile_bam), file(outfile_bai) into runIndelRealignerOutput_for_DepthOfCoverage, runIndelRealignerOutput_for_Multiple_Metrics 
	set indivID, dedup_bam into deleteRealignerTargetCreatorInput
 
 	script:
    outfile_bam = sampleID + ".clean.bam"            
    outfile_bai = sampleID + ".clean.bai"
  	
    """
    module load java/1.8.0_66

	java -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T IndelRealigner \
		-R ${REF} \
		-I ${dedup_bam} \
		-targetIntervals ${target_file} \
		-o ${outfile_bam}
	"""  
}




//#############################################################################################################
//
// Perform a several tasks to assess QC:
// 1) Depth of coverage over targets
// 2) Generate alignment stats, insert size stats, quality score distribution
// 3) Run FASTQC to assess read quality
//
//#############################################################################################################

process runDepthOfCoverage {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=10G"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/DepthOfCoverage"
	    
    input:
    set indivID, sampleID, bam, bai from runIndelRealignerOutput_for_DepthOfCoverage

    output:
    file("${prefix}*") into DepthOfCoverageOutput
    
    script:
    prefix = sampleID + "."
         
    """
    module load java/1.8.0_66

	java -XX:ParallelGCThreads=1 -Djava.io.tmpdir=tmp/ -Xmx10g -jar ${GATK} \
		-R ${REF} \
		-T DepthOfCoverage \
		-I ${bam} \
		--omitDepthOutputAtEachBase \
		-L ${TARGETS} \
		-ct 10 -ct 20 -ct 50 -ct 100 \
		-o ${sampleID}

	"""
}	



process runCollectMultipleMetrics {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=25G"
 	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics"
 	    
    input:
    set indivID, sampleID, bam, bai from runIndelRealignerOutput_for_Multiple_Metrics

    output:
    file("${prefix}*") into CollectMultipleMetricsOutput

    script:       
    prefix = sampleID + "."

    """
    module load java/1.8.0_66

	java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectMultipleMetrics \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics\
		INPUT=${bam} \
		REFERENCE_SEQUENCE=${REF} \
		ASSUME_SORTED=true \
		OUTPUT=${prefix} \
		TMP_DIR=tmp
	"""
}	


process runFastQC {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5G"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastQC/"
	    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastQC

    output:
    set file("*.zip"), file("*.html") into FastQCOutput
    	
    script:

    """
    module load fastqc/0.11.3
    fastqc -t 1 -o . ${fastqR1} ${fastqR2}
    """
}




//#############################################################################################################
//
// "Garbage Collection" processes designed to delete large BAM files once they are no longer needed to save on space
//
//#############################################################################################################

process runDeleteBWAInput {
    tag "${indivID}|${sampleID}|${rgID}|${libraryID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=1:00:00"
    
    input:
    set indivID, sampleID, libraryID, rgID, ubam, ubamxt from deleteBWAInput
    
    script:
    """
    rm -v ${ubam} ${ubamxt}
	"""
}

process runDeleteMarkDuplicatesInput {
    tag "${indivID}|${sampleID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=1:00:00"
    
    input:
    set indivID, sampleID, aligned_bam_list from deleteMarkDuplicatesInput
    
    script:
    cl = aligned_bam_list.getClass()
    cl_n = cl.getName()
    if(cl_n == "sun.nio.fs.UnixPath") {
      aligned_bai_list = aligned_bam_list.toString().replaceAll(".bam", ".bai")
    }
    else {
      aligned_bai_list = aligned_bam_list.join(" ").replaceAll(".bam", ".bai")
      aligned_bam_list = aligned_bam_list.join(" ")      
    }
    
    """
    rm -v ${aligned_bam_list} 
    rm -v ${aligned_bai_list} 
    """
}

process runDeleteRealignerTargetCreator {
    tag "${indivID}"
    executor = 'sge'
	clusterOptions = "-P ${PROJECT} -l h_rt=1:00:00"
    
    input:
    set indivID, dedup_bam_list from deleteRealignerTargetCreatorInput
    
    script:
    cl = dedup_bam_list.getClass()
    cl_n = cl.getName()
    if(cl_n == "sun.nio.fs.UnixPath") {
      dedup_bai_list = dedup_bam_list.toString().replaceAll(".bam", ".bai")
    }
    else {
      dedup_bai_list = dedup_bam_list.join(" ").replaceAll(".bam", ".bai")
      dedup_bam_list = dedup_bam_list.join(" ")
    }
  
    """
    rm -v ${dedup_bam_list} 
    rm -v ${dedup_bai_list} 
    """
}

