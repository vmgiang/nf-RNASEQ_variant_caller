process ADD_READ_GROUP {
    tag "$sampleId"
    label "mem_large"
    // publishDir "${params.results}/add_read_group/${sampleId}", mode : "copy"
    
    input:
        tuple val(sampleId), path(bam), path(bai)
    
    output:
        tuple val(sampleId), path('add_read_group.bam'), path('add_read_group.bam.bai')
    
    script:
    """
    # Add read group
    gatk AddOrReplaceReadGroups \
            -I $bam \
            -O add_read_group.bam \
            -RGID 4 \
            -RGLB lib1 \
            -RGPL illumina \
            -RGPU unit1 \
            -RGSM 20
    samtools index add_read_group.bam
    """
}


/*
 * Process 1: Split reads that contain Ns in their CIGAR string.
 *            Creates k+1 new reads (where k is the number of N cigar elements) 
 *            that correspond to the segments of the original read beside/between 
 *            the splicing events represented by the Ns in the original CIGAR.
 */

process RNASEQ_GATK_SPLITNCIGAR {
  tag "$replicateId"
  label 'mem_large'
  
  input: 
    path genome
    path index
    path genome_dict
    tuple val(replicateId), path(bam), path(index)

  output:
    tuple val(replicateId), path('split.bam'), path('split.bai')
  
  script:
  """
  # SplitNCigarReads and reassign mapping qualities
  gatk SplitNCigarReads \
            -R $genome \
            -I $bam \
            --refactor-cigar-string \
            -O split.bam
  """
}


/*
 * Process 2: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */

process RNASEQ_GATK_RECALIBRATE {
  tag "$replicateId"
  label "mem_large"

  input: 
    path genome
    path index
    path dict
    tuple val(replicateId), path(bam), path(index)
    tuple path(variants_file), path(variants_file_index)

  output:
    tuple \
      val(sampleId), \
      path("${replicateId}.final.uniq.bam"), \
      path("${replicateId}.final.uniq.bam.bai")
  
  script: 
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  # Indel Realignment and Base Recalibration
  gatk BaseRecalibrator \
          -R $genome \
          -I $bam \
          --known-sites $variants_file \
          --use-original-qualities \
          -O final.rnaseq.grp 

  gatk ApplyBQSR \
          -R $genome -I $bam \
          --use-original-qualities \
          --bqsr-recal-file final.rnaseq.grp \
          -O ${replicateId}.final.uniq.bam

  # Index BAM files
  samtools index ${replicateId}.final.uniq.bam
  """
}

/*
 * Process 3: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */

process RNASEQ_CALL_VARIANTS {
  tag "$sampleId"
  label "mem_xlarge"
  publishDir "${params.results}/call_variant/${sampleId}", mode : "copy"

  input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(bam), path(bai)
 
  output: 
    tuple val(sampleId), path('final.vcf')

  script:
  def bam_params = bam.collect{ "-I $it" }.join(' ')
  """
  # fix absolute path in dict file
  sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  
  # Variant calling
  gatk HaplotypeCaller \
          --native-pair-hmm-threads ${task.cpus} \
          --reference ${genome} \
          --output output.gatk.vcf.gz \
          ${bam_params} \
          --standard-min-confidence-threshold-for-calling 20.0 \
          --dont-use-soft-clipped-bases 

  # Variant filtering
  gatk VariantFiltration \
          -R ${genome} -V output.gatk.vcf.gz \
          --cluster-window-size 35 --cluster-size 3 \
          --filter-name FS --filter \"FS > 30.0\" \
          --filter-name QD --filter \"QD < 2.0\" \
          -O final.vcf
  """
}

process PREPARE_VCF_FILE {
  tag "$variantsFile.baseName"

  input: 
    path variantsFile

  output:
    tuple \
      path("${variantsFile}"), \
      path("${variantsFile}.tbi")
  
  script:  
  """
  tabix -p vcf $variantsFile
  """
}