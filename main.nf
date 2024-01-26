params.genomes = "/home/giangvm1/RNAseq/ASE/ref/Sequence/WholeGenomeFasta/genome.fa" // Here is reference genome
params.index = "/home/giangvm1/RNAseq/ASE/ref/Sequence/WholeGenomeFasta/genome.fa.fai" //index of reference genome
params.dict = "/home/giangvm1/RNAseq/ASE/ref/Sequence/WholeGenomeFasta/genome.dict" //dict of reference genome
params.results = "$baseDir/result2" //outdir
params.input = "/home/giangvm1/RNAseq/ASE/ASA_calling_star/VN_disease_33.txt" // input .csv file 
params.dbsnp="/home/giangvm1/RNAseq/ASE/script/resource_grch38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz" //dbsnp

include {
    RNASEQ_GATK_SPLITNCIGAR ;
    RNASEQ_GATK_RECALIBRATE ;
    RNASEQ_CALL_VARIANTS;
    PREPARE_VCF_FILE;
    ADD_READ_GROUP
} from "./modules.nf"

workflow {
    Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ",")
        .map {create_bam_bai_channel(it)}
        .set {bam_bai}
    
    PREPARE_VCF_FILE(params.dbsnp)

    ADD_READ_GROUP(bam_bai)
    RNASEQ_GATK_SPLITNCIGAR(params.genomes,
                            params.index,
                            params.dict,
                            ADD_READ_GROUP.out)

    RNASEQ_GATK_RECALIBRATE(params.genomes,
                            params.index,
                            params.dict,
                            RNASEQ_GATK_SPLITNCIGAR.out,
                            PREPARE_VCF_FILE.out)
    RNASEQ_CALL_VARIANTS(params.genomes,
                         params.index,
                         params.dict,
                         RNASEQ_GATK_RECALIBRATE.out)
}

def create_bam_bai_channel(LinkedHashMap row) {
    // create meta map

    id = row.sample

    // add path(s) of the bam file(s) to the meta map
    def bam_bai_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam}"
    }
    if (!file(row.bai).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> BAI file does not exist!\n${row.bai}"
    }

    bam_bai_meta = [ id, file(row.bam), file(row.bai)]

    return bam_bai_meta
}