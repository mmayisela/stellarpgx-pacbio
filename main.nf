#!/usr/bin/env nextflow


gene_name = params.gene
d_base = params.db_init
res_base = params.res_init
caller_base = params.caller_init
output_folder = params.out_dir

params.build='hg38'

if (params.build=='b37') {
    db = "${d_base}/${gene_name}/b37"
    res_dir = "${res_base}/${gene_name}/cyp_b37"
    caller_dir = "${caller_base}/${gene_name}/b37/bin"
    debug37 = "--minimum_extract_score_over_homref=0"
    debug38 = ""

    if (params.gene=='cyp2d6') {
        chrom = "22"
    	region_a1 = "22:42522000-42542000"
    	region_a2 = "042522000-042542000"
    	region_b1 = "22:42522300-42528400"
    	region_b2 = "042522300-042528400"

    } else if (params.gene=='cyp2a6') {
        chrom = "19"      
	region_a1 = "19:41339500-41389490"
	region_a2 = "041339500-041389490"
	region_b1 = "19:41348500-41358200"
	region_b2 = "041348500-041358200"

    } else if (params.gene=='cyp2b6') {	
        chrom = "19"
	region_a1 = "19:41487100-41534300"
	region_a2 = "041487100-041534300"
	region_b1 = "19:41494750-41522800"
	region_b2 = "041494750-041522800"

    } else if (params.gene=='cyp2c19') {
        chrom = "10"
	region_a1 = "10:96518000-96613000"
	region_a2 = "096518000-096613000"
	region_b1 = "10:96519000-96612750"
	region_b2 = "096519000-096612750"

    } else if (params.gene=='cyp2c9') {
        chrom = "10"
	region_a1 = "10:96695000-96749500"
	region_a2 = "096695000-096749500"
	region_b1 = "10:96695750-96748950"
	region_b2 = "096695750-096748950"

    } else if (params.gene=='cyp2c8') {
        chrom = "10"
	region_a1 = "10:96786500-96839250"
	region_a2 = "096786500-096839250"
	region_b1 = "10:96796750-96829600"
	region_b2 = "096796750-096829600"

    } else if (params.gene=='cyp2e1') {
        chrom = "10"
	region_a1 = "10:135330300-135352600"
	region_a2 = "135330300-135352600"
	region_b1 = "10:135339500-135351300"
	region_b2 = "135339500-135351300"

    } else if (params.gene=='cyp3a4') {
        chrom = "7"
	region_a1 = "7:99344600-99391800"
	region_a2 = "099344600-099391800"
	region_b1 = "7:99355450-99382600"
	region_b2 = "099355450-099382600"

    } else if (params.gene=='cyp3a5') {
        chrom = "7"
	region_a1 = "7:99235800-99287600"
	region_a2 = "099235800-099287600"
	region_b1 = "7:99245850-99278250"
	region_b2 = "099245850-099278250"

    } else if (params.gene=='cyp1a1') {
        chrom = "15"
	region_a1 = "15:75001900-75027900"
	region_a2 = "075001900-075027900"
	region_b1 = "15:75011600-75015320"
	region_b2 = "075011600-075015320"

    } else if (params.gene=='cyp1a2') {
        chrom = "15"
        region_a1 = "15:75031200-75058900"
        region_a2 = "075031200-075058900"
        region_b1 = "15:75038200-75047650"
        region_b2 = "075038200-075047650"

    } else if (params.gene=='cyp4f2') {
        chrom = "19"
	region_a1 = "19:15978850-16018850"
	region_a2 = "015978850-016018850"
	region_b1 = "19:15988850-16008850"
	region_b2 = "015988850-016008850"

    } else if (params.gene=='cypor') {
       chrom = "7"
       region_a1 = "7:75534473-75626173"
       region_a2 = "075534473-075626173"
       region_b1 = "7:75549473-75621173"
       region_b2 = "075549473-075621173"
    }


} else if (params.build=='hg19') {
    db = "${d_base}/${gene_name}/b37"
    res_dir = "${res_base}/${gene_name}/cyp_hg19"
    caller_dir = "${caller_base}/${gene_name}/b37/bin"
    debug37 = "--minimum_extract_score_over_homref=0"
    debug38 = ""

    if (params.gene=='cyp2d6') {
        chrom = "chr22"
    	region_a1 = "chr22:42522000-42542000"
    	region_a2 = "042522000-042542000"
    	region_b1 = "chr22:42522300-42528400"
    	region_b2 = "042522300-042528400"

    } else if (params.gene=='cyp2a6') {
        chrom = "chr19"
        region_a1 = "chr19:41339500-41389490"
        region_a2 = "041339500-041389490"
        region_b1 = "chr19:41348500-41358200"
        region_b2 = "041348500-041358200"


    } else if (params.gene=='cyp2b6') {
      	chrom = "chr19"
        region_a1 = "chr19:41487100-41534300"
        region_a2 = "041487100-041534300"
	region_b1 = "chr19:41494750-41522800"
	region_b2 = "041494750-041522800"

    } else if (params.gene=='cyp2c19') {
        chrom = "chr10"
         region_a1 = "chr10:96518000-96613000"
         region_a2 = "096518000-096613000"
         region_b1 = "chr10:96519000-96612750"
         region_b2 = "096519000-096612750"

    } else if (params.gene=='cyp2c9') {
        chrom = "chr10"
        region_a1 = "chr10:96695000-96749500"
        region_a2 = "096695000-096749500"
        region_b1 = "chr10:96695750-96748950"
        region_b2 = "096695750-096748950"

    } else if (params.gene=='cyp2c8') {
        chrom = "chr10"
	region_a1 = "chr10:96786500-96839250"
	region_a2 = "096786500-096839250"
	region_b1 = "chr10:96796750-96829600"
	region_b2 = "096796750-096829600"


    } else if (params.gene=='cyp2e1') {
        chrom = "chr10"
        region_a1 = "chr10:135330300-135352600"
	region_a2 = "135330300-135352600"
	region_b1 = "chr10:135339500-135351300"
	region_b2 = "135339500-135351300"

    } else if (params.gene=='cyp3a4') {
      	chrom = "chr7"
	region_a1 = "chr7:99344600-99391800"
	region_a2 = "099344600-099391800"
	region_b1 = "chr7:99355450-99382600"
	region_b2 = "099355450-099382600"


    } else if (params.gene=='cyp3a5') {
      	chrom = "chr7"
	region_a1 = "chr7:99235800-99287600"
	region_a2 = "099235800-099287600"
	region_b1 = "chr7:99245850-99278250"
	region_b2 = "099245850-099278250"

    } else if (params.gene=='cyp1a1') {
      	chrom = "chr15"
        region_a1 = "chr15:75001900-75027900"
	region_a2 = "075001900-075027900"
	region_b1 = "chr15:75011600-75015320"
	region_b2 = "075011600-075015320"


    } else if (params.gene=='cyp1a2') {
        chrom = "chr15"
	region_a1 = "chr15:75031200-75058900"
    	region_a2 = "075031200-075058900"
	region_b1 = "chr15:75038200-75047650"
	region_b2 = "075038200-075047650"

    } else if (params.gene=='cyp4f2') {
        chrom = "chr19"
        region_a1 = "chr19:15978850-16018850"
	region_a2 = "015978850-016018850"
	region_b1 = "chr19:15988850-16008850"
	region_b2 = "015988850-016008850"

    } else if (params.gene=='cypor') {
       chrom = "chr7"
       region_a1 = "chr7:75534473-75626173"
       region_a2 = "075534473-075626173"
       region_b1 = "chr7:75549473-75621173"
       region_b2 = "075549473-075621173"

    }


} else {
    db = "${d_base}/${gene_name}/hg38"
    res_dir = "${res_base}/${gene_name}/cyp_hg38"
    caller_dir = "${caller_base}/${gene_name}/hg38/bin"
    debug37 = ""
    debug38 = "--minimum_extract_score_over_homref=0"

    if (params.gene=='cyp2d6') {
        chrom = "chr22"
    	region_a1 = "chr22:42126000-42137500"
    	region_a2 = "042126000-042137500"
    	region_b1 = "chr22:42126300-42132400"
    	region_b2 = "042126300-042132400" 

    } else if (params.gene=='cyp2a6') {
        chrom = "chr19"
	region_a1 = "chr19:40833541-40887450"
	region_a2 = "040833541-040887450"
	region_b1 = "chr19:40842750-40852250"
	region_b2 = "040842750-040852250"
    
    } else if (params.gene=='cyp2b6') {
        chrom = "chr19"
	region_a1 = "chr19:40981280-41028400"
	region_a2 = "040981280-041028400"
	region_b1 = "chr19:40988800-41016900"
	region_b2 = "040988800-041016900"

    } else if (params.gene=='cyp2c19') {
        chrom = "chr10"
	region_a1 = "chr10:94752750-94865500"
	region_a2 = "094752750-094865500"
	region_b1 = "chr10:94759250-94853000"
	region_b2 = "094759250-094853000"

    } else if (params.gene=='cyp2c9') {
      	chrom = "chr10"
	region_a1 = "chr10:94935000-94990000"
	region_a2 = "094935000-094990000"
	region_b1 = "chr10:94936000-94989200"
	region_b2 = "094936000-094989200"

    } else if (params.gene=='cyp2c8') {
        chrom = "chr10"
	region_a1 = "chr10:95026750-95079500"
	region_a2 = "095026750-095079500"
	region_b1 = "chr10:95037050-95069800"
	region_b2 = "095037050-095069800"	

    } else if (params.gene=='cyp2e1') {
        chrom = "chr10"
	region_a1 = "chr10:133517350-133549100"
	region_a2 = "133517350-133549100"
	region_b1 = "chr10:133526050-133537800"
	region_b2 = "133526050-133537800"

    } else if (params.gene=='cyp3a4') {
        chrom = "chr7"
	region_a1 = "chr7:99746950-99794100"
	region_a2 = "099746950-099794100"
	region_b1 = "chr7:99757850-99784950"
	region_b2 = "099757850-099784950"

    } else if (params.gene=='cyp3a5') {	
        chrom = "chr7"
	region_a1 = "chr7:99638200-99690000"
	region_a2 = "099638200-099690000"
	region_b1 = "chr7:99648250-99680650"	
	region_b2 = "099648250-099680650"

    } else if (params.gene=='cyp1a1') {
        chrom = "chr15"
	region_a1 = "chr15:74709550-74735500"
	region_a2 = "074709550-074735500"
	region_b1 = "chr15:74719250-74725500"
	region_b2 = "074719250-074725500"

    } else if (params.gene=='cyp1a2') {
        chrom = "chr15"
	region_a1 = "chr15:74738850-74766600"
	region_a2 = "074738850-074766600"
	region_b1 = "chr15:74745850-74755280"
	region_b2 = "074745850-074755280"


    } else if (params.gene=='cyp4f2') {
        chrom = "chr19"
	region_a1 = "chr19:15868100-15908000"
	region_a2 = "015868100-015908000"
	region_b1 = "chr19:15878000-15898100"
	region_b2 = "015878000-015898100"

    } else if (params.gene=='cypor') {
       chrom = "chr7"
       region_a1 = "chr7:75905155-75996855"
       region_a2 = "075905155-075996855"
       region_b1 = "chr7:75910155-75991855"
       region_b2 = "075910155-075991855"

    }

}

params.format='binary'

if (params.format=='compressed') {
    ext = "cram"
    ind = "crai"
    cram_options = "--force_use_input_ref_for_cram_reading" 

} else { 
    ext = "bam"
    ind = "bai"
    cram_options = ""

}


align_file = Channel.fromFilePairs(params.in_bam, type: 'file') {  file -> file.name.replaceAll(/.${ext}|.${ind}$/,'') }

align_file.into { data1; data2; data3; data4; data5 }


ref_dir_val = new File("${params.ref_file}").getParent()
ref_genome = new File("${params.ref_file}").getName()



process call_snvs {
  label 'nanocaller'
  
  errorStrategy 'ignore'
  tag "${name}"
  cpus 10

  input:
     set val(name), file(bam) from data1
     path ref_dir from Channel.value("${ref_dir_val}")

  output:
     set val(name), path("${name}_files") into variants_ch
  script:

     """
     NanoCaller --bam ${bam[0]} --ref ${ref_dir}/${ref_genome} --cpu 10 --regions ${region_a1} --prefix ${name}_${gene_name} --sample ${name} --output ${name}_files

     """

}


variants_ch.into {variants_ch1; variants_ch2}


process get_depth {
//   maxForks 10
    label 'global'

    input:
    set val(name), file(bam) from data3
    path ref_dir from Channel.value("${ref_dir_val}")
    path res_dir

    output:
    set val(name), file("${name}_${gene_name}_ctrl.depth") into sv_ch3

    script:

    """   
    samtools bedcov --reference ${ref_dir}/${ref_genome} ${res_dir}/test3.bed ${name}.${ext} > ${name}_${gene_name}_ctrl.depth      

    """

}



// process format_snvs {
// //   maxForks 10


//     // publishDir "$output_folder/$gene_name/variants", pattern: '*vcf.gz', mode: 'copy', overwrite: 'true'
//     // publishDir "$output_folder/$gene_name/variants", pattern: '*vcf.gz.tbi', mode: 'copy', overwrite: 'true'

//     input:
//     set val(name), path(vcf_dir) from variants_ch
//     path caller_base

//     output:
//     set val(name), file("${name}_${gene_name}.vcf.gz"), file("${name}_${gene_name}.vcf.gz.tbi") into (var_norm1, var_norm2)

//     script:

//     """
//         bcftools norm -m - ${vcf} | bcftools view -e 'FILTER="RefCall" & VAF<0.25' > ${name}_${gene_name}_snvs.vcf
// 	python3 ${caller_base}/runtime/check_gt.py ${name}_${gene_name}_snvs.vcf | bcftools view -e 'GT="0/0"' > ${name}_${gene_name}.vcf
// 	bgzip ${name}_${gene_name}.vcf  
//         tabix ${name}_${gene_name}.vcf.gz

//     """

// }


// var_norm1.join(data4).set { prep_phase } 


// process phase_snvs {
//     label 'whatshap'    

// //    publishDir "$output_folder/phased_vcfs", mode: 'copy', overwrite: 'true'

//     input:
//     set val(name), file(unphased_vcf), file(unphased_tbi), file(bam) from prep_phase
//     path ref_dir from Channel.value("${ref_dir_val}")

//     output:
//     set val(name), file(phased_vcf), file(phased_tbi) into (phased_ch1, phased_ch2)

//     script:
//     phased_vcf = "${name}_${gene_name}_phased.vcf.gz"
//     phased_tbi = "${name}_${gene_name}_phased.vcf.gz.tbi"

//     """
//     whatshap phase \
//         --output $phased_vcf \
//         --reference  ${ref_dir}/${ref_genome} \
//         --chromosome $chrom \
//         $unphased_vcf \
//         ${bam[0]} --indels
//     tabix ${name}_${gene_name}_phased.vcf.gz
//     """
// }



process get_core_var {
//   maxForks 10
   
    errorStrategy 'ignore'
    tag "${name}"   

    input:
    set val(name), path("${name}_files") from variants_ch1
    path res_dir

    output:
    set val(name), path("${name}_int") into (core_vars1, core_vars2)

    script:
 
    """
    bcftools isec ${name}_files/${name}_${gene_name}.vcf.gz ${res_dir}/allele_def_var.vcf.gz -p ${name}_int -Oz
    bcftools norm -m - ${name}_int/0002.vcf.gz | bgzip -c > ${name}_int/${name}_${gene_name}_core.vcf.gz
    tabix ${name}_int/${name}_${gene_name}_core.vcf.gz

    """

}



process get_hap_snvs {
//   maxForks 10

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), path("${name}_int") from core_vars1

    output:
    set val(name), path("${name}_${gene_name}_hap1.dip"), path("${name}_${gene_name}_hap2.dip") into prep_ch

    script:
    """
    bcftools query -f'[%POS~%REF>%ALT\t%GT\n]' ${name}_int/${name}_${gene_name}_core.vcf.gz | grep -v '1|0' | cut -f1 > ${name}_${gene_name}_hap1.dip
    bcftools query -f'[%POS~%REF>%ALT\t%GT\n]' ${name}_int/${name}_${gene_name}_core.vcf.gz | grep -v '0|1' | cut -f1 > ${name}_${gene_name}_hap2.dip
    """

}




process query_dup_profile {
//   maxForks 10

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), path("${name}_files") from variants_ch2

    output:
    set val(name), file("${name}_${gene_name}_dup_phased_summary.txt") into dup_ch2

    script:

    """
    bcftools query -f'%POS~%REF>%ALT\t%QUAL\t[%GT\t%DP\t%FQ]\n' -i'GT="alt"' ${name}_files/${name}_${gene_name}.vcf.gz > ${name}_${gene_name}_dup_phased_summary.txt

    """

}


prep_ch.join(sv_ch3).set {fin_files1}
fin_files1.join(dup_ch2).set {fin_files}


process call_stars {
//   maxForks 10

    publishDir "$output_folder/$gene_name/alleles", mode: 'copy', overwrite: 'true'

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), path("${name}_${gene_name}_hap1.dip"), path("${name}_${gene_name}_hap2.dip"), file("${name}_${gene_name}_dp"), file("${name}_gene_dup_phased_summary.txt") from fin_files
    path db
    path caller_dir

    output:
    set val(name), file("${name}_${gene_name}.alleles") into star_ch

    script:
   
    """
    python3 ${caller_dir}/stellarpgx.py ${db}/alleles.dbs ${name}_${gene_name}_hap1.dip ${name}_${gene_name}_hap2.dip ${name}_${gene_name}_dp ${name}_gene_dup_phased_summary.txt ${db}/a_scores.dbs ${name} > ${name}_${gene_name}.alleles  

    """

}

