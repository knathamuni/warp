version 1.0

workflow CollectMetricsonFullVcf {
  input {
    File final_gather_vcf
    File final_gather_vcf_index
    String callset_name
    File dbsnp_vcf
    File dbsnp_vcf_index
    File eval_interval_list
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int large_disk
  }

    call CollectMetricsOnFullVcf {
      input:
        input_vcf = final_gather_vcf,
        input_vcf_index = final_gather_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size_gb = large_disk
    }

  }


task CollectMetricsOnFullVcf {

  input {
    File input_vcf
    File input_vcf_index
    String metrics_filename_prefix
    File dbsnp_vcf
    File dbsnp_vcf_index
    File interval_list
    File ref_dict
    Int disk_size_gb
    Int machine_mem_mb = 7500
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xms6000m -Xmx7000m" \
      CollectVariantCallingMetrics \
      --INPUT ~{input_vcf} \
      --DBSNP ~{dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ~{ref_dict} \
      --OUTPUT ~{metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ~{interval_list}
  >>>

  output {
    File detail_metrics_file = "~{metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "~{metrics_filename_prefix}.variant_calling_summary_metrics"
  }

  runtime {
    memory: "~{machine_mem_mb} MiB"
    cpu: 2
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size_gb + " HDD"
    preemptible: 1
    docker: gatk_docker
  }
}