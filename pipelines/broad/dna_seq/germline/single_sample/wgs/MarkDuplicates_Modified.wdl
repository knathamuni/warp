version 1.0

import "../../../../../../tasks/broad/BamProcessing.wdl" as Processing
import "../../../../../../tasks/broad/AggregatedBamQC.wdl" as AggregatedQC
import "../../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../../tasks/broad/BamToCram.wdl" as ToCram
import "../../../../../../tasks/broad/Utilities.wdl" as Utilities
import "../../../../../../pipelines/broad/dna_seq/germline/variant_calling/VariantCalling.wdl" as ToGvcf
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"
#import "../../tasks/broad/SplitLargeReadGroup.wdl" as SplitRG

workflow ReprocessFilesWorkflow {
    input {
        SampleAndUnmappedBams sample_and_unmapped_bams
        DNASeqSingleSampleReferences references
        DragmapReference? dragmap_reference
        VariantCallingScatterSettings scatter_settings
        PapiSettings papi_settings
        Array[File] input_files
        String sample_id
        String output_bam_basename
        String metrics_filename
        #Array[Float] sizes
        Int compression_level
        Int preemptible_tries
        String? read_name_regex
        Int memory_multiplier = 1
        Int additional_disk = 20
        File contamination_sites_ud
        File contamination_sites_bed
        File contamination_sites_mu
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String cross_check_fingerprints_by
        File haplotype_database_file
        Float lod_threshold
        String recalibrated_bam_basename
        Boolean hard_clip_reads = false
        Boolean unmap_contaminant_reads = true
        Boolean bin_base_qualities = true
        Boolean somatic = false
        Boolean perform_bqsr = true
        Boolean use_bwa_mem = true
        Boolean allow_empty_ref_alt = false
        #Array[File] sequence_grouping
        #Array[File] sequence_grouping_with_unmapped
        File? fingerprint_genotypes_file
        File? fingerprint_genotypes_index
        File wgs_coverage_interval_list
        Boolean? run_dragen_mode_variant_calling_
        Boolean use_spanning_event_genotyping_
        Boolean use_gatk3_haplotype_caller_ = false
        Boolean use_dragen_hard_filtering_
        Boolean provide_bam_output
        Int additional_disk_size = 20
        String gatk_docker = "broadinstitute/gatk:latest"
        String gatk_path = "/gatk/gatk"

    Float? sorting_collection_size_ratio
        Float gb_size_cutoff_for_preemptibles = 110.0


    }

    call Utilities.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
      ref_dict = references.reference_fasta.ref_dict,
      preemptible_tries = papi_settings.preemptible_tries
  }


    scatter (file in input_files) {
        Float cram_size = size(file, "GiB")
        call UnmarkDuplicates {
            input:
                input_bam = file,
                output_bam_basename = "output_~{basename(file)}",
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                total_input_size = cram_size,
                compression_level = compression_level,
                preemptible_tries = preemptible_tries,
                memory_multiplier = memory_multiplier,
                additional_disk = additional_disk
        }

        call FixSMTag {
            input:
                input_bam = UnmarkDuplicates.output_bam,
                output_bam_basename = "output_~{basename(file)}",
                sample_id = sample_id

        }
    }

    #Float agg_bam_size = size(FixSMTag.output_bam, "GiB")
    #Boolean data_too_large_for_preemptibles = agg_bam_size > gb_size_cutoff_for_preemptibles

    call Utilities.SumFloats as SumFloats {
        input:
            sizes = cram_size,
            preemptible_tries = papi_settings.preemptible_tries
  }

    Boolean data_too_large_for_preemptibles = SumFloats.total_size > gb_size_cutoff_for_preemptibles

    call MarkDuplicates {
        input:
            input_bams = FixSMTag.output_bam,
            output_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.unsorted.duplicates_marked",
            metrics_filename = sample_and_unmapped_bams.base_file_name + ".duplicate_metrics",
            total_input_size = SumFloats.total_size,
            compression_level = compression_level,
            preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
        }

    call Processing.SortSam as SortSampleBam {
            input:
                input_bam = MarkDuplicates.output_bam,
                output_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicate_marked.sorted",
                compression_level = compression_level,
                preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
          }
    #Float agg_bam_size = size(FixSMTag.output_bam, "GiB")
    #Boolean data_too_large_for_preemptibles = size(FixSMTag.output_bam, "GiB") > gb_size_cutoff_for_preemptibles

  if (defined(haplotype_database_file)) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [SortSampleBam.output_bam],
        input_bam_indexes = [SortSampleBam.output_bam_index],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = sample_and_unmapped_bams.base_file_name + ".crosscheck",
        total_input_size = SumFloats.total_size,
        lod_threshold = lod_threshold,
        cross_check_by = cross_check_fingerprints_by,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

    call Processing.CheckContamination as CheckContamination {
      input:
        input_bam = SortSampleBam.output_bam,
        input_bam_index = SortSampleBam.output_bam_index,
        contamination_sites_ud = contamination_sites_ud,
        contamination_sites_bed = contamination_sites_bed,
        contamination_sites_mu = contamination_sites_mu,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        output_prefix = sample_and_unmapped_bams.base_file_name + ".preBqsr",
        preemptible_tries = papi_settings.agg_preemptible_tries,
        contamination_underestimation_factor = 0.75
  }

  # We need disk to localize the sharded input and output due to the scatter for BQSR.
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space
  # to account for the fact that the data is not split evenly.
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel

  if (perform_bqsr) {
    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
      # Generate the recalibration model by interval
      call Processing.BaseRecalibrator as BaseRecalibrator {
        input:
          input_bam = SortSampleBam.output_bam,
          input_bam_index = SortSampleBam.output_bam_index,
          recalibration_report_filename = sample_and_unmapped_bams.base_file_name + ".recal_data.csv",
          sequence_group_interval = subgroup,
          dbsnp_vcf = references.dbsnp_vcf,
          dbsnp_vcf_index = references.dbsnp_vcf_index,
          known_indels_sites_vcfs = references.known_indels_sites_vcfs,
          known_indels_sites_indices = references.known_indels_sites_indices,
          ref_dict = references.reference_fasta.ref_dict,
          ref_fasta = references.reference_fasta.ref_fasta,
          ref_fasta_index = references.reference_fasta.ref_fasta_index,
          bqsr_scatter = bqsr_divisor,
          preemptible_tries = papi_settings.agg_preemptible_tries
      }
    }

    # Merge the recalibration reports resulting from by-interval recalibration
    # The reports are always the same size
    call Processing.GatherBqsrReports as GatherBqsrReports {
      input:
        input_bqsr_reports = BaseRecalibrator.recalibration_report,
        output_report_filename = sample_and_unmapped_bams.base_file_name + ".recal_data.csv",
        preemptible_tries = papi_settings.preemptible_tries
    }

    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
      # Apply the recalibration model by interval
      call Processing.ApplyBQSR as ApplyBQSR {
        input:
          input_bam = SortSampleBam.output_bam,
          input_bam_index = SortSampleBam.output_bam_index,
          output_bam_basename = recalibrated_bam_basename,
          recalibration_report = GatherBqsrReports.output_bqsr_report,
          sequence_group_interval = subgroup,
          ref_dict = references.reference_fasta.ref_dict,
          ref_fasta = references.reference_fasta.ref_fasta,
          ref_fasta_index = references.reference_fasta.ref_fasta_index,
          bqsr_scatter = bqsr_divisor,
          compression_level = compression_level,
          preemptible_tries = papi_settings.agg_preemptible_tries,
          bin_base_qualities = bin_base_qualities,
          somatic = somatic
      }
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call Processing.GatherSortedBamFiles as GatherBamFiles {
    input:
      input_bams = select_first([ApplyBQSR.recalibrated_bam, [SortSampleBam.output_bam]]),
      output_bam_basename = sample_and_unmapped_bams.base_file_name,
      total_input_size = SumFloats.total_size,
      compression_level = compression_level,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }
   call AggregatedQC.AggregatedBamQC as AggregatedBamQC {
    input:
      base_recalibrated_bam = GatherBamFiles.output_bam,
      base_recalibrated_bam_index = GatherBamFiles.output_bam_index,
      base_name = sample_and_unmapped_bams.base_file_name,
      sample_name = sample_and_unmapped_bams.sample_name,
      recalibrated_bam_base_name = recalibrated_bam_basename,
      haplotype_database_file = references.haplotype_database_file,
      references = references,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings
  }

  call ToCram.BamToCram as BamToCram {
    input:
      input_bam = GatherBamFiles.output_bam,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      duplication_metrics = MarkDuplicates.duplicate_metrics,
      chimerism_metrics = AggregatedBamQC.agg_alignment_summary_metrics,
      base_file_name = sample_and_unmapped_bams.base_file_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # QC the sample WGS metrics (stringent thresholds)
  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = sample_and_unmapped_bams.base_file_name + ".wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # QC the sample raw WGS metrics (common thresholds)
  call QC.CollectRawWgsMetrics as CollectRawWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = sample_and_unmapped_bams.base_file_name + ".raw_wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      run_dragen_mode_variant_calling = run_dragen_mode_variant_calling_,
      use_spanning_event_genotyping = use_spanning_event_genotyping_,
      calling_interval_list = references.calling_interval_list,
      evaluation_interval_list = references.evaluation_interval_list,
      haplotype_scatter_count = scatter_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = scatter_settings.break_bands_at_multiples_of,
      contamination = CheckContamination.contamination,
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      ref_str = references.reference_fasta.ref_str,
      dbsnp_vcf = references.dbsnp_vcf,
      dbsnp_vcf_index = references.dbsnp_vcf_index,
      base_file_name = sample_and_unmapped_bams.base_file_name,
      final_vcf_base_name = select_first([sample_and_unmapped_bams.final_gvcf_base_name]),
      agg_preemptible_tries = papi_settings.agg_preemptible_tries,
      use_gatk3_haplotype_caller = use_gatk3_haplotype_caller_,
      use_dragen_hard_filtering = use_dragen_hard_filtering_
  }

  if (provide_bam_output) {
    File provided_output_bam = GatherBamFiles.output_bam
    File provided_output_bam_index = GatherBamFiles.output_bam_index
  }



  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File? output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File output_bam = GatherBamFiles.output_bam
    File output_bam_index = GatherBamFiles.output_bam_index

    #Array[File] quality_yield_metrics = GatherBamFiles.quality_yield_metrics

   #Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = GatherBamFiles.unsorted_read_group_base_distribution_by_cycle_pdf
    #Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = GatherBamFiles.unsorted_read_group_base_distribution_by_cycle_metrics
    #Array[File] unsorted_read_group_insert_size_histogram_pdf = GatherBamFiles.unsorted_read_group_insert_size_histogram_pdf
    #Array[File] unsorted_read_group_insert_size_metrics = GatherBamFiles.unsorted_read_group_insert_size_metrics
    #Array[File] unsorted_read_group_quality_by_cycle_pdf = GatherBamFiles.unsorted_read_group_quality_by_cycle_pdf
   # Array[File] unsorted_read_group_quality_by_cycle_metrics = GatherBamFiles.unsorted_read_group_quality_by_cycle_metrics
   # Array[File] unsorted_read_group_quality_distribution_pdf = GatherBamFiles.unsorted_read_group_quality_distribution_pdf
   # Array[File] unsorted_read_group_quality_distribution_metrics = GatherBamFiles.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = AggregatedBamQC.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = AggregatedBamQC.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = AggregatedBamQC.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = AggregatedBamQC.read_group_gc_bias_summary_metrics

    File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = AggregatedBamQC.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = AggregatedBamQC.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = AggregatedBamQC.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = AggregatedBamQC.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = AggregatedBamQC.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
    File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

    File? fingerprint_summary_metrics = AggregatedBamQC.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = AggregatedBamQC.fingerprint_detail_metrics

    File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics

    File gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    File gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics


    File output_cram = BamToCram.output_cram
    File output_cram_index = BamToCram.output_cram_index
    File output_cram_md5 = BamToCram.output_cram_md5

    File validate_cram_file_report = BamToCram.validate_cram_file_report

    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
 }


task UnmarkDuplicates {
    input {
        File input_bam
        String output_bam_basename
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Float total_input_size
        Int compression_level
        Int preemptible_tries
        Int memory_multiplier = 1
        Int additional_disk = 20
    }
    Float md_disk_multiplier = 3
    #Float total_input_size = size(input_files, 'GiB')
    Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk
    Float memory_size = 7.5 * memory_multiplier
    Int java_memory_size = (ceil(memory_size) - 2)

    command {
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms5000m -Xmx5500m" \
            UnmarkDuplicates \
            -I ~{input_bam} \
            -O ~{output_bam_basename}.bam \
            -R ~{ref_fasta} \
            --sequence-dictionary ~{ref_dict}
    }
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        preemptible: preemptible_tries
        memory: "~{memory_size} GiB"
        disks: "local-disk " + disk_size + " HDD"
    }
    output {
        File output_bam = "~{output_bam_basename}.bam"
    }
}

task FixSMTag{
    input {
        File input_bam
        String sample_id
        String output_bam_basename
        #Runtime parameters
        Int disk_size = 80
        Int machine_mem_gb = 2
        Int preemptible_attempts = 3
    }

    Int command_mem_gb = machine_mem_gb - 1

    command {
        samtools view -H ~{input_bam} | sed 's/^\(@RG.*\tSM:\)\([^[:space:]]*\)/\1~{sample_id}/' | sed 's/^\(@RG.*\tLB:\)\([^[:space:]]*\)/\1~{sample_id}/' > new_header.sam
        samtools reheader new_header.sam ~{input_bam} > ~{output_bam_basename}.reheadered.bam
    }
    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
        disks: "local-disk " + disk_size + " HDD"
        memory: machine_mem_gb + " GB"
        preemptible: preemptible_attempts
    }
    output {
        File output_bam = "~{output_bam_basename}.reheadered.bam"
    }
}


task MarkDuplicates {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    Float total_input_size
    Int compression_level
    Int preemptible_tries

    String? read_name_regex
    Int memory_multiplier = 1
    Int additional_disk = 20

    Float? sorting_collection_size_ratio
  }

  Float md_disk_multiplier = 3
  Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk

  Float memory_size = 7.5 * memory_multiplier
  Int java_memory_size = (ceil(memory_size) - 2)


  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar /usr/picard/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="coordinate" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}