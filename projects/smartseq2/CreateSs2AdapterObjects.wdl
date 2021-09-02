version 1.0

import "../hca_mvp/tasks/AdapterTasks.wdl" as Tasks

workflow CreateSs2AdapterObjects {
  meta {
    description: "Creates json objects for indexing HCA 10x smartseq2 analysis data"
    allowNestedInputs: true
  }

  input {
    File? bam
    File? bai
    File? loom
    Int? ss2_index
    Array[String] process_input_ids # Array of space seperated strings...fastq for intermediate, intermediate looms for project level
    String input_id
    String project_id
    String version_timestamp
    String pipeline_type
    Boolean is_project_level
    String pipeline_version # parsed from metadata for intermediate, passed in for project level
    String reference_file_fasta # parsed from metadata for intermediate, passed in for project level
    File metadata
  }

  call Tasks.GetAnalysisFileMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      ss2_bam_file = select_first([bam]),
      ss2_bai_file = select_first([bai]),
      version_timestamp = version_timestamp,
      project_level = is_project_level,
      project_loom = loom
  }

  call Tasks.GetAnalysisProcessMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      references = reference_file_fasta,
      project_level = is_project_level,
      input_file = metadata,
      ss2_index = ss2_index
  }

  # pipeline_type is used for a dockstore URL here, so it needs to fit into this example:
  # "computational_method": "https://dockstore.org/workflows/github.com/broadinstitute/warp/Smartseq2_Multisample:MultiSampleSmartSeq2_v2.1.4",
  call Tasks.GetAnalysisProtocolMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = "Smartseq2_Multisample",
      version_timestamp = version_timestamp,
      project_level = is_project_level,
      pipeline_version = pipeline_version
  }

  if(defined(loom)) {
    call Tasks.GetCloudFileCreationDate as GetLoomFileCreationDate {
      input:
        file_path = select_first([loom])
    }

    call Tasks.GetFileDescriptor as GetLoomFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = select_first([loom]),
        file_path_string = select_first([loom]),
        input_uuid = input_id,
        creation_time = GetLoomFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }
  }

  # bam
  if (defined(bam)){
    call Tasks.GetCloudFileCreationDate  as GetBamFileCreationDate {
      input:
        file_path = select_first([bam])
    }

    call Tasks.GetFileDescriptor as GetBamFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = select_first([bam]),
        file_path_string = select_first([bam]),
        input_uuid = input_id,
        creation_time = GetBamFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }
  }

  # bai
  if (defined(bai)){
    call Tasks.GetCloudFileCreationDate  as GetBaiFileCreationDate {
    input:
      file_path = select_first([bai])
    }

    call Tasks.GetFileDescriptor as GetBaiFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = select_first([bai]),
        file_path_string = select_first([bai]),
        input_uuid = input_id,
        creation_time = GetBaiFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }
  }

# TODO: create one large links file for ss2
  # call Tasks.GetLinksFileMetadata {
  #   input:
  #     project_id = project_id,
  #     process_input_ids = process_input_ids, # for intermediate level use fastq_uuids from Terra, for project level use output_ids from intermediate files
  #     output_file_path = GetAnalysisFileMetadata.outputs_json,
  #     version_timestamp = version_timestamp,
  #     analysis_process_path = GetAnalysisProcessMetadata.analysis_process_outputs,
  #     analysis_protocol_path = GetAnalysisProtocolMetadata.analysis_protocol_outputs,
  #     project_level = is_project_level,
  #     file_name_string = input_id
  # }

  output {
    Array[File] analysis_file_outputs = GetAnalysisFileMetadata.analysis_file_outputs
    Array[File] analysis_process_outputs = GetAnalysisProcessMetadata.analysis_process_outputs
    Array[File] analysis_protocol_outputs = GetAnalysisProtocolMetadata.analysis_protocol_outputs
    # Array[File] links_outputs = GetLinksFileMetadata.links_outputs
    Array[File]? loom_file_descriptor_outputs = GetLoomFileDescriptor.file_descriptor_outputs
    Array[File]? bam_file_descriptor_outputs = GetBamFileDescriptor.file_descriptor_outputs
    Array[File]? bai_file_descriptor_outputs = GetBaiFileDescriptor.file_descriptor_outputs
  }
}

