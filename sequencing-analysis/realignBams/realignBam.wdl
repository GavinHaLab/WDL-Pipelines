version 1.0


## This workflow is used when the genome of an aligned bam file needs to be changed
## This process unmaps the bam, converts to a cram (in case you want to retain this unmapped file in an archive),
## and then reprocesses it for analysis later.  
#import "wdl-repo/processes/aligned-to-unmapped-cram/alignedToUnmappedCram.wdl" as unMap
#import "wdl-repo/processes/WGS-data-preprocessing/align_and_clean.wdl" as preProcess

import "https://fh-pi-paguirigan-a-eco.s3.amazonaws.com/wdl-repo/processes/aligned-to-unmapped-cram/alignedToUnmappedCram.wdl?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIARHEBKX7F3QKMBKPN/20220216/us-west-2/s3/aws4_request&X-Amz-Date=20220216T180604Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=119a9b182de97652031de68bb84f52f362f29c4568adcc098099fc1911a9e148" as unMap
import "https://fh-pi-paguirigan-a-eco.s3.amazonaws.com/wdl-repo/processes/WGS-data-preprocessing/align_and_clean.wdl?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIARHEBKX7F3QKMBKPN/20220216/us-west-2/s3/aws4_request&X-Amz-Date=20220216T180641Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Signature=526f50d6c48fc363b304fa29753d9670ab6f2f90e9ac51df3f8e2fb5679c22ee" as preProcess


workflow realignBam {
  meta {
        description: ""
        author: ""
        email:  ""
        allowNestedInputs: true
  }

  input {
    Array[bamSample] samplesToUnmap
    wgsReferenceData referenceDataSet
  }

  scatter (sample in samplesToUnmap) {
    call unMap.AlignedToUnMappedCram as unMap {
      input:
        samplesToUnmap = [sample]
    }
    ## https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/reprocessing/exome/ExomeReprocessing.wdl#L47
    preprocessSample unmappedData = object {
      sample_name: sample.sample_name,
      dataset_id: sample.dataset_id,
      unmappedCramBam: select_first(unMap.unmappedCram),
      unmappedCraiBai: select_first(unMap.unmappedCramIndex)
    }

    call preProcess.WGS_preprocess_for_variants as preProcess {
      input: 
        batchInputs = [unmappedData],
        referenceDataSet = referenceDataSet
    }
  }
  output {
    Array[Array[File]] analysisReadyBam = preProcess.analysisReadyBam
    Array[Array[File]] analysisReadyIndex = preProcess.analysisReadyIndex
  }
}






