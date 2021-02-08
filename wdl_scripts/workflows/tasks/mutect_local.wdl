version 1.0

struct RuntimeGATK {
	File gatk_path
	Int max_retries
	Int memory
	Int cpu
}

struct Reference {
	File ref_fasta
	File ref_fai
  File ref_dict
}


task Mutect2PON {
	input {
		Reference reference
		File normBam
		File normBamIdx
		String sampleName
		File germline	
		File germlineIdx
		
		# Runtime
		RuntimeGATK runtime_params	
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		# Note that as of May, 2019 -max-mnp-distance must be set to zero to avoid a bug in GenomicsDBImport. 
		gatk --java-options "-Xmx~{runtime_params.memory}m" Mutect2 \
			-R ~{reference.ref_fasta} \
			-I ~{normBam} \
			-tumor ~{sampleName} \
			-max-mnp-distance 0 \
			--germline-resource ~{germline} \
			-O ~{sampleName}_for_pon.vcf.gz
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File ponVCF = "~{sampleName}_for_pon.vcf.gz"
		File ponVCFIdx = "~{sampleName}_for_pon.vcf.gz.tbi"
	}
}


task PONDB {
	input {
		Array[File] ponVCFs
		Array[File] ponVCFIdxs
		Reference reference
		Array[String]+ intervals = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
																"chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
																"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
																"chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
		String panelDBName
		
		# Runtime
		RuntimeGATK runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" GenomicsDBImport \
			-R ~{reference.ref_fasta} \
			${sep=' ' prefix("-V ", ponVCFs)} \
			${sep=' ' prefix("-L ", intervals)} \
			--genomicsdb-workspace-path ~{panelDBName}
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File ponDB = "~{panelDBName}"
	}
}


task PON {
	input {
		Reference reference
		String ponDB
		String panelName
		
		# Runtime
		RuntimeGATK runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" CreateSomaticPanelOfNormals \
			-R ~{reference.ref_fasta} \
			-V gendb://~{ponDB} \
			-O ~{panelName}.vcf.gz
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File pon = "~{panelName}.vcf.gz"
	}
}
		

task Mutect2 {
	input {
		Reference reference
		File tumorBam
		File tumorBamIdx
		String sampleName
		File germline	
		File germlineIdx
		File pon
		File ponIdx
		
		# Runtime
		RuntimeGATK runtime_params	
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		# Note that as of May, 2019 -max-mnp-distance must be set to zero to avoid a bug in GenomicsDBImport. 
		gatk --java-options "-Xmx~{runtime_params.memory}m" Mutect2 \
			-R ~{reference.ref_fasta} \
			-I ~{tumorBam} \
			-tumor ~{sampleName} \
			-pon ~{pon} \
			--f1r2-tar-gz ~{sampleName}_f1r2.tar.gz \
			--germline-resource ~{germline} \
			-O ~{sampleName}.vcf.gz
			
		gatk --java-options "-Xmx~{runtime_params.memory}m" LearnReadOrientationModel \
			-I ~{sampleName}_f1r2.tar.gz \
			-O ~{sampleName}_read-orientation-model.tar.gz
			
		gatk --java-options "-Xmx~{runtime_params.memory}m" GetPileupSummaries \
			-I ~{tumorBam} \
			-V ~{germline} \
			-L ~{germline} \
			-O ~{sampleName}.pileups.table
			
		gatk --java-options "-Xmx~{runtime_params.memory}m" CalculateContamination \
			-I ~{sampleName}.pileups.table \
			--tumor-segmentation ~{sampleName}.segments.table \
			-O ~{sampleName}.contamination.table
		
		# FilterMutectCalls now requires a reference, which should be the same fasta file input to Mutect2.
		# Stats file is auto detected in folder
		gatk --java-options "-Xmx~{runtime_params.memory}m" FilterMutectCalls \
			-R ~{reference.ref_fasta} \
			-V ~{sampleName}.vcf.gz \
			--tumor-segmentation ~{sampleName}.segments.table \
			--contamination-table ~{sampleName}.contamination.table \
			-ob-priors ~{sampleName}_read-orientation-model.tar.gz \
			-O ~{sampleName}.filtered.vcf.gz
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File filteredVCF = "~{sampleName}.filtered.vcf.gz"
		File filteredVCFIdx = "~{sampleName}.filtered.vcf.gz.tbi"
	}
}
	
	
	
	
	
