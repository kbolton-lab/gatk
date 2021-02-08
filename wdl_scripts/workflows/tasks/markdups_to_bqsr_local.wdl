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


task markdupsIndiv {
	input {
		File bam
		String outputPrefix

		# runtime
		RuntimeGATK runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}

		gatk --java-options "-Xmx~{runtime_params.memory}m" MarkDuplicates \
			-I ~{bam} \
			-O ~{outputPrefix}.md.bam \
			-M ~{outputPrefix}.md.metrics
	}

	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File bamMD = "~{outputPrefix}.md.bam"
		File metrics = "~{outputPrefix}.md.metrics"
	}	
}


task pileups {
	input {
		File bam
		File bamIndex
		String outputPrefix
		Reference reference
		

		# Calculate pileups after if contamination
		Array[String]? intervals # need to change command if used
		File? variants_for_contamination
		File? variants_for_contamination_idx

		# runtime
		RuntimeGATK runtime_params
	}

	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" GetPileupSummaries \
			-R ~{reference.ref_fasta} \
			-I ~{bam} \
			-V ~{variants_for_contamination} \
			-L ~{variants_for_contamination} \
			-O ~{outputPrefix}.pileups.table
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File pileup = "~{outputPrefix}.pileups.table"
	}
}


task calcContam {
	input {
		File pileup
		String outputPrefix
		
		# runtime
		RuntimeGATK runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" CalculateContamination \
			-I ~{pileup} \
			-O ~{outputPrefix}.contamination.table \
			--tumor-segmentation ~{outputPrefix}.segments.table
			
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File contamination = "~{outputPrefix}.contamination.table"
		File segments = "~{outputPrefix}.segments.table"
	}
}


task bqsr {
	input {
		Reference reference
		File bamMD
		File bamMDIndex
		Array[File]+ known_sites
		Array[File]+ known_sites_idx
		Array[String] intervals = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
																"chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
																"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
																"chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
		String outputPrefix
		
		# runtime
		RuntimeGATK runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" BaseRecalibrator \
			-R ~{reference.ref_fasta} \
			-I ~{bamMD} \
			${sep=' ' prefix("--known-sites ", known_sites)} \
			-O ~{outputPrefix}.recal_data.table
	
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File recal_table = "~{outputPrefix}.recal_data.table"
	}
}


task apply_bqsr {
	input {
		Reference reference
		File bamMD
		File bamMDIndex
		String outputPrefix
		File bqsr_table
		
		# runtime
		RuntimeGATK runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_path}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" ApplyBQSR \
			-R ~{reference.ref_fasta} \
			-I ~{bamMD} \
			-bqsr ~{bqsr_table} \
			-O ~{outputPrefix}.bam
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File final = "~{outputPrefix}.bam"
	}
}
