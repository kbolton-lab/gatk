version 1.0

struct Runtime {
	File bwa_override
	File samtools_path
	File gatk_override
	File trimmomatic
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
		String bamPath
		String? pileupPath
		String outputPrefix

		# Calculate pileups after if contamination
		Array[String]? intervals
		Boolean contamination = false
		File? variants_for_contamination
		File? variants_for_contamination_idx

		# runtime
		Runtime runtime_params
	}
	
	command {
	
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_override}
		export SAMTOOLS=~{runtime_params.samtools_path}
		
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" MarkDuplicates \
			-I ~{bam} \
			-O ~{bamPath}/~{outputPrefix}.md.bam \
			-M ~{bamPath}/~{outputPrefix}.md.metrics
			
		$SAMTOOLS index ~{bamPath}/~{outputPrefix}.md.bam
			
		if [[ "~{contamination}" ]]; then
			mkdir -p ~{pileupPath}
			gatk --java-options "-Xmx~{runtime_params.memory}m" GetPileupSummaries \
				-I ~{bamPath}/~{outputPrefix}.md.bam \
				-V ~{variants_for_contamination} \
				-L ~{variants_for_contamination} \
				-O ~{pileupPath}/~{outputPrefix}.pileups.table
		fi
	}

	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File bamMD = "~{bamPath}/~{outputPrefix}.md.bam"
		File metrics = "~{bamPath}/~{outputPrefix}.md.metrics"
	}
		
}


task bqsr {
	input {
		Reference reference
		File bamMD
		String bamPath
		Array[File]+ known_sites
		Array[File]+ known_sites_idx
		Array[String]? intervals = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
																"chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
																"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
																"chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
		String outputPrefix
		
		# runtime
		Runtime runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_override}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" BaseRecalibrator \
			-R ~{reference.ref_fasta} \
			-I ~{bamMD} \
			${sep=' ' prefix("--known-sites ", known_sites)} \
			-O ~{bamPath}/~{outputPrefix}.recal_data.table
	
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File recal_table = "~{bamPath}/~{outputPrefix}.recal_data.table"
	}

}


task apply_bqsr {
	input {
		Reference reference
		File bamMD
		String bamPath
		String outputPrefix
		File bqsr_table
		
		# runtime
		Runtime runtime_params
	}
	
	command {
		set -e -o pipefail
		export GATK_LOCAL_JAR=~{runtime_params.gatk_override}
		
		gatk --java-options "-Xmx~{runtime_params.memory}m" ApplyBQSR \
			-R ~{reference.ref_fasta} \
			-I ~{bamMD} \
			-bqsr ~{bqsr_table} \
			-O ~{bamPath}/~{outputPrefix}.final.bam
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File final = "~{bamPath}/~{outputPrefix}.final.bam"
	}

}
