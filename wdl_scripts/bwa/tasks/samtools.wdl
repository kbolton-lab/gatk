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

task sort {
	input {
		File sam
		String bamPath	
		String outputPrefix
		# runtime
		Runtime runtime_params
	}
	
	command {
		set -e -o pipefail
		export SAMTOOLS=~{runtime_params.samtools_path}
		
		mkdir -p ~{bamPath}
		$SAMTOOLS sort ~{sam} -O BAM -o ~{bamPath}/~{outputPrefix}.bam
		$SAMTOOLS index ~{bamPath}/~{outputPrefix}.bam
		
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File outputBam = "~{bamPath}/~{outputPrefix}.bam"
	} 

}
	

task merge {
	input {
		Array[File] input_sams
		String bamPath	
		
		# runtime
		Runtime runtime_params
	}
	
	command {
		set -e -o pipefail
		export SAMTOOLS=~{runtime_params.samtools_path}

		$SAMTOOLS merge -f -r ~{bamPath}/merged.bam ${sep=' ' input_sams} 
		$SAMTOOLS index ~{bamPath}/merged.bam
		
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File mergedBam = "~{bamPath}/merged.bam"
	} 
	
}
