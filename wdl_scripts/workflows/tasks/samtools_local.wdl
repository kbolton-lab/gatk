version 1.0

struct RuntimeSamtools {
	File samtools_path
	Int max_retries
	Int memory
	Int cpu
}

struct bamPlusIndex {
	File bam
	File index
}

task sort {
	input {
		File sam
		String outputPrefix
		
		# runtime
		RuntimeSamtools runtime_params
	}
	
	command {
		set -e -o pipefail

		~{runtime_params.samtools_path} sort ~{sam} -O BAM -o ~{outputPrefix}.bam
		~{runtime_params.samtools_path} index ~{outputPrefix}.bam	
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		# we can zip these for input to Mutect2
		File bam = "~{outputPrefix}.bam"
		File bamIndex = "~{outputPrefix}.bam.bai"
	} 
}
	

task merge {
	input {
		Array[File] bams
		
		# runtime
		RuntimeSamtools runtime_params
	}
	
	command {
		set -e -o pipefail
		
		~{runtime_params.samtools_path} merge -f merged.bam ${sep=' ' bams} 
		~{runtime_params.samtools_path} index merged.bam
		
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File mergedBam = "merged.bam"
		File mergedBamIndex = "merged.bam.bai"
	} 
}


task index {
	input {
		File bam
		
		# runtime
		RuntimeSamtools runtime_params
	}
	
	String sample_basename = basename(bam)

	command {
		set -e -o pipefail

		# have to view it so index is in execution folder, else it goes to inputs folder
		~{runtime_params.samtools_path} view ~{bam} -O BAM -o ~{sample_basename}
		~{runtime_params.samtools_path} index ~{sample_basename}
	}
	
	runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
	
	output {
		File outBam = "~{sample_basename}"
		File indexedBam = "~{sample_basename}.bai"
	} 
}
