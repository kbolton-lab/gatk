version 1.0

struct BwaIndex {
	File fastaFile
	Array[File] indexFiles
}

struct InputFastq {
	File R1
	File? R2
}

struct RuntimeBWA {
	File bwa_path
	File trimmomatic_path
	Int max_retries
	Int memory
	Int cpu
}

task Mem {
	input {
		File read1
		File? read2
		Boolean trim
		BwaIndex bwaIndex
		String outputPrefix
		String? readgroup
		
		# runtime
		RuntimeBWA runtime_params
	}
	 

	command {
		set -e -o pipefail
		export BWA=~{runtime_params.bwa_path}
		export TRIMMOMATIC="java -jar ~{runtime_params.trimmomatic_path}"
		
		if [[ "~{trim}" ]]; then
			# when paired-end will need to do nested if 1 or 2
			$TRIMMOMATIC SE -threads ~{runtime_params.cpu} \
				~{read1} \
				/dev/stdout \
				SLIDINGWINDOW:4:20 \
				MINLEN:20 \
			| $BWA mem \
				~{"-t " + runtime_params.cpu} \
				~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
				~{bwaIndex.fastaFile} \
				/dev/stdin \
				> ~{outputPrefix}.sam \
				2> ~{outputPrefix}.log.bwamem
		else
			$BWA mem \
			~{"-t " + runtime_params.cpu} \
			~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
			~{bwaIndex.fastaFile} \
			~{read1} \
			~{read2} \
			> ~{outputPrefix}.sam \
			2> ~{outputPrefix}.log.bwamem
		fi
  }
  
  runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
  
	output {
		File outputSam = "~{outputPrefix}.sam"
	} 
	 
}


