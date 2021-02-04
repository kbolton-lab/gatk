version 1.0

struct BwaIndex {
	File fastaFile
	Array[File] indexFiles
}

struct InputFastq {
	File R1
	File? R2
}

struct Runtime {
   File bwa_override	
   File samtools_path
   File gatk_override
   File trimmomatic
   Int max_retries
   Int memory
   Int cpu
}

task Mem {
	input {
		#InputFastq inputFastq
		File read1
		File? read2
		Boolean trim
		BwaIndex bwaIndex
		String outputPath
		String outputPrefix
		String? readgroup
		
		# runtime
		Runtime runtime_params

	}


	command {
		set -e -o pipefail
		export BWA=~{runtime_params.bwa_override}
		export TRIMMOMATIC=~{runtime_params.trimmomatic}
		
		mkdir -p ~{outputPath}
		
		if [[ "~{trim}" ]]; then
			# when pair-end will need to do nested if 1 or 2
			java -jar $TRIMMOMATIC SE -threads ~{runtime_params.cpu} \
				~{read1} \
				/dev/stdout \
				SLIDINGWINDOW:4:20 \
				MINLEN:20 \
			| $BWA mem \
				~{"-t " + runtime_params.cpu} \
				~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
				~{bwaIndex.fastaFile} \
				/dev/stdin \
				> ~{outputPath}/~{outputPrefix}.sam \
				2> ~{outputPrefix}.log.bwamem
		else
			$BWA mem \
			~{"-t " + runtime_params.cpu} \
			~{"-R '" + readgroup}~{true="'" false="" defined(readgroup)} \
			~{bwaIndex.fastaFile} \
			~{read1} \
			~{read2} \
			> ~{outputPath}/~{outputPrefix}.sam \
			2> ~{outputPrefix}.log.bwamem
		fi
  }
  
  runtime {
		maxRetries: runtime_params.max_retries
		memory: runtime_params.memory + " MB"
		cpu: runtime_params.cpu
	}
  
	output {
		File outputSam = "~{outputPath}/~{outputPrefix}.sam"
	} 
	 
}


