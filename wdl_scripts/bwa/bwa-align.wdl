version 1.0

import "tasks/bwa.wdl" as bwa_task
import "tasks/samtools.wdl" as samtools
import "tasks/markdups_to_bqsr.wdl" as md_bqsr

struct Runtime {
   File bwa_override	
   File samtools_path
   File gatk_override
   File trimmomatic
   Int max_retries
   Int memory
   Int cpu
}

workflow AlignBwaMem {
	input {
		#InputFastq inputFastq
		Array[Map[String, String]]+ samples
		
		Boolean trim
		BwaIndex bwaIndex
		Reference reference
		String outputPath
		String bamPath
		String pileupPath
		
		# Calculate contamination using pileups, currently gnomad
		Boolean contamination = false
		File? variants_for_contamination
		File? variants_for_contamination_idx
		
		# BQSR
		Array[String]? intervals
		File gnomad # need at least 1 known-sites file
		File gnomad_idx
		File? mills
		File? mills_idx
		File? known_indels
		File? known_indels_idx
		
		
		# runtime
		File bwa_override
		File samtools_path
		File gatk_override
		File trimmomatic
		Int maxRetries = 2
		Int memory = 4
		Int threads = 2
	}
	
	Int max_retries_or_default = select_first([maxRetries, 2])

	Runtime standard_runtime = {"bwa_override": bwa_override,
															"samtools_path": samtools_path,
															"gatk_override": gatk_override,
															"trimmomatic" : trimmomatic,
															"max_retries": max_retries_or_default,
															"memory": memory * 1000,
															"cpu": threads}


	scatter(s in samples) {
	
		
		call bwa_task.Mem as bwa {
			input :
				#inputFastq = inputFastq,
				trim = trim,
				read1 = s.read1,
				#read2 = s.read2,
				bwaIndex = bwaIndex,
				outputPath = outputPath,
				outputPrefix = s.outputPrefix,
				readgroup = s.readgroup,
				runtime_params = standard_runtime	
		}
		
		
		call samtools.sort as samsort {
			input :
				sam = bwa.outputSam,
				bamPath = bamPath,
				outputPrefix = s.outputPrefix,
				runtime_params = standard_runtime
		}
		
		
		call md_bqsr.markdupsIndiv as md {
			input :
				bam = samsort.outputBam,
				bamPath = bamPath,
				pileupPath = pileupPath,
				outputPrefix = s.outputPrefix,
				contamination = contamination,
				variants_for_contamination = variants_for_contamination,
				variants_for_contamination_idx = variants_for_contamination_idx,		
				runtime_params = standard_runtime
		}
		

		call md_bqsr.bqsr as bqsr {
			input :
				reference = reference,
				bamMD = md.bamMD,
				bamPath = bamPath,
				outputPrefix = s.outputPrefix	,
				known_sites = select_all([gnomad, mills, known_indels]),
				known_sites_idx = select_all([gnomad_idx, mills_idx, known_indels_idx]),
				runtime_params = standard_runtime
		}
		
		
		call md_bqsr.apply_bqsr as apply_bqsr {
			input :
				reference = reference,
				bamMD = md.bamMD,
				bamPath = bamPath,
				bqsr_table = bqsr.recal_table,
				outputPrefix = s.outputPrefix,
				runtime_params = standard_runtime
		}
		
		
		call samtools.sort as samsort2 {
			input :
				sam = apply_bqsr.final,
				bamPath = bamPath,
				outputPrefix = s.outputPrefix + ".final",
				runtime_params = standard_runtime
		}
		
	}
	
	
	call samtools.merge {
		input : 
			input_sams = samsort.outputBam,
			bamPath = bamPath,
			runtime_params = standard_runtime	
	}
	
	output {
		File mergedBam = "~{bamPath}/merged.bam"
	}
	
}
