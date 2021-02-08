version 1.0

import "tasks/bwa_local.wdl" as bwa_task
import "tasks/samtools_local.wdl" as samtools
import "tasks/markdups_to_bqsr_local.wdl" as md_bqsr

struct bamPlusIndex {
	File bam
	File index
}

workflow AlignBwaMem {
	input {
		Array[Map[String, String]]+ samples
		Boolean trim
		BwaIndex bwaIndex
		Reference reference
		
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
		File bwa_path
		File trimmomatic_path
		File samtools_path
		File gatk_path
		Int maxRetries = 2
		Int memory = 4
		Int threads = 2
	}
	
	Int max_retries_or_default = select_first([maxRetries, 2])

	RuntimeBWA standard_runtime_bwa = {"bwa_path": bwa_path,
																		 "trimmomatic_path": trimmomatic_path,
																		 "max_retries": max_retries_or_default,
																		 "memory": memory * 1000,
																		 "cpu": threads}
	
	RuntimeSamtools standard_runtime_samtools = {"samtools_path": samtools_path,
																							 "max_retries": max_retries_or_default,
																							 "memory": memory * 1000,
																							 "cpu": threads}			
	
	RuntimeGATK standard_runtime_gatk = {"gatk_path": gatk_path,
																			 "max_retries": max_retries_or_default,
																			 "memory": memory * 1000,
																			 "cpu": threads}											
													


	scatter (s in samples) {
	
		call bwa_task.Mem as bwa {
			input :
				trim = trim,
				read1 = s.read1,
				#read2 = s.read2,
				bwaIndex = bwaIndex,
				outputPrefix = s.outputPrefix,
				readgroup = s.readgroup,
				runtime_params = standard_runtime_bwa	
		}
		
		
		call samtools.sort as samsort {
			input :
				sam = bwa.outputSam,
				outputPrefix = s.outputPrefix,
				runtime_params = standard_runtime_samtools
		}
		
		
		call md_bqsr.markdupsIndiv as md {
			input :
				bam = samsort.bam,
				outputPrefix = s.outputPrefix,
				runtime_params = standard_runtime_gatk
		}
		
		# sort and index after MarkDuplicates
		call samtools.index as samindex {
			input :
				bam = md.bamMD,
				runtime_params = standard_runtime_samtools
		}
		
		
		## Probably better to do this after Mutect2 call for tumor samples
		#if (contamination) {
		#	call md_bqsr.pileups as pu {
		#		input:
		#			bam = samindex.outBam,
		#			bamIndex = samindex.indexedBam,
		#			outputPrefix = s.outputPrefix,
		#			reference = reference,
		#			variants_for_contamination = variants_for_contamination,
		#			variants_for_contamination_idx = variants_for_contamination_idx,		
		#			runtime_params = standard_runtime_gatk
		#	}
			
		#	call md_bqsr.calcContam as cc {
		#		input:
		#			pileup = pu.pileup,
		#			outputPrefix = s.outputPrefix,
		#			runtime_params = standard_runtime_gatk
		#	}
		#}
		

		call md_bqsr.bqsr as bqsr {
			input :
				reference = reference,
				bamMD = samindex.outBam,
				bamMDIndex = samindex.indexedBam,
				outputPrefix = s.outputPrefix	,
				known_sites = select_all([gnomad, mills, known_indels]),
				known_sites_idx = select_all([gnomad_idx, mills_idx, known_indels_idx]),
				runtime_params = standard_runtime_gatk
		}
		
		
		call md_bqsr.apply_bqsr as apply_bqsr {
			input :
				reference = reference,
				bamMD = samindex.outBam,
				bamMDIndex = samindex.indexedBam,
				bqsr_table = bqsr.recal_table,
				outputPrefix = s.outputPrefix,
				runtime_params = standard_runtime_gatk
		}
		
		# sort and index after BQSR
		call samtools.sort as samsort2 {
			input :
				sam = apply_bqsr.final,
				outputPrefix = s.outputPrefix,
				runtime_params = standard_runtime_samtools
		}
		
	}
	
	# Merge all bams post BQSR
	call samtools.merge {
		input : 
			bams = samsort2.bam,
			runtime_params = standard_runtime_samtools	
	}
	
	
	output {
		Array[File] bam = samsort2.bam
		Array[File] bamIndex = samsort2.bamIndex
	}
	
}


