version 1.0

import "bwa-align-local.wdl" as align
import "tasks/mutect_local.wdl" as mutect2

workflow PON {
	input {
		Array[Map[String, String]]+ samples
		Boolean trim
		BwaIndex bwaIndex
		Reference reference
		String panelName
		
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
	
	RuntimeGATK standard_runtime_gatk = {"gatk_path": gatk_path,
																			 "max_retries": max_retries_or_default,
																			 "memory": memory * 1000,
																			 "cpu": threads}	

	call align.AlignBwaMem {
		input :
			samples = samples,
			trim = trim,
			bwaIndex = bwaIndex,
			reference = reference,
			
			# Calculate contamination using pileups, currently gnomad
			contamination = contamination,
			variants_for_contamination = variants_for_contamination,
			variants_for_contamination_idx = variants_for_contamination_idx,
			
			# BQSR
			gnomad = gnomad,
			gnomad_idx = gnomad_idx,
			mills = mills,
			mills_idx = mills_idx,
			known_indels = known_indels,
			known_indels_idx = known_indels_idx,
			
			# runtime
			bwa_path = bwa_path,
			trimmomatic_path = trimmomatic_path,
			samtools_path = samtools_path,
			gatk_path = gatk_path,
			maxRetries = maxRetries,
			memory = memory,
			threads = threads,
		}
		
		#output {
		#	Array[File] bam = AlignBwaMem.bam
		#	Array[File] bamIndex = AlignBwaMem.bamIndex
		#}
		
		scatter (bamIdx in zip(AlignBwaMem.bam, AlignBwaMem.bamIndex)) {
			call mutect2.Mutect2PON as M2PON {
				input:
					reference = reference,
					#normBam = bamIdx.bam,
					normBam = bamIdx.left,
					#normBamIdx = bamIdx.bamIndex,
					normBamIdx = bamIdx.right,
					#sampleName = basename(bamIdx.bam, ".bam"),
					sampleName = basename(bamIdx.left, ".bam"),
					germline = gnomad,
					germlineIdx = gnomad_idx,
					runtime_params = standard_runtime_gatk
			}
		}
		
	call mutect2.PONDB as PONDB {
		input:
			reference = reference,
			ponVCFs = M2PON.ponVCF,
			ponVCFIdxs = M2PON.ponVCFIdx,
			intervals = intervals,
			panelDBName = panelName,
			runtime_params = standard_runtime_gatk
	}
	
	call mutect2.PON as PON {
		input:
			reference = reference,
			ponDB = PONDB.ponDB,
			panelName = panelName,
			runtime_params = standard_runtime_gatk
	}
		
	output {
		File ponVCF = PON.pon
	}

}


