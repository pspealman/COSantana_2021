#Loading modules
module load biom/intel/2.1.5

#starting qiime
singularity shell core_2019.4.sif

#Joining reads
	#scrits_folder contains our custom python scripts
	scripts_folder=/scratch/ps163/Dr_Carolina/scripts/
	#fastq_folder contains forward and reverse reads:
	fastq_folder=/scratch/ps163/Dr_Carolina/UFBA_Carolina_Jan_2019-115734619/
	#vsearch_folder contains the results of qiime2 vsearch join
	vsearch_folder=/scratch/ps163/Dr_Carolina/qiime2/vsearch_join_fastq/
	#joined_folder contains the completed joined files
	joined_folder=/scratch/ps163/Dr_Carolina/qiime2/joined_fastq/
	
	
	
	#Import reads for P1-1:
	sample_file=P1_1-227388221
	new_sample_name=P1-1
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
			
	#Import reads for P1-2:
	sample_file=P1_2-227388216
	new_sample_name=P1-2
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
			
	#Import reads import reads for P1-3:
	sample_file=P1_3-227388220
	new_sample_name=P1-3
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
		
#Replicates Site 2
	#Import reads for P2-1:
	sample_file=P2_1-227388222
	new_sample_name=P2-1
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
			
	#Import reads for P2-2:
	sample_file=P2_2-227388214
	new_sample_name=P2-2
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
			
	#Import reads import reads for P2-3:
	sample_file=P2_3-227388217
	new_sample_name=P2-3
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
#
#Replicates Site 3
	#Import reads for P3-1:
	sample_file=P3_1-227388218
	new_sample_name=P3-1
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
			
	#Import reads for P3-2:
	sample_file=P3_2-2273882194
	new_sample_name=P3-2
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
			
	#Import reads import reads for P2-3:
	sample_file=P3_3-227388215
	new_sample_name=P3-3
	#	#This is the same each time
		echo "Starting on ${sample_file}..."
		
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${fastq_folder}/${sample_file} \
			--input-format CasavaOneEightSingleLanePerSampleDirFmt \
			--output-path ${vsearch_folder}/${new_sample_name}.qza
		
		qiime vsearch join-pairs \
			--i-demultiplexed-seqs ${vsearch_folder}/${new_sample_name}.qza \
			--o-joined-sequences ${vsearch_folder}/${new_sample_name}_joined.qza
			
		qiime tools export --input-path ${vsearch_folder}/${new_sample_name}_joined.qza --output-path ${vsearch_folder}/exported-artifact
		
		gunzip -f ${vsearch_folder}/exported-artifact/*.gz
		gunzip -f ${fastq_folder}/${sample_file}/*.gz
			
		python ${scripts_folder}/juntar.py \
			--fastq_1 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R1_001.fastq \
			--fastq_2 ${fastq_folder}/${sample_file}/${new_sample_name}_0_L001_R2_001.fastq \
			--vsearch ${vsearch_folder}/exported-artifact/${new_sample_name}_0_L001_R1_001.fastq \
			--output_file ${joined_folder}/${new_sample_name}_0_L001_R1_001.fastq
		
		echo "Completed ${new_sample_name}. Saved to ${joined_folder}"
