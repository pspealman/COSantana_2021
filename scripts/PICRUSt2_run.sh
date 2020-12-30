	

# source activate picrust2

# unzip rep-seqs-dada2.qza to get dna-sequences.fasta
# unzip table-dada2.qza to get feature-table.biom

picrust_dir=/scratch/ps163/Dr_Carolina/picrust2/

picrust2_pipeline.py \
	--stratified \
	--per_sequence_contrib \
	-p 4 \
	-s ${picrust_dir}dna-sequences.fasta \
	-i ${picrust_dir}feature-table.biom \
	-o picrust2_out_pipeline

# --max_nsti

picrust2_pipeline.py \
	--stratified \
	--per_sequence_contrib \
	-p 4 \
	--max_nsti 0.15 \
	-s ${picrust_dir}dna-sequences.fasta \
	-i ${picrust_dir}feature-table.biom \
	-o picrust2_015-nsti_out_pipeline
