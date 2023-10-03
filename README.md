# TEs_and_Transcript_emergence

Scripts for "TEs favour de novo transcripts emergence"

Part 1: De novo transcripts (Detection, filtering, characteristics)

	- 1. Extract Blast no match
		Input: Blast tabular output file, Transcriptome GTF file
		Output: List with all transcript that have no hit
	- 2. Filter annotated
		Input: List with all transcript that have no Blast hit, Transcriptome GTF file
		Output: List with filtered annotated transcripts (splicing + TPM)
	- 3. Filter de novo
		Input: List with all transcript that have no Blast hit, Transcriptome GTF file
		Output: List with filtered de novo transcripts (splicing + TPM)
	- 4. Calculate transcript properties
		Input: Transcriptome assembly fasta + gtf, Bedtools intersect gtf (with gene overlap info)
		Output: Csv file with transcript properties (exon/intron number, length etc.)

Part 2: Transcript homology (Find shared transcripts and non expressed homologues)

	- 1. Get shared transcripts info
		Input: Blast tabular output file (after blasting transcripts against each other), 
		       List of de novo transcript ids
		Output: File with populations to blast against, File with info on shared transcripts, File for 2.2
	- 2. Build orthogroups
		Input: File for 2.2 (from 2.1)
		Output: File with orthogroup info
	- 3. Blast transcripts against whole genomes
		Input: De novo transcript fasta unspliced, Genomes, file from 2.1 (populations to blast against)
		Output: Blast tabular output
	- 4. Get number of hit and no hit
		Input: Merged file from 2.1 (populations to blast against), blast output (2.3)
		Output: Info csv of which transcript had a sufficient blast hit
	- 5. Make homolog bedfile
		Input: Merged file from 2.1 (populations to blast against), blast output (2.3)
		Output: Bedfile for each population containing the homologs from that population

Part 3: Create the sequence files needed for the motif and TE search (up/downstream and intergenic)

	- 1. Turn GTF to bed
		Input: Gene/Transcript gtf file
		Output Gene/Transcript bedfile
	- 2. Create up/downstream and intergenic bedfiles
		Input: Gene/Transcript gtf file OR Bedtools sample output (intergenic)
		Output: Bedfile with up/downstream regions, Bedfile with intergenic regions

Part 4: TEs + Methylation

	- 1.Turn TE gtf to bedfile
		Input: TransposonUltimate gtf file, sequence_heads.txt
		Output: TE annotation as bedfile
	- 2. Bedtools intersect output csv
		Input: Bedtools intersect output (Transcript bed against TE bed)
		Output: Csv file
	- 3. Prepare TE overlap files for bedtools merge
		Input: Output files from 4.2 for each category
		Output: Bedfiles for bedtools merge (to get non-overlapping TE intervals)
	- 4. Calculate CpGoe
		Input: (Transcript) fasta file
		Output: Tsv file with info on CpGoe
	- 5. Compare TE Transcript Homology
		Input: Bedtools intersect output files of Transcripts - TE and Homologs -TE
		Output: Info file of which TEs are shared between transcript + homolog

Part 5: Motifs

	- 1. Search for motifs
		Input: Fasta file (e.g. with transcript upstream regions), Motif db (position frequency matrix as txt file)
		Output: Csv file with all identified motifs passing the treshold
	- 2. Get number of motifs
		Input: Output csv files from 5.1
		Output: Csv that sums up the number of motifs per transcript for two thresholds (complete upstream region)
	- 3. Get number of core
		Input: Output csv files from 5.1
		Output: Csv that sums up the number of core per transcript for two thresholds (region: -200/+100)
	- 4. Get number of individual motifs
		Input: Output csv files from 5.1
		Output: Csv with the count of each motif per transcript for two thresholds (complete upstream region)

Part 6: R scripts for plots, tables and stats

	- 1. De novo transcripts
		Input: Excel file with info on Blast + filtering steps, Merged csv from 1.4
		Output: Fig. 1A, supplemental figure with chromosome distribution
	- 2. Homology
		Input: TE info file (from 6.4), outputfile from 2.2 and outputfile from 2.4
		Output: Fig. 1E, supplemental table with blast hit info
	- 3. Get TE number
		Input: Bedtools intersect outputfiles merged
		Output: Plot Fig.3B + number TE
	- 4. Get relative TE overlap
		Input: Output files from bedtools merge followed by bedtools intersect to re-intersect with the transcripts/up-and downstream
		Output: Csv file with TE relative TE overlap
	- 5. Circular Plot CpG
		Input: Genome fai file, Transcript gtf file, TE annotation gtf, txt file with chromosome lengths, CpG info file (from 4.4)
		Output: Fig.2A
	- 5. Circular Plot TE overlap	
		Input: Genome fai file, Transcript gtf file, TE annotation gtf, txt file with chromosome lengths, Bedtools merge output file
		Output: Fig. 2B
	- 6. Compare TEs Transcript Homolog
		Input: TE info file (from 6.4), Outputfile from 4.5
		Output: Info on whether transcript + homolog share TEs
	- 7. Individual motifs plots and stats:
		Input: Outputfile from 5.4, 3 textfiles (Ids of TE transcripts/TE overlap transcripts and no TE overlap transcripts)
		Output: Csv with info on individual motifs that occur more in transcripts than noncoding 
	- 8. Merge motif and TE data
		Input: Outputfile from 5.2 (motifs) and 5.3 (core) , 3 textfiles (Ids of TE transcripts/TE overlap transcripts and no TE overlap transcripts)
		Output: Csv file for core and motifs (where TE info are added)
	- 9. Merge transcript and homolog data for stats
		Input: Csv files from 4.2 (for transcripts + homologs and their up/downstream), Bedtools merge output (for transcripts + homologs), Motif output csv fron 5.2 (motifs) and 5.3 (core) for transcripts and homologs
		Output: Merged csv file with info for transcripts and homologs
 	 - 10. Filter transcript homolog data
   		 Input: File from 6.9
   		 Output: Filtered file for 6.11
	- 11. Plots and stats
		Input: Merged csv file with info for transcripts and homologs (from 6.10), Csv file from 6.4, Csv file from 6.8, Merged outputfile from 1.4, CpG info files, transcript/te density files
		Output: Fig.3 A, C, D, E, Fig.4, Fig.5 + stats, supplemental figures
