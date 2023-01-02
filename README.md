# Extract-promoters-and-motif-scan


Scripts extract promoter regions and search for given motifs inside those promoters.

1. System requirements
	Windows
	All software dependencies and operating systems (including version numbers)
	Python 3.9.6
	Versions the software has been tested on
	Python 3.9.6
	Any required non-standard hardware
	None

2. Installation guide
	install biopython by running the following command in the command line: pip install biopython
	install pandas by running the following command in the command line: pip install pandas


3. Instructions for use
	Install python, biopython, and pandas
	
	Download :
		All files from this GitHub
		Tomato genome from "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL3.0.dna_sm.toplevel.fa.gz"
			after downloading unzip the file
		Arabidopsis genome from: "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz"
			after downloading unzip the file
			
	Make directories: 
		Extract promoter region\Genome_and_annotation
		Extract promoter region\input
		motif_check\input

	Place files in the correct directories:
	
		"Extract promoter region\extract_promoters_script.py"
		"Extract promoter region\Genome_and_annotation\Arabidospis_bed.bed"
		"Extract promoter region\Genome_and_annotation\Solanum_lycopersicum3.0_BED.bed"
		"Extract promoter region\Genome_and_annotation\Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa"
		"Extract promoter region\Genome_and_annotation\Solanum_lycopersicum.SL3.0.dna_sm.toplevel.fa.fasta"
		
		"motif_check\latters.xlsx"
		"motif_check\motif_check_script.py"
		"motif_check\motifs to check.xlsx"
		

	Place input files for screening (they can be CSV (tab separated) or excel (".xlsx") and they must have “GeneID” column with geneids from Arabidopsis or tomato (3.0 annotation)) in the directory "Extract promoter region\input". Only one organisam at the time. 

	Run script "extract_promoters_script.py" in "Extract promoter region" directory. Provide responses for inputs. 
	
	Results:
		All generated files will be in the directory "Extract promoter region\output" 
		# at this point, all desirable promoters are extracted 
	
	
	Copy files for screening from "Extract promoter region\output" to "motif_check\input"
	Provide all motif sequences 5'->3' in file "motif_check\motifs to check.xlsx" colomn "Motif seq" and name of the transcription factor in column "Name of TF"
	Run script "motif_check_script.py" in "motif_check" directory. Provide responses for inputs.  
	
	Results 
		In the output directory ("motif_check\output") one will find excel files (named as provided: Name of TF_Motif seq_percentage) which will contain the percentage of the genes with transcription factor binding site in their promoters
		In the output directory ("motif_check\output"), one will find directories (named as provided: Name of TF_Motif seq) which will contain files (.txt, tab separated) with extracted GeneID that contain at least one binding site for a given transcription factor
