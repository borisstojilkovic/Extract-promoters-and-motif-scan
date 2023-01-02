from Bio import SeqIO
import csv
import pandas as pd
import os
from Bio.Seq import Seq
import time
# Ask the user which species they want to work with
species=input("If you want to work with tomato type S, if you want to work with arabidopsis type A: ").upper()
# Set the path to the genome and BED file based on the species selected
if species== "S":
    GENOME = "Genome_and_annotation/Solanum_lycopersicum.SL3.0.dna_sm.toplevel.fa.fasta" # fasta file
    BED_File = "Genome_and_annotation/Solanum_lycopersicum3.0_BED.bed" #bed file name 

elif species== "A":
    GENOME = "Genome_and_annotation/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa" # fasta file
    BED_File = "Genome_and_annotation/Arabidospis_bed.bed" #bed file name 

# Read the BED file into a pandas DataFrame
bedDF = pd.read_csv(f"{BED_File}",sep="\t",low_memory=False) #open BED file as data frame  
# Ask the user for the length of the upstream and downstream regions  
length_of_upstream_region = int(input("How many bps you want from upstream region? "))
length_of_down_region = int(input(f"How many bps you want from downstream region? ")) 

#check if the folder output exsists if it doesnt make it
isExist = os.path.exists("output") 
if isExist!=True:
    os.mkdir("output")
    
start_time = time.time()
N=0
sequences_genome = SeqIO.parse(open(GENOME),'fasta')
# Create a dictionary to store the sequences from the genome
d={}
for sequence in sequences_genome:
    d[sequence.id]=sequence.seq
# Iterate over the files in the input folder
for file in os.listdir('input'): #check all files with gene ids
    ext=file.split(".")[-1] # Get the file extension
    # Read the file into a pandas DataFrame based on the file extension
    if ext=="tab":
        file_df = pd.read_csv(f"input/{file}",sep='\t',low_memory=False) #make df of input file with geneIDs
        file=file.replace(".tab","")
    elif ext=="xlsx": 
        file_df = pd.read_excel(f"input/{file}")
        file=file.replace(".xlsx","")
    geneids = file_df["GeneID"].tolist()# Convert the GeneID column to a list
    # Open the output file for writing
    with open(f"output/{length_of_upstream_region+length_of_down_region}bps({length_of_upstream_region}_{length_of_down_region})_{file}"[:200]+".Fasta", "w+") as f: #open result file to write 
        for ID in geneids:  # Iterate over the list of gene IDs
            try: # Extract the relevant information about the gene from the BED data frame
                chromosome_nb=bedDF.loc[bedDF["geneID"]==ID,"chromosome"].tolist()[0]
                start_pos=int(bedDF.loc[bedDF["geneID"]==ID,"start"].tolist()[0])
                end_pos=int(bedDF.loc[bedDF["geneID"]==ID,"end"].tolist()[0])
                strand_info=str(bedDF.loc[bedDF["geneID"]==ID,"strand"].tolist()[0])
                take_up_to=start_pos-length_of_upstream_region # Calculate the start position of the upstream region
                sequences_genome = SeqIO.parse(open(GENOME),'fasta') # Load the genome sequences from the FASTA file
                seq=d[str(chromosome_nb)] # Extract the genome sequence for the chromosome of the gene
                s=seq[take_up_to:(start_pos+length_of_down_region)] # Extract the sequence for the upstream region
                # If the upstream region extends beyond the start of the chromosome,
                # take only the part until the start of the chromosome
                if length_of_upstream_region>=start_pos:
                    take_up_to=0
                    #s=seq[take_up_to:start_pos]
                    s=seq[take_up_to:(start_pos+length_of_down_region)]
                if strand_info=="-":# If the gene is on the minus strand, reverse complement the sequence
                    #reverse complement
                    take_up_to=end_pos+length_of_upstream_region
                    s=seq[(end_pos-length_of_down_region):take_up_to]
                    # If the upstream region extends beyond the start of the chromosome,
                    # take only the part until the start of the chromosome
                    if length_of_upstream_region>=(len(seq))-end_pos:
                        take_up_to=int(len(seq))
                        s=seq[end_pos:take_up_to]
                    s=s.reverse_complement()#reverse complement
                #print(s)
                f.write(f">{ID}\n{s}\n")
            except IndexError:
                print(f"{ID} not found")
    print(f"{file} FINISHED!")
    N+=1
print(f"--- done {N} files in {round((time.time() - start_time),2)} s ---" )
    ###################

