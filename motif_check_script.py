from itertools import permutations
from Bio import SeqIO # import the SeqIO module from the BioPython library
import pandas as pd # import pandas for data manipulation
import os # import os for file and directory handling
import time # import time for measuring the run time of the script
import itertools # import itertools for combinations generation
import re # import re for string manipulation
import scipy.stats as stats
import sys
# create a list of DNA letters that are ambiguous and can be represented by "N" in the input motifs
listA=["N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V"]

# read in a list of the ambiguous DNA letters and their possible nucleotide replacements from an Excel file
lattersDF= pd.read_excel("latters.xlsx")
# start measuring the run time of the script
start_time = time.time()
# prompt the user for whether to search on the upper or both strands of DNA
RC= input("Do you wnt to search on the bottom (B), both bottom and upper (UB) or only on upper (U)?")
# read in the list of DNA motifs from an Excel file
motifsDF= pd.read_excel("motifs to check.xlsx")  
list_input=motifsDF["Motif seq"].tolist()
list_input = [s.upper() for s in list_input] # convert all motifs to uppercase
TF_name=motifsDF["Name of TF"].tolist() # create a list of the TF names from the "Name of TF" column of the dataframe
# initialize a counter for the current motif
L=0  
TF=0
# loop through each motif in the list
for t1 in list_input:  
    t1_lista=[*t1] # create a list of characters in the current motif
    intersection =[element for element in t1_lista if element in listA] # create a list of the ambiguous DNA letters in the current motif
    counterN=0 # initialize a counter for the number of ambiguous DNA letters in the current motif


    # check if the output directory for the current motif and strand exists. If it doesn't, create it.
    isExist = os.path.exists(f"output/{RC}{TF_name[TF]}_{t1}") 
    if isExist!=True:
        os.mkdir(f"output/{RC}{TF_name[TF]}_{t1}")
        
    # create an empty dataframe to store the positions of each occurrence of the current motif in each DNA sequence
    df_poition_of_motifs=pd.DataFrame()
    df_poition_of_motifs["a"]=range(1,100000)
    df_ratio=pd.DataFrame() # create an empty dataframe to store the ratio of the occurrences of the current motif to the length of each DNA sequence
    
    # loop through all fasta files in the "input" directory and its subdirectories
    for path, currentDirectory, files in os.walk("input"):
        for file in files:  #check all files with gene ids
            lista_startmotif=[] # create an empty list to store the starting positions of each occurrence of the current motif in the current DNA sequence
            basepairs_upstreams=int(file.split("bps")[0])# extract the number of basepairs upstream from the file name 
            # Load sequences from the current fasta file 
            sequences = SeqIO.parse(open(f"{os.path.join(path, file)}"),'fasta')
            file=file.replace(".fasta","").replace(".Fasta","").replace(".fa","") # remove the file extension from the file name
            new_neme=f"output/{RC}{TF_name[TF]}_{t1}/{t1}_{file}"[:200]+".txt" # create the output file name for the current DNA sequence and motif
            length_of_name=180 # limit the file name to 180 characters
            new_neme= (new_neme[:length_of_name]+".txt") if len(new_neme) > length_of_name else new_neme
            with open(new_neme, "w") as f: # open the output file for writing
                f.write(f"Motif\tseq\tis\ttimes_present\tX\tpresent\tin\tGeneID\n")  # write the header line to the output file
                total_genes=0 # initialize a counter for the total number of DNA sequences
                for seq in sequences: # loop through each DNA sequence in the current fasta file
                    total_genes+=1 # increment the total number of DNA sequences
                    if RC=="UB": # check if the search is being done on both strands or just the upper strand
                        s = seq.seq.upper() + " " +seq.seq.reverse_complement().upper()# # if both strands, concatenate the upper strand with the reverse complement of the lower strand
                    elif RC=="B":
                        s=seq.seq.reverse_complement().upper()
                    elif RC=="U": # if only the upper strand, just use the upper strand
                        s = seq.seq.upper()
                    else:
                        print("You entered wrong input strand try again")
                        sys.exit()    
                    times_present=0 # initialize a counter for the number of times the current motif is present in the current DNA sequence
                    #check if "N "is in the input motif t1 # check if the current motif contains any ambiguous DNA letters
                    list_L=[]

                    if intersection: # if it does, loop through each ambiguous DNA letter

                        #print (intersection)
                        for e in intersection: 
                            # create a list of possible nucleotide replacements for the current ambiguous DNA letter
                            globals()['latters%s' % e] = [latter for latter in lattersDF[e].to_list() if str(latter) != 'nan']
                            list_L.append(globals()['latters%s' % e]) # create a list of possible nucleotide replacements for the current ambiguous DNA letter

                            counterN = t1.count(e)+counterN
                            
                        #print (list_L)    
                        for element in itertools.product(*list_L):
                            #print (element) 
                            t = t1
                            el=0
                            for _ in element:
                                t=t.replace(intersection[el], _ ,1) # change 'N into tt with i latter from array N'
                                el+=1
                            #print(t)
                                times_present=s.count(t)+times_present #count how many times motif is present in seq ####can be UNINDENTED
                                lista_startmotif.extend([m.start() for m in re.finditer(str(t), str(s))])# append a list of the starting positions 
                    else:
                        times_present=s.count(t1) #count how many times motif is present in seq
                        lista_startmotif.extend([m.start() for m in re.finditer(str(t1), str(s))])# append a list of the starting positions # if you want to count overlaping f'(?={str(t1)})'
                    if times_present!=0:
                        f.write(f"Motif\t{t1}\tis\t{times_present}\tX\tpresent\tin\t{seq.id }\n") # write the results for the current DNA sequence and motif to the output file
            df=pd.read_csv(new_neme, sep="\t")
            number_of_BS=sum(df["times_present"])
            number_of_genes_with_motif=len(df["times_present"])      
            file_splited=file.split('_',1)
            try:
                new_file_name=f"output/{RC}{TF_name[TF]}_{t1}/{t1}_{file_splited[0]}_{number_of_genes_with_motif}({number_of_BS}bs)of{file_splited[1]}.txt"
                new_file_name = (new_file_name[:length_of_name]+".txt") if len(new_file_name) > length_of_name else new_file_name
                if (os.path.exists(new_file_name)):
                    os.remove(new_file_name)
                os.rename(new_neme,new_file_name)
            except IndexError: #if an index error apear
                print("file does not have _")
                new_file_name=f"output/{RC}{TF_name[TF]}_{t1}/{t1}_{number_of_genes_with_motif}({number_of_BS}bs)_{file}.txt"
                new_file_name= (new_file_name[:length_of_name]+".txt") if len(new_file_name) > length_of_name else new_file_name
                if (os.path.exists(new_file_name)):
                    os.remove(new_file_name)
                os.rename(new_neme,new_file_name)
            df_poition_of_motifs[f"{file}"]=pd.Series([x -basepairs_upstreams for x in lista_startmotif])#make a colomn with positions of found motif in promoter 
            dfs=df_poition_of_motifs[f"{file}"].value_counts()
            dfs.to_excel(f"output/{RC}{TF_name[TF]}_{t1}/{t1}_{number_of_genes_with_motif}({number_of_BS}bs)_{file}"[:200]+".xlsx")
            ######
            df_ratio[file]= pd.Series([number_of_genes_with_motif*100/total_genes,number_of_genes_with_motif,total_genes]) #calculate percentage
            ########
            print(f"{file} Finished!")
            L+=1 # increment the counter for the current motif
    df_poition_of_motifs=df_poition_of_motifs.drop('a', axis=1)
    df_poition_of_motifs.to_excel(f"output/{RC}{TF_name[TF]}_{t1}/{TF_name[TF]}_{t1}.xlsx", index=False) #save to excel
    #df_ratio=df_ratio.drop('a', axis=1)
    df_ratio=df_ratio.T
    df_ratio=df_ratio.rename(columns={0: 'Percentage', 1: f'Number of Promoters containing {TF_name[TF]} site', 2: "Total number of genes"})
    #############################
    #calculte statistics
    max = int(df_ratio["Total number of genes"].max()) #max_total_number_of_genes
    g_in_max=int(df_ratio.loc[df_ratio['Total number of genes'] == max, f'Number of Promoters containing {TF_name[TF]} site'].iloc[0]) #number of genes in max_total_number_of_genes with BS
    p_values=[]
    Odds_ratio=[]
    for index, row in df_ratio.iterrows():
        GFI=int(row[f'Number of Promoters containing {TF_name[TF]} site']) #gene of interest
        GOS=int(row['Total number of genes']) #total genes in sample usually DEGs
        # Create a contingency table with the counts of genes that have or do not have the binding site in their promoter regions
        # For example, if you have a set of 100 DEGs and 5000 genes in the genome, and the binding site is present in the promoter regions of 20 DEGs and 200 genes in the rest of the genome, the table would look like this:
        #contingency_table = pd.DataFrame({'Has binding site': [20, 200], 'Does not have binding site': [80, 4800]}, index=['DEGs', 'Non-DEGs'])
        # Perform Fisher's exact test on the contingency table
        contingency_table = pd.DataFrame({'Has binding site': [GFI, g_in_max], 'Does not have binding site': [GOS-GFI, max-g_in_max]}, index=['DEGs', 'Non-DEGs'])
        
        oddsratio, pvalue = stats.fisher_exact(contingency_table)
        p_values.append(pvalue)
        Odds_ratio.append(oddsratio)
    
    df_ratio["p_values"] = p_values
    df_ratio["oddsratio"]= Odds_ratio
    #############################
    df_ratio.to_excel(f"output/{RC}{TF_name[TF]}_{t1}_percentage.xlsx") #save to excel
    TF+=1 # increment the counter for the current TF
    print(f"done {t1}")
print(f"--- done {L} files in {round((time.time() - start_time),2)} s ---" ) # measure and print the run time of the script

