#Filter transcript homolog data
#Check that the transcripts within an orthroup have the same homolog
#If yes keep only one (from each population)

library(readr)
library(dplyr)
library(tidyverse)

#Load data

#Dataset with merged transcript/homolog data
df <- read_csv("Merged_data_de_novo_transcript_homolog_correct.csv") 

#Dataset with positions of transcripts/homologs
df2 <- read_csv("Positions_trans_hom.csv")

df2$Transcript_Id = df2$Id

#Merge the two datasets
final = merge(df, df2, by = c("Transcript_Id", "Population"))

#Exclude orthogroups with duplicated transcripts
#Split ID Alternative spliceforms have the same gene id and only differ in .1 .2 etc
#Transcirpt_ID == "Gene" ID of the transcript
#Rest = Population of transcript
all = subset(final,Transcription_status == "Transcript")
All_transcripts = all$Transcript_Id

data = final

data = data %>%
  separate(Transcript_Id, c("ID1", "Transcript_Id", "ID2", "Rest")) 

data$Pop_ID = paste(data$Transcript_Id, data$Population, data$Number_Orthogroup)

#Keep only one (transcript and max 6 homologs) for each alternative spliceform
data = data[!duplicated(data$Pop_ID), ] 
transcript_num = subset(data, Transcription_status == "Transcript")
all1 = subset(transcript_num,Transcription_status == "Transcript")
all1$Transcript_Id = paste(all1$ID1, all1$Transcript_Id, all1$ID2, sep = ".") #Recreate Transcript_ID column
all1$Transcript_Id = paste(all1$Transcript_Id, all1$Rest, sep = "::")
Remove_spliceforms = all1$Transcript_Id

#Remove the duplicated transcripts from the same population that are not alternative spliceforms
#Basically check if in an orthogroup there are two transcripts with the same population
#Rest is the population of the transcript meaning for each rest per population there should be max. 1 transcript/homolog
#If not OG is removed

filtered_df <- data %>%
  group_by(Number_Orthogroup, Population, Rest) %>%
  filter(n() > 1) %>%   
  ungroup()

filtered_df = filtered_df[order(filtered_df$Number_Orthogroup),] #just to check visually

#Remove orthogroups with duplicated transcripts that are NOT alternative spliceforms
data <- data %>%
  filter(!Number_Orthogroup %in% filtered_df$Number_Orthogroup)

#Merge data back together
data$Transcript_Id = paste(data$ID1, data$Transcript_Id, data$ID2, sep = ".") #Recreate Transcript_ID column
data$Transcript_Id = paste(data$Transcript_Id, data$Rest, sep = "::")
data$ID1 = NULL
data$ID2 = NULL
data$Rest = NULL

transcript_num = subset(data, Transcription_status == "Transcript")
all2 = subset(transcript_num,Transcription_status == "Transcript")
Remove_duplicates = all2$Transcript_Id

######################################################################################
#Split the data

data1 = subset(data, Number_seqs_in_orthogroup == 1) #Specific transcripts 

data2 = subset(data, Number_seqs_in_orthogroup > 1) #Shared/duplicated transcripts -> need more filtering

transcript_num = subset(data2, Transcription_status == "Transcript")
All_transcripts_multseq = subset(transcript_num,Transcription_status == "Transcript")
All_transcripts_multseq = All_transcripts_multseq$Transcript_Id

######################################################################################
#Check that number and populations of non expressed homologs match
#Meaning each shared transcript has a homolog in the same population
#If that doesnt match -> remove

counts <- data2 %>%
  group_by(Number_Orthogroup, Transcript_Id) %>%
  summarise(count = n())

okay_combinations <- counts %>%
  summarise(same_count = all(count == first(count)))

data2 <- data2 %>%
  inner_join(okay_combinations, by = c("Number_Orthogroup")) 

data2 = subset(data2, same_count == TRUE)

transcript_num = subset(data2, Transcription_status == "Transcript")

Remove_chrom = subset(transcript_num,Transcription_status == "Transcript")
Remove_diff_chrom = Remove_chrom$Transcript_Id

######################################################################################
#For all transcripts with multiple seqs in the orthogroup: Check that the 
#homologs start and end match in a 200 nt window

data2 = data2[order(data2$Seq_id),]

#Check that the start and end values match 
filtered_df <- data2 %>%
  group_by(Seq_id, Population) %>%
  summarise(
    Max_Start = max(transcript_start),
    Min_Start = min(transcript_start),
    Max_End = max(transcript_end),
    Min_End = min(transcript_end),
    Number_Orthogroup = Number_Orthogroup,
    Transcription_status = Transcription_status
  ) %>%
  ungroup() 

filtered_df$Start_filter = filtered_df$Max_Start - filtered_df$Min_Start
filtered_df$End_filter = filtered_df$Max_End - filtered_df$Min_End

filtered_df <- filtered_df %>%
  group_by(Seq_id, Population) %>%
  mutate(Max_Start_val = max(Start_filter)) %>%
  ungroup()

filtered_df <- filtered_df %>%
  group_by(Seq_id, Population) %>%
  mutate(Max_End_val = max(End_filter)) %>%
  ungroup()

filtered_df <- subset(filtered_df, Max_Start_val > 200 & Max_End_val > 200)

filtered_df = subset(filtered_df, Transcription_status != "Transcript" & Transcription_status != "Transcribed_TE")
#Exclude transcripts and transcribed TE (only want to check this based on Homologs)

#Merge the dfs back and include only the ones where this works. 
data2 <- data2 %>%
  filter(!Number_Orthogroup %in% filtered_df$Number_Orthogroup)

transcript_num = subset(data2, Transcription_status == "Transcript")
Removestartend = subset(transcript_num,Transcription_status == "Transcript")
Remove_startend = Removestartend$Transcript_Id

######################################################################################
#Now check that the homologs are all on the same chromosome
chrom_counts <- data2 %>%
  group_by(Seq_id, Population) %>%
  summarise(Unique_Chrom_Count = n_distinct(Chromosome),
            Number_Orthogroup = Number_Orthogroup)

sub_chrom = subset(chrom_counts, Unique_Chrom_Count > 1)

#Was 0 (as expected when positions match)

######################################################################################
#Remove orthogroups that contain transcribed TE and transcripts
type_counts <- data2 %>%
  group_by(Number_Orthogroup) %>%
  summarise(Unique_Type_Count = n_distinct(Transcription_status))

# Filter the original dataframe to retain only rows where all Chrom values are the same within each ID and Pop combination
data2 <- data2 %>%
  left_join(type_counts, by = c("Number_Orthogroup")) %>%
  filter(Unique_Type_Count <= 2) %>%
  select(-Unique_Type_Count)

data2 = data2[order(data2$Number_Orthogroup),]

transcript_num = subset(data2, Transcription_status == "Transcript")
RemoveTE = subset(transcript_num,Transcription_status == "Transcript")
Remove_TEs = RemoveTE$Transcript_Id

######################################################################################
data2$same_count = NULL
data2$same_seqs = NULL

#merge back
final_merge = rbind(data1, data2)
final_merge2 = subset(final_merge, Transcription_status == "Transcript" |Transcription_status == "Non_transcribed_Homolog")

#Only keep one homolog per orthogroup
filtered_df <- final_merge2 %>%
  group_by(Seq_id, Population) %>%
  filter(transcript_start == min(transcript_start)) %>%
  ungroup()

filtered_df = filtered_df[order(filtered_df$Number_Orthogroup),]

#If two homologs have the same start keep only one of them
filtered_df$Pop_ID = paste(filtered_df$Seq_id, filtered_df$Population)

filtered_df <- filtered_df[!duplicated(filtered_df$Pop_ID), ]
filtered_df$Pop_ID = NULL

#Check how many transcripts remain
x = subset(filtered_df, Transcription_status == "Transcript")


filtered_df2 = filtered_df[, c(1,2,4,5,6,7,8,9, 10, 11,12,13,14,15,16,17,18,19,20,21, 22,23,24,25)]

write.csv(filtered_df2, file = "Filtered_de_novo_transcripts_homologs.csv")

#hist(x$Number_populations_in_orthogroup)


#Now get the ids of the different transcripts removed in each step
Filter1_Removed_splieformes = setdiff(All_transcripts, Remove_spliceforms)
Filter2_Removed_duplicates = setdiff(Remove_spliceforms, Remove_duplicates)
Filter3_Removed_diff_chrom = setdiff(All_transcripts_multseq, Remove_diff_chrom)
Filter4_Remove_startend = setdiff(Remove_diff_chrom, Remove_startend)
Filter5_Remove_OGTE = setdiff(Remove_startend, Remove_TEs)

#Save them in a df
writeLines(Filter1_Removed_splieformes, "Step1_RemoveSpliceforms.txt")
writeLines(Filter2_Removed_duplicates, "Step2_RemoveDuplicates.txt")
writeLines(Filter3_Removed_diff_chrom, "Step3_RemoveDifferentHomologs.txt")
writeLines(Filter4_Remove_startend, "Step4_RemoveDifferentStartEnd.txt")
writeLines(Filter5_Remove_OGTE, "Step5_Remove_OG_TE.txt")
