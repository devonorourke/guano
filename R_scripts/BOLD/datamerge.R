#for combining dataframes containing complete BOLD datasets (pulled using the 'bold' R package)
#note that you need to have a directory containing only these files you want to combine - no other subdirectories or other files can be in this directory you set below

#set the directory to the parent folder to all files you want to combine:
setwd("/Users/devonorourke/Desktop/bold_data/copies_for_backup/")

#read in all the names of the files:
file_list <- list.files()

#combine the files
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.delim(file)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.delim(file)
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}


#you'll want to write this file into a text file to clear the memory from the system eventually...
###write.table(dataset, file = "/home/fosterLab/devonr/databases/BOLD/all_arthropod_fromRbind.txt")
###consider adding “row.names = false” to the output
