#***********************************************************
#*************************RNA SEQ **************************
#***********************************************************
set.seed(42)

#---------------------ALIGN READS----------------------

lapply(1:length(forward_files),function(x){
  
  #Update progress bar
  .GlobalEnv$prog_val<-prog_val+prog_inc
  setWinProgressBar(prog,label=paste("Aligning ",just_f_names[x],"...",sep=""),value=prog_val)
  
  #Align reads to index
  if(paired_end==FALSE){
    align(index = ann_file,
          readfile1 = forward_files[x],
          type=seq_val)
  } else if(paired_end==TRUE){
    align(index = ann_file,
          readfile1 = forward_files[x],
          readfile2 = reverse_files[x],
          type=seq_val)
  }
  
})
  
