#***********************************************************
#*************************RNA SEQ **************************
#***********************************************************

#Load design matrix
design_raw<-read.csv(paste(exp_directory,"/Metadata.csv",sep=""),header=TRUE,fileEncoding = 'UTF-8-BOM')
.GlobalEnv$design<-data.frame(Intercept=1,
                              GROUP=design_raw$GROUP)

set.seed(42)
if(length(levels(as.factor(design$GROUP)))>1){
  #---------------COUNTING FEATURES----------------------
  
  if(file.exists(paste(fastq_dir,"/rawfeaturecounts.csv",sep=""))==FALSE){
   
    #Get aligned .bam files in 1fastqfiles
    bam.files <- list.files(path = forward_dir, 
                            pattern = ".BAM$", 
                            full.names = TRUE)
    
    fc<-featureCounts(files=bam.files, 
                      annot.ext=ann_file,
                      isPairedEnd=paired_end,
                      isGTFAnnotationFile=TRUE,
                      GTF.featureType = ft_val,
                      GTF.attrType= att_val)
    
    countdata_raw<-as.data.frame(fc$counts)
    names(countdata_raw)<-str_replace_all(names(countdata_raw),".fastq.gz.subread.BAM","")
    
    setwd(fastq_dir)
    write.csv(countdata_raw,"rawfeaturecounts.csv",row.names=TRUE)
    
  } else{
    tk_messageBox(message="Feature counts already detected! Skipping counting.")
  }
} else{
  tk_messageBox(message="Cannot compare gene expression for one group.")
}
close(prog)
