#---------------------BUILD INDEX----------------------
set.seed(42)

#Get reference genome (index.file from parms)
setwd(annot_dir)
buildindex(basename=ann_file,
           reference=gen_file,
           gappedIndex=TRUE,
           indexSplit=TRUE)

