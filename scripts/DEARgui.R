
#**********************************************

#--------------------INSTALL NECESSARY PACKAGES

#**********************************************
tryCatch(find.package("utils"),error=function(e){
  install.packages(utils,repos="http://lib.stat.cmu.edu/R/CRAN/") 
})
library(utils)

start_val<-0
open_prog<-winProgressBar(title="Differential Expression Analysis in R",
                          label="Installing packages...",min=0,max=100,initial=start_val)

package.list<-c("this.path",
                "tcltk",
                "stringr",
                "BiocManager",
                "RColorBrewer",
                "gplots",
                "ggplot2",
                "dplyr")
n_pack<-length(package.list)
prog_inc<-25/n_pack
for(i in 1:length(package.list)){
  tryCatch(find.package(package.list[i]),
           error = function(e){
             .GlobalEnv$start_val<-start_val+prog_inc
             setWinProgressBar(open_prog,value=start_val,label=paste("Installing ",package.list[i],"...",sep=""))
             install.packages(package.list[i],repos="http://lib.stat.cmu.edu/R/CRAN/") 
           })
}

setWinProgressBar(open_prog,value=25,label="Loading this.path and tcltk...")
library(this.path)
library(tcltk)
setWinProgressBar(open_prog,value=35,label="Loading dplyr and stringr...")
library(dplyr)
library(stringr)
setWinProgressBar(open_prog,value=45,label="Loading BiocManager...")
library(BiocManager)
setWinProgressBar(open_prog,value=65,label="Loading ggplot2 and RColorBrewer...")
library(ggplot2)
library(RColorBrewer)
library(gplots)

#Install BiocManager packages
start_val<-70
setWinProgressBar(open_prog,value=70,label="Installing BiocManager packages...")
bioc.packages<-c("Rsubread", 
                 "edgeR",
                 "limma")
prog_inc<-(85-70)/length(bioc.packages)
for(i in 1:length(bioc.packages)){
  tryCatch(find.package(bioc.packages[i]),
           error=function(e) {
             .GlobalEnv$start_val<-start_val+prog_inc
             setWinProgressBar(open_prog,value=start_val,label=paste("Installing ",bioc.packages[i],"...",sep=""))
             BiocManager::install(bioc.packages[i])
           })
}

setWinProgressBar(open_prog,value=85,label="Loading Rsubread...")
library(Rsubread)
setWinProgressBar(open_prog,value=90,label="Loading edgeR...")
library(edgeR)
setWinProgressBar(open_prog,value=95,label="Loading limma...")
library(limma)
close(open_prog)

#Create fonts
title_font<-tkfont.create(size=14,weight="bold",family="Helvetica",underline=TRUE)
header_font<-tkfont.create(size=11,family="Helvetica",weight="bold")
normal_font<-tkfont.create(size=10,family="Helvetica")
underline_font<-tkfont.create(size=10,family="Helvetica",underline=TRUE)
file_display_font<-tkfont.create(size=6,family="Helvetica")
small_font<-tkfont.create(size=8,family="Helvetica")

#Variables
exp_var<-tclVar("test")
.GlobalEnv$ft_var<-tclVar("")
.GlobalEnv$att_var<-tclVar("")
.GlobalEnv$thresh_var<-tclVar("")
.GlobalEnv$sample_var<-tclVar("")
.GlobalEnv$seq_var<-tclVar("")
.GlobalEnv$reverse_reads<-tclVar(1)
.GlobalEnv$forward_files<-NULL
.GlobalEnv$reverse_files<-NULL
.GlobalEnv$gen_file<-NULL
.GlobalEnv$ann_file<-NULL
.GlobalEnv$paired_end<-NULL
.GlobalEnv$recount_ft_var<-tclVar(0)
.GlobalEnv$group_sel_choices<-c("Control","Experimental")

#GET DIRECTORY
directory<-gsub("/scripts","",this.dir())

#Create general subfolders
if(dir.exists(paste(directory,"/Files",sep=""))==FALSE){
  dir.create(paste(directory,"/Files",sep=""))
  dir.create(paste(directory,"/Files/genomes",sep=""))
  dir.create(paste(directory,"/Files/annotations",sep=""))
}
if(dir.exists(paste(directory,"/Experiments",sep=""))==FALSE){
  dir.create(paste(directory,"/Experiments",sep=""))
}

#Activate function
activate<-function(){
  tkconfigure(run_but,state="normal")
  tkconfigure(sel_f_fastq,state="normal")
  tkconfigure(sel_r_fastq,state="normal")
  tkconfigure(pairwise,state="normal")
  tkconfigure(load_gen,state="normal")
  tkconfigure(load_ann,state="normal")
  tkconfigure(meta_data,state="normal")
  tkconfigure(inst_files,state="normal")
  tkconfigure(view_plots,state="normal")
  tkconfigure(recount_ft,state="normal")
  tkconfigure(build_plot,state="normal")
}

#Read gtf function
read_gtf_function<-function(){
  
  tmp_prog<-winProgressBar(title="Reading Annotation",label="Reading annotation...",
                           min=0,max=100,initial=50,width=300)
  
  #Read GTF file
  gtf_names<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
  tryCatch(.GlobalEnv$gtf_file<-read.table(ann_file,header=FALSE,sep='\t',col.names = gtf_names,nrows=10000),
           error=function(e){
             tk_messageBox(message="There was an error reading your annotation file!")
             close(tmp_prog)
             tkconfigure(ann_sel,text="\n")
             stop()
           })
  
  #Get features
  gtf_features<-levels(as.factor(gtf_file$feature))
  tkconfigure(ft_type,values=gtf_features)
  
  #Get attributes
  gtf_attributes<-trimws(unlist(str_split(gtf_file$attribute[1],";")),which="left")
  gtf_attributes_cur<-str_split(gtf_attributes," ")
  just_attributes<-lapply(1:length(gtf_attributes_cur),function(x){
    tmp<-unlist(gtf_attributes_cur[x])
    out<-tmp[1]
  })
  just_attributes_cur<-unlist(just_attributes)
  if(length(which(just_attributes_cur==""))>0){
    just_attributes_cur<-just_attributes_cur[-which(just_attributes_cur=="")]
  }
  .GlobalEnv$just_attributes_cur<-just_attributes_cur
  tkconfigure(att_type,values=just_attributes_cur)
  
  close(tmp_prog)
}

#Open experiment
open_exp_func<-function(){
  
  name_exp_func<-function(){
    set_name_fun<-function(){
      experiment_name<-tclvalue(exp_var)
      if(str_detect(experiment_name,"[:alpha:]")==TRUE){
        
        .GlobalEnv$experiment_name_cur<-paste(str_sub(Sys.time(),end=10),"-",experiment_name,sep="")
        if(dir.exists(paste(directory,"/Experiments/",experiment_name_cur,sep=""))==TRUE){
          tk_messageBox(message="Warning: that experiment already exists!")
        } else{
          dir.create(paste(directory,"/Experiments/",experiment_name_cur,sep=""))
          dir.create(paste(directory,"/Experiments/",experiment_name_cur,"/fastq",sep=""))
          tkdestroy(exp_name)
        }
      } else{
        tk_messageBox(message="Please enter a valid experiment name.")
      }
    }
    
    #Define experiment name
    exp_name<-tktoplevel()
    tkwm.geometry(exp_name,"150x80+500+200")
    tkwm.title(exp_name,"DEAR-main V1.0.0")
    top_frame<-tkframe(exp_name,padx=10)
    tkgrid(top_frame)
    name_lbl<-tklabel(top_frame,text="Name your experiment:")
    tkgrid(name_lbl)
    ent_name<-tkentry(top_frame,textvariable=exp_var,width=20,justify="center")
    tkgrid(ent_name)
    exp_but<-tkbutton(top_frame,text="Enter",command=set_name_fun)
    tkgrid(exp_but,pady=5)
    tkwait.window(exp_name) 
  }
  
  #Create save dir function
  create_save_dirs<-function(){
    #Create save directories
    .GlobalEnv$fastq_dir<-paste(exp_directory,"/fastq",sep="")
    .GlobalEnv$forward_dir<-paste(fastq_dir,"/forward_reads",sep="")
    .GlobalEnv$reverse_dir<-paste(fastq_dir,"/reverse_reads",sep="")
    .GlobalEnv$gen_dir<-paste(directory,"/Files/genomes",sep="")
    .GlobalEnv$annot_dir<-paste(directory,"/Files/annotations",sep="")
    .GlobalEnv$plot_dir<-paste(exp_directory,"/plots",sep="")
  }
  
  #Load parameters
  load_parms<-function(){
    
      .GlobalEnv$parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
      if(is.na(parms$Reference_Genome)==FALSE){
        .GlobalEnv$gen_file<-parms$Reference_Genome 
        .GlobalEnv$gen_name<-paste(str_sub(gsub(paste(dirname(gen_file),"/",sep=""),"",gen_file),start=0,end=18),"...",sep="")
        tkconfigure(gen_sel,text=paste(gen_name,"\n",sep=""),font=small_font)
      } else{
        .GlobalEnv$gen_file<-NULL
        tkconfigure(gen_sel,text="\n",font=small_font)
      } 
      if(is.na(parms$Pairwise)==FALSE){
        .GlobalEnv$paired_end<-parms$Pairwise
        if(parms$Pairwise=="1"){
          .GlobalEnv$reverse_reads<-tclVar(1)
          tkconfigure(pairwise,variable=reverse_reads)
          tkconfigure(sel_r_fastq,state="normal")
        } else{
          .GlobalEnv$reverse_reads<-tclVar(0)
          tkconfigure(pairwise,variable=reverse_reads)
          tkconfigure(sel_r_fastq,state="disabled")
        }
      }
      if(is.na(parms$Forward_Files)==FALSE){
        .GlobalEnv$forward_files<-str_split(parms$Forward_Files,"\\+")[[1]]
        
        just_names<-lapply(1:length(forward_files),function(x){
          tmp<-forward_files[x]
          new_name<-gsub(paste(dirname(tmp),"/",sep=""),"",tmp)
        })
        .GlobalEnv$just_f_names<-unlist(just_names)
        
        if(length(just_f_names)<4){
          print_names<-just_f_names[1:3]
          print_names[which(is.na(print_names))]<-" "
          tkconfigure(for_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
        } else{
          print_names<-just_f_names[1:2]
          print_names<-c(print_names,"Additional files not shown...")
          tkconfigure(for_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
        }
        
      } else{
        .GlobalEnv$forward_files<-NULL
        tkconfigure(for_reads,text="\n\n",font=file_display_font)
      }
      if(is.na(parms$Reverse_Files)==FALSE){
        .GlobalEnv$reverse_files<-str_split(parms$Reverse_Files,"\\+")[[1]]
        
        just_names<-lapply(1:length(reverse_files),function(x){
          tmp<-reverse_files[x]
          new_name<-gsub(paste(dirname(tmp),"/",sep=""),"",tmp)
        })
        .GlobalEnv$just_r_names<-unlist(just_names)
        
        if(length(just_r_names)<4){
          print_names<-just_r_names[1:3]
          print_names[which(is.na(print_names))]<-" "
          tkconfigure(rev_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
        } else{
          print_names<-just_r_names[1:2]
          print_names<-c(print_names,"Additional files not shown...")
          tkconfigure(rev_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
        }
      } else{
        .GlobalEnv$reverse_files<-NULL
        tkconfigure(rev_reads,text="\n\n",font=file_display_font)
      }
      
      if(is.na(parms$SeqType)==FALSE){
        .GlobalEnv$seq_var<-tclVar(parms$SeqType)
      } else{
        .GlobalEnv$seq_var<-tclVar("rna")
      }
      tkconfigure(sq_type,textvariable=seq_var)
      
      if(is.na(parms$FtType)==FALSE){
        .GlobalEnv$ft_var<-tclVar(parms$FtType)
      } else{
        .GlobalEnv$ft_var<-tclVar("gene")
      }
      tkconfigure(ft_type,textvariable=ft_var)
      
      if(is.na(parms$AttType)==FALSE){
        .GlobalEnv$att_var<-tclVar(parms$AttType)
      } else{
        .GlobalEnv$att_var<-tclVar("gene_id")
      }
      tkconfigure(att_type,textvariable=att_var)
      
      if(is.na(parms$Threshold)==FALSE){
        .GlobalEnv$thresh_var<-tclVar(as.numeric(parms$Threshold))
      } else{
        .GlobalEnv$thresh_var<-tclVar(0.6)
      }
      tkconfigure(thresh_entry,textvariable=thresh_var)
      
      if(is.na(parms$Sample)==FALSE){
        .GlobalEnv$sample_var<-tclVar(as.numeric(parms$Sample))
      } else{
        .GlobalEnv$sample_var<-tclVar(1)
      }
      tkconfigure(sample_entry,textvariable=sample_var)
      if(is.na(parms$Genome_Annotation)==FALSE){
        .GlobalEnv$ann_file<-parms$Genome_Annotation 
        ann_name<-paste(str_sub(gsub(paste(dirname(ann_file),"/",sep=""),"",ann_file),start=0,end=18),"...",sep="")
        tkconfigure(ann_sel,text=paste(ann_name,"\n",sep=""),font=small_font)
        read_gtf_function()
      } else{
        .GlobalEnv$ann_file<-NULL
        tkconfigure(ann_sel,text="\n",font=small_font)
      }
    
    if(file.exists(paste(exp_directory,"/Metadata.csv",sep=""))==TRUE){
      .GlobalEnv$meta_file<-read.csv(paste(exp_directory,"/Metadata.csv",sep=""))
      tkconfigure(meta_disp,text="Read annotation detected!")
    } else{
      tkconfigure(meta_disp,text="No read annotation.")
    }
    if(is.na(parms$Groups)==TRUE||length(parms$Groups)==0){
      .GlobalEnv$group_sel_choices<-c("Control","Experimental")
    } else{
      .GlobalEnv$group_sel_choices<-str_split(parms$Groups,"\\+")[[1]]
    }
    #End load parms function
  }

  #Create new exp function
  create_new_exp_fun<-function(){
    name_exp_func()
    .GlobalEnv$exp_directory<-paste(directory,"/Experiments","/",experiment_name_cur,sep="")
    print(exp_directory)
    tkconfigure(cur_exp,text=paste("Experiment: ",str_sub(experiment_name_cur,1,18),"...",sep=""))
    
    #Create experiment subfolders
    if(dir.exists(paste(exp_directory,"/fastq",sep=""))==FALSE){
      dir.create(paste(exp_directory,"/fastq",sep=""))
    }
    if(dir.exists(paste(exp_directory,"/fastq/forward_reads",sep=""))==FALSE){
      dir.create(paste(exp_directory,"/fastq/forward_reads",sep=""))
    }
    if(dir.exists(paste(exp_directory,"/fastq/reverse_reads",sep=""))==FALSE){
      dir.create(paste(exp_directory,"/fastq/reverse_reads",sep=""))
    }
    if(dir.exists(paste(exp_directory,"/plots",sep=""))==FALSE){
      dir.create(paste(exp_directory,"/plots",sep=""))
    }
    
    #Create parameters file
    parms<-data.frame(Reference_Genome=NA,
                      Genome_Annotation=NA,
                      Forward_Files=NA,
                      Reverse_Files=NA,
                      Pairwise=1,
                      SeqType="rna",
                      FtType="gene",
                      AttType="gene_id",
                      Threshold=0.6,
                      Sample=1,
                      Groups=NA)
    write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
    tkconfigure(gen_sel,text="\n",font=small_font)
    tkconfigure(ann_sel,text="\n",font=small_font)
    tkconfigure(for_reads,text="\n\n",font=file_display_font) 
    tkconfigure(rev_reads,text="\n\n",font=file_display_font) 
    
    activate()
    load_parms()
    create_save_dirs()
    tkdestroy(exp_menu)
  }
  
  #Select existing experiment function
  sel_ex_fun<-function(){
    .GlobalEnv$experiment_name_cur<-tk_select.list(choices=list.dirs(paste(directory,'/Experiments',sep=""),full.names = FALSE,recursive = FALSE))
    if(experiment_name_cur==""){
      tk_messageBox(message="You did not select an experiment!")
      open_exp_func()
    } else{
      .GlobalEnv$exp_directory<-paste(directory,"/Experiments","/",experiment_name_cur,sep="")
      print(exp_directory)
      tkconfigure(cur_exp,text=paste("Experiment: ",str_sub(experiment_name_cur,1,18),"...",sep=""))
      activate()
      load_parms() 
    }
    create_save_dirs()
    tkdestroy(exp_menu)
  }
  
  exp_menu<-tktoplevel()
  tkwm.geometry(exp_menu,"250x140+200+100")
  tkwm.title(exp_menu,"DEAR-main")
  exp_frm<-tkframe(exp_menu)
  tkgrid(exp_frm,column=1,row=1,sticky="n",pady=10,padx=15)
  exp_op<-tklabel(exp_frm,text="Select an experiment option: ",font=header_font)
  tkgrid(exp_op)
  create_exp_but<-tkbutton(exp_frm,text="Create new experiment",command=create_new_exp_fun)
  tkgrid(create_exp_but,pady=10)
  sel_ex_but<-tkbutton(exp_frm,text="Select existing experiment",command=sel_ex_fun)
  tkgrid(sel_ex_but)
  tkwait.window(exp_menu)
  
  tkdestroy(exp_menu)
}

#Install files
install_files_func<-function(){
  
  #Install fastq function
    #http://ftp.sra.ebi.ac.uk/vol1/fastq/<accession-prefix>/<00-last-digit-of-full-accession>/<full-accession>
    #https://ena-docs.readthedocs.io/en/latest/faq/archive-generated-files.html?#archive-generated-files
  download_fastq<-function(){
      
      fastq_entry_file<-paste(exp_directory,"/fastq_entry.txt",sep="")
      if(file.exists(fastq_entry_file)==FALSE){
        #Prompt user to enter accession numbers
        tmp<-matrix(nrow=1,ncol=2)
        tmp[1,1]<-"SRR10820656"
        tmp[1,2]<-"SRR10820663"
        write.table(tmp,file=fastq_entry_file,row.names = FALSE,sep="\n",quote=FALSE,col.names = FALSE)
        shell.exec(fastq_entry_file)
      } else{
        shell.exec(fastq_entry_file)
      }
      cont<-tk_messageBox(message="Enter accession numbers, separated by a new line, and save your entry. Then click OK to continue.",type="okcancel")
      
      if(cont=="ok"){
        options(timeout = 10000)
        cat("\n",file=paste(exp_directory,"/fastq_entry.txt",sep=""),append=TRUE)
        .GlobalEnv$entered_fastq<-unlist(unname(read.delim2(paste(exp_directory,"/fastq_entry.txt",sep=""),header=FALSE,sep="\n")))

        #Remove fastq already in folder
        already_download<-str_sub(list.files(forward_dir,pattern=".fastq.gz"),1,11)
        if(length(already_download)>0){
          tk_messageBox(message=paste("Skipping: ",paste(already_download,collapse="\n"),sep="\n"))
          .GlobalEnv$entered_fastq<-entered_fastq[-which(entered_fastq%in%already_download)]
        }
        
        if(length(which(duplicated(entered_fastq)==TRUE))==0){
          .GlobalEnv$read_type<-tk_select.list(choices=c("Single-end","Paired-end"),title="Select read type.")
          if(read_type!=""){
            #Tmp bar
            .GlobalEnv$tmp_inst_inc<-80/(length(entered_fastq))
            .GlobalEnv$tmp_instal_val<-0
            tmp_instal<-winProgressBar(title="FASTQ Installation",
                                       label="Installing files...",
                                       min=0,max=100,width=300,
                                       initial=tmp_instal_val)
            
            lapply(entered_fastq,function(z){
              .GlobalEnv$tmp_instal_val<-tmp_instal_val+tmp_inst_inc
              setWinProgressBar(tmp_instal,value=tmp_instal_val,label=paste("Installing ",z,".fastq.gz...",sep=""))
              if(read_type=="Paired-end"){
                acc_prefix<-str_sub(z,0,6)
                last_dig<-paste("0",str_sub(z,-2),sep="")
                link<-paste("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/",acc_prefix,"/",last_dig,"/",z,sep="")
                
                file1<-paste(link,"/",z,"_1.fastq.gz",sep="")
                file2<-paste(link,"/",z,"_2.fastq.gz",sep="")
                lapply(c(file1,file2),function(x){
                  file_name<-gsub(paste(dirname(x),"/",sep=""),"",x)
                  if(str_detect(x,"_1.fastq.gz")==TRUE){
                    #If forward read
                    new_dir<-paste(forward_dir,"/",file_name,sep="")
                  } else if(str_detect(x,"_2.fastq.gz")==TRUE){
                    new_dir<-paste(reverse_dir,"/",file_name,sep="")
                  } else{
                    tk_messageBox(message="Expected forward/reverse format not detected!")
                    new_dir<-paste(fastq_dir,"/",file_name,sep="")
                  }
                  tryCatch(download.file(x,destfile=new_dir,quiet = TRUE,mode="wb"),error=function(e){
                    tk_messageBox(message="Error during download. Ensure that your file is accessible at the ENA.")
                    close(tmp_instal)
                    stop()
                  })
                })
              } else if(read_type=="Single-end"){
                #Format download link
                acc_prefix<-str_sub(z,0,6)
                last_dig<-paste("00",str_sub(z,-1),sep="")
                last_two_dig<-paste("0",str_sub(z,-2),sep="")
                link<-paste("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/",acc_prefix,"/",last_dig,"/",z,sep="")
                link_alt<-paste("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/",acc_prefix,"/",last_two_dig,"/",z,sep="")
                .GlobalEnv$file1<-paste(link,"/",z,".fastq.gz",sep="")
                .GlobalEnv$file1_alt<-paste(link_alt,"/",z,"_1.fastq.gz",sep="")
                
                file_name<-paste(z,".fastq.gz",sep="")
                new_dir<-paste(forward_dir,"/",file_name,sep="")
                
                tryCatch(download.file(file1,destfile=new_dir,quiet = TRUE,mode="wb"),error=function(e){
                  #Try alternative link with underscore
                  tryCatch(download.file(file1_alt,destfile=new_dir,quiet = TRUE,mode="wb"),error=function(w){
                    tk_messageBox(message="Error during download. Ensure that your file is accessible at the ENA.")
                    close(tmp_instal)
                    stop()
                  })
                })
              } 
            })
            
            close(tmp_instal) 
            tk_messageBox(message = "Download complete!")
          } else{
            tk_messageBox(message="You must select a read type! Download cancelled.")
          }
        } else{
          tk_messageBox(message="You have duplicate reads! Download cancelled.")
        }
      } else{
        
      }
      options(timeout = 60)
    }
      
  browse_fastq<-function(){
    browseURL("https://www.ncbi.nlm.nih.gov/sra")
  }
  
  manage_fastq<-function(){
    shell.exec(fastq_dir)
  }
  
  download_gen<-function(){
    
    ass_var<-tclVar("GCF_000005845.2")
    ass_var_2<-tclVar("ASM584v2")
    ass_entry_name<-tclVar("ecoli_k12")
    assembly_enr_func<-function(){
      .GlobalEnv$ref_seq<-tclvalue(ass_var)
      print(ref_seq)
      .GlobalEnv$ass_id<-tclvalue(ass_var_2)
      print(ass_id)
      .GlobalEnv$named_ass<-tclvalue(ass_entry_name)
      print(named_ass)
      
      if(ref_seq!=""&&ass_id!=""&&named_ass!=""){
        tkdestroy(assembly_entry)
        digits<-str_replace_all(ref_seq,"[^[:digit:]]","")
        first<-str_sub(digits,0,3)
        sec<-str_sub(digits,4,6)
        thir<-str_sub(digits,7,9)
        
        #Get FTP link to genome assembly
        #ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/308/975/GCF_000308975.1_ASM30897v2/GCF_000308975.1_ASM30897v2_genomic.fna.gz
        .GlobalEnv$base<-paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/",first,"/",sec,"/",thir,"/",ref_seq,"_",ass_id,"/",ref_seq,"_",ass_id,"_genomic.fna.gz",sep="")
        
        if(file.exists(paste(gen_dir,"/",named_ass,".fna.gz",sep=""))==FALSE){
          tmp_prog<-winProgressBar(title="Genome Download",
                                   label="Downloading genome...",
                                   min=0,max=100,width=300,initial = 50)
          
          options(timeout = 10000)
          tryCatch(download.file(base,destfile=paste(gen_dir,"/",named_ass,".fna.gz",sep=""),quiet = TRUE),
                   error=function(e){
                     tk_messageBox(message="Error downloading genome. Download cancelled.")
                     close(tmp_prog)
                     stop()
                   })
          close(tmp_prog)
          shell.exec(gen_dir)
          options(timeout = 60)
        } else{
          tk_messageBox(message="You already have a genome assembly with that name!")
        }
        
        #Get FTP link to genome annotation
        #ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/308/975/GCF_000308975.1_ASM30897v2/GCF_000308975.1_ASM30897v2_genomic.fna.gz
        .GlobalEnv$base2<-paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/",first,"/",sec,"/",thir,"/",ref_seq,"_",ass_id,"/",ref_seq,"_",ass_id,"_genomic.gtf.gz",sep="")
        if(file.exists(paste(annot_dir,"/",named_ass,"_genomic.gtf.gz",sep=""))==FALSE){
          tmp_prog<-winProgressBar(title="Annotation Download",
                                   label="Downloading genome annotation...",
                                   min=0,max=100,width=300,initial = 50)
          
          options(timeout = 10000)
          tryCatch(download.file(base2,destfile=paste(annot_dir,"/",named_ass,"_genomic.gtf.gz",sep=""),quiet = TRUE),
                   error=function(e){
                     tk_messageBox(message="Error downloading annotation. Download cancelled.")
                     close(tmp_prog)
                     stop()
                   })
          close(tmp_prog)
          shell.exec(annot_dir)
          options(timeout = 60)
        } else{
          tk_messageBox(message="You already have a genome annotation with that name!")
        }
        
      } else{
        tk_messageBox(message="You must fill out all fields!")
      }

    }
    
    assembly_entry<-tktoplevel()
    tkwm.geometry(assembly_entry,"160x160+400+100")
    tkwm.title(assembly_entry,"Assembly Entry")
    ass_frame<-tkframe(assembly_entry)
    tkgrid(ass_frame)
    ass_lbl<-tklabel(ass_frame,text="Assembly RefSeq:")
    tkgrid(ass_lbl)
    ass_entry<-tkentry(ass_frame,textvariable=ass_var,width=20)
    tkgrid(ass_entry)
    ass_lbl_2<-tklabel(ass_frame,text="Assembly name (no spaces):")
    tkgrid(ass_lbl_2)
    ass_entry_2<-tkentry(ass_frame,textvariable=ass_var_2)
    tkgrid(ass_entry_2)
    ass_name<-tklabel(ass_frame,text="Name the download file:")
    tkgrid(ass_name)
    ass_entry_3<-tkentry(ass_frame,textvariable=ass_entry_name,width=20)
    tkgrid(ass_entry_3)
    ass_but<-tkbutton(ass_frame,text="Enter",command=assembly_enr_func)
    tkgrid(ass_but,pady=5)
    
  }
  
  browse_gen<-function(){
    browseURL("https://www.ncbi.nlm.nih.gov/data-hub/genome/")
  }
  
  manage_gen<-function(){
    shell.exec(paste(directory,"/Files",sep=""))
  }
  
    down_fils<-tktoplevel()
    tkwm.title(down_fils,"Install files")
    tkwm.geometry(down_fils,"270x200+200+100")
    #Frame
    install_frame<-tkframe(down_fils)
    tkgrid(install_frame,padx=10)
    #Title
    down_title<-tklabel(install_frame,text="File Manager",font=header_font)
    tkgrid(down_title,pady=5,padx=20,column=1,columnspan=3,row=1)
    #Fastqlbl
    lbl_fastq<-tklabel(install_frame,text="SEQUENCES (.fastq.gz)",font=underline_font)
    tkgrid(lbl_fastq,pady=5,row=2,column=1,columnspan=3)
    #Download button
    fastq_download<-tkbutton(install_frame,text="Install",command=download_fastq,font=normal_font)
    tkgrid(fastq_download,pady=2,row=3,column=1,sticky="e",padx=3)
    #Browse button
    fastq_browse<-tkbutton(install_frame,text="Browse",command=browse_fastq,font=normal_font)
    tkgrid(fastq_browse,pady=2,row=3,column=2,padx=3)
    #Manage button
    fastq_manage<-tkbutton(install_frame,text="View",font=normal_font,command=manage_fastq)
    tkgrid(fastq_manage,pady=2,row=3,column=3,sticky="w",padx=3)
    #Gen lbl
    lbl_gen<-tklabel(install_frame,text="GENOMES (.fna) & ANNOTATIONS (.gtf)",font=underline_font)
    tkgrid(lbl_gen,pady=5,row=4,column=1,columnspan=3)
    #Download
    gen_download<-tkbutton(install_frame,text="Install",font=normal_font,command=download_gen)
    tkgrid(gen_download,pady=2,row=5,column=1,sticky="e",padx=3)
    #Browse
    gen_browse<-tkbutton(install_frame,text="Browse",font=normal_font,command=browse_gen)
    tkgrid(gen_browse,pady=2,row=5,column=2,padx=3)
    #Manage button
    gen_manage<-tkbutton(install_frame,text="View",font=normal_font,command=manage_gen)
    tkgrid(gen_manage,pady=2,row=5,column=3,sticky="w",padx=3)
    tkwait.window(down_fils)
    
}

#Get forward reads
forward_reads<-function(){
  .GlobalEnv$forward_files<-tk_choose.files(default=paste(forward_dir,"/SELECT_A_FILE",sep=""),multi=TRUE,
                                            filters=matrix(c(".fastq.gz",".fastq.gz",
                                            ".fastq",".fastq",
                                            "ALL","*"),3,2,byrow = TRUE),
                                            caption="Select forward reads.")
  
  if(length(forward_files)>0){
    just_names<-lapply(1:length(forward_files),function(x){
      tmp<-forward_files[x]
      new_name<-gsub(paste(dirname(tmp),"/",sep=""),"",tmp)
    })
    .GlobalEnv$just_f_names<-unlist(just_names)
    
    if(length(just_f_names)<4){
      print_names<-just_f_names[1:3]
      print_names[which(is.na(print_names))]<-" "
      tkconfigure(for_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
    } else{
      print_names<-just_f_names[1:2]
      print_names<-c(print_names,"Additional files not shown...")
      tkconfigure(for_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
    }
    
    #Update parms
    parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
    parms$Forward_Files<-paste(forward_files,collapse="+")
    write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
    
  } else{
    tkconfigure(for_reads,text="\n\n") 
  }
}

#Allow reverse reads
allow_reverse_reads<-function(){
  cur_reverse_reads<-tclvalue(reverse_reads)
  print(cur_reverse_reads)
  
  if(cur_reverse_reads=="1"){
    tkconfigure(sel_r_fastq,state="normal")
    
  } else if(cur_reverse_reads=="0"){
    tkconfigure(sel_r_fastq,state="disabled")
  }
  
  #Update parms
  tkconfigure(rev_reads,text="\n\n",font=file_display_font) 
  parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
  parms$Pairwise<-cur_reverse_reads
  parms$Reverse_Files<-NA
  write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
  
}

#Get reverse reads
get_reverse_reads<-function(){
  .GlobalEnv$reverse_files<-tk_choose.files(default=paste(reverse_dir,"/SELECT_A_FILE",sep=""),multi=TRUE,
                                            filters=matrix(c(".fastq.gz",".fastq.gz",
                                            ".fastq",".fastq",
                                            "ALL","*"),3,2,byrow = TRUE),
                                            caption="Select reverse reads.")
  
  if(length(reverse_files)>0){
    if(length(reverse_files)==length(forward_files)){
      just_names<-lapply(1:length(reverse_files),function(x){
        tmp<-reverse_files[x]
        new_name<-gsub(paste(dirname(tmp),"/",sep=""),"",tmp)
      })
      .GlobalEnv$just_r_names<-unlist(just_names)
      
      if(length(just_r_names)<4){
        print_names<-just_r_names[1:3]
        print_names[which(is.na(print_names))]<-" "
        tkconfigure(rev_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
      } else{
        print_names<-just_r_names[1:2]
        print_names<-c(print_names,"Additional files not shown...")
        tkconfigure(rev_reads,text=paste(print_names,collapse="\n"),font=file_display_font) 
      }
      
      #Update parms
      parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
      parms$Reverse_Files<-paste(reverse_files,collapse="+")
      write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
    } else{
      tk_messageBox(message="You must have the same number of forward and reverse reads!")
      tkconfigure(rev_reads,text="\n\n") 
    }
  } else{
    tkconfigure(rev_reads,text="\n\n") 
  }
}

#Load genome file
load_gen_file<-function(){
    .GlobalEnv$gen_file<-tk_choose.files(default=paste(gen_dir,"/SELECT_A_FILE",sep=""),multi=FALSE,filters=matrix(c("ALL","*",
                                                                                                                     ".fa.gz",".fa.gz",
                                                                                                                     ".fna.gz",".fna.gz",
                                                                                                                     ".fna",".fna",
                                                                                                                     ".fa",".fa"),5,2,byrow = TRUE),
                                         caption="Select a genome file.")
    if(length(gen_file)>0){
      
      .GlobalEnv$gen_name<-paste(str_sub(gsub(paste(dirname(gen_file),"/",sep=""),"",gen_file),start=0,end=18),"...",sep="")
      tkconfigure(gen_sel,text=paste(gen_name,"\n",sep=""),font=small_font)
      
      #Update parms
      parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
      parms$Reference_Genome<-gen_file
      write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
   
    } else{
      tk_messageBox(message="You didn't select a file!")
      tkconfigure(gen_sel,text="\n")
    }
  
  
}

#Load annotation file
load_ann_file<-function(){
  .GlobalEnv$ann_file<-tk_choose.files(default=paste(annot_dir,"/SELECT_A_FILE",sep=""),multi=FALSE,filters=matrix(c(".gtf.gz",".gtf.gz",
                                                                                        ".gtf",".gtf",
                                                                                        "ALL","*"),3,2,byrow = TRUE),
                                       caption="Select an annotation file.")
  if(length(ann_file)>0){
    
    .GlobalEnv$ann_name<-paste(str_sub(gsub(paste(dirname(ann_file),"/",sep=""),"",ann_file),start=0,end=18),"...",sep="")
    tkconfigure(ann_sel,text=paste(ann_name,"\n",sep=""),font=small_font)
    
    #Update parms
    parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
    parms$Genome_Annotation<-ann_file
    write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
    read_gtf_function()
  } else{
    tk_messageBox(message="You didn't select a file!")
    tkconfigure(ann_sel,text="\n")
  }
}

#Edit design matrix
edit_design<-function(){
  
  if(length(forward_files)>1){
    
    #Read existing file
    if(file.exists(paste(exp_directory,"/Metadata.csv",sep=""))==TRUE){
      get_meta<-read.csv(paste(exp_directory,"/Metadata.csv",sep=""))
      cur_classes<-get_meta$Classification
    }
    
    ann_reads_func<-function(){
      values<-lapply(1:length(forward_files),function(x){
        value<-tclvalue(eval(parse(text=paste("check",x,sep=""))))
      })
      .GlobalEnv$value_cur<-unlist(values)
      print(value_cur)
      
      if(" "%in%value_cur){
        tk_messageBox(message="You must annotate all reads!")
      } else{
        if(length(unique(value_cur))>1){
          #Create metadata
          meta_df<-data.frame(Reads=just_f_names,
                              Filepath=forward_files,
                              Classification=value_cur)
          convert<-data.frame(Classification=group_sel_choices,
                              GROUP=c(0,1))
          meta_df_cur<-left_join(meta_df,convert,by="Classification")
          .GlobalEnv$meta_file<-meta_df_cur
          
          #Write csv
          write.csv(meta_df_cur,paste(exp_directory,"/Metadata.csv",sep=""),row.names = FALSE)
          
          tkdestroy(annotate)
        } else{
          tk_messageBox(message="You must have at least 1 experimental and control condition.")
        }
      }
      
    }
    
    adjust_group_fun<-function(){
      
      control_var<-tclVar(group_sel_choices[1])
      exp_var<-tclVar(group_sel_choices[2])
      
      enter_fun<-function(){
        cont_val<-tclvalue(control_var)
        exp_val<-tclvalue(exp_var)
        .GlobalEnv$group_sel_choices<-c(cont_val,exp_val)
        print(group_sel_choices)
        
        #Save to parms
        parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
        parms$Groups<-paste(group_sel_choices,collapse="+")
        write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
        
        tkdestroy(adj_grp)
        tkdestroy(annotate)
        edit_design()
      }
      
      adj_grp<-tktoplevel()
      tkwm.geometry(adj_grp,"260x160+300+50")
      tkwm.title(adj_grp,"Adjust group names")
      adj_frame<-tkframe(adj_grp)
      tkgrid(adj_frame)
      grp_title<-tklabel(adj_frame,text="Group",font=header_font)
      tkgrid(grp_title,column=1,row=1,padx=10,pady=10)
      name_title<-tklabel(adj_frame,text="Name",font=header_font)
      tkgrid(name_title,column=2,row=1,padx=10,pady=10)
      control_lbl<-tklabel(adj_frame,text="Control",font=normal_font)
      tkgrid(control_lbl,column=1,row=2,padx=10,pady=5,sticky="w")
      exp_lbl<-tklabel(adj_frame,text="Experimental",font=normal_font)
      tkgrid(exp_lbl,column=1,row=3,padx=10,pady=5,sticky="w")
      control_entry<-tkentry(adj_frame,width=20,textvariable=control_var)
      tkgrid(control_entry,column=2,row=2,padx=10,pady=5)
      exp_entry<-tkentry(adj_frame,width=20,textvariable=exp_var)
      tkgrid(exp_entry,column=2,row=3,padx=10,pady=5)
      save_but<-tkbutton(adj_frame,text="Enter",command=enter_fun)
      tkgrid(save_but,column=1,columnspan=2,row=4,pady=5)
      
    }
    
    annotate<-tktoplevel()
    tkwm.geometry(annotate,paste("280x",70+32*length(just_f_names),"+200+100",sep=""))
    tkwm.title(annotate,"Annotate reads")
    top_frame<-tkframe(annotate)
    tkgrid(top_frame)
    top_lbl<-tklabel(top_frame,text="Annotate reads:",font=header_font)
    tkgrid(top_lbl,padx=50,row=1,column=1,columnspan=2)
    lapply(1:length(just_f_names),function(x){
      tmp<-just_f_names[x]
      if(str_count(tmp)>22){
        tmp<-str_sub(tmp,0,22)
      }
       
      #Label
      assign(paste("read",x,sep=""),
             tklabel(top_frame,text=tmp))
      tkgrid(eval(parse(text=paste("read",x,sep=""))),column=1,row=1+x,padx=10)
      #Combobox
      if(file.exists(paste(exp_directory,"/Metadata.csv",sep=""))==TRUE){
        assign(paste("group",x,sep=""),
               ttkcombobox(top_frame,values=group_sel_choices,textvariable=assign(paste("check",x,sep=""),tclVar(cur_classes[x]),envir = .GlobalEnv),width=12))
        
      } else{
        assign(paste("group",x,sep=""),
               ttkcombobox(top_frame,values=group_sel_choices,textvariable=assign(paste("check",x,sep=""),tclVar(" "),envir = .GlobalEnv),width=12))
        
      }
      tkgrid(eval(parse(text=paste("group",x,sep=""))),column=2,row=1+x,padx=10,pady=5)
    })
    #Enter button
    enter_ann<-tkbutton(top_frame,text="Save annotation",command=ann_reads_func)
    tkgrid(enter_ann,column=1,row=2+length(just_f_names),pady=5)
    #Change group label
    group_lbl<-tkbutton(top_frame,text="Adjust group name",command=adjust_group_fun)
    tkgrid(group_lbl,column=2,row=2+length(just_f_names),pady=5)
    tkwait.window(annotate)
    
  } else{
    tk_messageBox(message="You do not have enough reads! More than one reads are required.")
  }
  
  if(file.exists(paste(exp_directory,"/Metadata.csv",sep=""))==TRUE){
    tkconfigure(meta_disp,text="Read annotation detected!")
  }
  
}

#View plots
view_plot_fun<-function(){
  shell.exec(exp_directory)
}

#Run analysis
run_dear<-function(){
  continue<-FALSE
  scrape_parms()

  #Check for all inputs
  if(length(forward_files)>0){
    print(forward_files)
    
    if(length(forward_files)>=as.numeric(sample_val)){
      if(cur_reverse_reads=="1"){
        
        if(length(reverse_files)>0){
          
          if(length(reverse_files)==length(forward_files)){
            print(reverse_files)
            continue<-TRUE
            .GlobalEnv$paired_end<-TRUE
            
          } else{
            tk_messageBox(message="You must have the same number of forward and reverse reads!") 
          }
        } else{
          tk_messageBox(message="You have not selected any reverse reads!")
        }
        
      } else if(cur_reverse_reads=="0"){
        continue<-TRUE
        .GlobalEnv$paired_end<-FALSE
        print("No reverse reads")
        #No reverse reads
      }
    } else{
      tk_messageBox(message="Ensure that your Sample value is not larger than the number of entered reads.\n\nSample value restricts which genes are included in the analysis. Only genes with at least as many reads as the entered value will be included.")
    }
  } else{
    tk_messageBox(message="You have not selected any forward reads!")
    #No fastq files
  }
 
  if(continue==TRUE){
    if(length(gen_file)>0){
      if(length(ann_file)>0){
        if(file.exists(paste(exp_directory,"/Metadata.csv",sep=""))==TRUE){
          
          prog<-winProgressBar(title="Differential Expression Analysis in R",
                               label="Installing packages (one-time installation)...",
                               min=0,max=100,width=300)
          
          .GlobalEnv$re_count_stat<-tclvalue(recount_ft_var)
          if(re_count_stat=="1"&&file.exists(paste(fastq_dir,"/rawfeaturecounts.csv",sep=""))==TRUE){
            file.remove(paste(fastq_dir,"/rawfeaturecounts.csv",sep=""))
          }
          
          .GlobalEnv$prog_val<-33
          setWinProgressBar(prog,value=prog_val,label=paste("Building ",gen_name," index...",sep=""))
          #Build index
          #Check for index annotation
          ann_name<-gsub(paste(dirname(ann_file),"/",sep=""),"",ann_file)
          if(length(list.files(annot_dir,pattern=paste(ann_name,".log",sep="")))>0){
            tk_messageBox(message="Genome index already detected. Skipping indexing.")
          } else{
            tryCatch(source(paste(directory,"/scripts/1buildindex.R",sep="")),error=function(e){
              tk_messageBox(message="Error reading 1buildindex.R!")
              close(prog)
              stop()
            })
          }
          
          .GlobalEnv$prog<-prog
          .GlobalEnv$prog_inc<-(66-prog_val)/length(forward_files)
          #Align reads
          
          #Restrict forward files to those that haven't been analyzed
          bam_files<-list.files(path=forward_dir,pattern=".BAM$")
          tryCatch(detected_bams<-str_detect(bam_files,just_f_names),
                   warning=function(w){
                     tk_messageBox(message="There are too many .BAM files! Stopping analysis.")
                     close(prog)
                     stop()
                   })
          if(length(detected_bams)==0){
            #If none fastq aligned
            tryCatch(source(paste(directory,"/scripts/2alignreads.R",sep="")),error=function(e){
              tk_messageBox(message="Error reading 2alignreads.R!")
              close(prog)
              stop()
            }) 
          } else if(length(detected_bams)>0&&length(which(detected_bams==TRUE))==length(forward_files)){
            #If all fastq already aligned
            tk_messageBox(message=".BAM files already detected! Skipping read alignments.") 
          } else if(length(detected_bams)>0&&length(which(detected_bams==TRUE))<length(forward_files)){
            #If only some fastq already aligned
            tk_messageBox(message="Some .BAM files are missing! Only these will be generated!")
            .GlobalEnv$forward_files<-forward_files[which(detected_bams==FALSE)]
            .GlobalEnv$just_f_names<-just_f_names[which(detected_bams==FALSE)]
            
            if(cur_reverse_reads=="1"){
              .GlobalEnv$reverse_files<-reverse_files[which(detected_bams==FALSE)]
              .GlobalEnv$just_r_names<-just_r_names[which(detected_bams==FALSE)]
            }
            
            tryCatch(source(paste(directory,"/scripts/2alignreads.R",sep="")),error=function(e){
              tk_messageBox(message="Error reading 2alignreads.R!")
              close(prog)
              stop()
            })
          } 
          
          setWinProgressBar(prog,value=90,label="Counting features...")
          #DEA
          tryCatch(source(paste(directory,"/scripts/3deanalysis.R",sep="")),error=function(e){
            tk_messageBox(message="Error reading 3deanalysis.R!")
            close(prog)
          })
          close(prog)
        } else{
          #No read annotation
          tk_messageBox(message="You must annotate your reads!")
        }
      } else{
        #No ann file
        tk_messageBox(message="You have not selected an annotation file!")
      }
    } else{
     #No gen file
      tk_messageBox(message="You have not selected a genome!")
    }
  }
   
}

#Open plot builder
activate_plots<-function(){
  scrape_parms()
  tryCatch(source(paste(this.dir(),"/4plotbuilder.R",sep="")),
           error=function(e){
             tk_messageBox(message="Error reading 4plotbuilder.R!")
             close(calc_prog)
             stop()
           })
}

#Scrape parms
scrape_parms<-function(){
  #Update parms
  parms<-read.delim(paste(exp_directory,"/Parameters.txt",sep=""),sep=",")
  #Save feature type
  .GlobalEnv$ft_val<-tclvalue(ft_var)
  parms$FtType<-ft_val
  print(ft_val)
  #Save attributed type
  .GlobalEnv$att_val<-tclvalue(att_var)
  parms$AttType<-att_val
  print(att_val)
  #Save threshold
  .GlobalEnv$thresh_val<-tclvalue(thresh_var)
  parms$Threshold<-thresh_val
  print(thresh_val)
  #Save sample parm
  .GlobalEnv$sample_val<-tclvalue(sample_var)
  parms$Sample<-sample_val
  print(sample_val)
  #Sample sequence type 
  .GlobalEnv$seq_val<-tclvalue(seq_var)
  parms$SeqType<-seq_val
  print(seq_val)
  #Save reverse reads
  .GlobalEnv$cur_reverse_reads<-tclvalue(reverse_reads)
  parms$Pairwise<-cur_reverse_reads
  print(cur_reverse_reads)
  if(cur_reverse_reads=="1"){
    .GlobalEnv$paired_end<-TRUE
  } else{
    .GlobalEnv$paired_end<-FALSE
  }
  
  write.table(parms,file=paste(exp_directory,"/Parameters.txt",sep=""),sep=",",eol="\n",quote = FALSE,row.names = FALSE)
}

#Top level
top<-tktoplevel()
tkwm.geometry(top,"650x490+100+100")
tkwm.title(top,"DEAR-main V1.0.0")
#Frame
titleframe<-tkframe(top)
tkgrid(titleframe,row=1,column=1,columnspan=10,padx=150,pady=10)
#Title
title<-tklabel(titleframe,text="DEAR-main V1.0.0",font=title_font)
tkgrid(title,row=1,column=1,padx=30)
#Configure experiment
config_exp<-tkbutton(titleframe,text="Open Experiment",font=header_font,borderwidth=3,command=open_exp_func)
tkgrid(config_exp,pady=5,column=3,row=1)
#Current experiment label
cur_exp<-tklabel(titleframe,text="Experiment:",font=normal_font)
tkgrid(cur_exp,column=1,row=2,columnspan=2,sticky="w")
#Install fastq files
inst_files<-tkbutton(titleframe,text="Manage files",font=small_font,borderwidth=1,relief="raised",state="disabled",command=install_files_func)
tkgrid(inst_files,column=2,row=2,columnspan=3,sticky="w")
#View plots
view_plots<-tkbutton(titleframe,text="View files",font=small_font,borderwidth=1,relief="raised",state="disabled",command=view_plot_fun)
tkgrid(view_plots,column=3,row=2,columnspan=2,sticky="e")

#Left frame
left_frame<-tkframe(top,borderwidth=3,relief="raised")
tkgrid(left_frame,row=2,column=1,padx=20)
#File management label
file_man<-tklabel(left_frame,text="Load Sequence Files",font=header_font)
tkgrid(file_man)
#Fastq files label
fastq_lbl<-tklabel(left_frame,text="FASTQ FILES",font=underline_font)
tkgrid(fastq_lbl,pady=10)
#Pairwise reads
pairwise<-tkcheckbutton(left_frame,text="Paired end reads?",font=normal_font,command=allow_reverse_reads,variable=reverse_reads,state="disabled")
tkgrid(pairwise,pady=0)
#Select .fastq files
sel_f_fastq<-tkbutton(left_frame,text="Select foward reads (.fastq)",font=normal_font,command=forward_reads,state="disabled")
tkgrid(sel_f_fastq,pady=5)
#Select .fastq files
sel_r_fastq<-tkbutton(left_frame,text="Select reverse reads (.fastq)",font=normal_font,state="disabled",command=get_reverse_reads)
tkgrid(sel_r_fastq,pady=5,padx=10)
#Genome files
gen_files<-tklabel(left_frame,text="GENOME FILES",font=underline_font)
tkgrid(gen_files,pady=10)
#Load genome
load_gen<-tkbutton(left_frame,text="Load genome (.fa/.fna)",font=normal_font,command=load_gen_file,state="disabled")
tkgrid(load_gen,pady=5)
#Genome annotations
gen_ann<-tklabel(left_frame,text="GENOME ANNOTATIONS",font=underline_font)
tkgrid(gen_ann,pady=10)
#Load annot
load_ann<-tkbutton(left_frame,text="Load annotation (.gtf)",font=normal_font,command=load_ann_file,state="disabled")
tkgrid(load_ann,pady=5)

#Middle frame
mid_frame<-tkframe(top,borderwidth=3,relief="raised")
tkgrid(mid_frame,row=2,column=2,padx=5)
#Options label
file_man<-tklabel(mid_frame,text="Parameters",font=header_font,width=19)
tkgrid(file_man)
#Options
opt<-tklabel(mid_frame,text="OPTIONS",font=underline_font)
tkgrid(opt,pady=10)
#Sequence type
seq_lbl<-tklabel(mid_frame,text="Seq. type",font=normal_font)
tkgrid(seq_lbl,pady=5,sticky="w",padx=5,row=2)
sq_type<-ttkcombobox(mid_frame,values=c("rna","dna"),font=normal_font,width=7,textvariable=seq_var)
tkgrid(sq_type,pady=5,sticky="e",padx=5,row=2)
#Feature type
ft_lbl<-tklabel(mid_frame,text="Feature type",font=normal_font)
tkgrid(ft_lbl,pady=5,sticky="w",padx=5,row=3)
ft_type<-ttkcombobox(mid_frame,values="gene",font=normal_font,width=7,textvariable=ft_var)
tkgrid(ft_type,pady=5,sticky="e",padx=5,row=3)
#Attribute type
att_lbl<-tklabel(mid_frame,text="Attribute type",font=normal_font)
tkgrid(att_lbl,pady=5,sticky="w",padx=5,row=4)
att_type<-ttkcombobox(mid_frame,values="gene_id",font=normal_font,width=7,textvariable=att_var)
tkgrid(att_type,pady=5,sticky="e",padx=5,row=4)
#Threshold value
thresh_lbl<-tklabel(mid_frame,text="Threshold",font=normal_font)
tkgrid(thresh_lbl,pady=5,sticky="w",padx=5,row=5)
thresh_entry<-tkentry(mid_frame,textvariable=thresh_var,width=11)
tkgrid(thresh_entry,pady=5,sticky="e",padx=5,row=5)
#Sample value
sample_lbl<-tklabel(mid_frame,text="Sample value",font=normal_font)
tkgrid(sample_lbl,pady=5,sticky="w",padx=5,row=6)
sample_entry<-tkentry(mid_frame,textvariable=sample_var,width=11)
tkgrid(sample_entry,pady=5,sticky="e",padx=5,row=6)
#Configure matrices
meta_lbl<-tklabel(mid_frame,text="CONFIGURE INPUTS",font=underline_font)
tkgrid(meta_lbl,pady=10)
#Metadata
meta_data<-tkbutton(mid_frame,text="Annotate reads",font=normal_font,command=edit_design,state="disabled")
tkgrid(meta_data,pady=2)
#Re-count features
recount_ft<-tkcheckbutton(mid_frame,text="Re-count features?",font=normal_font,state="disabled",variable=recount_ft_var)
tkgrid(recount_ft,pady=2)

#RUN ANALYSIS
run<-tkframe(top,borderwidth=3,padx=170)
tkgrid(run,row=3,column=1,columnspan=2,pady=10,sticky="e")
run_but<-tkbutton(run,text="RUN DEAR",font=header_font,borderwidth=3,command=run_dear,state="disabled")
tkgrid(run_but)
#Build plots
build<-tkframe(top,borderwidth=3,padx=0)
tkgrid(build,row=3,column=2,columnspan=2,padx=0)
build_plot<-tkbutton(build,text="ANALYZE COUNTS",font=header_font,borderwidth=3,state="disabled",command=activate_plots)
tkgrid(build_plot)

#Preview selections
preview_frame<-tkframe(top,borderwidth=3,relief="raised")
tkgrid(preview_frame,column=3,row=2,padx=20)
#Preview label
prev_lbl<-tklabel(preview_frame,text="Preview selections",font=header_font)
tkgrid(prev_lbl,padx=10)
#Prev forward reads
for_lbl<-tklabel(preview_frame,text="Forward reads",font=normal_font)
tkgrid(for_lbl,pady=5)
for_reads<-tklabel(preview_frame,text="\n\n",font=file_display_font)
tkgrid(for_reads)
#Prev reverse reads
rev_lbl<-tklabel(preview_frame,text="Reverse reads",font=normal_font)
tkgrid(rev_lbl,pady=5)
rev_reads<-tklabel(preview_frame,text="\n\n",font=file_display_font)
tkgrid(rev_reads)
#Prev genome reads
sel_gen<-tklabel(preview_frame,text="Genome file",font=normal_font)
tkgrid(sel_gen,pady=5)
gen_sel<-tklabel(preview_frame,text="\n",font=small_font)
tkgrid(gen_sel)
#Prev annotation reads
sel_ann<-tklabel(preview_frame,text="Genome annotation",font=normal_font)
tkgrid(sel_ann,pady=5)
ann_sel<-tklabel(preview_frame,text="\n",font=small_font)
tkgrid(ann_sel)
#Read annotation
meta_disp<-tklabel(preview_frame,text="NO read annotation",font=normal_font)
tkgrid(meta_disp,pady=7)
#Wait for input
tkwait.window(top)
