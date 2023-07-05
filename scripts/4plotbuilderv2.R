

#--------------------RUN REGARDLESS

.GlobalEnv$design_path<-paste(exp_directory,"/Metadata.csv",sep="")

if(file.exists(design_path)==TRUE){
#Read design matrix
design_raw<-read.csv(paste(exp_directory,"/Metadata.csv",sep=""),header=TRUE,fileEncoding = 'UTF-8-BOM')
design_raw$Seq<-gsub(".fastq.gz","",design_raw$Reads)
.GlobalEnv$design_raw<-design_raw

#Get design matrix for significance fitting
.GlobalEnv$design<-data.frame(Intercept=1,
                              GROUP=design_raw$GROUP)

.GlobalEnv$seqs<-design_raw$Seq
.GlobalEnv$seqs_var<-tclVar(seqs)
.GlobalEnv$classes<-unique(design_raw$Classification)
.GlobalEnv$class_num<-length(classes)

.GlobalEnv$save_file<-paste(fastq_dir,"/analysisdata.Rdata",sep="")

if(file.exists(save_file)==FALSE){
  
  .GlobalEnv$thresh<-"auto"
  
  #--------------------FUNCTIONS
  
  #Start prog bar
  .GlobalEnv$an_prog<-winProgressBar(title="DEAR Analysis and Visualization",
                                     label="Loading analysis functions...",
                                     initial=10,min=0,max=100,width=300)
  
  #Read data fun, returns countdata
  .GlobalEnv$read_inputs<-function(){
    
    tryCatch(setWinProgressBar(an_prog,value=20,label="Reading raw feature counts..."),
             error=function(e)print("no prog"))
    
    #Load count data
    countdata_raw<-read.csv(paste(fastq_dir,"/rawfeaturecounts.csv",sep=""),header=TRUE,fileEncoding = 'UTF-8-BOM')
    #Read first col as rownames
    rownames(countdata_raw)<-countdata_raw$X
    countdata_raw<-countdata_raw[,-1,drop=FALSE]
    #GET COUNTDATA
    .GlobalEnv$countdata<-countdata_raw
    
    return(countdata)
  }
  
  #Filter lowly expressed genes
  .GlobalEnv$filter_genes<-function(countdata,thresh){
    
    tryCatch(setWinProgressBar(an_prog,value=30,label="Removing lowly expressed genes..."),
             error=function(e)print("no prog"))
    
    #Get CPM
    .GlobalEnv$cpm_data<-cpm(countdata)
    
    if(thresh=="auto"){
      #Determine thresh automatically, cpm that corresponds to counts of 10
      
      #For each sequence (column, get cpm that corresponds to 10)
      diff_to_10<-function(x){
        x_dat<-countdata[,x]
        vec<-abs(abs(x_dat)-10)
        ind<-which(vec==min(vec))[1]
        cor_cpm<-cpm_data[ind,x]
        return(cor_cpm)
      }
      counts_cor<-unlist(lapply(1:ncol(countdata),diff_to_10))
      
      #Get threshhold as mean
      use_thresh<-as.numeric(round(mean(counts_cor),1))
      
    } else{
      use_thresh<-as.numeric(thresh)
    }
    .GlobalEnv$use_thresh<-use_thresh
    
    #Get well expressed genes
    cpm_data_thresh<-cpm_data>use_thresh
    good_gene_inds<-which(apply(cpm_data_thresh,1,sum)>=1)
    
    print(use_thresh)
    print(nrow(cpm_data))
    print(nrow(countdata[good_gene_inds,]))
    
    return(countdata[good_gene_inds,])
    
  }
  
  #Process data fun
  .GlobalEnv$process_data<-function(countdata_cur,design){
    #Remove features with 0 reads
    .GlobalEnv$zero_count_features<-rownames(countdata_cur[rowSums(countdata_cur==0)>=1,])
    if(length(zero_count_features)>0){
      #Remove genes
      print(paste("removing ",length(zero_count_features)," zero  counts features"),sep="")
      countdata_cur<-countdata_cur[-which(rownames(countdata_cur)%in%zero_count_features),]
      .GlobalEnv$countdata_cur<-countdata_cur
    } 
    
    tryCatch(setWinProgressBar(an_prog,value=40,label="Converting to DGEList..."),
             error=function(e)print("no prog"))
    
    #Convert to a DGE obj
    dgeObj<-DGEList(countdata_cur)
    
    tryCatch(setWinProgressBar(an_prog,value=50,label="Normalizing"),
             error=function(e)print("no prog"))
    
    #Normalize
    dgeObj_norm<-calcNormFactors(dgeObj)
    
    tryCatch(setWinProgressBar(an_prog,value=60,label="Calculating between-group variance..."),
             error=function(e)print("no prog"))
    
    #Get between-group (total dataset) variation
    dgeObj_bwvar<-estimateCommonDisp(dgeObj_norm)
    
    tryCatch(setWinProgressBar(an_prog,value=70,label="Calculating within-group variance..."),
             error=function(e)print("no prog"))
    
    #Get within-group (within gene) variation
    dgeObj_wivar<-estimateGLMTrendedDisp(dgeObj_bwvar)
    dgeObj_tag<-estimateTagwiseDisp(dgeObj_wivar)
    
    tryCatch(setWinProgressBar(an_prog,value=80,label="Fitting linear model..."),
             error=function(e)print("no prog"))
    
    #Fit GLM
    fit<-glmFit(dgeObj_tag,design)
    
    #Conduct lilklihood ratio test for significance
    lrt<-glmLRT(fit)
    
    #Calculate FDR for all genes
    top_genes<-topTags(lrt,n=nrow(lrt$table))
    
    out_list<-list()
    out_list[[length(out_list)+1]]<-dgeObj
    out_list[[length(out_list)+1]]<-as.data.frame(top_genes)
    
    return(out_list)
  }
  
  #--------------------RUN ON OPEN
  
  #Read countdata
  countdata<-read_inputs()
  names(countdata)<-design_raw$Seq
  
  #Get well expressed genes
  .GlobalEnv$countdata_cur<-filter_genes(countdata,thresh)
  
  #Get significance
  .GlobalEnv$sig_list<-process_data(countdata_cur,design)
  .GlobalEnv$raw_dge<-sig_list[[1]]
  .GlobalEnv$annot_genes<-sig_list[[2]]
  
  #Order by feature for volcano plot
  .GlobalEnv$annot_genes_ord<-annot_genes[order(rownames(annot_genes)),]
  
  #Get all features for volcano plot
  .GlobalEnv$annot_fts<-rownames(annot_genes_ord)
  
  tryCatch(setWinProgressBar(an_prog,value=90,label="Formatting data..."),
           error=function(e)print("no prog"))
  
  #Aggregate data for reads/cpm and library distribution plot
  quality_control_list<-lapply(1:ncol(countdata_cur),function(x){
    tmp_counts<-countdata_cur[,x,drop=FALSE]
    seq=names(tmp_counts)
    features=rownames(tmp_counts)
    reads=tmp_counts[,1]
    cpm<-(reads/sum(reads))*1000000
    log2cpm<-log2(cpm)
    class<-design_raw[which(design_raw$Seq==seq),]$Classification
    
    df<-data.frame(Seq=seq,
                   Feature=features,
                   Reads=reads,
                   CPM=cpm,
                   Log2CPM=log2cpm,
                   Class=class)
    return(df)
  })
  quality_control_df<-bind_rows(quality_control_list)
  quality_control_df$ClassSeq<-paste(quality_control_df$Class,quality_control_df$Seq,sep="_")
  quality_control_df$FeatureClass<-paste(quality_control_df$Feature,quality_control_df$Class,sep="_")
  .GlobalEnv$quality_control_df<-quality_control_df
  
  #Get data for library sizes plot
  .GlobalEnv$lib_sizes<-data.frame(Seq=rownames(raw_dge$samples),
                                   Size=raw_dge$samples$lib.size,
                                   Class=design_raw[match(rownames(raw_dge$samples),design_raw$Seq),]$Classification)
  
  #Match data for heatmap plot
  #Feature, Class, Log2CPM, logFC, FDR 
  #Get LogFC and FDR
  .GlobalEnv$all_features<-unique(quality_control_df$Feature)
  
  heat_data<-data.frame(Feature=rep(all_features,class_num),
                        Class=unlist(lapply(classes,function(x){rep(x,times=nrow(annot_genes))})))
  heat_data$FeatureClass<-paste(heat_data$Feature,heat_data$Class,sep="_")
  
  #Get feature, class, Log2CPM
  #Average data over feature class
  na.rm.mean<-function(x){
    return(mean(x,na.rm=TRUE))
  }
  means<-tapply(quality_control_df$CPM,quality_control_df$FeatureClass,na.rm.mean)
  logmeans<-tapply(quality_control_df$Log2CPM,quality_control_df$FeatureClass,na.rm.mean)
  
  #Match to LogFC and FDR data
  heat_data$CPM<-means[match(heat_data$FeatureClass,names(means))]
  heat_data$Log2CPM<-logmeans
  
  annot_matches<-match(heat_data$Feature,rownames(annot_genes_ord))
  
  heat_data$LogFC<-annot_genes_ord[annot_matches,]$logFC
  heat_data$FDR<-annot_genes_ord[annot_matches,]$FDR
  
  #Sort by logFC
  heat_data_ord<-heat_data[order(abs(heat_data$LogFC),decreasing=TRUE),]
  .GlobalEnv$heat_data_ord<-heat_data_ord
  
  tryCatch(setWinProgressBar(an_prog,value=95,label="Saving analysis..."),
           error=function(e)print("no prog"))
  
  #Save dfs needed for analysis in Rdata
  save(list=c("use_thresh",
            "quality_control_df",
             "lib_sizes",
            "heat_data_ord",
            "all_features",
            "annot_genes_ord",
            "annot_fts"),file=save_file)
} else{
  
  #Start prog bar
  .GlobalEnv$an_prog<-winProgressBar(title="DEAR Analysis and Visualization",
                                     label="Loading analysis...",
                                     initial=50,min=0,max=100,width=300)
  
  load(save_file,envir=.GlobalEnv)
  
}



#Get variables for plotting from dfs

#Annotation vars for volcano plot
.GlobalEnv$annot_fts_var<-tclVar(annot_fts)
#All features for heatmap
.GlobalEnv$all_features_var<-tclVar(all_features)

#--------------------PLOT FUNCTIONS

.GlobalEnv$plot_type_var<-tclVar("Reads/CPM")

.GlobalEnv$show_raw_cpm_var<-tclVar("1")
.GlobalEnv$show_raw_cpm_val<-"1"
.GlobalEnv$displated_features<-all_features[1:10]
.GlobalEnv$cur_volc_point<-""

.GlobalEnv$plot_theme<-function(){
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17,face="bold"),
        strip.text = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12,face="bold"))
}

#Save the plot
.GlobalEnv$save_plot_function<-function(){
  setwd(plot_dir)
  ggsave(paste(str_replace_all(cur_plot_type,"/","_"),".png",sep=""),plot=plot,width=15,height=10)
  shell.exec(paste(plot_dir,"/",paste(str_replace_all(cur_plot_type,"/","_"),".png",sep=""),sep=""))
}

#Create the plot
.GlobalEnv$create_plot<-function(){
  
  print("rendering plot")
  
  if(cur_plot_type=="Reads/CPM"){
    
    .GlobalEnv$selected_seq_inds<-tkcurselection(seq_listbox)
    if(""%in%tclvalue(selected_seq_inds)){
      selected_seq_inds<-0
      tkselection.set(seq_listbox,0)
    }
    .GlobalEnv$selected_seqs<-seqs[as.numeric(selected_seq_inds)+1]
    print(selected_seqs)
    
    if(length(selected_seqs)>0){
      tmp_seqs<-quality_control_df[which(quality_control_df$Seq%in%selected_seqs),]
      .GlobalEnv$plot<-ggplot(tmp_seqs,aes(x=CPM,y=Reads))+
        geom_point()+
        xlab("Counts Per Million (CPM)")+
        scale_x_continuous(n.breaks=5,limits=c(0,use_thresh*2))+
        scale_y_continuous(limits=c(0,15),n.breaks=6)+
        geom_vline(xintercept=use_thresh,size=0.3,lty="dashed",col="red")+
        geom_hline(yintercept=10,size=0.3,lty="dashed",col="blue")+
        facet_grid(rows=vars(Seq))+
        theme_classic() 
      
      if(min(tmp_seqs$Reads)>=10){
        .GlobalEnv$plot<-plot+
          geom_text(label="No reads below 10 to remove",x=use_thresh,y=12)
      }
      
    }
    
  } else if(cur_plot_type=="Library Sizes"){
    
    #Generate library size plot
    .GlobalEnv$plot<-ggplot(lib_sizes,aes(x=Seq,y=Size,fill=Class))+
      geom_col(col="black")+
      theme_classic()+
      scale_fill_manual(values=c("gray60","gray40"))+
      ylab("Library Size")+
      xlab("")+
      scale_y_continuous(n.breaks=10)+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
    
  } else if(cur_plot_type=="Library Distribution"){
    
    #Generate library distribution plot
    .GlobalEnv$plot<-ggplot(quality_control_df,aes(x=Seq,y=Log2CPM))+
      geom_boxplot(size=0.1,outlier.size=0.1,outlier.alpha = 1,outlier.color="red",width=0.7)+
      xlab("")+
      theme_classic()+
      theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
    
  } else if(cur_plot_type=="Volcano Plot"){
    
    #Volcano plot
    .GlobalEnv$plot<-ggplot(annot_genes_ord,aes(x=logFC,y=-log10(FDR)))+
      geom_point(col=ifelse(annot_genes_ord$FDR<0.05,"red","black"),size=1)+
      scale_y_continuous(n.breaks=10)+
      scale_x_continuous(n.breaks=10)+
      geom_hline(yintercept=-log10(0.05),col="blue",lty="dashed")+
      xlab("Log(FC)")+
      ylab("-log10(False Discovery Rate (FDR))")+
      theme_classic()
    if(cur_volc_point!=""){
      
      sel_pt_dat<-annot_genes_ord[which(rownames(annot_genes_ord)==cur_volc_point),]
      sel_pt_label<-paste("Feature: ",rownames(sel_pt_dat),"\n",
                          "LogFC: ",round(sel_pt_dat$logFC,3),"\n",
                          "FDR: ",round(sel_pt_dat$FDR,6),sep="")
      
      .GlobalEnv$plot<-plot+
        geom_point(data=sel_pt_dat,col="purple",shape=21,fill="yellow",size=6)+
        geom_label(data=sel_pt_dat,fill="yellow",col="purple",size=6,aes(x=min(annot_genes_ord$logFC),y=0),label=sel_pt_label,vjust=-0.5,hjust=0,label.padding=unit(0.15,"in"))+
        geom_label(data=sel_pt_dat,label=cur_volc_point,col="purple",fill="yellow",hjust=-0.2,vjust=1.2)
      
    }
    
  } else if(cur_plot_type=="Heatmap"){

    #Plot top 10
    .GlobalEnv$tops<-heat_data_ord[which(heat_data_ord$Feature%in%displated_features),]
    
    .GlobalEnv$plot<-ggplot(tops,aes(x=Class,y=factor(Feature,levels=rev(displated_features)),fill=Log2CPM))+
      geom_tile()+
      scale_fill_gradient(low="yellow",high="red",name="Log2(CPM)",n.breaks=5)+
      ylab("Feature")+
      xlab("")+
      theme_classic()
    if(show_raw_cpm_val=="1"){
      .GlobalEnv$plot<-plot+
        geom_text(aes(label=round(Log2CPM,3)))
    }
    
  }
  
  plot_fun<-function(){
    
    .GlobalEnv$plot_to_plot<-plot+plot_theme()
    
    return(plot(plot_to_plot))
  }
  
  #Render the plot
  plot_frame<-tkframe(analyze_gui)
  tkgrid(plot_frame,column=2,row=1,sticky="w",rowspan=1000)
  
  plot_widg<-tkrplot(plot_frame,fun=plot_fun,hscale=2.6,vscale=2.2)
  tkgrid(plot_widg)
  
}

#Update graph parms fun
.GlobalEnv$update_graph_parms<-function(){
  
  .GlobalEnv$cur_plot_type<-tclvalue(plot_type_var)
  
  tryCatch(tkgrid.remove(graph_parm_frame),error=function(e)print("not rendered yet"))
  
  render_title<-function(){
    graph_parms_ttl<-tklabel(graph_parm_frame,text="Graph Parameters",font=underline_font)
    tkgrid(graph_parms_ttl,row=1,column=1)
  }
  
  .GlobalEnv$graph_parm_frame<-tkframe(analyze_gui)
  tkgrid(graph_parm_frame,column=1,row=3)
  
  
  if(cur_plot_type=="Reads/CPM"){
    
    render_title()
    
    .GlobalEnv$deselect_seqs<-function(){
      
      lapply(1:length(seqs),function(x){
        tkselection.clear(seq_listbox,x-1)
      })
  
    }
    
    .GlobalEnv$select_seqs<-function(){
      
      lapply(1:length(seqs),function(x){
        tkselection.set(seq_listbox,x-1)
      })
      create_plot()
      
    }
    
    #Select sequence
    #Scrollbar
    scroll<-tkscrollbar(graph_parm_frame,repeatinterval=1,command=function(...)tkyview(seq_listbox,...))
    tkgrid(scroll,row=2,sticky="nsw",padx=0,column=1)
    
    .GlobalEnv$seq_listbox<-tklistbox(graph_parm_frame,listvariable=seqs_var,width=20,height=6,selectmode="multiple",exportselection=FALSE,yscrollcommand=function(...)tkset(scroll,...))
    tkgrid(seq_listbox,row=2,column=1,padx=20)
    tkbind(seq_listbox,"<<ListboxSelect>>",create_plot)
    
    sel_but<-tkbutton(graph_parm_frame,text="Select All",command=select_seqs)
    tkgrid(sel_but,row=3,column=1,pady=5)
    
    desel_but<-tkbutton(graph_parm_frame,text="Deselect All",command=deselect_seqs)
    tkgrid(desel_but,row=4,column=1,pady=0)
    
    create_plot()
    
  } else if(cur_plot_type=="Library Sizes"){
    
    create_plot()
    
  } else if(cur_plot_type=="Library Distribution"){
    
    create_plot()
    
  } else if(cur_plot_type=="Heatmap"){
    
    render_title()
    
    check_cpms<-function(){
      
      .GlobalEnv$raw_cpm_val<-tclvalue(show_raw_cpm_var)
      if(raw_cpm_val=="0"){
        .GlobalEnv$show_raw_cpm_val<-"1"
      } else if(raw_cpm_val=="1"){
        .GlobalEnv$show_raw_cpm_val<-"0"
      }
      print(show_raw_cpm_val)
      create_plot()
      
    }
    
    get_displayed_fts<-function(){
      
      .GlobalEnv$top_ft<-as.numeric(tknearest(feature_list,1))
      .GlobalEnv$displayed_vec<-c(top_ft:(top_ft+9))+1
      .GlobalEnv$displated_features<-all_features[displayed_vec]
      create_plot()
      
    }
    
    select_a_feature<-function(){
      
      .GlobalEnv$cur_selected<-all_features[(as.numeric(tkcurselection(feature_list))+1)]
      .GlobalEnv$cur_selected_dat<-heat_data_ord[which(heat_data_ord$Feature==cur_selected),]
    
      message=paste("Feature: ",unique(cur_selected_dat$Feature),"\n",
                    "LogFC: ",round(unique(cur_selected_dat$LogFC),3),"\n",
                    "FDR: ",round(unique(cur_selected_dat$FDR),6),"\n",
                    "Group/Log2(CPM): ",paste(cur_selected_dat$Class,round(cur_selected_dat$Log2CPM,3),collapse=" "),sep="")
      tk_messageBox(message=message)
      
    }
    
    #Define scroll bar functions
    scroll_command<-function(...){
      tkset(scroll_fts,...)
      get_displayed_fts()
    }
    also_scroll_command<-function(...){
      tkyview(feature_list,...)
    }
    
    show_cpms<-tkcheckbutton(graph_parm_frame,text='Show Raw CPM',variable=show_raw_cpm_var)
    tkgrid(show_cpms,column=1,pady=5,row=1)
    tkbind(show_cpms,"<Button>",check_cpms)
    
    feature_ttl<-tklabel(graph_parm_frame,text="Features")
    tkgrid(feature_ttl,column=1,sticky="w",row=2)
    
    scroll_fts<-tkscrollbar(graph_parm_frame,repeatinterval=1,command=also_scroll_command)
    tkgrid(scroll_fts,column=1,row=3,sticky="nsw")
    
    feature_list<-tklistbox(graph_parm_frame,listvariable=all_features_var,height=10,width=20,selectmode="single",exportselection=FALSE,yscrollcommand=scroll_command)
    tkgrid(feature_list,column=1,sticky="w",padx=16,row=3)
    tkbind(feature_list,"<<ListboxSelect>>",select_a_feature)
    
    create_plot()
    
  } else if(cur_plot_type=="Volcano Plot"){
    
    select_volc_point<-function(){
      
      .GlobalEnv$cur_volc_point<-annot_fts[(as.numeric(tkcurselection(volc_feature_list))+1)]
      
      print(cur_volc_point)
      
      create_plot()
      
    }
    
    render_title()
    
    feature_ttl<-tklabel(graph_parm_frame,text="Features")
    tkgrid(feature_ttl,column=1,sticky="w",row=2)
    
    volcscroll<-tkscrollbar(graph_parm_frame,repeatinterval=1,command=function(...)tkyview(volc_feature_list,...))
    tkgrid(volcscroll,column=1,row=3,sticky="nsw")
    
    volc_feature_list<-tklistbox(graph_parm_frame,listvariable=annot_fts_var,height=10,width=20,selectmode="single",exportselection=FALSE,yscrollcommand=function(...)tkset(volcscroll,...))
    tkgrid(volc_feature_list,row=3,column=1,padx=16)
    tkbind(volc_feature_list,"<<ListboxSelect>>",select_volc_point)
    
    create_plot()
    
  }
  
}

tryCatch(close(an_prog),error=function(e)print("no prog"))

#--------------------RENDER GUI

analyze_gui<-tktoplevel()
tkwm.geometry(analyze_gui,"900x600+200+20")
tkwm.title(analyze_gui,"DEAR Plot Builder")

#Top frame
first_frame<-tkframe(analyze_gui)
tkgrid(first_frame,pady=10,padx=10,column=1,row=1)
  
  title_lbl<-tklabel(first_frame,text="DEAR Plot Builder",font=title_font)
  tkgrid(title_lbl)
  
#Plot lbl
plot_frame<-tkframe(analyze_gui)
tkgrid(plot_frame,column=2,row=2,sticky="n",rowspan=1000,padx=20,columnspan=100)

  plot_display<-tklabel(plot_frame,text="")
  tkgrid(plot_display)

#Type frame
type_frame<-tkframe(analyze_gui,borderwidth=3,relief="raised")
tkgrid(type_frame,pady=10,padx=10,column=1,row=2)

  #Select plot type
  type_lbl<-tklabel(type_frame,text="Select plot",justify="left",font=header_font)
  tkgrid(type_lbl,padx=35,row=1)
  
  type_sel<-ttkcombobox(type_frame,values=c("Reads/CPM","Library Sizes","Library Distribution","Heatmap","Volcano Plot"),textvariable=plot_type_var,width=15)
  tkgrid(type_sel,padx=35,pady=10,row=2)
  tkbind(type_sel,"<<ComboboxSelected>>",update_graph_parms)

#Generate plot
gen_frame<-tkframe(analyze_gui)
tkgrid(gen_frame,column=1,row=4,pady=10)

  #Save plot
  save_but<-tkbutton(gen_frame,text="Save plot",font=header_font,command=save_plot_function)
  tkgrid(save_but,pady=5,padx=15,row=1,column=2)

update_graph_parms()
  
tkwait.window(analyze_gui)
} else{
  
  tk_messageBox(message=paste("Design matrix file not found in ",exp_directory,".\n\nEnsure you have annotated your reads!",sep=""))
  
}


