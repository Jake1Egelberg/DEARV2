
if(file.exists(paste(fastq_dir,"/rawfeaturecounts.csv",sep=""))==TRUE){
  
  scrape_parms()
  
  #*****************************************************
  
  #******************CALCULATIONS AND PROCESSING*****************
  
  #******************************************************
  
  calc_prog<-winProgressBar(title="DEAR Plot Builder",
                            label="Calculating metrics...",
                            min=0,max=100,width=300,initial=50)
  
  #Load design matrix
  design_raw<-read.csv(paste(exp_directory,"/Metadata.csv",sep=""),header=TRUE,fileEncoding = 'UTF-8-BOM')
  design_raw$Seq<-gsub(".fastq.gz","",design_raw$Reads)
  
  #Name sequences
  .GlobalEnv$raw_classes<-levels(as.factor(design_raw$Classification))
  .GlobalEnv$class_seq<-paste(design_raw[which(design_raw$Classification==raw_classes[1]),]$Classification,1:length(design_raw[which(design_raw$Classification==raw_classes[1]),]$Classification),sep=" ")
  .GlobalEnv$class_seq<-c(class_seq,paste(design_raw[which(design_raw$Classification==raw_classes[2]),]$Classification,1:length(design_raw[which(design_raw$Classification==raw_classes[2]),]$Classification),sep=" "))
  .GlobalEnv$classes<-design_raw$Classification
  
  #Load count data
  countdata_raw<-read.csv(paste(fastq_dir,"/rawfeaturecounts.csv",sep=""),header=TRUE,fileEncoding = 'UTF-8-BOM')
  rownames(countdata_raw)<-countdata_raw$X
  countdata_raw<-countdata_raw[,-1,drop=FALSE]
  
  #GET COUNTDATA
  .GlobalEnv$countdata<-countdata_raw
  #Get sequence names
  seqs<-names(countdata_raw)
  
  #Remove features with 0 reads
  .GlobalEnv$zero_count_features<-rownames(countdata[rowSums(countdata==0)>=1,])
  if(length(zero_count_features)>0){
    #Remove genes
    .GlobalEnv$countdata<-countdata[-which(rownames(countdata)%in%zero_count_features),]
  }
  
  #FORMAT FOR GGPLOT
  mycpm<-lapply(1:length(seqs),function(x){
    tmp_counts<-countdata[,x,drop=FALSE]
    df<-data.frame(Seq=names(tmp_counts),Feature=rownames(tmp_counts),Reads=tmp_counts[,1])
    df$CPM<-(df$Reads/sum(df$Reads))*1000000
    df$Log2CPM<-log2(df$CPM)
    df$Class<-classes[x]
    df$ClassSeq<-class_seq[x]
    out<-df
  })
  .GlobalEnv$mycpm_cur<-bind_rows(mycpm)
  
  #READS/CPM (get data for <30 reads)
  .GlobalEnv$reads_plot_data<-mycpm_cur[which(mycpm_cur$Reads<quantile(mycpm_cur$Reads,0.1)[[1]]),]
  
  #HEATMAP (get average log2(cpm) within each group, calculate var and FC b/w groups)
  class_avg<-lapply(classes,function(x){
    tmp<-mycpm_cur[which(mycpm_cur$Class==x),]
    avgs<-tapply(tmp$Log2CPM,tmp$Feature,mean)
    out<-data.frame(Feature=names(avgs),Avg=avgs,Class=x)
  })
  class_avg_cur<-bind_rows(class_avg)
  rownames(class_avg_cur)<-NULL
  #Get var and fc(Mean(Log2CPM)) b/w groups
  fold_change_func<-function(x){
    y<-unique(x)
    out<-y[1]/y[2]
  }
  vars<-tapply(class_avg_cur$Avg,class_avg_cur$Feature,var)
  fcs<-tapply(class_avg_cur$Avg,class_avg_cur$Feature,fold_change_func)
  #Create df matching feature with var(mean(log2(cpm))) and fc
  vars_df<-data.frame(Feature=names(vars),
                      Var=vars,
                      FC=fcs)
  .GlobalEnv$vars_df<-vars_df[order(vars,decreasing=TRUE),]
  
  #Join avg of groups to mycpm_cur
  class_avg_cur$Class_Ft<-paste(class_avg_cur$Class,"_",class_avg_cur$Feature,sep="")
  mycpm_cur$Class_Ft<-paste(mycpm_cur$Class,"_",mycpm_cur$Feature,sep="")
  class_avg_cur<-class_avg_cur[,-c(1,3)]
  class_avg_cur<-unique(class_avg_cur)
  .GlobalEnv$mycpm_cur<-left_join(mycpm_cur,class_avg_cur,by="Class_Ft")
  
  #VOLCANO PLOT
  calc_tests<-lapply(unique(mycpm_cur$Feature),function(x){
    y<-mycpm_cur[which(mycpm_cur$Feature==x),]
    #T-test
    test<-t.test(y$Log2CPM~y$Class,alternative="two.sided",var.equal=FALSE,paired=FALSE)
    #Get one-sided raw p value
    p_val<-(test$p.value)/2
    df<-data.frame(Feature=x,P_Raw=p_val)
  })
  raw_sig_fts<-bind_rows(calc_tests)
  
  #Calculate false discovery rate
  raw_sig_fts$FDR<-p.adjust(raw_sig_fts$P_Raw,method="fdr")
  
  #Join to data
  .GlobalEnv$vars_df<-left_join(vars_df,raw_sig_fts,by="Feature")
  rownames(vars_df)<-vars_df$Feature
  
  #Calc metrics
  .GlobalEnv$vars_df$Log2AbsFC<-log2(abs(vars_df$FC))
  .GlobalEnv$vars_df$LogFDR<-log(vars_df$FDR,base=10)*-1
  vars_df$Sig<-NA
  if(length(which(vars_df$FDR<0.05))>0){
    .GlobalEnv$vars_df[which(vars_df$FDR<0.05),]$Sig<-"Y" 
  }
  .GlobalEnv$vars_df[which(vars_df$FDR>=0.05),]$Sig<-"N"
  .GlobalEnv$vars_df<-vars_df[order(vars_df$LogFDR,decreasing=TRUE),]
  
  #Get sig fts
  sig_fts<-vars_df[which(vars_df$Sig=="Y"),]

  #plot(vars_df$Log2AbsFC,vars_df$LogFDR,col=ifelse(vars_df$Sig=="Y","red","black"))
  #abline(h=-log(0.05),col="red")
  
  #Close progress bar
  close(calc_prog)
  
  #*****************************************************
  
  #******************GRAPH GENERATION*****************
  
  #******************************************************
  
  #Define variables
  #PLOT THEME
  .GlobalEnv$plot_theme<-function(){
    theme(axis.text.x = element_text(size=3,angle=45,hjust=1.1),
          axis.text.y = element_text(size=3),
          axis.title.x = element_text(size=4),
          axis.title.y = element_text(size=4),
          strip.text.y = element_text(size=3),
          panel.grid.major = element_line(size=0.3),
          panel.grid.minor=element_blank(),
          panel.border = element_rect(size=0.3),
          legend.text = element_text(size=3),
          legend.title = element_text(size=4),
          legend.key.width = unit(0.3,"cm"),
          legend.key.height=unit(0.3,"cm"),
          legend.margin = margin(t=0,unit="cm"))
  }
  #Number of parameters for parm frame
  parm_num<-1:10
  plot_type_var<-tclVar("Reads/CPM")
  #Reads/CPM
  scale_thresh<-tclVar(thresh_val)
  .GlobalEnv$saved_thresh<-thresh_val
  .GlobalEnv$thresh_val_start<-thresh_val
  new_thresh_val<-thresh_val
  max_thresh<-ceiling(quantile(mycpm_cur$CPM,0.1)[[1]])
  scale_resolution<-max_thresh/25
  max_x<-ceiling(quantile(mycpm_cur$CPM,0.1)[[1]])
  x_res<-floor(max_x/25)
  scale_x<-tclVar(max_x/2)
  max_y<-quantile(mycpm_cur$Reads,0.1)[[1]]
  y_res<-floor(max_y/25)
  scale_y<-tclVar(max_y/2)
  #Selected seqs for display in listbox
  all_seqs_raw<-levels(as.factor(mycpm_cur$Seq))
  all_seqs<-tclVar(levels(as.factor(mycpm_cur$Seq)))
  .GlobalEnv$selected_seqs<-levels(as.factor(mycpm_cur$Seq))
  lens<-0:(length(levels(as.factor(mycpm_cur$Seq)))-1)
  #If apply thresh
  switch_seq<-tclVar(1)
  #Show outliers boxplot
  show_outs<-tclVar(1)
  #Starting index for heatmap display
  .GlobalEnv$top_ft<-0
  .GlobalEnv$text_size<-0.9
  heat_show_txt_var<-tclVar(1)
  heat_show_var<-tclVar(1)
  .GlobalEnv$lim_max<-max(mycpm_cur$Avg)*1.05
  .GlobalEnv$lim_min<-min(mycpm_cur$Avg)/1.05
  #Volcano plit
  volc_app_thresh<-tclVar(1)
  .GlobalEnv$selected_sig_ft<-NA
  #Apply threshold across plots 
  .GlobalEnv$apply_threshold_var<-tclVar(0)
  
  #Generate parms frame
  generate_parms_frame<-function(){
    #Plot parameters frame
    .GlobalEnv$plot_parms_frame<-tkframe(analyze_gui,borderwidth=3,relief="raised")
    tkgrid(plot_parms_frame,pady=10,padx=10,column=1,row=3)
    parms_head<-tklabel(plot_parms_frame,text="Plot parameters",font=header_font)
    tkgrid(parms_head,padx=35,pady=5,column=1,columnspan=2)
  }
  
  #Save plot function
  save_plot_function<-function(){
    setwd(plot_dir)
    ggsave(paste(str_replace_all(plot_type_val,"/","_"),".png",sep=""),plot=plot,width=600,height=500,units = "px")
    shell.exec(paste(plot_dir,"/",paste(str_replace_all(plot_type_val,"/","_"),".png",sep=""),sep=""))
  }
  
  #Reset graph parms function
  reset_graph_parms<-function(){
    #Reset parms frame
    tkgrid.remove(plot_parms_frame)
    generate_parms_frame()
    
    #Remove display
    tkconfigure(plot_display,image="")
    
  }
  
  #Update graph parms function
  update_graph_parms<-function(){
    
    #Get selected plot type
    .GlobalEnv$plot_type_val<-tclvalue(plot_type_var)
    print(plot_type_val)
    
    select_listbox<-function(){
      selected_nums<-tkcurselection(sel_seqs)
      .GlobalEnv$selected_seqs<-all_seqs_raw[as.numeric(selected_nums)+1]
      if(length(selected_seqs)==0){
        .GlobalEnv$selected_seqs<-all_seqs_raw
      }
      print(selected_seqs)
      build_plot_fun()
    }
    
    select_all_func<-function(){
      lapply(lens,function(x){
        tkselection.set(sel_seqs,x)
      })
      select_listbox()
    }
    
    deselect_all_func<-function(){
      lapply(lens,function(x){
        tkselection.clear(sel_seqs,x)
      })
      select_listbox()
    }
    
    generate_listbox<-function(col_num,row_num){
      #Select sequences lbl
      .GlobalEnv$seq_entry<-tklabel(plot_parms_frame,text="Select sequences",font=underline_font)
      tkgrid(seq_entry,row=row_num,column=col_num,columnspan=2)
      #Scrollbar
      .GlobalEnv$scroll<-tkscrollbar(plot_parms_frame,repeatinterval=1,command=function(...)tkyview(seq_list,...))
      tkgrid(scroll,column=col_num,row=row_num+1,sticky="nsw",padx=10)
      #Listbox
      .GlobalEnv$sel_seqs<-tklistbox(plot_parms_frame,listvariable=all_seqs,width=20,height=6,selectmode="multiple",exportselection=FALSE,yscrollcommand=function(...)tkset(scroll,...))
      tkgrid(sel_seqs,column=col_num,row=row_num+1,columnspan=2,pady=10)
      tkbind(sel_seqs,"<<ListboxSelect>>",select_listbox)
      #Select all
      .GlobalEnv$seq_all<-tkbutton(plot_parms_frame,text="Select all",command=select_all_func)
      tkgrid(seq_all,row=row_num+2,column=col_num,pady=3)
      #Deselectall
      .GlobalEnv$desel_all<-tkbutton(plot_parms_frame,text="Deselect all",command=deselect_all_func)
      tkgrid(desel_all,row=row_num+2,column=col_num+1,pady=3)
    }
    
    #If selected reads/cpm plot
    if(plot_type_val=="Reads/CPM"){
      
      #Reset graph parms
      reset_graph_parms()
      
      #Build graph
      build_plot_fun()
      .GlobalEnv$thresh_val_start<-tclvalue(scale_thresh)
      
      adjust_plot<-function(...){
        build_plot_fun()
      }
      
      save_threshhold<-function(){
        .GlobalEnv$saved_thresh<-tclvalue(scale_thresh)
        tk_messageBox(message=paste("You have selected ",saved_thresh,sep=""))
      }
      
      #Set x axis limits
      parm1<-tklabel(plot_parms_frame,text="X-axis max:",font=normal_font)
      tkgrid(parm1,pady=5,sticky="w",padx=10,column=1)
      x_axis_entry<-tkscale(plot_parms_frame,width=15,orient="horizontal",from=0,to=max_x,resolution=x_res,length=80,width=10,sliderlength=20,variable=scale_x,command=adjust_plot)
      tkgrid(x_axis_entry,row=1,column=2,sticky="w",padx=3)
      tkset(x_axis_entry,tclvalue(scale_x))
      #Set y axis limits
      parm2<-tklabel(plot_parms_frame,text="Y-axis max:",font=normal_font)
      tkgrid(parm2,plot_parms_frame,pady=5,sticky="w",padx=10,column=1)
      y_axis_entry<-tkscale(plot_parms_frame,width=15,orient="horizontal",from=0,to=max_y,resolution=y_res,length=80,width=10,sliderlength=20,variable=scale_y,command=adjust_plot)
      tkgrid(y_axis_entry,row=2,column=2,sticky="w",padx=3)
      tkset(y_axis_entry,tclvalue(scale_y))
      #Configure parm to set thresh
      parm3<-tklabel(plot_parms_frame,text="Threshold",font=normal_font)
      tkgrid(parm3,pady=5,sticky="w",padx=10,column=1)
      thresh_entry<-tkscale(plot_parms_frame,width=15,orient="horizontal",from=0,to=max_thresh,resolution=scale_resolution,length=80,width=10,sliderlength=20,variable=scale_thresh,command=adjust_plot)
      tkgrid(thresh_entry,row=3,column=2,sticky="w",padx=3)
      tkset(thresh_entry,thresh_val_start)
      #Listbox
      generate_listbox(1,4)
      #Save threshold
      save_thresh<-tkbutton(plot_parms_frame,width=15,text="Save threshold",command=save_threshhold,font=small_font)
      tkgrid(save_thresh,row=7,column=1,columnspan=2,pady=5)
      
    } else if(plot_type_val=="Library Sizes"){
      
      #Reset graph parms
      reset_graph_parms()
      
      #Build graph
      build_plot_fun()
      
      #Generate listbox and scroll
      generate_listbox(1,1)
      #Apply threshhold
      switch_seq<-tkcheckbutton(plot_parms_frame,text="Apply thresh",variable=apply_threshold_var,command=build_plot_fun)
      tkgrid(switch_seq,row=4,column=1,columnspan=2,pady=10)
      
    } else if(plot_type_val=="Library Distribution"){
      
      #Reset graph parms
      reset_graph_parms()
      
      #Build graph
      build_plot_fun()
      
      #Generate listbox and scroll
      generate_listbox(1,1)
      
      #Apply threshhold
      switch_seq<-tkcheckbutton(plot_parms_frame,text="Apply thresh",variable=apply_threshold_var,command=build_plot_fun)
      tkgrid(switch_seq,row=4,column=1,columnspan=2,pady=5)
      #Show outliers
      show_out<-tkcheckbutton(plot_parms_frame,text="Show outliers",variable=show_outs,command=build_plot_fun)
      tkgrid(show_out,row=5,column=1,columnspan=2,pady=5)
      
    } else if(plot_type_val=="Heatmap"){
      
        #Select all sequences
        .GlobalEnv$selected_seqs<-levels(as.factor(mycpm_cur$Seq))
      
        #Reset graph parms
        reset_graph_parms()
      
        #Build graph
        build_plot_fun()
        
        feature_sel_fun<-function(){
          build_plot_fun()
          
          tmp_sel<-tkcurselection(ft_list)
          .GlobalEnv$tmp_ft_sel<-vars_df$Feature[as.numeric(tmp_sel)+1]
          if(length(tmp_ft_sel)>0){
            tmp_df<-mycpm_cur[which(mycpm_cur$Feature==tmp_ft_sel),]
            tmp_df_cur<-data.frame(Sequence=tmp_df$Seq,
                                   Feature=tmp_df$Feature,
                                   Counts=tmp_df$Reads,
                                   CPM=round(tmp_df$CPM,digits=2),
                                   Log2CPM=round(tmp_df$Log2CPM,digits=2),
                                   Group=tmp_df$Class,
                                   MeanLog2CPM=round(tmp_df$Avg,digits=2),
                                   FC=round(unique(tmp_df$Avg)[1]/unique(tmp_df$Avg)[2],digits=2))
            fix(tmp_df_cur) 
          }
        }
        
        #Define scroll bar
        scroll_command<-function(...){
          tkset(scroll_bar,...)
          .GlobalEnv$top_ft<-as.numeric(tknearest(ft_list,1))
          build_plot_fun()
        }
        
        #Define order of features in listbox
        vars_ordered_tcl<-tclVar(paste(vars_df$Feature,", FC= ",round(vars_df$FC,2),sep=""))
        
        heat_show_txt<-tkcheckbutton(plot_parms_frame,text="Apply threshold?",variable=apply_threshold_var,command=build_plot_fun)
        tkgrid(heat_show_txt,column=1,row=1,sticky="w",padx=10)
        heat_txt<-tkcheckbutton(plot_parms_frame,text="Show Mean(Log2(CPM))?",variable=heat_show_var,command=build_plot_fun)
        tkgrid(heat_txt,column=1,row=2,sticky="w",padx=10)
        parm2<-tklabel(plot_parms_frame,text="Most variable",font=normal_font)
        tkgrid(parm2,pady=5,sticky="w",padx=10,column=1)
        scroll_bar<-tkscrollbar(plot_parms_frame,repeatinterval=1,command=function(...)tkyview(ft_list,...))
        tkgrid(scroll_bar,column=1,row=3,sticky="nsw")
        ft_list<-tklistbox(plot_parms_frame,listvariable=vars_ordered_tcl,height=15,width=18,yscrollcommand=scroll_command)
        tkgrid(ft_list,row=3,column=1,sticky="w",padx=20)
        tkbind(ft_list,"<<ListboxSelect>>",feature_sel_fun)
        
    } else if(plot_type_val=="Volcano Plot"){
      sig_ft_var<-tclVar(rownames(vars_df))
      .GlobalEnv$use_curs<-FALSE
      
      #Reset graph parms
      reset_graph_parms()
      
      sel_sig_ft<-function(){
        if(use_curs==FALSE){
          .GlobalEnv$selected_sig_ft<-rownames(vars_df)[as.numeric(tkcurselection(volc_sigs))+1]
        } else{
          .GlobalEnv$selected_sig_ft<-rownames(vars_df_cur)[as.numeric(tkcurselection(volc_sigs))+1]
        }
        print(selected_sig_ft)
        tkconfigure(det1_det,text=round(vars_df[which(vars_df$Feature==selected_sig_ft),]$Log2AbsFC,digits=3))
        tkconfigure(det2_det,text=round(vars_df[which(rownames(vars_df)==selected_sig_ft),]$LogFDR,digits=3))
        build_plot_fun()
      }
      
      volcano_apply_thresh<-tkcheckbutton(plot_parms_frame,text="Apply threshold?",variable=apply_threshold_var,command=build_plot_fun)
      tkgrid(volcano_apply_thresh,column=1,row=1,sticky="w",padx=10)
      volc_scroll<-tkscrollbar(plot_parms_frame,repeatinterval=1,command=function(...)tkyview(volc_sigs,...))
      tkgrid(volc_scroll,column=1,row=2,sticky="nsw")
      .GlobalEnv$volc_sigs<-tklistbox(plot_parms_frame,listvariable=sig_ft_var,height=15,width=18,yscrollcommand=function(...)tkset(volc_scroll,...))
      tkgrid(volc_sigs,row=2,column=1,sticky="w",padx=20,pady=5)
      tkbind(volc_sigs,"<<ListboxSelect>>",sel_sig_ft)
      display_info<-tklabel(plot_parms_frame,text="Feature details",font=normal_font)
      tkgrid(display_info,row=3,column=1,columnspan=2)
      det1<-tklabel(plot_parms_frame,text="Fold change =",font=normal_font)
      tkgrid(det1,row=4,column=1,pady=2,sticky="w",padx=5)
      det1_det<-tklabel(plot_parms_frame,text="",font=normal_font)
      tkgrid(det1_det,row=4,column=1,pady=2,sticky="e",padx=40)
      det2<-tklabel(plot_parms_frame,text="FDR = ",font=normal_font)
      tkgrid(det2,row=5,column=1,pady=2,sticky="w",padx=50)
      det2_det<-tklabel(plot_parms_frame,text="",font=normal_font)
      tkgrid(det2_det,row=5,column=1,pady=2,sticky="e",padx=40)
      
      #Build graph
      build_plot_fun()
      
    }
    
    
    #End different parm layouts
  }
  
  
  #-----------------------BUILD PLOT----------------------
  
  #Build plot function
  build_plot_fun<-function(){
    #Get selected plot type
    .GlobalEnv$plot_type_val<-tclvalue(plot_type_var)
    print(plot_type_val)
    
    #Get if apply thresh
    .GlobalEnv$apply_threshold_val<-tclvalue(apply_threshold_var)
    print(apply_threshold_val)
    
    #Graphing...
    if(plot_type_val=="Reads/CPM"){
      
      #Select seqs
      reads_plot_data<-reads_plot_data[which(reads_plot_data$Seq%in%selected_seqs),]
      
      #Retrieve new thresh entries
      new_thresh_val<-as.numeric(tclvalue(scale_thresh))
      .GlobalEnv$new_x_axis<-as.numeric(tclvalue(scale_x))
      .GlobalEnv$new_y_axis<-as.numeric(tclvalue(scale_y))
      
      #Graph
      .GlobalEnv$plot<-ggplot(reads_plot_data,aes(x=CPM,y=Reads))+
        geom_point(size=0.3)+
        facet_grid(rows=vars(Seq))+
        xlab("Counts Per Million (CPM)")+
        ylab("Raw Counts")+
        scale_x_continuous(n.breaks=5,limits=c(0,new_x_axis))+
        scale_y_continuous(limits=c(0,new_y_axis),n.breaks=6)+
        geom_vline(xintercept=new_thresh_val,size=0.3,lty="dashed",col="red")+
        geom_hline(yintercept=10,size=0.3,lty="dashed",col="blue")+
        theme_bw()+
        plot_theme()+
        theme(axis.text.x = element_text(angle=0,hjust=0))
      
    } else if(plot_type_val=="Library Sizes"){
      
      #Preprocessing
      plot_data<-mycpm_cur
      if(apply_threshold_val=="1"){
        plot_data<-mycpm_cur[which(mycpm_cur$CPM>as.numeric(saved_thresh)),] 
      }
      sum<-as.data.frame(tapply(plot_data$Reads,plot_data$Seq,sum))
      names(sum)<-"Sum"
      sum$Seq<-rownames(sum)
      sum$Class<-classes
      #Select seqs
      sum<-sum[which(sum$Seq%in%selected_seqs),]
      .GlobalEnv$sum<-sum
      
      #Graph
      .GlobalEnv$plot<-ggplot(sum,aes(x=Seq,y=Sum,fill=Class))+
        geom_col(width=0.7,col="black",size=0.1)+
        xlab("Sequence Accession")+
        ylab("Total Counts")+
        scale_y_continuous(n.breaks=10)+
        scale_fill_grey(start=0.5,end=0.8,name="Group")+
        theme_bw()+
        plot_theme()
      
    } else if(plot_type_val=="Library Distribution"){
      
      #Preprocessing
      plot_data<-mycpm_cur
      if(apply_threshold_val=="1"){
        plot_data<-mycpm_cur[which(mycpm_cur$CPM>as.numeric(saved_thresh)),] 
      }
      plot_data<-plot_data[which(plot_data$Seq%in%selected_seqs),]
      .GlobalEnv$show_outs_val<-tclvalue(show_outs)
      if(show_outs_val=="1"){
        out_alpha=1
      } else if(show_outs_val=="0"){
        out_alpha=0
      }
      
      #Graph
      .GlobalEnv$plot<-ggplot(plot_data,aes(x=Seq,y=Log2CPM,fill=Class))+
        geom_boxplot(size=0.1,outlier.size=0.1,outlier.alpha = out_alpha,outlier.color="red",width=0.7)+
        scale_y_continuous(n.breaks=10,limits=c(0,NA))+
        theme_bw()+
        scale_fill_grey(start=0.5,end=0.8,name="Group")+
        ylab("Log2(CPM)")+
        xlab("Sequence Accession")+
        plot_theme()
      
    } else if(plot_type_val=="Heatmap"){
      
      #Get highly variable genes  
      var_subset<-vars_df$Feature[(top_ft+1):(top_ft+15)]
      plot_data<-mycpm_cur[which(mycpm_cur$Feature%in%var_subset),]
      #Manage labels
      heat_show_val<-tclvalue(heat_show_var)
      if(heat_show_val=="1"){
        plot_data$Label<-round(plot_data$Avg,2)
      } else{
        plot_data$Label<-""
      }
      #Manage thresh application
      plot_data$alpha_val<-1
      if(apply_threshold_val=="1"){
        #Get genes below thresh
        below<-unique(mycpm_cur[which(mycpm_cur$CPM<saved_thresh),]$Feature)
        if(length(intersect(below,plot_data$Feature))>0){
          plot_data[which(plot_data$Feature%in%below),]$alpha_val<-0
          plot_data[which(plot_data$Feature%in%below),]$Label<-""
        }
      }
      
      #Graph
      .GlobalEnv$plot<-ggplot(plot_data,aes(x=Class,y=factor(Feature,levels=rev(var_subset)),fill=Avg))+
        geom_tile(alpha=plot_data$alpha_val)+
        geom_text(aes(label=Label),size=text_size)+
        scale_fill_gradient(low="yellow",high="red",name="Mean Log2(CPM)",limits=c(lim_min,lim_max),n.breaks=10)+
        ylab("Feature")+
        xlab("Group")+
        theme_bw()+
        plot_theme()
      
      
    } else if(plot_type_val=="Volcano Plot"){
      
      if(apply_threshold_val=="1"){
        #Get genes below thresh
        .GlobalEnv$below<-unique(mycpm_cur[which(mycpm_cur$CPM<saved_thresh),]$Feature)
        
        #Reset selected gene if below thresh
        if(is.na(selected_sig_ft)==FALSE&&selected_sig_ft%in%below){
          selected_sig_ft<-NA
        }
        
        if(length(below)>0){
          #Remove genes below thresh
          vars_df_cur<-vars_df[-which(rownames(vars_df)%in%below),] 
          tkconfigure(volc_sigs,listvariable=tclVar(rownames(vars_df_cur)))
          .GlobalEnv$use_curs<-TRUE
          
        } else{
          vars_df_cur<-vars_df
        }
      } else{
        vars_df_cur<-vars_df
        tkconfigure(volc_sigs,listvariable=tclVar(rownames(vars_df_cur)))
      }
      
      #Set default display options
      vars_df_cur$Col="black"
      if(length(which(vars_df_cur$Sig=="Y"))>0){
        vars_df_cur[which(vars_df_cur$Sig=="Y"),]$Col<-"red" 
      }
      vars_df_cur$Size<-0.2
      vars_df_cur$Text<-0
      
      #Get selected sig ft
      if(length(selected_sig_ft)>0){
        if(is.na(selected_sig_ft)==FALSE){
          vars_df_cur[which(rownames(vars_df_cur)==selected_sig_ft),]$Col<-"darkgreen"
          vars_df_cur[which(rownames(vars_df_cur)==selected_sig_ft),]$Text<-1.5
          .GlobalEnv$sel_sig_dat<-vars_df_cur[which(rownames(vars_df_cur)==selected_sig_ft),]
        }  else{
          selected_sig_ft<-NA
        }
      } else{
        selected_sig_ft<-NA
      }
      

      .GlobalEnv$vars_df_cur<-vars_df_cur
      
      .GlobalEnv$plot<-ggplot(data=vars_df_cur,aes(x=Log2AbsFC,y=LogFDR))+
        geom_point(size=vars_df_cur$Size,col=vars_df_cur$Col,alpha=1)+
        geom_hline(yintercept = -log(0.05,base=10),size=0.1,col="red")+
        xlab("Log2(CPM) Fold Change")+
        ylab("-log10(FDR)")+
        scale_y_continuous(n.breaks=10,limits=c(NA,1.1*max(vars_df_cur$LogFDR)))+
        theme_bw()+
        plot_theme()+
        theme(axis.text.x = element_text(angle=0,hjust=0))
      if(is.na(selected_sig_ft)==FALSE){
        .GlobalEnv$plot<-plot+geom_text(inherit.aes=FALSE,data=sel_sig_dat,aes(x=Log2AbsFC,y=LogFDR),label=selected_sig_ft,size=sel_sig_dat$Text,col=sel_sig_dat$Col,vjust=-1,fontface="bold")
      }
      
      
    }
    
    setwd(directory)
    ggsave("tmp.png",plot=plot,width=600,height=500,units = "px")
    tst2<-tkimage.create("photo","test",
                         file=paste(directory,"/tmp.png",sep=""))
    tkconfigure(plot_display,image=tst2)
    file.remove(paste(directory,"/tmp.png",sep=""))
    
  }
  
  #Gui
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
  
  #Generate parms frame
  generate_parms_frame()
  
  #Generate plot
  gen_frame<-tkframe(analyze_gui)
  tkgrid(gen_frame,column=1,row=4)
  gen_but<-tkbutton(gen_frame,text="Build plot",font=header_font,command=build_plot_fun)
  tkgrid(gen_but,pady=5,padx=15,sticky="e",row=1,column=1)
  #Save plot
  save_but<-tkbutton(gen_frame,text="Save plot",font=header_font,command=save_plot_function)
  tkgrid(save_but,pady=5,padx=15,sticky="w",row=1,column=2)
  #Update graph parms
  update_graph_parms()
  tkwait.window(analyze_gui)
  
  
  
} else{
  tk_messageBox(message="Feature counts not detected!")
}

