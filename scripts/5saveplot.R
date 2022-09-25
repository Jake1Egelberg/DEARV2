
#Get selected plot type
.GlobalEnv$plot_type_val<-tclvalue(plot_type_var)
print(plot_type_val)

#Get if apply thresh
.GlobalEnv$switch_seq_val<-tclvalue(switch_seq)
print(switch_seq_val)

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
  
  #Select seqs
  sum<-sum[which(sum$Seq%in%selected_seqs),]
  
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
  if(switch_seq_val=="1"){
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
  
  #Get genes below thresh
  below<-unique(mycpm_cur[which(mycpm_cur$CPM<saved_thresh),]$Feature)
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
  heat_show_txt_val<-tclvalue(heat_show_txt_var)
  plot_data$alpha_val<-1
  if(heat_show_txt_val=="1"){
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
  
  .GlobalEnv$volc_app_thresh_val<-tclvalue(volc_app_thresh)
  
  if(volc_app_thresh_val=="1"){
    #Get genes below thresh
    .GlobalEnv$below<-unique(mycpm_cur[which(mycpm_cur$CPM<saved_thresh),]$Feature)
    if(length(below)>0){
      #Remove genes below thresh
      tops_cur<-tops[-which(rownames(tops)%in%below),] 
      if(length(which(rownames(sig_fts)%in%below))>0){
        sig_fts_cur<-sig_fts[-which(rownames(sig_fts)%in%below),]
        new_names<-rownames(sig_fts)
        new_names[which(rownames(sig_fts)%in%below)]<-paste(new_names[which(rownames(sig_fts)%in%below)]," BELOW",sep="")
        tkconfigure(volc_sigs,listvariable=tclVar(new_names))
      }
      
    } else{
      tops_cur<-tops
    }
  } else{
    tops_cur<-tops
    tkconfigure(volc_sigs,listvariable=tclVar(rownames(tops)))
  }
  
  tops_cur$Size<-0.2
  #Get selected sig ft
  if(is.na(selected_sig_ft)==FALSE){
    tops_cur[which(rownames(tops_cur)==selected_sig_ft),]$Size<-1
  } 
  
  #Ensure that -log10(FDR=0) is not Inf
  if(length(which(-log10(tops_cur$FDR)==Inf))>0){
    tops_cur[which(-log10(tops_cur$FDR)==Inf),]$FDR<-0.01
  }
  
  .GlobalEnv$plot<-ggplot(data=tops_cur,aes(x=logFC,y=-log10(FDR)))+
    geom_point(size=tops_cur$Size,col=ifelse(tops_cur$FDR<0.05,"red","black"))+
    geom_hline(yintercept = -log10(0.05),size=0.1,col="red")+
    xlab("Log2(CPM) Fold Change")+
    ylab("-log10(False Discovery Rate)")+
    scale_y_continuous(n.breaks=10,limits=c(NA,1.1*max(-log10(tops_cur$FDR))))+
    theme_bw()+
    plot_theme()+
    theme(axis.text.x = element_text(angle=0,hjust=0))
  
}

setwd(directory)
ggsave("tmp.png",plot=plot,width=600,height=500,units = "px")
tst2<-tkimage.create("photo","test",
                     file=paste(directory,"/tmp.png",sep=""))
tkconfigure(plot_display,image=tst2)
file.remove(paste(directory,"/tmp.png",sep=""))
