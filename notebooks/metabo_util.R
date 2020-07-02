#-------------------------------
# R funtions for metabolomics analysis.
#-------------------------------
# Last update: 20181216, HuangShi

p <- c("reshape","ggplot2","pheatmap","pROC","combinat","plyr","vegan","optparse", "igraph")
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE, repos =c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
    suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
suppressWarnings(suppressMessages(invisible(lapply(p, usePackage))))

#clustering for "types"
#-------------------------------
## KL/JS divergence measure for relative-abundance(density/frequency) data
#-------------------------------
JSD<-function(object, eps=10^-4, overlap=TRUE,...)
{
    if(!is.numeric(object))
        stop("object must be a numeric matrix\n")
    
    z <- matrix(NA, nrow=ncol(object), ncol=ncol(object))
    colnames(z) <- rownames(z) <- colnames(object)
    
    w <- object < eps
    if (any(w)) object[w] <- eps
## If you takes as input a matrix of density values 
## with one row per observation and one column per 
## distribution, add following statement below.
 # object <- sweep(object, 2, colSums(object) , "/")
    
    for(k in seq_len(ncol(object)-1)){
      for(l in 2:ncol(object)){
        ok <- (object[, k] > eps) & (object[, l] > eps)
        if (!overlap | any(ok)) {
          m=0.5*(object[,k]+object[,l])
          z[k,l] <- sqrt(0.5*sum(object[,k] *(log(object[,k]) - log(m)))+0.5*sum(object[,l] *(log(object[,l]) - log(m))))
          z[l,k] <- sqrt(0.5*sum(object[,l] *(log(object[,l]) - log(m)))+0.5*sum(object[,k] *(log(object[,k]) - log(m))))
        }
      }
    }
    diag(z)<-0
    z
}

#-------------------------------
PAM.best<-function(matrix,dist.mat){
#-------------------------------
# matrix: 
#         row.names	Sample_id
#         col.names	Varibles
# For example, data should be organized like this:
# Sample_id	V1	V2	V3	etc...
# sample_0001	15	6	25
# sample_0002	7	9	32
# etc...
#-------------------------------
 if(!is.numeric(matrix))
        stop("matrix must be a numeric matrix\n")
 if(!is.numeric(dist.mat) && class(dist.mat)=="dist")
        stop("dist.mat must be numeric distance matrix\n")
#-------------------------------
# nc - number_of_clusters
#-------------------------------
min_nc=2
if(nrow(matrix)>20){
max_nc=20} else {
max_nc=nrow(matrix)-1}
res <- array(0,c(max_nc-min_nc+1, 2))
res[,1] <- min_nc:max_nc
siavgs <- array(0,c(max_nc-min_nc+1, 2))
siavgs[,1] <- min_nc:max_nc
clusters <- NULL
for (nc in min_nc:max_nc)
{
cl <- pam(dist.mat, nc, diss=TRUE)
res[nc-min_nc+1,2] <- CH <- index.G1(matrix,cl$cluster,d=dist.mat,centrotypes="medoids")
siavgs[nc-1,2]<-cl$silinfo$avg.width
clusters <- rbind(clusters, cl$cluster)
}
CH_nc<-(min_nc:max_nc)[which.max(res[,2])]
Si_nc<-(min_nc:max_nc)[which.max(siavgs[,2])]
print(paste("max CH for",CH_nc,"clusters=",max(res[,2])))
print(paste("max Si for",Si_nc,"clusters=",max(siavgs[,2])))

CH_cluster<-clusters[which.max(res[,2]),]
Si_cluster<-clusters[which.max(siavgs[,2]),]
objectList      <- list()
    objectList$min_nc <- min_nc
    objectList$max_nc <- max_nc
    objectList$CH     <- CH_cluster
    objectList$Si     <- Si_cluster
    objectList$CH_nc  <- CH_nc
    objectList$Si_nc  <- Si_nc
    objectList$res    <- res
    objectList$siavgs <- siavgs
    return(objectList)
}
#--------------------------------------------------
log.mat<-function(mat,base=2){
  mat[mat==0]<-0.000001
  log.mat<-log(mat,base)
  return(log.mat)
}
#--------------------------------------------------
rangeScaling <- function(v) {
themin <- min(v)
themax <- max(v) 
v <- (v- themin) / (themax - themin)
return(v)
}
#--------------------------------------------------
paretoScaling <- function(v) {
themean <- mean(v)
thesd <- sd(v) 
v <- (v- themean) / sqrt(thesd)
return(v)
}
#--------------------------------------------------
PlotCorrHeatMap<-function(cor.method, colors, data){
    main <- xlab <- ylab <- NULL;
    #print (dim(data));
    #write.csv(data,file="test.csv"); # debug
    
    # Decide size of margins and length of labels
    labelLen<-max(nchar(colnames(data)))
    margins<-c(12, 12)

    if(ncol(data) > 500){
        filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
        rk <- rank(-filter.val, ties.method='random');
        data <- as.data.frame(data[,rk < 500]);
        
        print("Data is reduced to 500 vars ..");
    }

    #colnames(data)<-substr(colnames(data), 1, labelLen);
     corr.mat<-cor(data, method=cor.method);

    # set up parameter for heatmap
    suppressMessages(require(RColorBrewer));
    suppressMessages(require(gplots));
    if(colors=="jet"){
        colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(256) 
    }else if(colors=="gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(colors == "heat"){
        colors <- heat.colors(256);
    }else if(colors == "topo"){
        colors <- topo.colors(256);
    }else{
        colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256));
    }
    
    heatmap<-heatmap.2(corr.mat,
             Rowv=TRUE,
             Colv=TRUE,
            #dendrogram= c("none"),
             distfun = dist,
             hclustfun = hclust,
             xlab = xlab,
             ylab = ylab,
             key=TRUE,
             keysize=0.8, # size of the key in the chart
             trace="none",
             density.info=c("none"),
             margins=margins,
            #col=brewer.pal(10,"PiYG")
             main = main,
             col=colors
             )
    objectList      <- list()
    objectList$heatmap  <- heatmap
    objectList$corr.mat <- corr.mat
    return(objectList)
}
#--------------------------------------------------
CleanData <-function(bdata, removeNA=T, removeNeg=T){
    if(sum(bdata==Inf)>0){
        inx <- bdata == Inf;
        bdata[inx] <- NA;
        bdata[inx] <- max(bdata, na.rm=T)*2
    }
    if(sum(bdata==-Inf)>0){
        inx <- bdata == -Inf;
        bdata[inx] <- NA;
        bdata[inx] <- min(bdata, na.rm=T)/2
    }
    if(removeNA){
        if(sum(is.na(bdata))>0){
            bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
        }
    }
    if(removeNeg){
        if(sum(bdata<=0) > 0){
            inx <- bdata <= 0;
            bdata[inx] <- NA;
            bdata[inx] <- min(bdata, na.rm=T)/2
        }
    }
    bdata;
}

#--------------------------------------------------Last update: 20170115
DistBoxplot<-function(dm,group,dm_name='',group_name='',IndividualID=NULL,outpath='./'){
    if(  ncol(dm)!=nrow(dm) & any(is.na(dm))==TRUE)
    stop('The distance matrix is not squared')
    if( length(unique(group))==1)
    stop('At least two levels for a given sample category in your metadata file are required.')
    if( length(group)!=nrow(dm))
    stop('The number of rows in metadata and distance matrix are not equal')
    if(nlevels(group)>length(group)*0.9)
    stop('The number of levels in a certain category can not exceed 90% of total number of samples')
    #--------------------------------
    dm<-dm[order(rownames(dm)),order(colnames(dm))]
    dm<-dm[order(group),order(group)]
    group_ordered<-group[order(group)]
    fac<-factor(group_ordered)
    names(fac)<-rownames(dm)
    
    if(!is.null(IndividualID)){
        Ind<-IndividualID[order(group)]
        names(Ind)<-rownames(dm)
    }
    #--------------------------------
    #install.packages('combinat')
    #require(combinat)
    require(plyr)
    require(pheatmap)
    require(ggplot2)
    #--------------------------------
    if(is.null(IndividualID)){
        colnames(dm)<-rownames(dm)<-paste(rownames(dm),fac,sep="____") }else{
        colnames(dm)<-rownames(dm)<-paste(rownames(dm),fac,Ind,sep="____") }
    
    dm[lower.tri(dm)]<-NA
    melt_dm<-melt(data.matrix(dm))
    melt_dm<-melt_dm[!is.na(melt_dm$value),]
    melt_dm<-melt_dm[which(melt_dm[, 1]!=melt_dm[, 2]),]
    
    Row_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[, 1]),"____",fixed=TRUE)))
    Col_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[, 2]),"____",fixed=TRUE)))
    VS<-paste(Row_Info[,2],"_VS._",Col_Info[,2],sep="")
    dm_value<-data.frame(VS,Row_Info,Col_Info,d=melt_dm$value)
    
    if(is.null(IndividualID)){
    colnames(dm_value)<-c("GroupPair","Sample_1","Group_1","Sample_2","Group_2","Dist")
    DistType<-as.factor(dm_value$Group_1==dm_value$Group_2)
    DistType<-factor(DistType,levels=levels(DistType),labels=c("AllBetween","AllWithin"))
    dm_value<-data.frame(DistType,dm_value)
    }else{
    colnames(dm_value)<-c("GroupPair","Sample_1","Group_1","Subject_1","Sample_2","Group_2","Subject_2","Dist")
    Ind_dm_value<-dm_value[which(dm_value$Subject_1==dm_value$Subject_2),]
    }

    #--------------------------------Output distance table and boxplot
    if(!is.null(IndividualID)){
        group_name<-paste(group_name,"IndividualID",sep="4")
        filepath<-sprintf('%s%s%s%s%s',outpath,dm_name,'.',group_name,'.Ind_dm_values.xls')
        sink(filepath);write.table(Ind_dm_value,quote=FALSE,sep='\t',row.names = FALSE);sink()
        
        #-----Output distance boxplot
        
        plot<-qplot(x=GroupPair, y=Dist,  data=Ind_dm_value, geom='boxplot',main='',xlab="Group pair", ylab=paste(dm_name,'_Distance',sep='')) + coord_flip() + theme_bw()
     suppressMessages(ggsave(filename=paste(outpath,dm_name,'.',group_name,'.boxplot.ggplot.pdf',sep=''),plot=plot, height=ifelse(nlevels(fac)>2,nlevels(fac),2)))
    
    invisible(Ind_dm_value)
        
    }else{
        
        #-----Objects to return
        if(length(dm_value$Dist)>4){
            p_t<-with(dm_value,t.test(Dist~DistType))$p.value
            p_w<-with(dm_value,wilcox.test(Dist~DistType))$p.value}else{
                p_t<-NA
                p_w<-NA}
            
            #-----Output distance-values table
            filepath<-sprintf('%s%s%s%s%s',outpath,dm_name,'.',group_name,'.dm_values.xls')
            sink(filepath);write.table(dm_value,quote=FALSE,sep='\t',row.names=FALSE);sink()
            
            #-----Output distance boxplot
            if(nlevels(group)<30){
            plot<-qplot(x=GroupPair, y=Dist, data=dm_value, geom='boxplot',position='dodge',main='',xlab="Group pair",ylab=paste(dm_name,'_Distance',sep='')) + coord_flip() +  theme_bw()
            suppressMessages(ggsave(filename=paste(outpath,dm_name,'.',group_name,'.boxplot.ggplot.pdf',sep=''),plot=plot))
            #-----Output statistical results
            p<-pairwise.wilcox.test(dm_value$Dist,factor(dm_value$GroupPair))$p.value
            sink(paste(outpath,dm_name,'.',group_name,'.dm_stats_P_values(Wilcoxon-test).xls',sep=''));cat('\t');write.table(p,quote=FALSE,sep='\t');sink()
            if(length(unique(as.vector(p)[!is.na(p)]))>1)
            pheatmap(p,cluster_rows = FALSE, cluster_cols = FALSE,display_numbers = T,main=group_name,filename=paste(outpath,dm_name,'.',group_name,'.dm_stats_P_values(Wilcoxon-test).pdf',sep=''))
            }
            
            objectList      <- list()
            objectList$dm_value     <- dm_value
            objectList$p_t    <- p_t
            objectList$p_w    <- p_w
            
            invisible(objectList)
    }
}

#-------------------------------Last update: 20181219
BetweenGroup.test <-function(data, group, p.adj.method="bonferroni", Enr_level=NA, q_cutoff=0.2, paired=FALSE){
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  
  n_group<-nlevels(group)
  Enr_level<-ifelse(is.na(Enr_level), levels(group)[1], Enr_level)
  if(!is.numeric(n_group) | n_group==1)
    stop("group must be a numeric and up to two levels\n")
  if(n_group==2){
    test_output<-matrix(NA,ncol=4,nrow=ncol(data))
    rownames(test_output)<-colnames(data)
    colnames(test_output)<- c("T.test_p","Wilcoxon.test_p","T.test_p.adj","Wilcoxon.test_p.adj")
    for(i in 1:ncol(data))
    {
      if(var(data[, i])==0){
        test_output[i,]<-rep(1, 4)
      }else{
        test_output[i,1]<-t.test(data[,i]~group,paired=paired)$p.value
        test_output[i,2]<-wilcox.test(data[,i]~group, paired=paired, conf.int=TRUE, exact=FALSE, correct=FALSE)$p.value
      }
    }
    test_output[,3]<-p.adjust(test_output[,1], method = p.adj.method, n = ncol(data))
    test_output[,4]<-p.adjust(test_output[,2], method = p.adj.method, n = ncol(data))
    
    func_bygroup_list<-c("mean", "sd", "median", "log10_median", "OccRate")
    OccRate<-function(x) sum(x!=0)/length(x)
    log10_median<-function(x, base=10) log.mat(median(x), base=base)
    tmp<-apply(data,2,function(x) tapply(x, group, function(x) c(mean(x), sd(x), median(x), log10_median(x), OccRate(x) ) ));
    stats_bygroup_df<-t(sapply(tmp, unlist));  colnames(stats_bygroup_df)<-unlist(lapply(levels(group), function(x) paste(func_bygroup_list, x, sep="_")))
    #-------------------------------
    auroc <- function(response, predictor) {
      require(ROCR)
      pred <- prediction(predictor, response)
      auroc  <- ROCR::performance(pred, "auc")@y.values[[1]]
      return(auroc)
    }
    func_all_list<-c("mean_all", "var_all", "sd_all", "OccRate_all", "AUC")
    stats_all_df<-data.frame(t(apply(data,2,function(x) c(mean(x), var(x), sd(x), OccRate(x), auroc(group,x)) )))
    colnames(stats_all_df)<-func_all_list
    #-------------------------------Enrichment
    require("plyr")
    IfSig<-as.factor(ifelse(test_output[, "Wilcoxon.test_p.adj"]< q_cutoff, "Sig", "NotSig"))
    logMeanAbd<-log.mat(t(apply(data,2,function(x) tapply(x, group, mean))))
    mean_logfc<-logMeanAbd[, 1]-logMeanAbd[, 2]
    Enr0<-factor(ifelse(mean_logfc>0, paste(Enr_level, "enriched", sep="_"), paste(Enr_level, "depleted", sep="_")))
    IfSigEnr=interaction(IfSig, Enr0)
    Enr<-mapvalues(IfSigEnr,c(paste("NotSig.",Enr_level,"_depleted",sep=""), 
                              paste("NotSig.",Enr_level,"_enriched",sep=""),
                              paste("Sig.",Enr_level,"_depleted",sep=""),
                              paste("Sig.",Enr_level,"_enriched",sep="")
    ), 
    c("Neutral", "Neutral", 
      paste(Enr_level,"_depleted",sep=""),
      paste(Enr_level,"_enriched",sep="")
    ))
    test_output<-data.frame(mean_logfc=mean_logfc, test_output, IfSig, IfSigEnr, Enr)
    #-------------------------------
    output1<-data.frame(stats_all_df, stats_bygroup_df, test_output)
    return(data.frame(output1))
  }else{
    output2<-matrix(NA,ncol=n_group+5,nrow=ncol(data))
    rownames(output2)<-colnames(data)
    colnames.output2<-array(NA)
    for(j in 1:ncol(output2)){
      if(j<=n_group){
        colnames.output2[j]<-c(paste("mean_",levels(group)[j],sep=""))
      }else{
        colnames.output2[(n_group+1):(n_group+5)]<-c("Var.test","Oneway-test_p","Kruskal.test_p",
                                                     "Oneway-test_p.adj","Kruskal.test_p.adj")
      }
    }
    colnames(output2)<-colnames.output2
    for(i in 1:ncol(data))
    {
      for(j in 1:n_group)
      {
        output2[i,j]<-mean(data[which(group==levels(group)[j]),i])
      }
      output2[i,(n_group+1)]<-bartlett.test(data[,i]~group)$p.value
      if(output2[i,(n_group+1)]<0.01)
        output2[i,(n_group+2)]<-oneway.test(data[,i]~group)$p.value
      else
        output2[i,(n_group+2)]<-oneway.test(data[,i]~group, var.equal=T)$p.value
      output2[i,(n_group+3)]<-kruskal.test(data[,i]~group)$p.value
      output2[i,(n_group+4)]<-NA
      output2[i,(n_group+5)]<-NA
    }
    output2[ ,(n_group+4)]<-p.adjust(output2[,(n_group+2)], method = p.adj.method, n = ncol(data))
    output2[ ,(n_group+5)]<-p.adjust(output2[,(n_group+3)], method = p.adj.method, n = ncol(data))
    return(data.frame(output2))
  }
}


#--------------------------------------------------

vectorize_dm<-function(dm, group=NULL, duplicate=TRUE)
{
  if(ncol(dm)!=nrow(dm) & any(is.na(dm))==TRUE)
    stop('The distance matrix is not squared')
  dm<-data.matrix(dm)
  require("reshape2")
  if(!is.null(group)){
    if( length(unique(group))==1)
      stop('At least two levels for a given sample category in your metadata file are required.')
    if( length(group)!=nrow(dm))
      stop('The number of rows in metadata and distance matrix are not equal')
    if(is.factor(group) & nlevels(group)>length(group)*0.9)
      stop('The number of levels in a certain category can not exceed 90% of total number of samples')
    
    colnames(dm)<-rownames(dm)<-paste(rownames(dm), group, sep="____")
    if(duplicate){ 
      melt_dm<-subset(melt(dm), value!=0 & value!=1) 
    }else{ 
      dm[lower.tri(dm)]<-NA
      melt_dm<-subset(melt(dm), ! is.na(value) & value!=0 & value!=1) 
    }
    
    Row_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[,1]),"____",fixed=TRUE)))
    Col_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[,2]),"____",fixed=TRUE)))
    VS<-paste(Row_Info[,2],"_VS._",Col_Info[,2],sep="")
    dm_value<-data.frame(VS,Row_Info,Col_Info,d=melt_dm$value)
    
    colnames(dm_value)<-c("GroupPair","Sample_1","Group_1","Sample_2","Group_2","value")
    if(is.factor(group)){
      DistType<-as.factor(dm_value$Group_1==dm_value$Group_2)
      DistType<-factor(DistType,levels=levels(DistType),labels=c("AllBetween","AllWithin"))
      dm_value<-data.frame(DistType,dm_value)}
  }else{
    if (duplicate) { 
      dm_value<-subset(melt(dm), value!=0 & value!=1)  
    }else{ 
      dm[lower.tri(dm)] <- NA
      dm_value <- subset(melt(dm), ! is.na(value) & value!=0 & value!=1)  }
    colnames(dm_value)<-c("Sample_1","Sample_2","value")
  }
  
  dm_value
  
}

getPermuteMatrix <- function(perm, N,  strata = NULL) {
  ## 'perm' is either a single number, a how() structure or a
  ## permutation matrix
  if (length(perm) == 1) {
    perm <- how(nperm = perm) 
  }
  ## apply 'strata', but only if possible: ignore silently other cases
  if (!missing(strata) && !is.null(strata)) {
    if (inherits(perm, "how") && is.null(getBlocks(perm)))
      setBlocks(perm) <- strata
  }
  ## now 'perm' is either a how() or a matrix
  if (inherits(perm, "how"))
    perm <- shuffleSet(N, control = perm)
  ## now 'perm' is a matrix (or always was). If it is a plain
  ## matrix, set minimal attributes for printing. This is a dirty
  ## kluge: should be handled more cleanly.
  if (is.null(attr(perm, "control")))
    attr(perm, "control") <-
      structure(list(within=list(type="supplied matrix"),
                     nperm = nrow(perm)), class = "how")
  perm
}
#-------------------------------
log_mat<-function(mat, base=2, peudozero=1e-6){
  mat[mat==0] <- peudozero
  log_mat <- log(mat, base)
  return(log_mat)
}
#-------------------------------
PAM.best<-function(matrix,dm){
  if(!is.numeric(matrix))
    stop("matrix must be a numeric matrix\n")
  if(!is.numeric(dist.mat) && class(dm)=="dist")
    stop("dist.mat must be numeric distance matrix\n")
  require("cluster")
  require("fpc")
  require("clusterSim")
  
  min_nc=2
  if(nrow(matrix)>20){
    max_nc=20} else {
      max_nc=nrow(matrix)-1}
  res <- array(0,c(max_nc-min_nc+1, 2))
  res[,1] <- min_nc:max_nc
  siavgs <- array(0,c(max_nc-min_nc+1, 2))
  siavgs[,1] <- min_nc:max_nc
  clusters <- NULL
  for (nc in min_nc:max_nc)
  {
    cl <- pam(dm, nc, diss=TRUE)
    res[nc-min_nc+1,2] <- CH <- index.G1(matrix,cl$cluster,d=dm,centrotypes="medoids")
    siavgs[nc-1,2]<-cl$silinfo$avg.width
    clusters <- rbind(clusters, cl$cluster)
  }
  CH_nc<-(min_nc:max_nc)[which.max(res[,2])]
  Si_nc<-(min_nc:max_nc)[which.max(siavgs[,2])]
  print(paste("max CH for",CH_nc,"clusters=",max(res[,2])))
  print(paste("max Si for",Si_nc,"clusters=",max(siavgs[,2])))
  
  CH_cluster<-clusters[which.max(res[,2]),]
  Si_cluster<-clusters[which.max(siavgs[,2]),]
  objectList      <- list()
  objectList$min_nc <- min_nc
  objectList$max_nc <- max_nc
  objectList$CH     <- CH_cluster
  objectList$Si     <- Si_cluster
  objectList$CH_nc  <- CH_nc
  objectList$Si_nc  <- Si_nc
  objectList$res    <- res
  objectList$siavgs <- siavgs
  return(objectList)
}
# ggplot2 color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


Trim_v_corr_mat<-function(v_Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6){
  if (Pos_Edge & Neg_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold | value < (-Threshold))
  }else if (Pos_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold & value < 1)
  }else if (Neg_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > -1 & value < (-Threshold))
  }else stop('Please check the pos/neg edges parameters')
  v_Corr_mat 
}   

Trim0_v_corr_mat<-function(v_Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6){
  if (Pos_Edge & Neg_Edge){
    v_Corr_mat[!with(v_Corr_mat, value > Threshold | value < (-Threshold)), "value"]<-0
    
  }else if (Pos_Edge){
    v_Corr_mat[!with(v_Corr_mat, value > Threshold & value < 1), "value"]<-0
  }else if (Neg_Edge){
    v_Corr_mat[!with(v_Corr_mat, value > -1 & value < (-Threshold)), "value"]<-0
  }else stop('Please check the pos/neg edges parameters')
  v_Corr_mat 
}


Creat_graph1<-function(Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6 ){
  if (Pos_Edge & Neg_Edge){
    g <- graph.adjacency(abs(Corr_mat) > Threshold & abs(Corr_mat)<1,mode="undirected")
  }else if (Pos_Edge){
    g <- graph.adjacency(Corr_mat > Threshold & Corr_mat < 1,mode="undirected")
  }else if (Neg_Edge){
    g <- graph.adjacency(Corr_mat < -1 * Threshold & Corr_mat > -1 ,mode="undirected")
  }else stop('Please check the pos/neg edges parameters')
  
  Pos_Edge_N <- 0
  Neg_Edge_N <- 0
  
  for(i in 2:nrow(Corr_mat)){
    for(j in 1:(i-1)){
      if ((Neg_Edge) & (Corr_mat[i,j] < -1 * Threshold) & (Corr_mat[i,j] > -1)){
        E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$color <- 'green'
        E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$weight <- round(abs(Corr_mat[i,j]),3)
        Neg_Edge_N <- Neg_Edge_N + 1
      }
      if ((Pos_Edge) & (Corr_mat[i,j] > Threshold) & (Corr_mat[i,j] < 1)){
        E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$color <- 'red'
        E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$weight <- round(abs(Corr_mat[i,j]),3)
        Pos_Edge_N <- Pos_Edge_N + 1
      }       
    }
  }
  g
}

Creat_graph0<-function(Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6 ){
  if (Pos_Edge & Neg_Edge){
    g <- graph.adjacency(abs(Corr_mat) > Threshold & abs(Corr_mat)<1,mode="undirected")
  }else if (Pos_Edge){
    g <- graph.adjacency(Corr_mat > Threshold & Corr_mat < 1,mode="undirected")
  }else if (Neg_Edge){
    g <- graph.adjacency(Corr_mat < -1 * Threshold & Corr_mat > -1 ,mode="undirected")
  }else stop('Please check the pos/neg edges parameters')
  
  v_Corr_mat<-vectorize_dm(Corr_mat, duplicate=FALSE)
  v_Corr_mat_trimmed<-Trim_v_corr_mat(v_Corr_mat, Pos_Edge=Pos_Edge, Neg_Edge=Neg_Edge, Threshold=Threshold)
  
  lapply(1:nrow(v_Corr_mat_trimmed), function(x){ 
    if(v_Corr_mat_trimmed[x, 3]<0){E(g, path=c(v_Corr_mat_trimmed[x, 1], v_Corr_mat_trimmed[x, 2]))$color<- "green"}else{
      E(g, path=c(v_Corr_mat_trimmed[x, 1], v_Corr_mat_trimmed[x, 2]))$color<- "red"
    } 
    E(g, path=c(v_Corr_mat_trimmed[x, 1], v_Corr_mat_trimmed[x, 2]))$weight<- abs(v_Corr_mat_trimmed[x, 3]) 
  })
  g
}

Creat_graph<-function(Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.9){
  v_Corr_mat<-vectorize_dm(Corr_mat, duplicate=FALSE)
  v_Corr_mat_trimmed<-Trim_v_corr_mat(v_Corr_mat, Pos_Edge=Pos_Edge, Neg_Edge=Neg_Edge, Threshold=Threshold)
  g<-graph.edgelist(as.matrix(v_Corr_mat_trimmed[, c(1,2)]), directed=FALSE)
  E(g)$weight<- round(abs(v_Corr_mat_trimmed[, 3]), 6)
  v_Corr_mat_trimmed[, "color"]<-ifelse(v_Corr_mat_trimmed[, 3]<0, "green", "red")
  E(g)$color<- v_Corr_mat_trimmed[, "color"]
  g
} 

map <- function(x, range = c(0,1), from.range=NA) {
  if(any(is.na(from.range))) from.range <- range(x, na.rm=TRUE)
  
  ## check if all values are the same
  if(!diff(from.range)) return(
    matrix(mean(range), ncol=ncol(x), nrow=nrow(x), 
           dimnames = dimnames(x)))
  
  ## map to [0,1]
  x <- (x-from.range[1])
  x <- x/diff(from.range)
  ## handle single values
  if(diff(from.range) == 0) x <- 0 
  
  ## map from [0,1] to [range]
  if (range[1]>range[2]) x <- 1-x
  x <- x*(abs(diff(range))) + min(range)
  
  x[x<min(range) | x>max(range)] <- NA
  
  x
}

Plot_network_graph<-function(Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.5, Edge_width="Corr", node_size="degree", node_color="", layout="layout.fruchterman.reingold", outdir="./out.pdf", width=10, height=10){
  v_Corr_mat<-vectorize_dm(Corr_mat, duplicate=FALSE)
  v_Corr_mat_trimmed<-Trim_v_corr_mat(v_Corr_mat, Pos_Edge=Pos_Edge, Neg_Edge=Neg_Edge, Threshold=Threshold)
  g<-graph.edgelist(as.matrix(v_Corr_mat_trimmed[, c(1,2)]), directed=FALSE)
  v_Corr_mat_trimmed[, "color"]<-ifelse(v_Corr_mat_trimmed[, 3]<0, gg_color_hue(2)[1], gg_color_hue(2)[2])
  E(g)$color<- v_Corr_mat_trimmed[, "color"]
  
  Pos_Edge_N <- length(which(v_Corr_mat_trimmed[, "value"]>0))
  Neg_Edge_N <- length(which(v_Corr_mat_trimmed[, "value"]<0))
  
  Node_num<-length(V(g))
  Edge_num<-length(E(g))
  Ave_degree<-round(mean(igraph::degree(g)), 4)
  Degree_assortativity<-round(assortativity_degree(g), 4)
  Ave_path_length<-round(average.path.length(g), 4)
  Density<-round(graph.density(g), 4)
  Diameter<-diameter(g)
  Radius<-radius(g)
  Cluster_num<-clusters(g)$no
  Modularity<-modularity(g, clusters(g)$membership)
  Transitivity<-transitivity(g)
  Closeness_centralization<-round(centralization.closeness(g)$centralization, 4)
  Betweenness_centralization<-round(centralization.betweenness(g)$centralization, 4)
  Degree_centralization<-round(centralization.degree(g)$centralization, 4)
  Cluster_membership<-clusters(g)$membership
  # node size setup
  if(node_size=="degree" & var(igraph::degree(g))!=0) {
    vertex_size<-map(igraph::degree(g),c(1,20))
  }else if(node_size=="Betweenness_centralization" & var(igraph::centralization.betweenness(g)$res)!=0) {
    vertex_size<-map(igraph::centralization.betweenness(g)$res,c(1,10))
  }else{
    vertex_size<-rep(1, Node_num)
  }
  # node color setup
  if(node_color=="degree" & var(igraph::degree(g))!=0){
    vertex_color<-try(map(igraph::degree(g), c(1,10))) 
  }else if (node_color=="Cluster_membership"){
    vertex_color<-Cluster_membership
  }else if (node_color==""){
    vertex_color<-rep("grey80", Node_num)
  }
  # edge width setup
  Corr<-round(abs(v_Corr_mat_trimmed[, 3]), 6)
  if(Edge_width=="Corr" && var(Corr)!=0) {
    E(g)$weight<- Corr
    edge_width<-try(map(Corr,c(1, 5))) 
    }else if(Edge_width=="Betweenness_centralization" && var(igraph::centralization.betweenness(g)$res)!=0) {
      E(g)$weight<- igraph::centralization.betweenness(g)$res
      edge_width<-try(map(igraph::centralization.betweenness(g)$res, c(1,5))) 
      }else{
      edge_width=rep(2, Edge_num)
  }
  
  if(!is.null(outdir)){
    pdf(outdir, width=width, height=height)
    op<-par(mar=c(2,2,6,2))
    str = paste(
      "  Threshold: ",Threshold,
      "  Nodes: ", Node_num,
      "  Edges: ", Edge_num,
      "  # of Pos edge: ",Pos_Edge_N,  
      "  # of Neg edge: ",Neg_Edge_N,
      "\nAverage_degree: ", Ave_degree,
      "  Degree_assortativity: ", Degree_assortativity,
      "  Clusters: ", Cluster_num, 
      "  Modularity: ", Modularity,
      "  Average_path_length: ",Ave_path_length,
      "\nDensity: ", Density,
      "  Diameter: ", Diameter,
      "  Radius: ",Radius,  
      "  Transitivity: ",Transitivity,
      "\nCloseness_centralization:", Closeness_centralization,
      "  Betweenness_centralization:", Betweenness_centralization,
      "  Degree_centralization:", Degree_centralization
    )
    l<-get(layout)(g)
    plot(g, main=str, edge.label=E(g)$weight, 
         edge.label.cex = 0.25, 
         vertex.frame.color="white", 
         vertex.size=vertex_size, 
         vertex.color=vertex_color, 
         edge.width=edge_width,
         vertex.label.cex=0.5) 
    if(Ave_degree!=0){
      legend("bottomright", c("neg","pos"), lty=c(1,1), col=gg_color_hue(2) )
    }
    par(op)
    invisible(dev.off())
  }
  objectList      <- list()
  objectList$Node_num                      <- Node_num
  objectList$Edge_num                      <- Edge_num 
  objectList$Pos_Edge_num                  <- Pos_Edge_N 
  objectList$Neg_Edge_num                  <- Neg_Edge_N 
  objectList$Ave_degree                    <- Ave_degree
  objectList$Degree_assortativity          <- Degree_assortativity 
  objectList$Ave_path_length               <- Ave_path_length
  objectList$Density                       <- Density
  objectList$Diameter                      <- Diameter
  objectList$Radius                        <- Radius
  objectList$Cluster_num                   <- Cluster_num
  objectList$Modularity                    <- Modularity
  objectList$Transitivity                  <- Transitivity
  objectList$Closeness_centralization      <- Closeness_centralization
  objectList$Betweenness_centralization    <- Betweenness_centralization
  objectList$Degree_centralization         <- Degree_centralization
  objectList$g                             <- g
  objectList$Cluster_membership            <- Cluster_membership
  
  invisible(objectList)
  
}                              
