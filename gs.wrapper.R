gs.wrapper = function(setlist, backgroundset, 
                      filename=NULL, 
                      functiondir=NULL,resourcedir=NULL,
                      gs.range=c(25,2500), ecut=0.05, ocut=5, include.identifiers=FALSE){
  # function to run basic pathway annotation on an input list of gene sets, using predefined annotation sets.  
  # please adjust locations of functions and annotation files for your system
  # backgroundset is a data.frame with fields ID, Entrez.ID and Anno.Symbol  
  # Entrez.ID must be HUMAN IDs
  # setlist is a list containing vectors of identifiers in backgroundset ID
  #  each element is named with the display name for the gene set 
  
  if(is.null(functiondir)) functiondir=getwd()
  if(is.null(resourcedir)) resourcedir=getwd()
  
  # imports
  
  source(file.path(functiondir,'gs.pval.R'))
  tmp=suppressPackageStartupMessages(require(AnnotationDbi))
  tmp=suppressPackageStartupMessages(require(data.table))
  
  # check arguments
  # setlist should be list of nonzero-length vectors
  if(!is.list(setlist) || sum(lengths(setlist)>0) != length(setlist) ){
    stop("setlist must be a list of nonzero-length vectors")
  }
  
  # load annotations
  # need to set up DB to enable non-hardcoded method!
  # annotations with EntrezGeneIDs
  load(file.path(resourcedir,"hsa.kegg.gs.RData")); kegg=hsa.kegg.gs
  kegg = lapply(hsa.kegg.gs,function(x) as.numeric(x))
  load(file.path(resourcedir,"msigdb.v5.0.entrez.gmt.RData"))
  msigdb5 = lapply(msigdb5,function(x) as.numeric(x))
  # all 3 GOs in one list; genes as UniProt IDs in pipe-delimited field 
  gotry =   try(load(paste0(resourcedir,"GO_processed.gmt.RData")))
  if(class(gotry)=="try-error") { #This is slow, so load processed GO if possible
    load(file.path(resourcedir,"GO.gmt.RData"))
    # pull out individual GOs and translate their UniProt accessions to EntrezIDs
    GOBP= lapply(GO_gmt[[1]],function(x){
      as.numeric(intraIDMapper(sub('^[^|]+\\|([^|]+)\\|.*$','\\1',x), species="HOMSA", srcIDType="UNIPROT", destIDType="EG", keepMultGeneMatches=FALSE)) }
    )
    GOMF= lapply(GO_gmt[[2]],function(x){
      as.numeric(intraIDMapper(sub('^[^|]+\\|([^|]+)\\|.*$','\\1',x), species="HOMSA", srcIDType="UNIPROT", destIDType="EG", keepMultGeneMatches=FALSE)) }
    )
    GOCC= lapply(GO_gmt[[3]],function(x){
      as.numeric(intraIDMapper(sub('^[^|]+\\|([^|]+)\\|.*$','\\1',x), species="HOMSA", srcIDType="UNIPROT", destIDType="EG", keepMultGeneMatches=FALSE)) }
    )
    rm(GO_gmt)
    save(GOBP,GOMF,GOCC,file=file.path(resourcedir,"GO_processed.gmt.RData"))
  }
  # annotations with geneSymbols; not a robust identifier type
  load(file.path(resourcedir,"pathwayCommons.gmt.RData"))
  load(file.path(resourcedir,"NCI-Nature_Curated.gmt.RData"));nci=nci_gmt; rm(nci_gmt)
  load(file.path(resourcedir,"BioCarta.gmt.RData"));BioCarta=BioCarta_gmt;rm(BioCarta_gmt)
  load(file.path(resourcedir,"Reactome.gmt.RData"));Reactome=Reactome_gmt;rm(Reactome_gmt)
  annolist = c('kegg','msigdb5','GOBP','GOMF','GOCC','pC','nci','BioCarta','Reactome')
  idtypes = paste('ids',c(rep('ncbi',5),rep('gene',4)),sep='.')
  use.gsr = c(F,T,T,T,T,T,F,F,F)
  # define universes for each gene set source
  uni.kegg = unique(unlist(kegg,F,F))
  uni.pC = unique(unlist(pC,F,F)); uni.pC = toupper(uni.pC)
  uni.msigdb5 = unique(unlist(msigdb5,F,F))
  uni.GOBP = unique(unlist(GOBP,F,F))
  uni.GOMF = unique(unlist(GOMF,F,F))
  uni.GOCC = unique(unlist(GOCC,F,F))
  uni.nci = unique(unlist(nci,F,F)); uni.nci=toupper(uni.nci)
  uni.BioCarta = unique(unlist(BioCarta,F,F)); uni.BioCarta = toupper(uni.BioCarta)
  uni.Reactome = unique(unlist(Reactome,F,F)); uni.Reactome = toupper(uni.Reactome)
  
  # annotate up- and down-regulated selected genes in files
  # expand this to include universe size, sig size and collapsed identifiers if requested
  out.dt = data.table(rn="",E.Value="",P.Value="",Set.Size="",Overlap="",Sig.Size="",Uni.Size="",Uni.Size.Intersect.Data="",Source="",Condx="")
  if(include.identifiers){
    out.dt = cbind(out.dt,ID.list="",alt.ID.list="")
  }
  out.dt = out.dt[0,]
  # loop over anno sources 
  for( src in annolist ){
    a = which(annolist==src)
    if( use.gsr[a] ){ min.size=gs.range[1]; max.size=gs.range[2]
    }else{ min.size=ocut; max.size=Inf }
    if(idtypes[a]=="ids.ncbi"){ tmp = backgroundset$Entrez.ID
    }else{ tmp = backgroundset$Anno.Symbol }
    for( i in 1:length(setlist) ){
      signame = names(setlist)[i]; sig = setlist[[i]]
      if(include.identifiers){
        if(idtypes[a]=="ids.ncbi"){ # provide gene symbols for human readability
          GS = gs.pval(gs=get(src), rnames=backgroundset$Entrez.ID, is.sig = backgroundset$ID %in% sig, 
                       uni.size=length(unique(get(paste('uni',src,sep='.')))),
                       min.size=min.size, max.size=max.size, include.identifiers=TRUE, 
                       alt.rnames=backgroundset$Anno.Symbol)
          uni.size.int = sum(backgroundset$Entrez.ID%in%unique(get(paste('uni',src,sep='.'))))
        } else {
          GS = gs.pval(gs=get(src), rnames=backgroundset$Anno.Symbol, is.sig = backgroundset$ID %in% sig, 
                       uni.size=length(unique(get(paste('uni',src,sep='.')))),
                       min.size=min.size, max.size=max.size, include.identifiers=TRUE, 
                       alt.rnames=backgroundset$Entrez.ID)
          uni.size.int = sum(backgroundset$Anno.Symbol%in%unique(get(paste('uni',src,sep='.'))))
        }
      } else {
        GS = gs.pval(gs=get(src), rnames=tmp, is.sig = backgroundset$ID %in% sig, 
                     uni.size=length(unique(get(paste('uni',src,sep='.')))),
                     min.size=min.size, max.size=max.size, include.identifiers=FALSE)
        uni.size.int = sum(tmp%in%unique(get(paste('uni',src,sep='.'))))
      }
      sig.n = sum(GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut,na.rm=T)
      print(paste(signame, src,"uni.size",attr(GS,"Uni.Size"),"intersect size",uni.size.int, "sig.n", sig.n))
      if(sig.n){
        # convert slice of GS to data.frame to preserve structure if n==1
        out.dt = rbindlist(list(out.dt,data.table(GS,keep.rownames=T)[GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut,]),fill=T)
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Sig.Size"] = attr(GS,"Sig.Size")
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Uni.Size"] = attr(GS,"Uni.Size")
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Uni.Size.Intersect.Data"] = uni.size.int
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Source"] = src
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Condx"] = signame
        if(include.identifiers){
          out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"ID.list"] = attr(GS,"ID.list")[GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut]
          out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"alt.ID.list"] = attr(GS,"alt.ID.list")[GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut]
        }
      }
    }
  }
  names(out.dt)[names(out.dt)=="rn"] = "Set.Name"
  if(!is.null(filename)){
    write.csv(out.dt,file=filename,row.names=F,quote=F)
  }
  return(out.dt)
  
}

