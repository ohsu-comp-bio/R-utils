gs.wrapper = function(setlist, backgroundset, include.identifiers=FALSE, 
                      filename=NULL, resourceset=NULL, 
                      functiondir=NULL, resourcedir=NULL, return.OM=FALSE,
                      gs.range=c(25,2500), ecut=0.05, ocut=5, anno.uni=TRUE,  
                      annolist = c('kegg','msigdb5','GOBP','GOMF','GOCC',
                                   'pC','nci','BioCarta','Reactome'), 
                      idtypes = c(rep('Entrez.ID',5),rep('Anno.Symbol',4)),
                      use.gsr = c(F,T,T,T,T,T,F,F,F) ){
  # function to run basic pathway annotation on an input list of gene sets, using predefined annotation sets.  
  # gs.range limits the sizes of gene sets used for annotation
  # ecut sets the maximum E-value for annotations to return
  # ocut sets the minimum identifier number for annotations to return
  # return.OM is a flag indicating whether occurrence matrix is to be returned
  #   occurrence matrix == logical matrix, genes x sets
  # anno.uni is a flag: if TRUE, use the annotation source's universe size; 
  #   or if FALSE, use the overlap of backgroundset & anno universes
  # functiondir = optional directory with functions to be loaded (default wd)
  # resourcedir = optional directory with genesets to be loaded (default wd)
  # please adjust locations of functions and annotation files for your system
  # if resourceset is NULL (default) a default set of annotation sources
  #   is loaded, with all arguments ready as defaults
  # otherwise:
  # resourceset is a vector of filenames containing annotation
  #   sources as lists of named elements, each of which contains a vector
  #   of identifiers. Each file is an RData file. The filename should be
  #   the same as that of the list varname inside the file.
  # idtypes is a vector of the same length as resourceset that contains
  #   backgroundset colnames corresponding to each resourceset's id type
  # annolist is a vector of annotation source variable names as strings
  # use.gsr is a logical vector setting whether to use gs.range to restrict
  #   annotation source by gene set size
  # backgroundset is a data.table with fields ID, Entrez.ID and Anno.Symbol  
  # Entrez.ID must be HUMAN IDs to use default annotation sources
  # any other identifier types required by anno sources are in addnl columns
  # setlist is a list containing vectors of identifiers in backgroundset ID
  #  each element is named with the display name for the gene set 
  # Note that setlist MUST include a names() attribute.
  
  if(is.null(functiondir)) functiondir=getwd()
  if(is.null(resourcedir)) resourcedir=getwd()
  
  # imports
  source(file.path(functiondir,'gs.pval.R'))
  tmp=suppressPackageStartupMessages(require(data.table))
  
  # check arguments
  # setlist should be list of nonzero-length vectors
  if(!is.list(setlist) || sum(lengths(setlist)>0) != length(setlist) ){
    stop("setlist must be a list of nonzero-length vectors")
  }
  # get ready to save occurrence matrices if requested
  if( return.OM ){ OM_ls = NULL }
  
  # load annotations
  if( is.null(resourceset) ){
    message("Using default annotation sets")
    # annotations with EntrezGeneIDs
    load(file.path(resourcedir,"hsa.kegg.gs.RData")); kegg=hsa.kegg.gs
    kegg = lapply(hsa.kegg.gs,function(x) as.numeric(x))
    load(file.path(resourcedir,"msigdb.v5.0.entrez.gmt.RData"))
    msigdb5 = lapply(msigdb5,function(x) as.numeric(x))
    load(paste0(resourcedir,"GO_processed.gmt.RData"))
    # annotations with geneSymbols; not a robust identifier type
    load(file.path(resourcedir,"pathwayCommons.gmt.RData"))
    load(file.path(resourcedir,"NCI-Nature_Curated.gmt.RData"));nci=nci_gmt; rm(nci_gmt)
    load(file.path(resourcedir,"BioCarta.gmt.RData"));BioCarta=BioCarta_gmt;rm(BioCarta_gmt)
    load(file.path(resourcedir,"Reactome.gmt.RData"));Reactome=Reactome_gmt;rm(Reactome_gmt)
    # define universes for each gene set source
    uni = list(uni.kegg = unique(unlist(kegg,F,F)),
      uni.pC = unique(unlist(pC,F,F)),
      uni.pC = toupper(uni.pC),
      uni.msigdb5 = unique(unlist(msigdb5,F,F)),
      uni.GOBP = unique(unlist(GOBP,F,F)),
      uni.GOMF = unique(unlist(GOMF,F,F)),
      uni.GOCC = unique(unlist(GOCC,F,F)),
      uni.nci = toupper(unique(unlist(nci,F,F))),
      uni.BioCarta = toupper(unique(unlist(BioCarta,F,F))),
      uni.Reactome = toupper(unique(unlist(Reactome,F,F))) 
    )
  } else { #non-default load
    uni = NULL
    for( i in 1:length(resourceset) ){
      resourcename = sub('.RData','',resourceset[i])
      message("loading ",file.path(resourcedir,resourceset[i]))
      stopifnot( file.exists(file.path(resourcedir,resourceset[i])) )
      load(file.path(resourcedir,resourceset[i]))
      uni[[ paste("uni",resourcename,sep='.') ]] = 
             unique(unlist(get(resourcename)))
    }
    for( src in annolist ){
      message(src," : ",length(uni[[paste('uni',src,sep='.')]]))
    }
  }  
  # annotate up- and down-regulated selected genes in files
  out.dt = data.table(rn="",E.Value="",P.Value="",Set.Size="",Overlap="",Sig.Size="",Uni.Size="",Source="",Condx="",FoldEnrich="")
  if(include.identifiers){
    out.dt = cbind(out.dt,ID.list="",alt.ID.list="")
  }
  out.dt = out.dt[0,]
  # loop over anno sources 
  for( src in annolist ){
    a = which(annolist==src)
    # set anno parameters
    # range of gene set sizes
    if( use.gsr[a] ){ min.size=gs.range[1]; max.size=gs.range[2]
    }else{ min.size=ocut; max.size=Inf }
    # ID type(s) to use
    tmp = backgroundset[,get(idtypes[a])]
    alt.idtype = "Anno.Symbol" # default == Symbol for human readability
    if( idtypes[a]=="Anno.Symbol" ){ alt.idtype = "Entrez.ID" }
    # universe size to use
    if( anno.uni ){
      uni.size=length(uni[[paste('uni',src,sep='.')]])
    } else { uni.size = NULL }
    # annotate each set
    for( i in 1:length(setlist) ){
      signame = names(setlist)[i]; sig = setlist[[i]]
      if(include.identifiers){
        anno = gs.pval(gs=get(src), rnames=tmp, 
            is.sig = backgroundset$ID %in% sig, uni.size=uni.size,
            min.size=min.size, max.size=max.size, include.identifiers=TRUE, 
            alt.rnames=backgroundset$Anno.Symbol,return.OM=return.OM)
      } else {
        anno = gs.pval(gs=get(src), rnames=tmp, 
            is.sig = backgroundset$ID %in% sig, uni.size=uni.size, 
            min.size=min.size, max.size=max.size, include.identifiers=FALSE,
            return.OM=return.OM)
      }
      GS = anno$GS
      if( return.OM ){ 
        OM_ls = c(OM_ls, anno['OM'])
        names(OM_ls)[length(OM_ls)] = signame 
      }
      sig.n = sum(GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut,na.rm=T)
      print(paste(signame, src,"uni.size",attr(GS,"Uni.Size"),"sig.n", sig.n))
      if(sig.n){
        # convert slice of GS to data.frame to preserve structure if n==1
        out.dt = rbindlist(list(out.dt,data.table(GS,keep.rownames=T)[GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut,]),fill=T)
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Sig.Size"] = attr(GS,"Sig.Size")
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Uni.Size"] = attr(GS,"Uni.Size")
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Source"] = src
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Condx"] = signame
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"FoldEnrich"] = 
      ( as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Overlap"]) /
        as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Sig.Size"]) ) /
      ( as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Set.Size"]) /
        as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Uni.Size"]) ) 
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
  return(list(out.dt=out.dt, OM_ls=OM_ls))
  
}

