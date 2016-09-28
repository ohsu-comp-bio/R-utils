GMT_from_BioPax = function(biopaxstring){
  # extracts gene sets from OWL file as gene lists in GMT format
  # currently set up for BioPax2 only
  # returns list in GMT format
  require('rBiopaxParser')
  # process argument
  if( is.character(biopaxstring) ){ # assume file name
    biopax = readBiopax(biopaxstring)
  } else { biopax = biopaxstring }
  # check structure
  if( is.null(names(biopax)) | sum(names(biopax)=="dt")==0 | 
      !"biopax" %in% class(biopax) ){
    stop(paste(biopaxstring,"appears not to be a BioPax object"))
  }
  pw_component_list = listPathwayComponents(biopax, listInstances(biopax,class="pathway")$id, returnIDonly = T)
  if (length(pw_component_list) == 0) {
      stop("Pathways seem to have no pathway components")
  }
  if( sum(biopax$dt$class %chin% c("control", "catalysis", "modulation", "templatereactionregulation")) == 0 ){
      stop("Pathways have no regulatory components" )
  }
  # find pathways and pathway components
  p2pc = subset(biopax$dt,class=='pathway' & property=='PATHWAY-COMPONENTS',c('id','property_attr_value'))
  p2pc$property_attr_value = sub('^#','',p2pc$property_attr_value)
  # find next ID downstream of relevant pathway components
  pc2x = subset(biopax$dt,property %chin% toupper(c('controller','controlled','left','right','product')) & id %chin% p2pc$property_attr_value,c('class','id','property','property_attr_value'))
  setnames(pc2x,'property_attr_value','next_id')
  setnames(pc2x,'id','property_attr_value')
  pc2x$next_id = sub('^#','',pc2x$next_id)
  # track next IDs back to pathways 
  p2next = merge(p2pc,pc2x,by='property_attr_value',all.x=T,allow.cartesian=T)
  setnames(p2next,'property_attr_value','pathway_component')
  # at this point, some IDs are one step from a physical entity, others farther
  # retrieve and merge next layer..
  pc2pe = subset(biopax$dt, property=='PHYSICAL-ENTITY' & id %chin% pc2x$next_id, c('class','id','property','property_attr_value'))
  pc2n = subset(biopax$dt, property %chin% toupper(c('controller','controlled','left','right','product')) & id %chin% pc2x$next_id, c('class','id','property','property_attr_value'))
  setnames(pc2n,old=c('id','property_attr_value','class','property'),new=c('next_id','id3','next_class','next_property'))
  setnames(pc2pe,'id','next_id')
  # second round of mapping for next-layer IDs mapping to non-physical entities
  pc2n$id3 = sub('^#','',pc2n$id3)
  pc2pe2= subset(biopax$dt, property=='PHYSICAL-ENTITY' & id %chin% pc2n$id3, c('class','id','property','property_attr_value'))
  # bring in first-layer IDs for mapping back to pathways
  setnames(pc2pe2,'id','id3')
  pc2pe2n = merge(pc2n,subset(pc2pe2,select=c('id3','class','property','property_attr_value')),by='id3',all.y=TRUE)
  # check retrieval on remaining non-physical mappings
  pc2n2= subset(biopax$dt, property %chin% toupper(c('controller','controlled','left','right','product')) & id %chin% pc2n$property_attr_value, c('class','id','property','property_attr_value')) # expect 0 rows
  if( length(pc2n2$id)>0 ){
    warning(paste(length(unique(pc2n2$id)),"recursive mappings not pursued"))}
  # combine harvested mappings; adding mappings doesn't depend on new entities
  pc2pe.2x = rbindlist(list(pc2pe,pc2pe2n),use.names=T,fill=T)
  pc2pe.2x$property_attr_value = sub('^#','',pc2pe.2x$property_attr_value)
  pc2pe.2x = pc2pe.2x[!data.table:::duplicated.data.table(pc2pe.2x,by=NULL)]
  # pull mapped entities
  mymk.pe = biopax$dt$property_attr == "rdf:datatype" & biopax$dt$property == "NAME" & biopax$dt$class !='complex' & ! sub('^#','',biopax$dt$property_attr_value) %chin% pc2pe.2x$property_attr_value & biopax$dt$id %chin% pc2pe.2x$property_attr_value
  # pull mapped complexes
  pc2c = subset(biopax$dt,property_attr == "rdf:resource" & property %chin% c("COMPONENTS","PHYSICAL-ENTITY") & ! sub('^#','',property_attr_value) %chin% pc2pe.2x$property_attr_value & id %chin% pc2pe.2x$property_attr_value,c('id','property_attr_value'))
  pc2c = pc2c[!duplicated(pc2c)]
  setnames(pc2c,old=names(pc2c),new=c('property_attr_value','id4'))
  pc2c$id4 = sub('^#','',pc2c$id4)
  pc2c2pe = subset(biopax$dt, property=='PHYSICAL-ENTITY' & id %chin% pc2c$id4, c('class','id','property','property_attr_value'))
  setnames(pc2c2pe,old=c('id','property_attr_value'),new=c('id4','component_prop_attr_val'))
  pc2c2pe$component_prop_attr_val = sub('^#','',pc2c2pe$component_prop_attr_val)
  # bring in previous IDs for mapping back to pathways
  pc2pe2c = merge(pc2c,pc2c2pe,by='id4',all=TRUE)
  # pull entities mapped within complexes
  mymk.cpe= biopax$dt$property_attr == "rdf:datatype" & biopax$dt$property == "NAME" & biopax$dt$class !='complex' & ! sub('^#','',biopax$dt$property_attr_value) %chin% pc2pe2c$component_prop_attr_val & biopax$dt$id %chin% pc2pe2c$component_prop_attr_val
  # entity to name translation
  pe2name = subset(biopax$dt,mymk.pe|mymk.cpe,c('id','property_value'))
  setnames(pe2name,old=names(pe2name),new=c("property_attr_value","name"))
  # bring final mappings and entity names back to pathways
  pc2pe2name1 = merge(subset(pc2pe.2x,!pc2pe.2x$property_attr_value %chin% pc2pe2c$property_attr_value,c('property_attr_value','next_id','id3')),pe2name,all.x=T,by='property_attr_value')
  setnames(pe2name,'property_attr_value','component_prop_attr_val')
  pc2pe2name2 = merge(merge(subset(pc2pe.2x,pc2pe.2x$property_attr_value %chin% pc2pe2c$property_attr_value,c('property_attr_value','next_id','id3')),subset(pc2pe2c,select=c('property_attr_value','id4','component_prop_attr_val')),all.x=T,by='property_attr_value',allow.cartesian=T),pe2name,all.x=T,by='component_prop_attr_val')
  pc2pe2name = rbindlist(list(pc2pe2name1,subset(pc2pe2name2,select=names(pc2pe2name1))),use.names=T,fill=T)
  pc2pe2name = pc2pe2name[!data.table:::duplicated.data.table(pc2pe2name,by=NULL)]
  # bring mappings and names back to pathway IDs
  p2name = merge(p2next,pc2pe2name,all.x=T,by='next_id',allow.cartesian=T)
  p2name = p2name[!is.na(p2name$name)]
  # convert to GMT format
  pw_list = listInstances(biopax,class="pathway")
  mymk = pw_list$id %chin% unique(p2name$id)
  if( sum(mymk) < length(pw_list$id) ){
    warning(paste(sum(mymk),"of",length(pw_list$id),"pathways gave mapped named physical entities"))
    pw_list = pw_list[mymk,]
  }
  pw_gmt = vector(mode="list",length=length(pw_list$id))
  names(pw_gmt) = pw_list$name
  for(i in 1:length(pw_list$id) ){
    pw_gmt[[i]] = unique(p2name$name[p2name$id==pw_list$id[i]])
  }

  return(pw_gmt)
}


GMT_from_GOowl = function(GOowlstring,GOannostring,trim_regex=NULL){
  # TODO: add alternate reading, and argument protection
  # reads in GO ontology in GOowl format and as parsed in GO.db
  #  and reads in species-specific gene mappings from GOanno file
  # optionally trims gene identifiers using trim_regex
  # returns list of GMT-format lists @ for GO BP, MF and CC
  require("XML")
  require("data.table")
  # import processed GO tree
  require("GO.db")
  GOtags = list(P = 'BP', F = 'MF', C = 'CC')
  # read go OWL
  tmp = xmlRoot(xmlTreeParse(GOowlstring))
  tmp = tmp$children[xmlSApply(tmp,xmlName)=='Class']
  # fill GO ontology data table
  GOtab = data.table(rowcount=1:length(tmp),id="",label="",hasOBONamespace="", subClassOf="",key="rowcount")
  for(i in 1:length(tmp)){
    tmp2 = xmlApply(tmp[[i]],xmlValue)
    GOtab[i,'deprecated'] = length(xmlElementsByTagName(tmp[[i]],'deprecated'))>0
    # parse rest if not deprecated (data may not be available)
    if( !GOtab[i,get('deprecated')] ){
      GOtab[i,c("id","label","hasOBONamespace")] = tmp2[c("id","label","hasOBONamespace")]
      GOtab[i,'subClassOf'] = paste(sub('^.*\\/(GO_[0-9]+).*$','\\1',xmlElementsByTagName(tmp[[i]],'subClassOf')),sep='',collapse=';')
    }
  }
  GOtab = subset(GOtab,subset=deprecated==F,select=c("id","label","hasOBONamespace","subClassOf"))
  rm(tmp,tmp2)
  # read in tab-delimited GO-gene annotation in GAF format
  #  rows with NOT in column 4 refer to FALSE mappings to exclude
  tmp = fread(GOannostring,skip='\t')
  tmp = tmp[!grepl('NOT',tmp$V4)]; setkeyv(tmp,c('V5','V9','V2','V1'))
  tmp = subset(tmp,select=paste('V',c(1,2,3,5,9),sep=''))
  tmp = tmp[!data.table:::duplicated.data.table(tmp,by=NULL)]
  # join data tables to form GMT
  GOtype = unique(GOtab$hasOBONamespace)
  names(GOtype) = toupper(sub('^[^_]+_(.).*$','\\1',GOtype))#matches GOtab
  GO_gmt = vector(mode="list",length=length(GOtype))
  names(GO_gmt) = paste("GO",GOtags,sep='')
  for(j in 1:length(GOtype)){
    GOdown = as.list(get(paste("GO",GOtags[names(GOtype)[j]],"OFFSPRING",sep='')))
    tmp2 = GOtab[hasOBONamespace==GOtype[j]]
    setnames(tmp2,old=c('id','hasOBONamespace'),new=c('V5','V9'))
    tmp2$V9 = names(GOtype)[j]
    tmp2 = merge(tmp2,tmp,by=c('V5','V9'))
    # group ID columns
    tmp3 = as.data.table(cbind(paste(tmp2$V9,tmp2$V5,tmp2$label,sep=":"),paste(tmp2$V1,tmp2$V2,tmp2$V3,sep="|"),tmp2$V5))
    lfun = function(x){if(length(x)==0){I(list(''))}else{I(list(x))}}
    GOj = unique(tmp3$V1); GOjID = sub('^.:(GO:[0-9]+):.*$','\\1',GOj)
    GO_gmt[[j]] = vector(mode='list',length=length(GOj))
    # stuff list with genes assoc w each GO term and its offspring
    for(i in 1:length(GOj)){
      GOset = c(GOjID[i],GOdown[[ GOjID[i] ]])
      GO_gmt[[j]][i] = lfun( unique(tmp3$V2[tmp3$V3 %in% GOset]) )
    }
    names(GO_gmt[[j]]) = GOj
  }
  if( !is.null(trim_regex) ){
    for(j in 1:length(GOtype)){
      GO_gmt[[j]] = lapply(GO_gmt[[j]],function(x){sub(trim_regex,'\\1',x)})
    }
  }
  return(GO_gmt)
}
