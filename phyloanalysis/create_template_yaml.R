require("yaml")
create_template_yaml=function(cfg,patient,sex,yaml_output){
  dat=read.table(cfg,head=T,stringsAsFactors = FALSE,colClasses=c("character"),comment.char = "#")
  tmp=dat[which(dat$TYPE=="C"),]
  dat=rbind(dat[which(dat$TYPE!="C"),],tmp[order(tmp$LABEL),])
  
  header=readLines(con = cfg)
  header=gsub("#","",header[grep("^#",header)]) ## could just do gsub("=",": ",header)
  keyvalues=lapply(header,function(x) unlist(strsplit(x,"=")))
  header=sapply(keyvalues, function(element){
    if(length(element)!=2){
      stop("Error in header...")
    }
    sprintf("%s: %s",trimws(element[1]),trimws(element[2]))
  })
  
  body=unlist(strsplit(as.yaml(dat,column.major = FALSE),split="\n"))
  body=sprintf("  %s",body)
  #browser()
  header=c(sprintf(c("PATIENT: %s","SEX: %s","PDID: %s","COMMENT: %s"),c(patient,sex,patient,"None")),c(header,"SAMPLES:"))
  cat(paste(header,collapse="\n"))
  yaml_output=gsub("cfg$","yaml",cfg)
  cat("writing to template YAML file: ",yaml_output,"\n")
  if(file.exists(yaml_output)){
    stop("Output yaml file exists")
  }
  writeLines(c(header,body),con=yaml_output)
  yaml_output
}


args = commandArgs(trailingOnly=TRUE)
if(length(args)<3){
  stop("Rscript create_template_yaml <cfg> <patient id> <sex>")
}
create_template_yaml(args[1],args[2],args[3])


