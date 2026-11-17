# Build documentation from definitions only: never source experiment/tool scripts.
# This keeps roxygen useful without converting the research repository to a package.
source("R/load.R");load_mvgc();require_packages("roxygen2")
args<-commandArgs(TRUE);strict<-"--require-complete"%in%args
run_examples<-"--examples"%in%args
destination<-"docs/reference";dir.create(destination,recursive=TRUE,showWarnings=FALSE)
files<-sort(unique(c(list.files("R",pattern="\\.R$",recursive=TRUE,full.names=TRUE),
 list.files("config",pattern="\\.R$",full.names=TRUE),
 unlist(lapply(c("tools","experiments","figures","tables"),function(d)
  list.files(d,pattern="\\.R$",recursive=TRUE,full.names=TRUE))))))
definitions<-list()
for(path in files) {
 lines<-readLines(path,warn=FALSE);parsed<-parse(path,keep.source=TRUE);locations<-attr(parsed,"srcref")
 for(i in seq_along(parsed)) {
  expr<-parsed[[i]]
  if(!is.call(expr)||(!identical(expr[[1]],as.name("<-"))&&!identical(expr[[1]],as.name("=")))||length(expr)!=3L||
   !is.symbol(expr[[2]])||!is.call(expr[[3]])||!identical(expr[[3]][[1]],as.name("function")))next
  name<-as.character(expr[[2]]);location<-locations[[i]]
  first<-location[1];last<-location[3];before<-first-1L
  while(before>=1L&&grepl("^\\s*(#|$)",lines[before]))before<-before-1L
  preceding<-if(first>before+1L)lines[seq.int(before+1L,first-1L)]else character()
  doc<-preceding[grepl("^\\s*#'",preceding)]
  parameters<-names(expr[[3]][[2]])
  param_tags<-grep("^\\s*#'\\s*@param\\s+",doc,value=TRUE)
  declared<-unlist(strsplit(sub("^\\s*#'\\s*@param\\s+([^[:space:]]+).*","\\1",param_tags),",",fixed=TRUE))
  missing<-setdiff(parameters,declared);extra<-setdiff(declared,parameters)
  definitions[[length(definitions)+1L]]<-list(name=name,source=path,line=first,
   expression=expr,documentation=doc,parameters=parameters,
   missing=missing,extra=extra,has_return=any(grepl("@returns?\\b",doc)),
   has_examples=any(grepl("@examples?(If)?\\b",doc)))
 }
}
names_all<-vapply(definitions,`[[`,"","name")
rows<-list();example_env<-new.env(parent=globalenv());load_mvgc(example_env)
# Function assignments have no workflow side effects. Do not evaluate other AST nodes.
for(d in definitions)if(startsWith(d$source,"R/")||startsWith(d$source,"config/"))eval(d$expression,example_env)
for(d in definitions) {
 documented<-length(d$documentation)>0L
 complete<-documented&&!length(d$missing)&&!length(d$extra)&&d$has_return
 id<-if(sum(names_all==d$name)>1L)paste0(gsub("[^A-Za-z0-9]","_",d$source),"__",d$name)else d$name
 rd_path<-html_path<-"";examples_status<-"not_requested"
 if(documented) {
  input<-paste(c(d$documentation,deparse(d$expression,width.cutoff=500L)),collapse="\n")
  topics<-roxygen2::roc_proc_text(roxygen2::rd_roclet(),input)
  if(length(topics)!=1L)stop("Expected one documentation topic for ",d$source,":",d$name)
  rd_path<-file.path(destination,paste0(id,".Rd"));html_path<-file.path(destination,paste0(id,".html"))
  text<-format(topics[[1]])
  text<-sub("% Please edit documentation in [^\r\n]*",paste0("% Please edit documentation in ",d$source),text)
  rd_lines<-strsplit(text,"\n",fixed=TRUE)[[1]]
  if(any(grepl("(^|[^\\\\])%",rd_lines)&!grepl("^[[:space:]]*%",rd_lines)))
    stop("Unescaped percent would truncate generated Rd for ",d$source,":",d$name,call.=FALSE)
  writeLines(text,rd_path,useBytes=TRUE)
  # Rd comments begin with a literal percent sign. A parse warning can silently
  # truncate a contract or its braces; neither a generated file nor a green
  # parameter count is sufficient evidence that the help page is readable.
  rd<-withCallingHandlers(tools::parse_Rd(rd_path),warning=function(w)
    stop("Malformed generated Rd for ",d$source,":",d$name,": ",conditionMessage(w),call.=FALSE))
  withCallingHandlers(tools::checkRd(rd),warning=function(w)
    stop("Invalid generated Rd for ",d$source,":",d$name,": ",conditionMessage(w),call.=FALSE))
  tools::Rd2HTML(rd,out=html_path,texmath="ascii")
  if(run_examples&&d$has_examples) {
   script<-tempfile(fileext=".R");tools::Rd2ex(rd,out=script)
   # Examples execute in a separate environment, never in the live data job.
   env<-new.env(parent=example_env)
   sys.source(script,env);examples_status<-"PASS"
  }else if(!d$has_examples)examples_status<-"missing"
 }
 rows[[length(rows)+1L]]<-data.frame(function_name=d$name,source=d$source,line=d$line,
  documented=documented,complete=complete,missing_params=paste(d$missing,collapse=","),
  extra_params=paste(d$extra,collapse=","),return_documented=d$has_return,
  examples=examples_status,rd=rd_path,html=html_path,stringsAsFactors=FALSE)
}
coverage<-do.call(rbind,rows)
index <- coverage[c("function_name","source","documented","html")]
write.csv(index,file.path(destination,"topics.csv"),row.names=FALSE,eol="\n")
private <- Sys.getenv("MVGC_PRIVATE_VERIFICATION_DIR", "")
if (nzchar(private)) { dir.create(private,recursive=TRUE,showWarnings=FALSE); write.csv(coverage,file.path(private,"roxygen_coverage.csv"),row.names=FALSE) }
#' Escape text for the local function-reference HTML index
#' @param x Character text for an HTML text node, not a URL or attribute value.
#' @return Character text with ampersand and angle brackets escaped. Quotes are
#'   deliberately unchanged because the result is inserted only as visible text.
escape<-function(x){x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)}
links<-vapply(which(coverage$documented),function(i)paste0("<li><a href='",basename(coverage$html[i]),"'>",
 escape(coverage$function_name[i]),"</a> â€” ",escape(coverage$source[i]),"</li>"),"")
writeLines(c("<!doctype html><html lang='en'><meta charset='utf-8'><title>Article 1 function reference</title>",
 "<h1>Article 1 function reference</h1>",
 "<ul>",links,"</ul></html>"),file.path(destination,"index.html"))
cat(sum(coverage$documented),"documented;",sum(coverage$complete),"complete;",nrow(coverage),"functions inventoried\n")
if(strict&&any(!coverage$complete))stop("Function documentation remains incomplete; see private documentation verification output")

# Keep generated public artifacts byte-stable across Windows and Unix checkouts.
for (path in list.files(destination, pattern="\\.(Rd|html|csv)$", full.names=TRUE)) {
  lines <- readLines(path, warn=FALSE, encoding="UTF-8")
  writeBin(charToRaw(paste0(paste(lines,collapse="\n"),"\n")),path)
}
