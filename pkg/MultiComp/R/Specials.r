
Specials <- function(f,spec,split2="+",...) {
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
   if (is.null(pos)) return(NULL)
   x <- rownames(attributes(tt)$factors)[pos]
   st <- gsub(" ","",x)
   res <- unlist(strsplit(st,"[()]"))[2]
   if (is.null(split2)) return(res)
   unlist(strsplit(res,"+",fixed=TRUE))
}


