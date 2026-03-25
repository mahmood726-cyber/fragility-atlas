env <- new.env()
rda_files <- list.files("C:/Models/Pairwise70/data", pattern="^CD001431_", full.names=TRUE)
cat("RDA files:", rda_files, "\n")
load(rda_files[1], envir=env)
df <- get(ls(env)[1], envir=env)
analyses <- aggregate(Study ~ Analysis.group + Analysis.number, data=df, FUN=length)
names(analyses)[3] <- "k"
analyses <- analyses[order(-analyses$k),]
print(head(analyses, 5))
cat("\nR selects: group=", analyses[1, "Analysis.group"], " num=", analyses[1, "Analysis.number"], " k=", analyses[1, "k"], "\n")
