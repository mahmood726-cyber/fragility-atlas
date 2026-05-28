resolve_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE))
  }
  normalizePath("..", winslash = "/", mustWork = TRUE)
}

first_existing_path <- function(paths) {
  for (path in paths) {
    if (!is.na(path) && nzchar(path) && file.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

PROJECT_ROOT <- resolve_project_root()
pairwise_dir <- first_existing_path(c(
  Sys.getenv("PAIRWISE70_DATA_DIR", unset = ""),
  file.path(dirname(PROJECT_ROOT), "Projects", "mahmood789", "Pairwise70", "data"),
  file.path(dirname(PROJECT_ROOT), "Projects", "Pairwise70", "data"),
  file.path(dirname(PROJECT_ROOT), "Models", "Pairwise70", "data"),
  file.path(dirname(PROJECT_ROOT), "Pairwise70", "data")
))
if (is.null(pairwise_dir)) {
  stop("Pairwise70 data directory not found. Set PAIRWISE70_DATA_DIR to inspect CD001431.")
}

env <- new.env()
rda_files <- list.files(pairwise_dir, pattern="^CD001431_", full.names=TRUE)
cat("RDA files:", rda_files, "\n")
load(rda_files[1], envir=env)
df <- get(ls(env)[1], envir=env)
analyses <- aggregate(Study ~ Analysis.group + Analysis.number, data=df, FUN=length)
names(analyses)[3] <- "k"
analyses <- analyses[order(-analyses$k),]
print(head(analyses, 5))
cat("\nR selects: group=", analyses[1, "Analysis.group"], " num=", analyses[1, "Analysis.number"], " k=", analyses[1, "k"], "\n")
