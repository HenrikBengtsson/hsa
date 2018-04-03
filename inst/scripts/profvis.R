library("profvis")
`[.profvis` <- function(x, i, ...) { prof <- x$x$message$prof; prof <- prof[i, ...]; x$x$message$prof <- prof; x }

source("hsa/cstruct1.R")

outpath <- file.path("tests", ".checks", "profvis")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
outprefix <- file.path(outpath, "hap1_full_all_22_100000")
cat(sprintf("Prefix of output files: %s\n", sQuote(outprefix)))

# Input data
contacts <- file.path("inst", "exdata", "contact_hap1_full_all_22_100000.txt")

## Parameters
data <- read.table(file = contacts, header = FALSE)
data <- as.matrix(data)
lsmap0 <- list(data)

set.seed(12345)

options(digits = 7L, scipen = 0L) ## Reproducible print() output

## AD HOC: Trick cstruct1.R code to write files with 12 digits
## (still plenty) instead of 15 digits for easier 'diff' comparisons
write.table <- function(x, ...) {
  for (kk in seq_along(x)) {
    if (is.double(x[[kk]])) x[[kk]] <- round(x[[kk]], digits = 12L)
  }
  utils::write.table(x, ...)
}

## Benchmarking history:
## commit 0880732: ~270 secs
## commit ccc77de: ~265 secs
pv <- profvis::profvis({
  out <- fmain(lsmap0 = lsmap0, lscov0 = 0L, outfile = outprefix,
               Maxiter = 1L, submaxiter = 100L, lambda = 25, Leapfrog = 20,
               epslon = 0.003, mkfix = 0, rho = 0, mk = 1L)
})
print(out)
cat("\n\n")

cat("Files produced:\n")
file1 <- sprintf("%s.txt", outprefix)
file2 <- sprintf("%sbeta.txt", outprefix)
cat(sprintf("- file 1: %s (%d bytes; %s)\n",
            sQuote(basename(file1)), file.size(file1), tools::md5sum(file1)))
cat(sprintf("- file 2: %s (%d bytes; %s)\n",
            sQuote(basename(file2)), file.size(file2), tools::md5sum(file2)))
cat("\n\n")

file3 <- sprintf("%s.profvis.rds", outprefix)
saveRDS(pv, file = file3)

cat("\n\n")
print(sessionInfo())

## Render profvis results
print(head(pv, n = 100e3))



