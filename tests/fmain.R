source("hsa/cstruct1.R")

outpath <- file.path("tests", ".checks")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
outprefix <- file.path(outpath, "hap1_full_all_22_100000")
cat(sprintf("Prefix of output files: %s\n", sQuote(outprefix)))

# Input data
path <- file.path("inst", "exdata")
contacts <- file.path(path, "contact_hap1_full_all_22_100000.txt")

## Parameters
mak <- 1L
Iscovfile <- FALSE
K <- 1L
lscov0 <- 0L

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
dt <- system.time({
  out <- fmain(lsmap0 = lsmap0, lscov0 = lscov0, outfile = outprefix,
               Maxiter = 3L, submaxiter = 100L, lambda = 25, Leapfrog = 20,
               epslon = 0.003, mkfix = 0, rho = 0, mk = mak)
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

cat("\n\n")
print(sessionInfo())

cat("\nProcessing time:\n")
print(dt)

cat("\nValidation:\n")
data0_1 <- read.table(file.path(path, basename(file1)), header = FALSE)
data1_1 <- read.table(file1, header = FALSE)
stopifnot(all.equal(data1_1, data0_1))

data0_2 <- read.table(file.path(path, basename(file2)), header = FALSE)
data1_2 <- read.table(file2, header = FALSE)
stopifnot(all.equal(data1_2, data0_2))
