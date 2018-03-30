source("hsa/cstruct1.R")

contacts <- file.path("inst", "exdata", "contact_hap1_full_all_22_100000.txt")
outprefix <- tempfile(pattern = "cstruct1_fmain_")
message("Prefix of output files: ", sQuote(outprefix))

## Parameters
mak <- 1L
Iscovfile <- FALSE
K <- 1L
lscov0 <- 0L

data <- read.table(file = contacts, header = FALSE)
data <- as.matrix(data)
lsmap0 <- list(data)

set.seed(12345)

## Benchmarking history:
## commit 0880732: ~ 270 secs = 4.5 mins
dt <- system.time({
  out <- fmain(lsmap0 = lsmap0, lscov0 = lscov0, outfile = outprefix,
               Maxiter = 3L, submaxiter = 100L, lambda = 25, Leapfrog = 20,
               epslon = 0.003, mkfix = 0, rho = 0, mk = mak)
})
print(dt)
saveRDS(out, file = sprintf("%s.out.rds", outfile))
