# Regenerate the R reference outputs used by tests/test_r_reference.py.
# Run inside the dev container from this folder:  Rscript generate.R
suppressMessages({library(vegan); library(RLdbRDA)})
num <- function(out) {
  for (c in c("Fa", "r2", "r2adj", "N", "pval")) out[[c]] <- as.numeric(out[[c]])
  out
}
save <- function(out, file) write.table(num(out), file, sep = "\t", quote = FALSE, col.names = NA)

spec <- read.delim("../../../data/varespec.tsv", row.names = 1)
chem <- read.delim("../../../data/varechem.tsv", row.names = 1)
for (s in c(1, 42)) { set.seed(s); save(rldbrda(spec, chem), sprintf("varespec_seed%d.tsv", s)) }

for (case in c("mixed", "collinear")) {
  df <- read.delim(sprintf("%s_data.tsv", case), row.names = 1)
  meta <- read.delim(sprintf("%s_meta.tsv", case), row.names = 1)
  for (method in c("bray", "jaccard", "euclidean")) {
    set.seed(7)
    out <- tryCatch(suppressMessages(rldbrda(df, meta, method = method, p_cutoff = 0.5)),
                    error = function(e) NULL)
    while ("meta" %in% search()) detach("meta")  # RLdbRDA leaves it attached on error
    if (is.null(out)) { message("R failed for ", case, "/", method); next }
    save(out, sprintf("%s_%s_seed7.tsv", case, method))
  }
}

x <- as.matrix(read.delim("small_counts.tsv", row.names = 1))
methods <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower",
             "morisita", "horn", "mountford", "jaccard", "raup", "binomial", "chao",
             "altGower", "cao", "mahalanobis", "clark", "chisq", "chord", "hellinger")
d <- sapply(methods, function(m) as.vector(suppressWarnings(vegdist(x, m))))
write.table(signif(d, 17), "small_vegdist.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
