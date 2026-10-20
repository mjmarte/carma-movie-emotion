# Runs every script in src/ in numbered order, stopping at the first error. Each script reads config.json.
scripts <- sort(list.files("src", pattern = "\\.(R|py)$", full.names = TRUE))
for (s in scripts) {
  cat("==", s, "\n")
  status <- if (grepl("\\.py$", s)) system2("python3", s) else system2("Rscript", s)
  if (status != 0) stop("failed: ", s)
}
