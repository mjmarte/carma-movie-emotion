# Reads the local configuration (config.json at the repository root, not tracked) that maps each input to its file on this
# machine. Copy config.example.json to config.json and set the paths. Scripts refer to inputs only by key.
suppressPackageStartupMessages(library(jsonlite))
.cfg_file <- Sys.getenv("CARMA_CONFIG", "config.json")
if (!file.exists(.cfg_file)) stop("config.json not found: copy config.example.json to config.json and set the paths")
CFG <- jsonlite::read_json(.cfg_file)
input_file <- function(key) {
  p <- CFG$inputs[[key]]
  if (is.null(p) || !file.exists(p)) stop("input '", key, "' is not set in config.json or the file does not exist")
  p
}
output_path <- function(...) {
  p <- file.path(CFG$output_dir, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}
# Participant-level rules kept in config.json, so no participant identifiers appear in the code: groups and participants
# outside the analytic sample, and recorded IDs to correct.
RULES <- CFG$participant_rules
