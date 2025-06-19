library(ggplot2)
library(animint2)
library(PeakSegJoint)
library(data.table)

# Load data
data_path <- "\\PSJ.RData"
if (!file.exists(data_path)) {
  stop("PSJ.RData not found at ", data_path)
}
load(data_path)

# Check required objects
required_objects <- c("coverage", "regions.by.problem", "modelSelection.by.problem", 
                      "problems", "problem.labels", "filled.regions", "error.total.chunk", 
                      "error.total.all", "peaks.by.problem")
missing_objects <- setdiff(required_objects, names(PSJ))
if (length(missing_objects) > 0) {
  stop("Missing required objects in PSJ.RData: ", paste(missing_objects, collapse = ", "))
}

# Validate key data frames
validate_data <- function(df, name, required_cols = NULL) {
  if (is.null(df) || nrow(df) == 0) stop(name, " is empty")
  cat(name, " rows:", nrow(df), "\n")
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop("Missing columns in ", name, ": ", paste(missing_cols, collapse = ", "))
    }
  }
}
validate_data(PSJ$coverage, "PSJ$coverage", c("base", "count", "sample.id"))
validate_data(PSJ$problems, "PSJ$problems", c("problemStart", "problemEnd", "problem.i", "problem.name"))
validate_data(do.call(rbind, PSJ$regions.by.problem), "PSJ$regions.by.problem", 
              c("chromStart", "chromEnd", "problem.name", "bases.per.problem"))
validate_data(PSJ$error.total.chunk, "PSJ$error.total.chunk", c("bases.per.problem", "errors", "regions"))
validate_data(PSJ$filled.regions, "PSJ$filled.regions", c("chromStart", "chromEnd", "annotation"))
validate_data(do.call(rbind, PSJ$peaks.by.problem), "PSJ$peaks.by.problem", 
              c("chromStart", "chromEnd", "peaks", "problem.name"))

# Inspect problem.name consistency
cat("Unique problem.name in prob.regions:", length(unique(PSJ$problems$problem.name)), "\n")
cat("Unique problem.name in regions:", length(unique(do.call(rbind, PSJ$regions.by.problem)$problem.name)), "\n")

# Prepare data
res.error <- PSJ$error.total.chunk
ann.colors <- c(
  noPeaks = "#f6f4bf",
  peakStart = "#ffafaf",
  peakEnd = "#ff4c4c",
  peaks = "#a445ee"
)

# Problem regions
all.regions <- do.call(rbind, PSJ$regions.by.problem)
prob.regions.names <- c("bases.per.problem", "problem.i", "problem.name", "chromStart", "chromEnd")
prob.regions <- unique(data.frame(all.regions)[, prob.regions.names])
prob.regions$sample.id <- "problems"
validate_data(prob.regions, "prob.regions", prob.regions.names)

# Model selection
all.modelSelection <- do.call(rbind, PSJ$modelSelection.by.problem)
validate_data(all.modelSelection, "all.modelSelection", 
              c("min.log.lambda", "max.log.lambda", "peaks", "problem.name"))
modelSelection.errors <- all.modelSelection[!is.na(all.modelSelection$errors), ]
penalty.range <- with(all.modelSelection, c(min(max.log.lambda), max(min.log.lambda)))
penalty.mid <- mean(penalty.range)
cat("penalty.range:", penalty.range, "\n")

# Coverage counts
coverage.counts <- table(PSJ$coverage$sample.id)
facet.rows <- length(coverage.counts) + 1
cat("facet.rows:", facet.rows, "\n")

# Error resolution
dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval - 1) / 2 + 1
res.error$min.bases.per.problem <- res.error$bases.per.problem / dval2
res.error$max.bases.per.problem <- res.error$bases.per.problem * dval2

# Model selection labels
modelSelection.labels <- unique(with(all.modelSelection, {
  data.frame(
    problem.name = problem.name,
    bases.per.problem = bases.per.problem,
    problemStart = problemStart,
    problemEnd = problemEnd,
    min.log.lambda = penalty.mid,
    peaks = max(peaks) + 0.5
  )
}))
validate_data(modelSelection.labels, "modelSelection.labels", c("problem.name", "min.log.lambda", "peaks"))

# Sample peaks
sample.peaks <- do.call(rbind, PSJ$peaks.by.problem)
prob.peaks.names <- c("bases.per.problem", "problem.i", "problem.name", "peaks", "chromStart", "chromEnd")
problem.peaks <- unique(data.frame(sample.peaks)[, prob.peaks.names])
problem.peaks$sample.id <- "problems"
validate_data(sample.peaks, "sample.peaks", prob.peaks.names)

# Dynamic x-axis limits
xlim <- range(c(PSJ$coverage$base, PSJ$problems$problemStart, PSJ$problems$problemEnd, 
                all.regions$chromStart, all.regions$chromEnd, sample.peaks$chromStart, 
                sample.peaks$chromEnd), na.rm = TRUE) / 1e3
cat("X-axis limits (kb):", xlim, "\n")

# Validate or set default first selector
if (!exists("first", where = PSJ) || !is.list(PSJ$first) || 
    is.null(PSJ$first$problem.name) || is.null(PSJ$first$bases.per.problem)) {
  cat("PSJ$first is missing or incomplete, setting default\n")
  valid_problem <- as.character(prob.regions$problem.name[1])
  valid_bases <- as.numeric(prob.regions$bases.per.problem[1])
  if (is.na(valid_problem) || is.na(valid_bases)) {
    stop("Cannot set default PSJ$first: no valid problem.name or bases.per.problem in prob.regions")
  }
  PSJ$first <- list(
    problem.name = valid_problem,
    bases.per.problem = valid_bases
  )
}
cat("PSJ$first selectors:\n")
print(PSJ$first)
if (length(PSJ$first$problem.name) == 0 || !PSJ$first$problem.name %in% prob.regions$problem.name) {
  stop("PSJ$first$problem.name (", PSJ$first$problem.name, ") is invalid or not found in prob.regions$problem.name")
}
if (length(PSJ$first$bases.per.problem) == 0 || !PSJ$first$bases.per.problem %in% prob.regions$bases.per.problem) {
  stop("PSJ$first$bases.per.problem (", PSJ$first$bases.per.problem, ") is invalid or not found in prob.regions$bases.per.problem")
}

# Construct visualization
cat("Constructing data viz with .variable .value\n")
print(system.time({
  viz <- list(
    coverage = ggplot() +
      geom_line(aes(base / 1e3, count), data = PSJ$coverage, color = "grey50") +
      geom_tallrect(aes(xmin = chromStart / 1e3, xmax = chromEnd / 1e3, fill = annotation),
                    alpha = 0.5, color = "grey", data = PSJ$filled.regions) +
      geom_segment(aes(chromStart / 1e3, problem.i, xend = chromEnd / 1e3, yend = problem.i),
                   showSelected = "bases.per.problem",
                   clickSelects = "problem.name",
                   data = prob.regions) +
      geom_segment(aes(problemStart / 1e3, problem.i, xend = problemEnd / 1e3, yend = problem.i),
                   showSelected = "bases.per.problem",
                   clickSelects = "problem.name",
                   size = 5, data = PSJ$problems) +
      geom_text(aes(chromStart / 1e3, problem.i,
                    label = sprintf("%d problems mean size %.1f kb", problems, mean.bases / 1e3)),
                showSelected = "bases.per.problem",
                data = PSJ$problem.labels, hjust = 0) +
      geom_tallrect(aes(xmin = chromStart / 1e3, xmax = chromEnd / 1e3, linetype = status),
                    showSelected = "problem.name",
                    showSelected2 = "bases.per.problem",
                    data = all.regions, fill = NA, color = "black") +
      geom_segment(aes(chromStart / 1e3, 0, xend = chromEnd / 1e3, yend = 0),
                   clickSelects = "problem.name",
                   showSelected = "problem.name",
                   showSelected2 = "bases.per.problem",
                   data = sample.peaks, size = 7, color = "deepskyblue") +
      geom_segment(aes(chromStart / 1e3, problem.i, xend = chromEnd / 1e3, yend = problem.i),
                   clickSelects = "problem.name",
                   showSelected = "problem.name",
                   showSelected2 = "bases.per.problem",
                   data = problem.peaks, size = 7, color = "deepskyblue") +
      scale_y_continuous("aligned read coverage", breaks = function(limits) floor(limits[2])) +
      scale_linetype_manual("error type",
                            limits = c("correct", "false negative", "false positive"),
                            values = c(correct = 0, "false negative" = 3, "false positive" = 1)) +
      scale_x_continuous("position on chr11 (kilo bases = kb)") +
      scale_fill_manual(values = ann.colors) +
      coord_cartesian(xlim = xlim) +
      theme_bw() +
      theme_animint(width = 1500, height = facet.rows * 100) +
      theme(panel.margin = unit(0, "cm")) +
      facet_grid(sample.id ~ ., labeller = function(df) {
        df$sample.id <- sub("McGill0", "", sub(" ", "\n", df$sample.id))
        df
      }, scales = "free"),
    resError = ggplot() +
      geom_tallrect(aes(xmin = min.bases.per.problem, xmax = max.bases.per.problem),
                    clickSelects = "bases.per.problem",
                    alpha = 0.5, data = res.error) +
      geom_line(aes(bases.per.problem, errors / regions * 100, color = chunks, size = chunks),
                data = data.frame(res.error, chunks = "this")) +
      geom_line(aes(bases.per.problem, errors / regions * 100, color = chunks, size = chunks),
                data = data.frame(PSJ$error.total.all, chunks = "all")) +
      scale_x_log10() +
      ggtitle("select problem size") +
      ylab("minimum percent incorrect regions") +
      theme_bw(),
    modelSelection = ggplot() +
      geom_segment(aes(min.log.lambda, peaks, xend = max.log.lambda, yend = peaks),
                   showSelected = "problem.name",
                   showSelected2 = "bases.per.problem",
                   data = data.frame(all.modelSelection, what = "peaks"), size = 5) +
      geom_text(aes(min.log.lambda, peaks,
                    label = sprintf("%.1f kb in problem %s", (problemEnd - problemStart) / 1e3, problem.name)),
                showSelected = "problem.name",
                showSelected2 = "bases.per.problem",
                data = data.frame(modelSelection.labels, what = "peaks")) +
      geom_segment(aes(min.log.lambda, as.integer(errors), xend = max.log.lambda, yend = as.integer(errors)),
                   showSelected = "problem.name",
                   showSelected2 = "bases.per.problem",
                   data = data.frame(modelSelection.errors, what = "errors"), size = 5) +
      geom_tallrect(aes(xmin = min.log.lambda, xmax = max.log.lambda),
                    clickSelects = "problem.name",
                    showSelected = "problem.name",
                    showSelected2 = "bases.per.problem",
                    data = all.modelSelection, alpha = 0.5) +
      ggtitle("select number of samples with 1 peak") +
      ylab("") +
      facet_grid(what ~ ., scales = "free") +
      theme_bw(),
    title = "Animint compiler with .variable .value aesthetics",
    first = PSJ$first
  )
}))

# Compile and export visualization
out_dir <- "\\PSJ"
cat("Compiling data viz to", out_dir, "\n")
print(system.time({
  animint2pages(viz, out_dir)
}))

# Post-compilation checks
tsv_files <- list.files(out_dir, pattern = "*.tsv", full.names = TRUE)
cat("Number of TSV files generated:", length(tsv_files), "\n")
if (length(tsv_files) > 0) {
  cat("Sample TSV file sizes:\n")
  for (f in head(tsv_files, 3)) {
    cat(basename(f), ": ", file.size(f), " bytes\n")
  }
}