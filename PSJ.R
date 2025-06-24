library(ggplot2)
library(animint2)
library(PeakSegJoint)
library(data.table)

# Load data
data_path <- "PSJ.RData"
if (!file.exists(data_path)) stop("PSJ.RData not found at ", data_path)
load(data_path)

# Check required objects
required_objects <- c("coverage", "regions.by.problem", "modelSelection.by.problem", 
                     "problems", "problem.labels", "filled.regions", "error.total.chunk", 
                     "error.total.all", "peaks.by.problem")
missing_objects <- setdiff(required_objects, names(PSJ))
if (length(missing_objects) > 0) stop("Missing required objects in PSJ.RData: ", paste(missing_objects, collapse = ", "))

# Validate data frames
validate_data <- function(df, name, required_cols = NULL) {
  if (is.null(df) || nrow(df) == 0) stop(name, " is empty")
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) stop("Missing columns in ", name, ": ", paste(missing_cols, collapse = ", "))
  }
}
validate_data(PSJ$coverage, "PSJ$coverage", c("base", "count", "sample.id"))
validate_data(PSJ$problems, "PSJ$problems", c("problemStart", "problemEnd", "problem.i", "problem.name"))
validate_data(do.call(rbind, PSJ$regions.by.problem), "PSJ$regions.by.problem", c("chromStart", "chromEnd", "problem.name", "bases.per.problem"))
validate_data(PSJ$error.total.chunk, "PSJ$error.total.chunk", c("bases.per.problem", "errors", "regions"))
validate_data(PSJ$filled.regions, "PSJ$filled.regions", c("chromStart", "chromEnd", "annotation"))
validate_data(do.call(rbind, PSJ$peaks.by.problem), "PSJ$peaks.by.problem", c("chromStart", "chromEnd", "peaks", "problem.name"))

# Prepare data
res.error <- PSJ$error.total.chunk
ann.colors <- c(noPeaks="#f6f4bf", peakStart="#ffafaf", peakEnd="#ff4c4c", peaks="#a445ee")

# Problem regions
all.regions <- do.call(rbind, PSJ$regions.by.problem)
prob.regions <- unique(data.frame(all.regions)[, c("bases.per.problem", "problem.i", "problem.name", "chromStart", "chromEnd")])
prob.regions$sample.id <- "problems"

# Model selection
all.modelSelection <- do.call(rbind, PSJ$modelSelection.by.problem)
modelSelection.errors <- all.modelSelection[!is.na(all.modelSelection$errors), ]
penalty.range <- with(all.modelSelection, c(min(max.log.lambda), max(min.log.lambda)))
penalty.mid <- mean(penalty.range)

modelSelection.labels <- unique(with(all.modelSelection, {
  data.frame(problem.name, bases.per.problem, problemStart, problemEnd,
             min.log.lambda = penalty.mid, peaks = max(peaks) + 0.5)
}))

# Sample peaks
sample.peaks <- do.call(rbind, PSJ$peaks.by.problem)
problem.peaks <- unique(data.frame(sample.peaks)[, c("bases.per.problem", "problem.i", "problem.name", "peaks", "chromStart", "chromEnd")])
problem.peaks$sample.id <- "problems"

# Combine problem line data into coverage
problems.coverage <- data.frame(
  base = prob.regions$chromStart,
  count = NA,
  sample.id = "problems"
)

# Combine with original coverage
combined.coverage <- rbind(
  PSJ$coverage[, c("base", "count", "sample.id")],
  problems.coverage
)

# Ensure 'problems' is at top
sample.levels <- c("problems", setdiff(unique(as.character(PSJ$coverage$sample.id)), "problems"))
combined.coverage$sample.id <- factor(combined.coverage$sample.id, levels = sample.levels)
prob.regions$sample.id <- factor(prob.regions$sample.id, levels = sample.levels)
sample.peaks$sample.id <- factor(sample.peaks$sample.id, levels = sample.levels)
problem.peaks$sample.id <- factor(problem.peaks$sample.id, levels = sample.levels)

PSJ$coverage <- combined.coverage

# Error resolution window
dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval - 1) / 2 + 1
res.error$min.bases.per.problem <- res.error$bases.per.problem / dval2
res.error$max.bases.per.problem <- res.error$bases.per.problem * dval2

# X-axis limits
xlim <- range(c(PSJ$coverage$base, PSJ$problems$problemStart, PSJ$problems$problemEnd,
              all.regions$chromStart, all.regions$chromEnd, sample.peaks$chromStart, sample.peaks$chromEnd),
              na.rm = TRUE) / 1e3

# First selector
default.problem <- prob.regions$problem.name[1]
default.bases <- prob.regions$bases.per.problem[1]
PSJ$first <- list(problem.name = default.problem, bases.per.problem = default.bases)

# Function to create peak variable names
peakvar <- function(position) {
  paste0(gsub("[-:]", ".", position), "peaks")
}

# Build the visualization with proper parameter syntax

viz <- list(
  coverage = ggplot() +
    geom_tallrect(
      aes(xmin = chromStart / 1e3, xmax = chromEnd / 1e3, fill = annotation),
      data = PSJ$filled.regions,
      alpha = 0.5,
      color = "grey"
    ) +
    geom_line(
      aes(base / 1e3, count),
      data = PSJ$coverage,
      color = "grey50"
    ) +
    geom_segment(
      aes(chromStart / 1e3, 1, xend = chromEnd / 1e3, yend = 1),
      data = transform(prob.regions, problem.i = 1),
      size = 1,
      color = "black",
      showSelected = "bases.per.problem",
      clickSelects = "problem.name"
    ) +
    geom_segment(
      aes(problemStart / 1e3, 1, xend = problemEnd / 1e3, yend = 1),
      data = transform(PSJ$problems, problem.i = 1),
      size = 5,
      color = "black",
      showSelected = "bases.per.problem",
      clickSelects = "problem.name"
    ) +
    geom_text(
      aes(chromStart / 1e3, 1,
          label = sprintf("%d problems mean size %.1f kb", 
                         nrow(PSJ$problems), 
                         mean(PSJ$problems$problemEnd - PSJ$problems$problemStart) / 1e3)),
      data = transform(prob.regions, problem.i = 1),
      hjust = 0,
      vjust = -0.5,
      showSelected = "bases.per.problem"
    ) +
    geom_tallrect(
      aes(xmin = chromStart / 1e3, xmax = chromEnd / 1e3, linetype = status),
      data = all.regions,
      fill = NA,
      color = "black",
      showSelected = c("problem.name", "peaks"),  # Fixed this line
      showSelected2 = "bases.per.problem"
    ) +
    geom_segment(
      aes(chromStart / 1e3, 0, xend = chromEnd / 1e3, yend = 0),
      data = sample.peaks,
      size = 7,
      color = "deepskyblue",
      clickSelects = "problem.name",
      showSelected = c("problem.name", "peaks"),  
      showSelected2 = "bases.per.problem"
    ) +
    geom_segment(
      aes(chromStart / 1e3, problem.i, xend = chromEnd / 1e3, yend = problem.i),
      data = problem.peaks,
      size = 7,
      color = "deepskyblue",
      clickSelects = "problem.name",
      showSelected = c("problem.name", "peaks"), 
      showSelected2 = "bases.per.problem"
    ) +
      scale_y_continuous("aligned read coverage", breaks = function(limits) floor(limits[2]))+
    scale_linetype_manual("error type",
                         limits = c("correct", "false negative", "false positive"),
                         values = c(correct = 0, "false negative" = 3, "false positive" = 1)) +
    scale_x_continuous("position on chr11 (kilo bases = kb)") +
    scale_fill_manual(values = ann.colors) +
    coord_cartesian(xlim = c(118167.406, 118238.833)) +
    theme_bw() +
    theme_animint(width = 1500, height = length(sample.levels) * 100) +
    theme(panel.margin = grid::unit(0, "cm")) +
    facet_grid(sample.id ~ ., scales = "free"),

  resError = ggplot() +
    ggtitle("select problem size") +
    ylab("minimum percent incorrect regions") +
    geom_tallrect(
      aes(xmin = min.bases.per.problem, xmax = max.bases.per.problem),
      data = res.error,
      alpha = 0.5,
      clickSelects = "bases.per.problem"
    ) +
    scale_x_log10() +
    geom_line(
      aes(bases.per.problem, errors / regions * 100, color = chunks, size = chunks),
      data = data.frame(res.error, chunks = "this")
    ) +
    geom_line(
      aes(bases.per.problem, errors / regions * 100, color = chunks, size = chunks),
      data = data.frame(PSJ$error.total.all, chunks = "all")
    ),

  modelSelection = ggplot() +
    geom_segment(
      aes(min.log.lambda, peaks, xend = max.log.lambda, yend = peaks),
      data = data.frame(all.modelSelection, what = "peaks"),
      size = 5,
      showSelected = "problem.name",
      showSelected2 = "bases.per.problem"
    ) +
    geom_text(
      aes(min.log.lambda, peaks,
          label = sprintf("%.1f kb in problem %s", 
                         (problemEnd - problemStart) / 1e3, 
                         problem.name)),
      data = data.frame(modelSelection.labels, what = "peaks"),
      showSelected = "problem.name",
      showSelected2 = "bases.per.problem"
    ) +
    geom_segment(
      aes(min.log.lambda, as.integer(errors), xend = max.log.lambda, yend = as.integer(errors)),
      data = data.frame(modelSelection.errors, what = "errors"),
      size = 5,
      showSelected = "problem.name",
      showSelected2 = "bases.per.problem"
    ) +
    ggtitle("select number of samples with 1 peak") +
    ylab("") +
    geom_tallrect(
      aes(xmin = min.log.lambda, xmax = max.log.lambda),
      data = all.modelSelection,
      alpha = 0.5,
      clickSelects = "peaks",
      showSelected = "problem.name",
      showSelected2 = "bases.per.problem"
    ) +
    facet_grid(what ~ ., scales = "free"),

  title = "PeakSegJoint Interactive Visualization",
  first = PSJ$first,
  selector.types = list(
    problem.name = "single",
    bases.per.problem = "single"
  ),
  duration = list(
    problem.name = 1000,
    bases.per.problem = 1000
  )
)

# Compile the visualization
animint2dir(viz, out.dir = "PSJ_interactive_fixed")
