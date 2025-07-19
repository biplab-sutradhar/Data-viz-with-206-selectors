library(ggplot2)
library(animint2)
library(PeakSegJoint)
library(data.table)
library(grid)

# Assuming works_with_R is a custom function; commented out as it may not be available
# works_with_R("3.2.2", ggplot2="2.1.0", PeakSegJoint="2016.3.2", "tdhock/animint@03735869af84629d269556442345b2ea506ab42a")

load("PSJ.RData")

res.error <- PSJ$error.total.chunk

ann.colors <- c(noPeaks="#f6f4bf", peakStart="#ffafaf", peakEnd="#ff4c4c", peaks="#a445ee")

# Data preparation remains identical
all.regions <- do.call(rbind, PSJ$regions.by.problem)
prob.regions.names <- c("bases.per.problem", "problem.i", "problem.name", "chromStart", "chromEnd")
prob.regions <- unique(data.frame(all.regions)[, prob.regions.names])
prob.regions$sample.id <- "problems"

all.modelSelection <- do.call(rbind, PSJ$modelSelection.by.problem)
modelSelection.errors <- all.modelSelection[!is.na(all.modelSelection$errors), ]
penalty.range <- with(all.modelSelection, c(min(max.log.lambda), max(min.log.lambda)))
penalty.mid <- mean(penalty.range)

coverage.counts <- table(PSJ$coverage$sample.id)
facet.rows <- length(coverage.counts) + 1
dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval - 1) / 2 + 1
res.error$min.bases.per.problem <- res.error$bases.per.problem / dval2
res.error$max.bases.per.problem <- res.error$bases.per.problem * dval2

modelSelection.labels <- unique(with(all.modelSelection, {
  data.frame(problem.name = problem.name, bases.per.problem = bases.per.problem,
             problemStart = problemStart, problemEnd = problemEnd,
             min.log.lambda = penalty.mid, peaks = max(peaks) + 0.05)
}))

# Set sample.id factor levels to ensure "problems" is at the top
sample.levels <- c("problems", setdiff(unique(as.character(PSJ$coverage$sample.id)), "problems"))
PSJ$coverage$sample.id <- factor(PSJ$coverage$sample.id, levels = sample.levels)
prob.regions$sample.id <- factor(prob.regions$sample.id, levels = sample.levels)
PSJ$filled.regions$sample.id <- factor(PSJ$filled.regions$sample.id, levels = sample.levels)
PSJ$problems$sample.id <- factor(PSJ$problems$sample.id, levels = sample.levels)
PSJ$problem.labels$sample.id <- factor(PSJ$problem.labels$sample.id, levels = sample.levels)

# Visualization with for loops
cat("constructing data viz with for loops\n")
print(system.time({
  viz.for <- list(
    coverage = ggplot() +
      geom_segment(aes(chromStart/1e3, problem.i, xend = chromEnd/1e3, yend = problem.i + 0.05),
                   data = prob.regions, size = 1, color = "black",
                   showSelected = "bases.per.problem", clickSelects = "problem.name") +
      ggtitle("select problem") +
      geom_text(aes(chromStart/1e3, problem.i,
                    label = sprintf("%d problems mean size %.1f kb", problems, mean.bases/1e3)),
                data = PSJ$problem.labels, hjust = 0, showSelected = "bases.per.problem") +
      geom_segment(aes(problemStart/1e3, problem.i, xend = problemEnd/1e3, yend = problem.i),
                   data = PSJ$problems, size = 5, color = "black",
                   showSelected = "bases.per.problem", clickSelects = "problem.name") +
      scale_y_continuous("aligned read coverage", breaks = function(limits) floor(limits[2])) +
      scale_linetype_manual("error type", limits = c("correct", "false negative", "false positive"),
                            values = c(correct = 0, "false negative" = 3, "false positive" = 1)) +
      scale_x_continuous("position on chr11 (kilo bases = kb)") +
      coord_cartesian(xlim = c(118167.406, 118238.833)) +
      geom_tallrect(aes(xmin = chromStart/1e3, xmax = chromEnd/1e3, fill = annotation),
                    data = PSJ$filled.regions, alpha = 0.5, color = "grey") +
      scale_fill_manual(values = ann.colors) +
      theme_bw() +
      theme_animint(width = 1500, height = facet.rows * 100) +
      theme(
        panel.margin = grid::unit(0, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        plot.title = element_text(margin = margin(b = 20))
      ) +
      facet_grid(sample.id ~ ., labeller = function(df) {
        df$sample.id <- sub("McGill0", "", sub(" ", "\n", df$sample.id))
        df
      }, scales = "free") +
      geom_line(aes(base/1e3, count), data = PSJ$coverage, color = "grey50"),

    resError = ggplot() +
      ggtitle("select problem size") +
      ylab("minimum percent incorrect regions") +
      geom_tallrect(aes(xmin = min.bases.per.problem, xmax = max.bases.per.problem),
                    data = res.error, alpha = 0.5, clickSelects = "bases.per.problem") +
      scale_x_log10() +
      geom_line(aes(bases.per.problem, errors/regions * 100, color = chunks, size = chunks),
                data = data.frame(res.error, chunks = "this")) +
      geom_line(aes(bases.per.problem, errors/regions * 100, color = chunks, size = chunks),
                data = data.frame(PSJ$error.total.all, chunks = "all")) +
      theme(
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        plot.margin = unit(c(2, 1, 1, 1), "cm")
      ),

    modelSelection = ggplot() +
      geom_segment(aes(min.log.lambda, peaks, xend = max.log.lambda, yend = peaks),
                   data = data.frame(all.modelSelection, what = factor("peaks", levels = c("peaks", "errors"))), size = 5,
                   showSelected = c("problem.name", "bases.per.problem")) +
      geom_text(aes(min.log.lambda, peaks,
                    label = sprintf("%.1f kb in problem %s", (problemEnd - problemStart)/1e3, problem.name)),
                data = data.frame(modelSelection.labels, what = factor("peaks", levels = c("peaks", "errors"))),
                showSelected = c("problem.name", "bases.per.problem")) +
      geom_segment(aes(min.log.lambda, as.integer(errors), xend = max.log.lambda, yend = as.integer(errors)),
                   data = data.frame(modelSelection.errors, what = factor("errors", levels = c("peaks", "errors"))), size = 5,
                   showSelected = c("problem.name", "bases.per.problem")) +
      ggtitle("select number of samples with 1 peak") +
      ylab("") +
      facet_grid(what ~ ., scales = "free", labeller = as_labeller(c(peaks="Number of peaks", errors="Label errors"))) +
      theme(
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank()
      ),

    title = "Animint compiler with for loops",
    first = PSJ$first
  )

  for (problem.dot in names(PSJ$modelSelection.by.problem)) {
    regions.dt <- PSJ$regions.by.problem[[problem.dot]]
    regions.dt[[problem.dot]] <- regions.dt$peaks
    if (!is.null(regions.dt)) {
      regions.dt$sample.id <- factor(regions.dt$sample.id, levels = sample.levels)
      viz.for$coverage <- viz.for$coverage +
        geom_tallrect(aes(xmin = chromStart/1e3, xmax = chromEnd/1e3, linetype = status),
                      data = data.frame(regions.dt), fill = NA, color = "black",
                      showSelected = c(problem.dot, "bases.per.problem"))
    }
    if (problem.dot %in% names(PSJ$peaks.by.problem)) {
      peaks <- PSJ$peaks.by.problem[[problem.dot]]
      peaks[[problem.dot]] <- peaks$peaks
      prob.peaks.names <- c("bases.per.problem", "problem.i", "problem.name", "chromStart", "chromEnd", problem.dot)
      prob.peaks <- unique(data.frame(peaks)[, prob.peaks.names])
      prob.peaks$sample.id <- factor("problems", levels = sample.levels)
      viz.for$coverage <- viz.for$coverage +
        geom_segment(aes(chromStart/1e3, problem.i, xend = chromEnd/1e3, yend = problem.i + 0.05),
                     data = prob.peaks, size = 7, color = "deepskyblue",
                     clickSelects = "problem.name", showSelected = c(problem.dot, "bases.per.problem")) +
        geom_segment(aes(chromStart/1e3, 0, xend = chromEnd/1e3, yend = 0),
                     data = peaks, size = 7, color = "deepskyblue",
                     clickSelects = "problem.name", showSelected = c(problem.dot, "bases.per.problem"))
    }
    modelSelection.dt <- PSJ$modelSelection.by.problem[[problem.dot]]
    modelSelection.dt[[problem.dot]] <- modelSelection.dt$peaks
    viz.for$modelSelection <- viz.for$modelSelection +
      geom_tallrect(aes(xmin = min.log.lambda, xmax = max.log.lambda),
                    data = modelSelection.dt, alpha = 0.5,
                    clickSelects = problem.dot, showSelected = c("problem.name", "bases.per.problem"))
  }
}))

cat("compiling data viz with for loops\n")
print(system.time({
  animint2dir(viz.for, out.dir = "PSJ-for-loops")
}))

# Visualization with .variable .value aesthetics
sample.peaks <- do.call(rbind, PSJ$peaks.by.problem)
prob.peaks.names <- c("bases.per.problem", "problem.i", "problem.name", "peaks", "chromStart", "chromEnd")
problem.peaks <- unique(data.frame(sample.peaks)[, prob.peaks.names])
problem.peaks$sample.id <- "problems"

peakvar <- function(position) {
  paste0(gsub("[-:]", ".", position), "peaks")
}

cat("constructing data viz with .variable .value\n")
print(system.time({
  viz <- list(
    coverage = ggplot() +
      geom_segment(aes(chromStart/1e3, problem.i, xend = chromEnd/1e3, yend = problem.i + 0.05),
                   data = prob.regions, size = 1, color = "black",
                   showSelected = "bases.per.problem", clickSelects = "problem.name") +
      ggtitle("select problem") +
      geom_text(aes(chromStart/1e3, problem.i,
                    label = sprintf("%d problems mean size %.1f kb", problems, mean.bases/1e3)),
                data = PSJ$problem.labels, hjust = 0, showSelected = "bases.per.problem") +
      geom_segment(aes(problemStart/1e3, problem.i, xend = problemEnd/1e3, yend = problem.i),
                   data = PSJ$problems, size = 5, color = "black",
                   showSelected = "bases.per.problem", clickSelects = "problem.name") +
      scale_y_continuous("aligned read coverage", breaks = function(limits) floor(limits[2])) +
      scale_linetype_manual("error type", limits = c("correct", "false negative", "false positive"),
                            values = c(correct = 0, "false negative" = 3, "false positive" = 1)) +
      scale_x_continuous("position on chr11 (kilo bases = kb)") +
      coord_cartesian(xlim = c(118167.406, 118238.833)) +
      geom_tallrect(aes(xmin = chromStart/1e3, xmax = chromEnd/1e3, fill = annotation),
                    data = PSJ$filled.regions, alpha = 0.5, color = "grey") +
      scale_fill_manual(values = ann.colors) +
      theme_bw() +
      theme_animint(width = 1500, height = facet.rows * 100) +
      theme(
        panel.margin = grid::unit(0, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        plot.title = element_text(margin = margin(b = 20))
      ) +
      facet_grid(sample.id ~ ., labeller = function(df) {
        df$sample.id <- sub("McGill0", "", sub(" ", "\n", df$sample.id))
        df
      }, scales = "free") +
      geom_line(aes(base/1e3, count), data = PSJ$coverage, color = "grey50") +
      geom_tallrect(aes(xmin = chromStart/1e3, xmax = chromEnd/1e3, linetype = status),
                    data = all.regions, fill = NA, color = "black",
                    showSelected.variable = peakvar(problem.name), showSelected.value = "peaks",
                    showSelected = "bases.per.problem") +
      geom_segment(aes(chromStart/1e3, 0, xend = chromEnd/1e3, yend = 0),
                   data = sample.peaks, size = 7, color = "deepskyblue",
                   clickSelects = "problem.name", showSelected.variable = peakvar(problem.name),
                   showSelected.value = "peaks", showSelected = "bases.per.problem") +
      geom_segment(aes(chromStart/1e3, problem.i, xend = chromEnd/1e3, yend = problem.i + 0.05),
                   data = problem.peaks, size = 7, color = "deepskyblue",
                   clickSelects = "problem.name", showSelected.variable = peakvar(problem.name),
                   showSelected.value = "peaks", showSelected = "bases.per.problem"),

    resError = ggplot() +
      ggtitle("select problem size") +
      ylab("minimum percent incorrect regions") +
      geom_tallrect(aes(xmin = min.bases.per.problem, xmax = max.bases.per.problem),
                    data = res.error, alpha = 0.5, clickSelects = "bases.per.problem") +
      scale_x_log10() +
      geom_line(aes(bases.per.problem, errors/regions * 100, color = chunks, size = chunks),
                data = data.frame(res.error, chunks = "this")) +
      geom_line(aes(bases.per.problem, errors/regions * 100, color = chunks, size = chunks),
                data = data.frame(PSJ$error.total.all, chunks = "all")) +
      theme(
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        plot.margin = unit(c(2, 1, 1, 1), "cm")
      ),

    modelSelection = ggplot() +
      geom_segment(aes(min.log.lambda, peaks, xend = max.log.lambda, yend = peaks),
                   data = data.frame(all.modelSelection, what = "peaks"), size = 5,
                   showSelected = c("problem.name", "bases.per.problem")) +
      geom_text(aes(min.log.lambda, peaks,
                    label = sprintf("%.1f kb in problem %s", (problemEnd - problemStart)/1e3, problem.name)),
                data = data.frame(modelSelection.labels, what = "peaks"),
                showSelected = c("problem.name", "bases.per.problem")) +
      geom_segment(aes(min.log.lambda, as.integer(errors), xend = max.log.lambda, yend = as.integer(errors)),
                   data = data.frame(modelSelection.errors, what = "errors"), size = 5,
                   showSelected = c("problem.name", "bases.per.problem")) +
      ggtitle("select number of samples with 1 peak") +
      ylab("") +
      facet_grid(what ~ ., scales = "free", labeller = as_labeller(c(peaks="Number of peaks", errors="Label errors"))) +
      theme(
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank()
      ),

    title = "Animint compiler with .variable .value aesthetics",
    first = PSJ$first
  )
}))

cat("compiling data viz with .variable .value\n")
print(system.time({
  animint2dir(viz, out.dir = "PSJ-var-value")
}))

# Visualization with lookup functions
not.problem.name <- names(PSJ$first)[names(PSJ$first) != "problem.name"]

problems.by.res <- split(PSJ$problems, PSJ$problems$bases.per.problem)

# get_ functions
get_tallrect <- function(...) {
  L <- list(...)
  res.problems <- problems.by.res[[L$bases.per.problem]]
  problem.name.vec <- paste(res.problems$problem.name)
  peaks.by.problem <- list()
  for (problem.name in problem.name.vec) {
    selector.name <- peakvar(problem.name)
    prob.peaks <- PSJ$peaks.by.problem[[selector.name]]
    peaks.by.problem[[problem.name]] <- subset(prob.peaks, peaks == L[[selector.name]])
  }
  show.peaks <- do.call(rbind, peaks.by.problem)
  error.regions <- PeakErrorSamples(show.peaks, PSJ$filled.regions)
  error.regions
}

problem.name <- "chr11:118174946-118177139"  # Example; adjust if needed
bases.per.problem <- "6516"  # Example; adjust if needed
get_segment <- function(problem.name, bases.per.problem, ...) {
  L <- list(...)
  res.problems <- problems.by.res[[bases.per.problem]]
  is.selected <- res.problems$problem.name == problem.name
  other.problems <- res.problems[!is.selected, ]
  peaks.by.problem <- list()
  for (other.name in other.problems$problem.name) {
    selector.name <- peakvar(other.name)
    prob.peaks <- PSJ$peaks.by.problem[[selector.name]]
    peaks.by.problem[[other.name]] <- subset(prob.peaks, peaks == L[[selector.name]])
  }
  other.peaks <- do.call(rbind, peaks.by.problem)
  
  selector.name <- peakvar(problem.name)
  prob.peaks <- PSJ$peaks.by.problem[[selector.name]]
  modelSelection <- PSJ$modelSelection.by.problem[[selector.name]]
  modelSelection$errors <- NA
  for (model.i in 1:nrow(modelSelection)) {
    model.row <- modelSelection[model.i, ]
    model.peaks <- subset(prob.peaks, peaks == model.row$peaks)
    all.peaks <- rbind(other.peaks, model.peaks)
    error.regions <- PeakErrorSamples(all.peaks, PSJ$filled.regions)
    modelSelection$errors[model.i] <- with(error.regions, sum(fp + fn))
  }
  modelSelection
}

# Cache functions
cache_tallrect <- function(bases.per.problem) {
  res.problems <- problems.by.res[[bases.per.problem]]
  peakvar(res.problems$problem.name)
}
cache_segment <- function(bases.per.problem, ...) {
  res.problems <- problems.by.res[[bases.per.problem]]
  is.selected <- res.problems$problem.name == problem.name
  other.problems <- res.problems[!is.selected, ]
  peakvar(other.problems$problem.name)
}

cat("constructing data viz with lookup function\n")
print(system.time({
  viz <- list(
    coverage = ggplot() +
      geom_segment(aes(chromStart/1e3, problem.i, xend = chromEnd/1e3, yend = problem.i + 0.05),
                   data = prob.regions, size = 1, color = "black",
                   showSelected = "bases.per.problem", clickSelects = "problem.name") +
      ggtitle("select problem") +
      geom_text(aes(chromStart/1e3, problem.i,
                    label = sprintf("%d problems mean size %.1f kb", problems, mean.bases/1e3)),
                data = PSJ$problem.labels, hjust = 0, showSelected = "bases.per.problem") +
      geom_segment(aes(problemStart/1e3, problem.i, xend = problemEnd/1e3, yend = problem.i),
                   data = PSJ$problems, size = 5, color = "black",
                   showSelected = "bases.per.problem", clickSelects = "problem.name") +
      scale_y_continuous("aligned read coverage", breaks = function(limits) floor(limits[2])) +
      scale_linetype_manual("error type", limits = c("correct", "false negative", "false positive"),
                            values = c(correct = 0, "false negative" = 3, "false positive" = 1)) +
      scale_x_continuous("position on chr11 (kilo bases = kb)") +
      coord_cartesian(xlim = c(118167.406, 118238.833)) +
      geom_tallrect(aes(xmin = chromStart/1e3, xmax = chromEnd/1e3, fill = annotation),
                    data = PSJ$filled.regions, alpha = 0.5, color = "grey") +
      scale_fill_manual(values = ann.colors) +
      theme_bw() +
      theme_animint(width = 1500, height = facet.rows * 100) +
      theme(
        panel.margin = grid::unit(0, "cm"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        plot.title = element_text(margin = margin(b = 20))
      ) +
      facet_grid(sample.id ~ ., labeller = function(df) {
        df$sample.id <- sub("McGill0", "", sub(" ", "\n", df$sample.id))
        df
      }, scales = "free") +
      geom_line(aes(base/1e3, count), data = PSJ$coverage, color = "grey50") +
      geom_tallrect(aes(xmin = chromStart/1e3, xmax = chromEnd/1e3, linetype = status),
                    updateWhenChanged = not.problem.name, cache = cache_tallrect,
                    data = data.frame(get_tallrect), fill = NA, color = "black",
                    showSelected.variable = peakvar(problem.name), showSelected.value = "peaks",
                    showSelected = "bases.per.problem") +
      geom_segment(aes(chromStart/1e3, 0, xend = chromEnd/1e3, yend = 0),
                   data = sample.peaks, size = 7, color = "deepskyblue",
                   clickSelects = "problem.name", showSelected.variable = peakvar(problem.name),
                   showSelected.value = "peaks", showSelected = "bases.per.problem") +
      geom_segment(aes(chromStart/1e3, problem.i, xend = chromEnd/1e3, yend = problem.i + 0.05),
                   data = problem.peaks, size = 7, color = "deepskyblue",
                   clickSelects = "problem.name", showSelected.variable = peakvar(problem.name),
                   showSelected.value = "peaks", showSelected = "bases.per.problem"),

    resError = ggplot() +
      ggtitle("select problem size") +
      ylab("minimum percent incorrect regions") +
      geom_tallrect(aes(xmin = min.bases.per.problem, xmax = max.bases.per.problem),
                    data = res.error, alpha = 0.5, clickSelects = "bases.per.problem") +
      scale_x_log10() +
      geom_line(aes(bases.per.problem, errors/regions * 100, color = chunks, size = chunks),
                data = data.frame(res.error, chunks = "this")) +
      geom_line(aes(bases.per.problem, errors/regions * 100, color = chunks, size = chunks),
                data = data.frame(PSJ$error.total.all, chunks = "all")) +
      theme(
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank(),
        plot.margin = unit(c(2, 1, 1, 1), "cm")
      ),

    modelSelection = ggplot() +
      ggtitle("select number of samples with 1 peak") +
      xlab("model complexity penalty log(lambda)") +
      ylab("") +
      facet_grid(what ~ ., scales = "free",
                 labeller = as_labeller(c(peaks="Number of peaks", errors="Label errors"))) +
      geom_tallrect(
        aes(xmin = min.log.lambda, xmax = max.log.lambda),
        data = all.modelSelection, alpha = 0.5,
        clickSelects = "peaks",
        showSelected = c("problem.name", "bases.per.problem")
      ) +
      geom_segment(
        aes(min.log.lambda, peaks, xend = max.log.lambda, yend = peaks),
        data = data.frame(all.modelSelection, what = factor("peaks", levels = c("peaks", "errors"))), 
        size = 5,
        showSelected = c("problem.name", "bases.per.problem")
      ) +
      geom_text(
        aes(min.log.lambda, peaks,
            label = sprintf("%.1f kb in problem %s", (problemEnd - problemStart) / 1e3, problem.name)),
        data = data.frame(modelSelection.labels, what = factor("peaks", levels = c("peaks", "errors"))), 
        hjust = 0.5,
        showSelected = c("problem.name", "bases.per.problem")
      ) +
      geom_segment(
        aes(min.log.lambda, as.integer(errors), xend = max.log.lambda, yend = as.integer(errors)),
        data = data.frame(modelSelection.errors, what = factor("errors", levels = c("peaks", "errors"))), 
        size = 5,
        showSelected = c("problem.name", "bases.per.problem")
      ) +
      theme(
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_blank()
      ),

    title = "Animint compiler with lookup function",
    first = PSJ$first
  )
}))

cat("compiling data viz with lookup function\n")
print(system.time({
  animint2dir(viz, out.dir = "PSJ-lookup")
}))