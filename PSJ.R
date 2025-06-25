library(ggplot2)
library(animint2)
library(PeakSegJoint)
library(data.table)
library(grid)

load("PSJ.RData")

required_objects <- c(
  "coverage", "regions.by.problem", "modelSelection.by.problem",
  "problems", "problem.labels", "filled.regions", "error.total.chunk",
  "error.total.all", "peaks.by.problem"
)
missing_objects <- setdiff(required_objects, names(PSJ))
if (length(missing_objects) > 0) {
  stop("The following required data objects are missing from PSJ.RData: ", paste(missing_objects, collapse=", "))
}

all.regions <- do.call(rbind, PSJ$regions.by.problem)
all.modelSel <- do.call(rbind, PSJ$modelSelection.by.problem)
sample.peaks <- do.call(rbind, PSJ$peaks.by.problem)
res.error <- PSJ$error.total.chunk

sample.levels <- c("problems", setdiff(unique(as.character(PSJ$coverage$sample.id)), "problems"))

prob.regions <- unique(data.frame(all.regions)[, c("bases.per.problem", "problem.i", "problem.name", "chromStart", "chromEnd")])
prob.regions$sample.id <- "problems"

problem.peaks <- unique(data.frame(sample.peaks)[, c("bases.per.problem", "problem.i", "problem.name", "peaks", "chromStart", "chromEnd")])
problem.peaks$sample.id <- "problems"

PSJ$problems$sample.id <- "problems"
PSJ$problem.labels$sample.id <- "problems"

problems.coverage <- data.frame(base = prob.regions$chromStart, count = NA, sample.id = "problems")
combined.coverage <- rbind(PSJ$coverage[, c("base", "count", "sample.id")], problems.coverage)

combined.coverage$sample.id <- factor(combined.coverage$sample.id, levels = sample.levels)
prob.regions$sample.id      <- factor(prob.regions$sample.id,      levels = sample.levels)
sample.peaks$sample.id      <- factor(sample.peaks$sample.id,      levels = sample.levels)
problem.peaks$sample.id     <- factor(problem.peaks$sample.id,     levels = sample.levels)
all.regions$sample.id       <- factor(all.regions$sample.id,       levels = sample.levels)
PSJ$filled.regions$sample.id<- factor(PSJ$filled.regions$sample.id,levels = sample.levels)
PSJ$problems$sample.id      <- factor(PSJ$problems$sample.id,      levels = sample.levels)
PSJ$problem.labels$sample.id<- factor(PSJ$problem.labels$sample.id,levels = sample.levels)

PSJ$coverage <- combined.coverage

modelSel.errs <- subset(all.modelSel, !is.na(errors))

pen.range <- with(all.modelSel, c(min(max.log.lambda), max(min.log.lambda)))
pen.mid <- mean(pen.range)
modelSel.lbls <- unique(with(all.modelSel, data.frame(
  problem.name, bases.per.problem, problemStart, problemEnd,
  min.log.lambda = pen.mid, peaks = max(peaks) + 0.5
)))

PSJ$problems <- PSJ$problems[order(PSJ$problems$problemStart), ]
PSJ$problems$problem.i <- seq_along(PSJ$problems$problemStart) * 0.05

prob.regions <- merge(prob.regions, PSJ$problems[, c("problem.name", "problem.i")], by = "problem.name", suffixes = c("", ".new"))
prob.regions$problem.i <- prob.regions$problem.i.new
prob.regions$problem.i.new <- NULL

problem.peaks <- merge(problem.peaks, PSJ$problems[, c("problem.name", "problem.i")], by = "problem.name", suffixes = c("", ".new"))
problem.peaks$problem.i <- problem.peaks$problem.i.new
problem.peaks$problem.i.new <- NULL

dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval - 1) / 2 + 1
res.error$min.bases.per.problem <- res.error$bases.per.problem / dval2
res.error$max.bases.per.problem <- res.error$bases.per.problem * dval2

PSJ$first <- list(
  problem.name = prob.regions$problem.name[1],
  bases.per.problem = prob.regions$bases.per.problem[1]
)

ann.colors <- c(noPeaks = "#f6f4bf", peakStart = "#ffafaf", peakEnd = "#ff4c4c", peaks = "#a445ee")

viz <- list(
  title = "PeakSegJoint Interactive Visualization",
  
  coverage = ggplot() +
    facet_grid(sample.id ~ ., scales = "free_y") +
    geom_tallrect(
      aes(xmin = chromStart / 1e3, xmax = chromEnd / 1e3, fill = annotation),
      data = PSJ$filled.regions, alpha = 0.5
    ) +
    geom_line(
      aes(base / 1e3, count),
      data = PSJ$coverage, color = "grey50"
    ) +
    geom_segment(
      aes(x = chromStart / 1e3, xend = chromEnd / 1e3, y = problem.i, yend = problem.i + 0.05),
      data = prob.regions,
      size = 1, color = "black",
      showSelected = "bases.per.problem", clickSelects = "problem.name",
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = problemStart / 1e3, xend = problemEnd / 1e3, y = problem.i, yend = problem.i),
      data = PSJ$problems,
      size = 5, color = "black",
      showSelected = "bases.per.problem", clickSelects = "problem.name",
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(chromStart / 1e3, problem.i -50,
          label = sprintf("%d problems mean size %.1f kb", problems, mean.bases / 1e3)),
      data = PSJ$problem.labels,
      hjust = 0,
      size = 10,
      color = "black",
      showSelected = "bases.per.problem"
    ) +
    geom_tallrect(
      aes(xmin = chromStart / 1e3, xmax = chromEnd / 1e3, linetype = status),
      data = all.regions, fill = NA, color = "black",
      showSelected = c("problem.name", "peaks"),
      showSelected2 = "bases.per.problem",
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = chromStart / 1e3, xend = chromEnd / 1e3, y = 1, yend = 1),
      data = subset(sample.peaks, sample.id != "problems"),
      size = 7, color = "deepskyblue",
      clickSelects = "problem.name",
      showSelected = c("problem.name", "peaks"),
      showSelected2 = "bases.per.problem",
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = chromStart / 1e3, xend = chromEnd / 1e3, y = problem.i, yend = problem.i),
      data = problem.peaks,
      size = 7, color = "deepskyblue",
      clickSelects = "problem.name",
      showSelected = c("problem.name", "peaks"),
      showSelected2 = "bases.per.problem",
      inherit.aes = FALSE
    ) +
    scale_y_continuous("aligned read coverage", breaks = function(l) floor(l[2])) +
    scale_linetype_manual(
      "error type",
      limits = c("correct", "false negative", "false positive"),
      values = c(correct = 0, "false negative" = 3, "false positive" = 1)
    ) +
    scale_x_continuous("position on chr11 (kilo bases = kb)") +
    scale_fill_manual(values = ann.colors) +
    coord_cartesian(xlim = c(118167.406, 118238.833)) +
    theme_bw() +
    theme_animint(width = 1500, height = length(sample.levels) * 120) +
    theme(panel.margin = grid::unit(0, "cm")),

   resError = ggplot() +
    ggtitle("select problem size") +
    ylab("minimum percent incorrect regions") +
    geom_tallrect(
      aes(xmin=min.bases.per.problem, xmax=max.bases.per.problem),
      data=res.error, alpha=0.5, clickSelects="bases.per.problem"
    ) +
    scale_x_log10() +
    geom_line(
      aes(bases.per.problem, errors/regions*100, color=chunks, size=chunks),
      data=data.frame(res.error, chunks="this")
    ) +
    geom_line(
      aes(bases.per.problem, errors/regions*100, color=chunks, size=chunks),
      data=data.frame(PSJ$error.total.all, chunks="all")
    ),

  modelSelection = ggplot() +
    ggtitle("select number of samples with 1 peak") +
    xlab("model complexity penalty log(lambda)") +
    ylab("") +
    facet_grid(what ~ ., scales = "free", 
               labeller = as_labeller(c(peaks="Number of peaks", errors="Label errors"))) +
    geom_tallrect(
      aes(xmin = min.log.lambda, xmax = max.log.lambda),
      data = all.modelSel, alpha = 0.5,
      clickSelects = "peaks",
      showSelected = "problem.name", showSelected2 = "bases.per.problem"
    ) +
    geom_segment(
      aes(min.log.lambda, peaks, xend = max.log.lambda, yend = peaks),
      data = data.frame(all.modelSel, what = "peaks"),
      size = 5,
      showSelected = "problem.name", showSelected2 = "bases.per.problem"
    ) +
    geom_text(
      aes(min.log.lambda, peaks,
          label = sprintf("%.1f kb in problem %s", (problemEnd - problemStart) / 1e3, problem.name)),
      data = data.frame(modelSel.lbls, what = "peaks"),
      hjust=0,
      showSelected = "problem.name", showSelected2 = "bases.per.problem"
    ) +
    geom_segment(
      aes(min.log.lambda, as.integer(errors), xend = max.log.lambda, yend = as.integer(errors)),
      data = data.frame(modelSel.errs, what = "errors"),
      size = 5,
      showSelected = "problem.name", showSelected2 = "bases.per.problem"
    ),
    
  first = PSJ$first,
  selector.types = list(problem.name = "single", bases.per.problem = "single"),
  duration = list(problem.name = 1000, bases.per.problem = 1000)
)

animint2dir(viz, out.dir = "PSJ_interactive_corrected_final")
