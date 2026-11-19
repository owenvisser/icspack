# Figure 2: Original DFS=0 versus DFS>0 periodontal distributions,
# stratified by the patient's total observed tooth count.
# MEC-weighted tooth percentages sum to 100% within each DFS group.
stopifnot(requireNamespace("here", quietly = TRUE))
source(here::here("R", "common.R"))
d <- load_teeth()
fig_dir <- here::here("figures", "nhanes")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
filename <- file.path(fig_dir, "F2_NHANES_DFS_only_CAL_PD_pattern")
levels6 <- c("<1", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "5+")
status_levels <- c("Not present", "Present")
count_levels <- c("10-15 teeth", "16-21 teeth", "22-28 teeth")

# Use all observed teeth when determining patient tooth count, before
# outcome-specific filtering. Retain the original Figure 2 eligibility rule.
if (anyDuplicated(d[c("SEQN", "tooth_number")])) {
  stop("Duplicate patient/tooth records")
}
tooth_counts <- table(d$SEQN)

make_dist <- function(outcome, group) {
  required <- c("SEQN", outcome, "DFS", "WTMEC4YR")
  ok <- stats::complete.cases(d[, required, drop = FALSE]) &
    is.finite(d[[outcome]]) & is.finite(d$DFS) &
    is.finite(d$WTMEC4YR) & d$WTMEC4YR > 0
  z <- d[ok, required, drop = FALSE]
  complete_counts <- table(z$SEQN)
  z <- z[as.character(z$SEQN) %in%
           names(complete_counts)[complete_counts >= 2L], , drop = FALSE]
  z$tooth_count <- as.integer(tooth_counts[as.character(z$SEQN)])
  z$tooth_group <- cut(z$tooth_count, breaks = c(9, 15, 21, 28),
                       labels = count_levels)
  z <- z[z$tooth_group == group, , drop = FALSE]
  if (!nrow(z)) stop("No eligible teeth for ", outcome, ", ", group)
  z$status <- factor(ifelse(z$DFS > 0, "Present", "Not present"),
                     levels = status_levels)
  z$category <- cut(z[[outcome]], breaks = c(-Inf, 1, 2, 3, 4, 5, Inf),
                    right = FALSE, labels = levels6)
  sums <- stats::xtabs(WTMEC4YR ~ status + category, data = z,
                       drop.unused.levels = FALSE)
  if (any(rowSums(sums) == 0)) stop("Empty DFS group: ", outcome, ", ", group)
  props <- sweep(sums, 1L, rowSums(sums), "/")
  stopifnot(all(abs(rowSums(props) - 1) < 1e-10))
  rows <- expand.grid(dfs_status = status_levels, category = levels6,
                      stringsAsFactors = FALSE)
  rows$proportion <- mapply(function(s, c) props[s, c],
                            rows$dfs_status, rows$category)
  rows$periodontal <- outcome
  rows$tooth_count_group <- group
  rows$participants <- length(unique(z$SEQN))
  rows$teeth <- nrow(z)
  rows$teeth_in_dfs_group <- as.integer(table(z$status)[rows$dfs_status])
  list(props = props, summary = rows)
}

panels <- list()
summaries <- list()
for (outcome in c("CAL", "PD")) {
  for (group in count_levels) {
    key <- paste(outcome, group, sep = "_")
    panels[[key]] <- make_dist(outcome, group)
    summaries[[key]] <- panels[[key]]$summary
  }
}
utils::write.csv(do.call(rbind, summaries),
  file.path(fig_dir, "F2_NHANES_DFS_by_tooth_count.csv"), row.names = FALSE)

# Original black-and-white hatches: crosshatch (DFS=0), diagonal (DFS>0).
hatched_rect <- function(left, bottom, right, top, crossed) {
  graphics::rect(left, bottom, right, top, col = "white", border = NA)
  graphics::rect(left, bottom, right, top, density = 15, angle = 45,
                 col = "black", border = NA)
  if (crossed) graphics::rect(left, bottom, right, top, density = 15,
                              angle = -45, col = "black", border = NA)
  graphics::rect(left, bottom, right, top, border = "black", lwd = 0.65)
}

max_share <- max(vapply(panels, function(x) max(x$props), numeric(1)))
ylim_top <- max(0.62, max_share * 1.12)
ticks <- pretty(c(0, ylim_top), n = 5)
ticks <- ticks[ticks >= 0 & ticks <= ylim_top]

plot_panel <- function(outcome, group, show_y) {
  p <- panels[[paste(outcome, group, sep = "_")]]$props
  graphics::par(mar = c(4.7, if (show_y) 4.8 else 3.1, 2.1, 0.55),
                mgp = c(2.45, 0.65, 0))
  mids <- graphics::barplot(p, beside = TRUE, width = 1,
                           space = c(0.12, 0.65),
                           names.arg = levels6, col = "white", border = NA,
                           axes = FALSE, main = group,
                           ylim = c(0, ylim_top), cex.names = 0.81,
                           cex.main = 1.02,
                           ylab = if (show_y) "Percent within DFS group" else "",
                           xlab = paste("Max", outcome, "(mm)"))
  for (j in seq_len(ncol(p))) {
    for (i in seq_len(nrow(p))) {
      hatched_rect(mids[i, j] - 0.5, 0, mids[i, j] + 0.5, p[i, j],
                   crossed = (i == 1L))
    }
  }
  graphics::axis(2, at = ticks, labels = sprintf("%g%%", 100 * ticks),
                 las = 1, cex.axis = 0.82)
  graphics::box(bty = "l")
}

plot_legend <- function() {
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::text(0.32, 0.5, "DFS:", adj = c(1, 0.5), cex = 0.98)
  hatched_rect(0.35, 0.28, 0.38, 0.72, crossed = TRUE)
  graphics::text(0.39, 0.5, "Not present", adj = c(0, 0.5), cex = 0.94)
  hatched_rect(0.56, 0.28, 0.59, 0.72, crossed = FALSE)
  graphics::text(0.60, 0.5, "Present", adj = c(0, 0.5), cex = 0.94)
}

plot_f2 <- function() {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::layout(1); graphics::par(oldpar) }, add = TRUE)
  graphics::layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7, 7),
                          nrow = 3, byrow = TRUE),
                   heights = c(1, 1, 0.14))
  for (outcome in c("CAL", "PD")) {
    for (j in seq_along(count_levels)) {
      plot_panel(outcome, count_levels[j], show_y = (j == 1L))
    }
  }
  plot_legend()
}

grDevices::pdf(paste0(filename, ".pdf"), width = 15, height = 8.3,
               useDingbats = FALSE)
plot_f2()
grDevices::dev.off()

grDevices::png(paste0(filename, ".png"), width = 4500, height = 2490,
               res = 300)
plot_f2()
grDevices::dev.off()
message("Saved Figure 2, stratified by patient tooth count:\n",
        "  ", filename, ".pdf\n",
        "  ", filename, ".png\n",
        "  ", file.path(fig_dir, "F2_NHANES_DFS_by_tooth_count.csv"))
