
# Figure F2 ONLY: black-and-white DFS presence by CAL and PD categories.
# Not present = crosshatch; present = diagonal hatch. No CI plots.

stopifnot(requireNamespace("here", quietly = TRUE))

source(here::here("R", "common.R"))

d <- load_teeth()

fig_dir <- here::here("figures", "nhanes")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

levels6 <- c("(-Inf,1)", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,Inf)")
status_levels <- c("Not present", "Present")

make_dist <- function(x) {
  required <- c("SEQN", x, "DFS", "WTMEC4YR")

  ok <- stats::complete.cases(d[, required, drop = FALSE]) &
    is.finite(as.numeric(d[[x]])) &
    is.finite(as.numeric(d$DFS)) &
    is.finite(d$WTMEC4YR) &
    d$WTMEC4YR > 0

  z <- d[ok, required, drop = FALSE]

  counts <- table(z$SEQN)

  z <- z[
    as.character(z$SEQN) %in% names(counts)[counts >= 2L],
    ,
    drop = FALSE
  ]

  if (!nrow(z)) {
    stop("No eligible complete tooth pairs for ", x, "/DFS")
  }

  category <- cut(
    as.numeric(z[[x]]),
    breaks = c(-Inf, 1, 2, 3, 4, 5, Inf),
    right = FALSE,
    labels = levels6
  )

  status <- factor(
    ifelse(z$DFS > 0, "Present", "Not present"),
    levels = status_levels
  )

  sums <- matrix(
    0,
    nrow = 2L,
    ncol = length(levels6),
    dimnames = list(status_levels, levels6)
  )

  for (s in status_levels) {
    for (c in levels6) {
      take <- !is.na(category) & status == s & category == c
      sums[s, c] <- sum(z$WTMEC4YR[take])
    }
  }

  if (any(rowSums(sums) <= 0)) {
    stop("Empty DFS status group for ", x)
  }

  props <- sums / rowSums(sums)

  rows <- expand.grid(
    dfs_status = status_levels,
    category = levels6,
    stringsAsFactors = FALSE
  )

  rows$periodontal <- x

  rows$proportion <- mapply(
    function(st, ca) props[st, ca],
    rows$dfs_status,
    rows$category
  )

  rows$complete_tooth_rows <- nrow(z)
  rows$participants <- length(unique(z$SEQN))

  list(props = props, summary = rows)
}

cal <- make_dist("CAL")
pd <- make_dist("PD")

# Drawing two sets of line hatches on the same white rectangle gives true
# cross-hatching. 'barplot(density=..., angle=c(45,-45))' only gives opposite
# diagonal patterns, not crossed lines.

hatched_rect <- function(left, bottom, right, top, crossed) {
  graphics::rect(left, bottom, right, top, col = "white", border = NA)

  graphics::rect(
    left,
    bottom,
    right,
    top,
    density = 16,
    angle = 45,
    col = "black",
    border = NA
  )

  if (crossed) {
    graphics::rect(
      left,
      bottom,
      right,
      top,
      density = 16,
      angle = -45,
      col = "black",
      border = NA
    )
  }

  graphics::rect(left, bottom, right, top, border = "black", lwd = 0.7)
}

plot_f2 <- function() {
  oldpar <- graphics::par(no.readonly = TRUE)

  on.exit(
    {
      graphics::layout(1)
      graphics::par(oldpar)
    },
    add = TRUE
  )

  graphics::layout(
    matrix(c(1, 2, 3, 3), ncol = 2, byrow = TRUE),
    heights = c(8.5, 1.3)
  )

  ylim_top <- max(0.60, max(cal$props, pd$props) * 1.10)

  ticks <- pretty(c(0, ylim_top), n = 6)
  ticks <- ticks[ticks >= 0 & ticks <= ylim_top]

  for (x in c("CAL", "PD")) {
    props <- if (x == "CAL") cal$props else pd$props

    graphics::par(
      mar = c(4.6, if (x == "CAL") 5.0 else 2.6, 2.3, 0.8),
      mgp = c(2.6, 0.7, 0)
    )

    mids <- graphics::barplot(
      props,
      beside = TRUE,
      width = 1,
      space = c(0.12, 0.75),
      names.arg = levels6,
      col = "white",
      border = NA,
      axes = FALSE,
      main = x,
      ylim = c(0, ylim_top),
      ylab = if (x == "CAL") "Proportion within DFS group" else "",
      cex.names = 0.8,
      cex.main = 1.15
    )

    for (j in seq_len(ncol(props))) {
      for (i in seq_len(nrow(props))) {
        # First row: absent = crossed lines; second row: present = diagonal lines.
        hatched_rect(
          mids[i, j] - 0.5,
          0,
          mids[i, j] + 0.5,
          props[i, j],
          crossed = (i == 1L)
        )
      }
    }

    graphics::axis(
      2,
      at = ticks,
      labels = sprintf("%g%%", 100 * ticks),
      las = 1,
      cex.axis = 0.85
    )

    graphics::box(bty = "l")
  }

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

  prefix <- "DFS:"
  a <- "Not present"
  b <- "Present"

  swatch <- 0.045
  gap <- 0.014
  item_gap <- 0.075

  p_w <- graphics::strwidth(prefix, cex = 0.95)
  a_w <- graphics::strwidth(a, cex = 0.95)
  b_w <- graphics::strwidth(b, cex = 0.95)

  total <- p_w + item_gap + swatch + gap + a_w + item_gap + swatch + gap + b_w
  xpos <- (1 - total) / 2

  graphics::text(xpos, 0.55, prefix, adj = c(0, 0.5), cex = 0.95)

  xpos <- xpos + p_w + item_gap

  hatched_rect(xpos, 0.37, xpos + swatch, 0.73, crossed = TRUE)

  xpos <- xpos + swatch + gap
  graphics::text(xpos, 0.55, a, adj = c(0, 0.5), cex = 0.95)

  xpos <- xpos + a_w + item_gap

  hatched_rect(xpos, 0.37, xpos + swatch, 0.73, crossed = FALSE)

  graphics::text(xpos + swatch + gap, 0.55, b, adj = c(0, 0.5), cex = 0.95)
}

filename <- file.path(fig_dir, "F2_NHANES_DFS_only_CAL_PD_pattern")

grDevices::pdf(
  paste0(filename, ".pdf"),
  width = 11,
  height = 5.4,
  useDingbats = FALSE
)

plot_f2()
grDevices::dev.off()

grDevices::png(
  paste0(filename, ".png"),
  width = 3300,
  height = 1620,
  res = 300
)

plot_f2()
grDevices::dev.off()

message("Saved only the black-and-white NHANES F2 figure")