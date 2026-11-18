
# Produce ONLY T1 and the supplementary numerical tables for F1 and S1.

stopifnot(requireNamespace("here", quietly = TRUE))

source(here::here("R", "common.R"))
source(here::here("R", "latex.R"))
source(here::here("R", "simulation.R"))

cfg <- simulation_config()
folder <- here::here("results", "simulation", cfg$name)

scenario <- utils::read.csv(file.path(folder, "scenario_summary.csv"))
pooled <- utils::read.csv(file.path(folder, "pooled_MAE_by_K.csv"))

methods <- c("NW", "CW", "EQ", "ER")
settings <- simulation_settings(cfg)

tab <- here::here("manuscript", "generated", cfg$name)

get_s <- function(rho, eta, K, method) {
  z <- scenario[
    scenario$rho_uv == rho &
      scenario$eta == eta &
      scenario$K == K &
      scenario$method == method,
    ,
    drop = FALSE
  ]

  if (nrow(z) != 1L) {
    stop("Missing/duplicated simulation result")
  }

  z
}

# T1: all methods' mean estimate (mean estimated delta SE), at K=20.

lines <- c(
  tex_header(
    sprintf(
      "Mean Pearson correlation (mean delta-method cluster-sandwich standard error) for each simulation setting at $K=20$. Each entry averages %d simulation replicates, with %d potential observations in each of %d clusters. EQ and ER use averaged cluster influence contributions; ER is exhaustive only when all split pairs are evaluated, otherwise it is a capped reference.",
      cfg$reps,
      cfg$nmax,
      cfg$M
    ),
    "tab:sim_estse20",
    "ccrrrr"
  ),
  tex_row(c("$\\rho_{uv}$", "$\\eta$", "$\\rho_0$", methods)),
  tex_cmd("midrule")
)

for (eta in cfg$eta) {
  if (eta != cfg$eta[1]) {
    lines <- c(lines, tex_cmd("addlinespace[3pt]"))
  }

  for (rho in cfg$rho_uv) {
    cells <- c(
      sprintf("%.2f", rho),
      as.character(eta),
      sprintf("%.3f", (rho + cfg$sigma2 / 2) / (1 + cfg$sigma2))
    )

    cells <- c(cells, vapply(
      methods,
      function(m) {
        z <- get_s(rho, eta, 20L, m)
        fmt_est_se(z$mean_estimate, z$mean_se)
      },
      character(1)
    ))

    lines <- c(lines, tex_row(cells))
  }
}

write_tex(
  c(lines, tex_footer()),
  file.path(tab, "T1_simulation_estimate_meanSE.tex")
)

# Supplement: one baseline table for NW/CW/ER (no repeated K values).

lines <- c(
  tex_header(
    "Mean absolute error relative to the complete-data correlation $\\rho_0$ for NW, CW and ER over the nine settings. These methods do not depend on $K$. ER denotes the exhaustive or capped numerical reference.",
    "tab:supp_mae_reference",
    "ccrrr",
    colsep = "7pt"
  ),
  tex_row(c("$\\rho_{uv}$", "$\\eta$", "NW", "CW", "ER")),
  tex_cmd("midrule")
)

for (eta in cfg$eta) {
  if (eta != cfg$eta[1]) {
    lines <- c(lines, tex_cmd("addlinespace[3pt]"))
  }

  for (rho in cfg$rho_uv) {
    vals <- vapply(
      c("NW", "CW", "ER"),
      function(m) {
        sprintf("%.4f", get_s(rho, eta, 20L, m)$mae)
      },
      character(1)
    )

    lines <- c(lines, tex_row(c(sprintf("%.2f", rho), eta, vals)))
  }
}

write_tex(
  c(lines, tex_footer()),
  file.path(tab, "S_MAE_reference_values.tex")
)

# Supplement: EQ MAE for all K, 3 correlation columns, 3 retention blocks.

lines <- c(
  tex_header(
    "Mean absolute error of EQ relative to $\\rho_0$ at each rank-grid resolution $K$, across the nine generating settings.",
    "tab:supp_mae_eq",
    "crrr",
    colsep = "7pt"
  ),
  tex_row(c(
    "$K$",
    "$\\rho_{uv}=0.00$",
    "$\\rho_{uv}=0.25$",
    "$\\rho_{uv}=0.50$"
  )),
  tex_cmd("midrule")
)

for (eta in cfg$eta) {
  if (eta != cfg$eta[1]) {
    lines <- c(lines, tex_cmd("addlinespace[4pt]"))
  }

  lines <- c(lines, tex_row(paste0(
    tex_cmd("multicolumn{4}{l}{"),
    tex_cmd("textit{$"),
    tex_cmd("eta"),
    "=", eta, "$}}"
  )))

  for (K in cfg$K) {
    cells <- c(K, vapply(
      cfg$rho_uv,
      function(rho) {
        sprintf("%.4f", get_s(rho, eta, K, "EQ")$mae)
      },
      character(1)
    ))

    lines <- c(lines, tex_row(cells))
  }
}

write_tex(
  c(lines, tex_footer()),
  file.path(tab, "S_MAE_EQ_values.tex")
)

# Supplement: numerical values for pooled F1, by K and method.

lines <- c(
  tex_header(
    "Mean absolute error relative to $\\rho_0$, pooled equally across the nine simulation settings, underlying Figure~\\ref{fig:mae_pooled}.",
    "tab:supp_mae_pooled",
    "crrrr",
    colsep = "7pt"
  ),
  tex_row(c("$K$", methods)),
  tex_cmd("midrule")
)

for (K in cfg$K) {
  vals <- vapply(
    methods,
    function(m) {
      z <- pooled[
        pooled$K == K & pooled$method == m,
        ,
        drop = FALSE
      ]

      if (nrow(z) != 1L) {
        stop("Pooled MAE row missing")
      }

      sprintf("%.4f", z$pooled_mae)
    },
    character(1)
  )

  lines <- c(lines, tex_row(c(K, vals)))
}

write_tex(
  c(lines, tex_footer()),
  file.path(tab, "S_MAE_pooled_values.tex")
)
