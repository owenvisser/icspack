
# Produce ONLY T2/T3 (DFS) and the two supplementary five-outcome tables.

stopifnot(requireNamespace("here", quietly = TRUE))

source(here::here("R", "common.R"))
source(here::here("R", "latex.R"))

d <- utils::read.csv(
  here::here("results", "nhanes", "correlations_delta.csv")
)

methods <- c("NW", "CW", "ER")
ages <- c("All", "30-49", "50-64", "65+")

folder <- here::here("manuscript", "generated", "nhanes")

value <- function(age, x, y, method) {
  z <- d[
    d$age_group == age &
      d$periodontal == x &
      d$caries == y &
      d$method == method,
    ,
    drop = FALSE
  ]

  if (nrow(z) != 1L) {
    stop("NHANES result missing or duplicated")
  }

  fmt_est_se(z$estimate, z$standard_error)
}

for (x in c("CAL", "PD")) {
  for (part in c("main", "supplement")) {
    y <- if (part == "main") "DFS" else c("DS", "DSI", "FS", "FSI", "DFSI")
    appendix <- identical(part, "supplement")

    caption <- if (!appendix) {
      sprintf(
        "MEC survey-weighted Pearson correlations between tooth-level %s and DFS, overall and by age. Entries are correlations (delta-method cluster-sandwich standard errors). NW, CW and ER are defined in the Methods. Standard errors treat patients as independent clusters, without stratified-PSU design adjustments.",
        x
      )
    } else {
      sprintf(
        "MEC survey-weighted Pearson correlations between tooth-level %s and additional caries measures, overall and by age. Entries are correlations (delta-method cluster-sandwich standard errors) under NW, CW and ER. Standard errors treat patients as independent clusters.",
        x
      )
    }

    label <- if (appendix) {
      paste0("tab:supp_nhanes_", tolower(x))
    } else {
      paste0("tab:nhanes_", tolower(x))
    }

    lines <- c(
      tex_header(
        caption,
        label,
        paste0("ll", paste(rep("r", length(y)), collapse = "")),
        colsep = if (appendix) "2pt" else "4pt",
        size = if (appendix) "scriptsize" else "small",
        environment = if (appendix) "table*" else "table"
      ),
      tex_row(c("Ages", "Weight", y)),
      tex_cmd("midrule")
    )

    for (age in ages) {
      if (age != ages[1]) {
        lines <- c(lines, tex_cmd("addlinespace[4pt]"))
      }

      for (j in seq_along(methods)) {
        method <- methods[j]

        agecell <- if (j == 1L) {
          tag <- if (age == "65+") {
            "$65+$"
          } else if (age == "30-49") {
            "30--49"
          } else if (age == "50-64") {
            "50--64"
          } else {
            "All"
          }

          tex_cmd(paste0("multirow[t]{3}{*}{", tag, "}"))
        } else {
          ""
        }

        lines <- c(lines, tex_row(c(
          agecell,
          method,
          vapply(
            y,
            function(outcome) {
              value(age, x, outcome, method)
            },
            character(1)
          )
        )))
      }
    }

    file <- if (appendix) {
      paste0("S_NHANES_", x, "_other_caries.tex")
    } else {
      paste0(if (x == "CAL") "T2" else "T3", "_NHANES_", x, "_DFS.tex")
    }

    write_tex(
      c(lines, tex_footer(environment = if (appendix) "table*" else "table")),
      file.path(folder, file)
    )
  }
}
