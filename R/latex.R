# Minimal LaTeX writer used only for the paper's specified tables.
tex_slash <- intToUtf8(92L)
tex_cmd <- function(s) paste0(tex_slash, s)
tex_row <- function(cells) {
  paste0(
    paste(cells, collapse = " & "), " ",
    tex_slash, tex_slash
  )
}
tex_header <- function(caption, label, format, colsep = "3pt", size = "small",
                       environment = "table") {
  c(
    tex_cmd(paste0("begin{", environment, "}[!t]")),
    tex_cmd(paste0("caption{", caption, "}")),
    tex_cmd(paste0("label{", label, "}")),
    paste0(
      tex_cmd("centering"), tex_cmd(size),
      tex_cmd(paste0("setlength{", tex_cmd("tabcolsep"), "}{", colsep, "}"))
    ),
    tex_cmd(paste0("begin{tabular}{", format, "}")), tex_cmd("toprule")
  )
}
tex_footer <- function(environment = "table") {
  c(
    tex_cmd("bottomrule"), tex_cmd("end{tabular}"),
    tex_cmd(paste0("end{", environment, "}"))
  )
}
