data(extraversion, package = "diffIRT")

Extra <- data.frame(
  sbj = seq_len(nrow(extraversion)),
  setNames(extraversion[, 1:10],  paste0("Item", 1:10, "_resp")),
  setNames(extraversion[, 11:20], paste0("Item", 1:10, "_rt"))
)

Extra <- tidyr::pivot_longer(
  Extra,
  cols = -sbj,
  names_to = c("item", ".value"),
  names_pattern = "(Item\\d+)_(resp|rt)"
)

Extra$item <- factor(Extra$item)
