# Purpose: Shared plotting helpers; exact matching implementations only.
# Helpers are rebound to the calling stage scope to retain original lexical lookup.
# 仅合并实现完全相同的函数；作用域保留在各独立阶段。

shared_save_figure_pair_028 <- function(plot, stem, width, height, dpi = 600L) {
    dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
    ggplot2::ggsave(paste0(stem, ".tiff"), plot = plot, width = width, height = height, units = "in", dpi = dpi, compression = "lzw")
    ggplot2::ggsave(paste0(stem, "_preview.png"), plot = plot, width = width, height = height, units = "in", dpi = 180)
}
