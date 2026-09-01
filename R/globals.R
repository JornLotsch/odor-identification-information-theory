# ==============================================================================
# globals.R
#
# Shared helper definitions and plotting utilities for the information-theoretic
# analysis of response-key structure in clinical olfactory tests.
#
# Global packages, variables and helper functions
# ============================================================================== 

# ---- Packages ----------------------------------------------------------------

REQUIRED_PACKAGES <- c(
  "ggplot2",
  "ggthemes",
  "pbmcapply",
  "patchwork",
  "ComplexHeatmap",
  "cABCanalysis",
  "grid"
)

invisible(lapply(REQUIRED_PACKAGES, library, character.only = TRUE))

# ---- Parameters --------------------------------------------------------------

SEED <- 42
TRIALS <- 100
TRAINING_PARTITION_SIZE <- 0.67
VALIDATION_PARTITION_SIZE <- 0.8
# No project-specific external helper sources are required here.
# The repository currently relies on the active analytical scripts and their
# shared helper definitions in this file alone.

# Plotting templates
# ======================== #
# Colors
# ======================== #

# Extended colorblind palette
cb_palette <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7"
)

# Extended colorblind palette with IoT as cornsilk shades
cornsilk1_palette <- c(
  "cornsilk1", "cornsilk2", "cornsilk3", "cornsilk4", "grey85",
  "lightgoldenrod1", "lightgoldenrod2", "lightgoldenrod3", "lightgoldenrod4",
  "gold1", "gold2", "gold3", "gold4",
  "lemonchiffon1", "lemonchiffon2", "lemonchiffon3", "lemonchiffon4"
)

cornsilk2_palette <- c(
  "cornsilk1", "cornsilk2", "cornsilk3", "cornsilk4", "grey85",
  "gold1", "gold2", "gold3", "gold4",
  "khaki1", "khaki2", "khaki3", "khaki4",
  "#DDD6B2", "#D2CC9F", "#C8C28C", "#BEB879"
)


breaks <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9)
nyt_colors <- c(
  "ghostwhite",
  "#fbfbfb",
  "#e6f0fa",
  "#c9def9",
  "#add0fa",
  "#7bb8fa",
  "dodgerblue2",
  "#041a58"
)



actual_palette <- cornsilk2_palette

#' Custom publication-quality ggplot2 theme
#'
#' @return ggplot2 theme object
theme_plot <- function() {
  ggplot2::theme_minimal(base_family = "Libre Franklin") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "plain", size = 12, color = "#222222",
        hjust = 0, margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(face = "plain", size = 10, color = "#444444"),
      axis.text = ggplot2::element_text(face = "plain", size = 10, color = "#444444"),
      plot.caption = ggplot2::element_text(
        size = 8, color = "#888888",
        hjust = 0, margin = ggplot2::margin(t = 10)
      ),
      panel.grid.major.y = ggplot2::element_line(
        color = "#dddddd", linetype = "dashed", size = 0.3
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks = ggplot2::element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks.length = grid::unit(5, "pt"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.direction = "vertical",
      plot.margin = ggplot2::margin(20, 20, 20, 20),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12, color = "#222222")
    )
}

