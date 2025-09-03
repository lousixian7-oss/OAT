library(forestploter)
library(readr)

data = read_csv("res.csv")
data$' ' <- paste(rep(" ", 15), collapse = " ") 

data$'OR(95% CI)' <- ifelse(is.na(data$or), "",
                            sprintf("%.2f (%.2f to %.2f)",
                                    data$or,
                                    data$or_lci95, data$or_uci95))

data$P = ifelse(is.na(data$P), "", data$P)

data$Exposure = ifelse(is.na(data$Exposure), "", data$Exposure)
data$`No.of SNP` = ifelse(is.na(data$`No.of SNP`), "", data$`No.of SNP`)

# Forest plot theme  You can adjust various styles here
tm <- forest_theme(base_size = 10, # Base size
                   base_family = "serif", # Font (use windowsFonts() to check available fonts)
                   # Shape, line type, color, and width of confidence interval points
                   ci_pch = 16,
                   ci_col = "#4575b4", # #762a83
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Add short vertical lines at both ends of CI
                   
                   # Reference line width, type, and color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   
                   # Fill and border color of summary diamonds
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   
                   # Footnote size, font, and color
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "blue",
)

plot <- forest(
  # Required columns (default unchanged). Adjust here to display the data you need.
  # Check the data frame in the environment panel and fill in the correct column indices
  data[, c(1,2,3,4,9,8)], 
  # OR values and their corresponding confidence intervals
  est = data$or, 
  lower = data$or_lci95,
  upper = data$or_uci95,
  ci_column = 5,  # The empty column defined earlier is used to draw the forest plot. Modify as needed.
  ref_line = 1,  # Reference line. Change to 0 if using beta values
  # Range of the interval display
  xlim = c(0.80, 1.30), 
  # Tick marks for the x-axis. Remove the # before ticks_at to enable
  ticks_at = c(0.80,0.90, 1, 1.10,1.20,1.30), 
  theme = tm,
)
# Plotting. You can change p1.pdf to another filename
# You can also adjust width and height to control the figure size
pdf("p1.pdf", family = "GB1", width = 10, height = 8)
plot
dev.off()
