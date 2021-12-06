#!/usr/bin/env Rscript



data = "Long_results_RF.txt"
color_variable = "variable"
manual_color_vector = "#34495E,#BAC2C6,#F7DC6F"
yvariable = "value"
xvariable = "IS_Family"
melted = TRUE
title = "NULL"
x_label = "IS Family"
y_label = "Presence in 5% most important variables (%)"
color_variable_order = "NULL"
xvariable_order = "NULL"
y_add = 0
group_variable = "NULL"
add_point = FALSE
yaxis_scale_mode = "NULL"
facet_variable = "NULL"
facet_variable_order = "NULL"
stat = "identity"
bar_mode = "dodge"
facet_nrow = NULL
facet_ncol = NULL
error_bar_variable = "NULL"
base_font_size = 11
legend.position = "right"
xtics = TRUE
xtics_angle = 0
ytics = TRUE
facet_scales = "free"
coordinate_flip = FALSE
add_text = FALSE
add_bar_link = FALSE
font_path = "NULL"
width = 20
height = 9
outputprefix = "Figure_2_github"
outputpictype = "pdf"
saveppt = FALSE
savehtml = FALSE

library(ImageGP)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(tidyr)
library(htmlwidgets)
library(plotly)
if (data == "") {
  script = sub(".*=", "", commandArgs()[4])
  #print(script)
  system(paste(script, "-h"))
  stop("At least -f is required!")
}



if (outputprefix == "") {
  outputprefix = data
}

filename = paste0(outputprefix, '.bar.', outputpictype)



color_variable_order = sp_string2vector(color_variable_order)
xvariable_order = sp_string2vector(xvariable_order)
facet_variable_order = sp_string2vector(facet_variable_order)
manual_color_vector = sp_string2vector(manual_color_vector)

cat(sp_current_time(), "Starting...\n")

sp_barplot(
  data = data,
  color_variable = color_variable,
  yvariable = yvariable,
  xvariable = xvariable,
  melted = melted,
  title = title,
  x_label = x_label,
  y_label = y_label,
  color_variable_order = color_variable_order,
  xvariable_order = xvariable_order,
  y_add = y_add,
  group_variable = group_variable,
  add_point = add_point,
  yaxis_scale_mode = yaxis_scale_mode,
  facet_variable = facet_variable,
  stat = stat,
  bar_mode = bar_mode,
  facet_variable_order = facet_variable_order,
  facet_nrow = facet_nrow,
  facet_ncol = facet_ncol,
  error_bar_variable = error_bar_variable,
  base_font_size = base_font_size,
  legend.position = legend.position,
  xtics = xtics,
  xtics_angle = xtics_angle,
  ytics = ytics,
  manual_color_vector = manual_color_vector,
  facet_scales = facet_scales,
  coordinate_flip = coordinate_flip,
  add_text = add_text,
  add_bar_link = add_bar_link,
  font_path = font_path,
  height = height,
  width = width,
  filename = filename,
  saveppt = saveppt,
  savehtml = savehtml
)

