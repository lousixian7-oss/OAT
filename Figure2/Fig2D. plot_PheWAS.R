# Install packages ----
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}
if(!"dplyr" %in% installed.packages()){install.packages('dplyr')}
if(!"tidyverse" %in% installed.packages()){install.packages('tidyverse')}
if(!"scales" %in% installed.packages()){install.packages('scales')}


# Load packages ----
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(scales)

# Database used ----
# https://azphewas.com/


# Plotting function ----
plot_PheWAS <- function(
    data = data,
    type = 'binary', # Data type: binary = Binary traits, cont = Continuous traits
    title = 'ATP2A1', # Figure title
    pic_type = 'png', # Figure type: png/tiff/pdf
    dir_save = 'result' # Output directory
){
  # Create result folder
  if(!dir.exists(dir_save)){dir.create(dir_save, recursive = T)}
  
  # Figure title
  title_name <- ifelse(type == 'binary',paste0('Binary traits PheWAS association with ',title),
                       paste0('Continuous traits PheWAS association with ',title))
  
  # Data formatting
  traits <- read.xlsx('traits.xlsx')
  merge_df <- merge(traits, data, by = 'Phenotypic.category')
  merge_df <- merge_df[order(merge_df$class,decreasing = F), ]
  merge_df$class <- factor(merge_df$class, levels = unique(merge_df$class))
  
  # Risk type
  if(type == 'binary'){
    merge_df <- merge_df %>%
      dplyr::mutate(risk = case_when(Odds.ratio > 1 ~ 'Odds.ratio > 1',
                                     Odds.ratio < 1 ~ 'Odds.ratio < 1'))
    shape_set <- c('Odds.ratio > 1' = 24, 'Odds.ratio < 1' = 25)
  }
  if(type == 'cont'){
    merge_df <- merge_df %>%
      dplyr::mutate(risk = case_when(Effect.size > 0 ~ 'Effect.size > 0',
                                     Effect.size < 0 ~ 'Effect.size < 0'))
    shape_set <- c('Effect.size > 0' = 24, 'Effect.size < 0' = 25)
  }
  
  # Identify traits beyond threshold
  merge_df <- merge_df %>%
    dplyr::mutate(sig_label = ifelse(P.value < 1e-08, merge_df$Phenotype,NA))
  
  # Save processed data
  write.xlsx(merge_df, paste0(dir_save,'/', title_name,"_",type,".xlsx"))
  
  
  # X-axis grouping
  merge_df$index <- 1:nrow(merge_df)
  axis_set <- merge_df %>% 
    group_by(class) %>% 
    summarise(center=(max(index)+min(index))/2)
  
  # Point colors
  color_set <- c('#BFD641', '#B01066','#BFD641', '#B01066','#BFD641', '#B01066',
                 '#BFD641', '#B01066','#BFD641', '#B01066','#BFD641', '#B01066',
                 '#BFD641', '#B01066','#BFD641', '#B01066','#BFD641', '#B01066',
                 '#BFD641', '#B01066','#BFD641', '#B01066','#BFD641', '#B01066')
  
  # Plot, format 1
  p <- ggplot(merge_df) +
    geom_point(aes(x = index, 
                   y = -log10(P.value), 
                   color = class, 
                   fill = class,
                   shape = risk),
               alpha = 0.7
    ) +
    geom_point(data = merge_df %>%
                 dplyr::filter(P.value < 1e-08),
               aes(x = index, y = -log10(P.value), 
                   color = class, 
                   fill = class,
                   shape = risk),
               show.legend = F, 
               color = "#000000") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, 
                                 angle = 90, 
                                 hjust = 1, 
                                 vjust = .5, 
                                 face = 'bold.italic'),
      axis.title.y = element_text(size = 20, face = 'bold.italic'),
      axis.text.y = element_text(size = 18),
      legend.title = element_blank(),
      legend.text = element_text(size = 15, face = 'italic'),
      legend.key.height = unit(1, 'cm'),
      panel.border = element_rect(linewidth = 1.5),
      panel.grid = element_blank(),
      plot.title = element_text(size = 22, hjust = .5, face = 'bold.italic')
    ) +
    geom_hline(yintercept = -log10(1e-08), linetype = 4, color = 'gray') +
    labs(title = title_name, x = '') +
    guides(color = 'none', 
           fill = 'none',
           shape = guide_legend(override.aes = list(size = 5))) +
    scale_x_continuous(
      labels = axis_set$class,
      expand = expansion(0.01, 0),
      breaks = axis_set$center,
      limits = c(min(merge_df$index),max(merge_df$index))
    ) +
    scale_shape_manual(values = shape_set) +
    scale_color_manual(values = pals::glasbey()) +
    scale_fill_manual(values = pals::glasbey()) 
  
  # Save figure
  ggsave(plot = p,
         device = pic_type,
         width = 12,
         height = 6,
         dpi = 300,
         paste0(dir_save,'/',title,"_",type,'_1.',pic_type))
  
  
  
  # Plot, format 2
  p <- ggplot(merge_df) +
    geom_point(aes(x = index, 
                   y = -log10(P.value), 
                   color = class, 
                   fill = class,
                   shape = risk),
               alpha = 0.7
    ) +
    geom_point(data = merge_df %>%
                 dplyr::filter(P.value < 1e-08),
               aes(x = index, y = -log10(P.value), 
                   color = class, 
                   fill = class,
                   shape = risk),
               show.legend = F, 
               color = "#000000") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 20, face = 'bold.italic'),
      axis.text.y = element_text(size = 18),
      legend.title = element_blank(),
      legend.text = element_text(size = 15, face = 'italic'),
      legend.key.height = unit(1, 'cm'),
      panel.border = element_rect(linewidth = 1.5),
      panel.grid = element_blank(),
      plot.title = element_text(size = 22, hjust = .5, face = 'bold.italic')
    ) +
    geom_hline(yintercept = -log10(1e-08), linetype = 4, color = 'gray') +
    labs(title = title_name, x = '') +
    guides(color = guide_legend(override.aes = list(size = 5), ncol = 2), 
           fill = 'none',
           shape = guide_legend(override.aes = list(size = 5))) +
    scale_x_continuous(
      labels = axis_set$class,
      expand = expansion(0.01, 0),
      breaks = axis_set$center,
      limits = c(min(merge_df$index),max(merge_df$index))
    ) +
    scale_shape_manual(values = shape_set) +
    scale_color_manual(values = pals::glasbey()) +
    scale_fill_manual(values = pals::glasbey()) 
  
  # Save figure
  ggsave(plot = p,
         device = pic_type,
         width = 15,
         height = 6,
         dpi = 300,
         paste0(dir_save,'/',title,"_",type,'_2.',pic_type))
}


# Load data ----
data <- read.csv("OAT_binary_0830T123455Z_AZ_PheWAS_Portal.csv")


# Call plotting function ----
plot_PheWAS(
  data = data,
  type = 'binary', # Data type: binary/cont; binary = Binary traits, cont = Continuous traits
  title = 'OAT', # Figure title
  pic_type = 'png', # Figure type: png/tiff/pdf
  dir_save = './plot' # Output directory
)
