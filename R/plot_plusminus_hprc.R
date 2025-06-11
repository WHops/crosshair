#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)
library(viridis)

# # Function to display usage information
# usage <- function() {
#   cat("Usage: Rscript plot_paf.R <paf_file>\n")
#   cat("Arguments:\n")
#   cat("  paf_file  Path to the PAF file\n")
#   quit(status = 1)
# }

# Function to read and process the PAF file
read_paf <- function(file) {
  lines <- readLines(file)
  non_header_lines <- lines[!grepl("^#", lines)]
  split_lines <- strsplit(non_header_lines, "\t")
  paf_mat <- do.call(rbind, lapply(split_lines, function(x) {
    length(x) <- 17
    return(x)
  }))
  
  paf_df <- as.data.frame(paf_mat)
  colnames(paf_df) <- c("query", "query_len", "query_start", "query_end", 
                        "strand", "target", "target_len", "target_start", 
                        "target_end", "num_matching_bases", "alignment_block_len", 
                        "mapq", "NM_string", '14', '15','16','TP')
  
  paf_df <- paf_df %>%
    mutate(
      query_len = as.numeric(query_len),
      query_start = as.numeric(query_start),
      query_end = as.numeric(query_end),
      target_len = as.numeric(target_len),
      target_start = as.numeric(target_start),
      target_end = as.numeric(target_end),
      num_matching_bases = as.numeric(num_matching_bases),
      alignment_block_len = as.numeric(alignment_block_len),
      queryname = str_extract(query, "^[^:]+"),
      querystart = as.numeric(str_extract(query, "(?<=:)[0-9]+(?=-)")),
      queryend = as.numeric(str_extract(query, "(?<=-)[0-9]+(?=/)")),
      querylength = as.numeric(str_extract(query, "(?<=/)[0-9]+")),
      nm = as.numeric(str_extract(NM_string, "(?<=NM:i:)[0-9]+")),
    )
  
  return(paf_df)
}

determine_wellmapping_utigs <- function(paf_data, min_fract_perfect, min_n_perfect) {
  valid_queries <- paf_data %>%
    group_by(queryname) %>%
    summarise(total = n(), zero_nm = sum(nm == 0)) %>%
    filter(zero_nm >= total * min_fract_perfect) %>%
    filter(zero_nm >= min_n_perfect) %>%
    pull(queryname)
  
  return(valid_queries)
}

# Function to highlight entries based on query start and end positions
annotate_utig_starts_ends <- function(paf_data) {
  paf_data <- paf_data %>%
    mutate(
      highlight_start = if_else(querystart == 1, "green", NA_character_),
      highlight_end = if_else(abs(querylength - queryend) < 100, "red", NA_character_)
    )
  return(paf_data)
}

# Function to plot the PAF data
plot_paf <- function(paf_data, xlabel='Target assembly', breakpoints=NULL, colorlim = 100) {
  start_end_width <- median(paf_data$query_len, na.rm = TRUE)
  
  p <- ggplot(paf_data, aes(x = target_start, y = strand)) +
    geom_tile(aes(width = queryend - querystart, fill = nm), color = NA) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, colorlim),
      oob = scales::squish,
      breaks = c(0, 1, 2, 3),
      labels = c("0", "1", "2", "3")
    ) +
    geom_tile(data = subset(paf_data, !is.na(highlight_start)), aes(width = start_end_width), fill = "green", color = NA) +
    geom_tile(data = subset(paf_data, !is.na(highlight_end)), aes(width = start_end_width), fill = "red", color = NA) +
    labs(
      x = paste0("Position along sequence ", xlabel),
      y = "Parental Utig",
      fill = "SNPs per window (1000 bp)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          panel.grid = element_blank())
  
  # Add vertical line for breakpoint if provided
  if (!is.null(breakpoints)) {
    for (breakpoint in breakpoints){
      p <- p + geom_vline(xintercept = breakpoint, linetype = "dashed", color = "red")
    }
  }
  
  return(p)
}


# Function to find edges using convolution
find_edges_convolution <- function(vector, window_size = 10) {
  half_window <- floor(window_size / 2)
  
  # Create the step function kernel
  kernel <- c(rep(-1, half_window), rep(1, half_window))
  
  # Apply convolution
  conv_result <- convolve(vector, rev(kernel), type = "filter")
  
  return(conv_result)
}

# Function to detect extrema
detect_extrema <- function(conv_result, window_size, detection_threshold) {
  max_value <- max(abs(conv_result))
  
  if (max_value > detection_threshold) {
    max_index <- which(abs(conv_result) == max_value)[1]  # Take the first occurrence
    return(max_index)
  } else {
    return(NULL)
  }
}

round_custom <- function(x, f){
  return(round(x/f) * f)
}







#paf_file = '../hg38_HG00621.2_10000.paf'
#paf_file = '../hprc_pafs/test/refto_hg38_HG005.2_y.fa_10000.paf'
#paf_file = '../pafs/HG01978.2_HG00621.2_1000.paf'
paf_file = '../hprc_pafs/to_hg38/refto_hg38_HG00621.2_y.fa_10000.paf'
paf_file = '../hprc_pafs/to_hg38/refto_hg38_HG01258.1_y.fa_10000.paf'
paf_file = '../hprc_pafs/refto_hg38_HG01928.1_y.fa_10000.paf'
paf_file = '../hprc_pafs/refto_hg38_HG005.2_y.fa_50000.paf'

#paf_file = '../hg38_HG00621.2_10000.paf'
#paf_file = '../HG01978.2_HG00621.2_1000.paf'
seqname_x = 'HPRC_testplot'
chunklen = 10000
crossover_detection_window_size <- 100
limit = 100
genetrack_link = '../../../data/genomes/hg38/15q13_golgas/15q13_GOLGA.bed'
sdtrack_link = '../../../data/segdups/hg38/simple_bed/15q13_ABC.bed'
chrx_startcoord = 25000000

paffolder = '../hprc_pafs/to_hg38/noinvs/'

fullthing = c()
fullthing_df = c()

n_files = 0
#Iterate through each paf file
for (paf_file in list.files(paffolder, pattern = ".paf$", full.names = TRUE)) {
  n_files = n_files + 1

  if (n_files > 20){
    next
  }
  # Part 1 
  extremum_detection_threshold <- (crossover_detection_window_size / 2) * 0.6
  paf_data <- read_paf(paf_file)

  # For every query, keep only the best alignment (num_matchin_bases maximized)
  # paf_data <- paf_data %>%
  #   group_by(query) %>%                      # Group by the query column
  #   slice_max(num_matching_bases, n = 1) %>% # Keep the row with the maximum value of num_matching_bases
  #   ungroup()
  paf_data = paf_data[paf_data$TP == 'tp:A:P' & paf_data$mapq > 0,]
    
  valid_queries <- determine_wellmapping_utigs(paf_data, 0.0, 0) #0.5, 5)#0.5, 5)
  paf_data <- paf_data %>% filter(queryname %in% valid_queries)
  paf_data <- annotate_utig_starts_ends(paf_data)
  paf_data = paf_data[order(paf_data$num_matching_bases, decreasing=F),]
  paf_data$threshold <- "no"
  #plot_paf(paf_data, seqname_x, colorlim = 10)

  # Part 2
  df = paf_data
  
  if (sum(paf_data[paf_data$target_start > 8000000,]$strand == '-') > sum(paf_data[paf_data$target_start > 8000000,]$strand == '+')){
    
    #print("flip")
    # flip the strand
    df = transform(df, strand = ifelse(strand == "+", "-", "+"))
    
  }
  
  df$target_start = round_custom(df$target_start, chunklen)
  df$nm =  df$num_matching_bases #chunklen -

  df <- df %>%
    group_by(strand, target_start) %>%
    mutate(nm = sum(nm)) %>%
    ungroup()

  df <- df %>%
    tidyr::complete(strand, target_start, fill = list(nm = 0))

  df = df %>% group_by(strand, target_start) %>% slice(1)

  df <- df %>%
    mutate(strand = recode(strand, '+' = 'plus', '-' = 'minus'))

  df$samples = paf_file
  fullthing_df = rbind(fullthing_df, df)
  df_merged <- df[,c('target_start','strand','nm')] %>%
    tidyr::pivot_wider(names_from = strand, values_from = nm, names_prefix = "nm_")# %>%
  df_merged$plus_aln = df_merged$nm_plus < df_merged$nm_minus 

  df_merged$diff = df_merged$nm_minus - df_merged$nm_plus
  df_merged[df_merged$diff > limit,]$diff = limit
  df_merged[df_merged$diff < -limit,]$diff = -limit

  df_merged$sample = paf_file
  
  fullthing = rbind(fullthing, df_merged)
  
}
# Part 3
genetrack = read.table(genetrack_link, header=F)
sdtrack = read.table(sdtrack_link, header=F)

colnames(genetrack) = c('chr','start','end','genename')
colnames(sdtrack) = c('chr','start','end','strand', 'SDname')

genetrack[c("start", "end")] <- genetrack[c("start", "end")] - chrx_startcoord
sdtrack[c("start", "end")] <- sdtrack[c("start", "end")] - chrx_startcoord

tip_off = 0#.99
ft = fullthing

unique_samples <- ft %>%
  dplyr::distinct(sample) %>%
  dplyr::mutate(y_pos = as.numeric(factor(sample)) * (chunklen*(1-tip_off)*2)) %>%
  group_by(sample) %>% slice(1)

unique_samples$sample_name <- sub(".*hg38_(.*?)_y\\.fa.*", "\\1", unique_samples$sample)
fullthing_df$sample_name <- sub(".*hg38_(.*?)_y\\.fa.*", "\\1", fullthing_df$samples)

# Define the order of samples and their corresponding genotypes
# sample_order <- data.frame(
#   sample_name = c("HG005.2", "HG00735.1", "HG00741.2", "HG01123.1", "HG01928.2", "HG01952.1", "HG00621.2", "HG01258.1", "HG01928.1", "HG01978.2", "HG02148.2"),
#   genotype = c("Inv_Del", "Inv_Del", "Inv_Del", "Inv_Del", "Inv_Del", "Inv_Del", "Inv_largest", "Inv_largest", "Inv_largest", "Inv_largest", "Inv_largest")
# )

# sample_order <- data.frame(
#   sample_name = c("HG00621.2", "HG01258.1", "HG01928.1", "HG01978.2", "HG02148.2"),
#   genotype = c("Inv_largest", "Inv_largest", "Inv_largest", "Inv_largest", "Inv_largest")
# )
# 
# sample_order <- data.frame(
#   sample_name = c("HG00621.2"),#, "HG01258.1", "HG01928.1", "HG01978.2", "HG02148.2"),
#   genotype = c("Inv_largest"))#, "Inv_largest", "Inv_largest", "Inv_largest", "Inv_largest")
# #)

# Merge the order into your existing data frames
# unique_samples <- unique_samples %>%
#   left_join(sample_order, by = "sample_name") %>%
#   arrange(genotype, factor(sample_name, levels = sample_order$sample_name))
# 
# fullthing_df <- fullthing_df %>%
#   left_join(sample_order, by = "sample_name") %>%
#   arrange(genotype, factor(sample_name, levels = sample_order$sample_name))
# 
# # Update the factor levels for sample_name to reflect this order
# unique_samples$sample_name <- factor(unique_samples$sample_name, levels = sample_order$sample_name)
# fullthing_df$sample_name <- factor(fullthing_df$sample_name, levels = sample_order$sample_name)
# 


q = ggplot()  +
  # Add the segments (lines) for the lollipop sticks
  geom_segment(data = fullthing_df, 
               aes(x = target_start, 
                   xend = target_start,  
                   y = as.numeric(factor(sample_name))*2 * (chunklen*(1-tip_off)*2),
                   yend = as.numeric(factor(sample_name))*2 * (chunklen*(1-tip_off)*2) + (as.numeric(as.factor(strand)) - 1.5) * 2 * pmax(0,(nm - (chunklen*tip_off))),
                   color = strand),
               
              size=1) +
  geom_text(data = unique_samples, 
            aes(x = 9000000, 
                y = as.numeric(factor(sample_name))*2 * (chunklen*(1-tip_off)*2),
                label = sample_name), 
            hjust = 1,
            vjust = 1,
            size = 5) +
  theme_minimal() +
  labs(x='[hg38] chr15:25000000+') +
  #ylim(0,5000) +
# + 
  # Overlay the difference as a line plot
  #geom_line(data = df_diff, aes(x = target_start, y = diff), color = "black", size = 1) +
  
  xlim(c(4000000, 9000000)) + 
  geom_rect(data=genetrack, aes(xmin=start, xmax=end, ymin= (-chunklen/100) -1000, ymax= (-(chunklen/100)+5000)), fill="blue", alpha=0.5) +
  geom_rect(data=sdtrack, aes(xmin=start, xmax=end, ymin= (-chunklen/100) -5000, ymax= (-(chunklen/100)-1000)), fill="black", alpha=0.5) +
  geom_rect(data=genetrack, aes(xmin=start, xmax=end, ymin= (-chunklen/100) +500000, ymax= (-(chunklen/100)-6000)), fill="blue", alpha=0.2)
  #geom_rect(data=sdtrack, aes(xmin=start, xmax=end, ymin= (-chunklen/100) +500000, ymax= (-(chunklen/100)-60)), fill="black", alpha=0.2)
q

