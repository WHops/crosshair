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
    length(x) <- 13
    return(x)
  }))
  
  paf_df <- as.data.frame(paf_mat)
  colnames(paf_df) <- c("query", "query_len", "query_start", "query_end", 
                        "strand", "target", "target_len", "target_start", 
                        "target_end", "num_matching_bases", "alignment_block_len", 
                        "mapq", "NM_string")
  
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
      nm = as.numeric(str_extract(NM_string, "(?<=NM:i:)[0-9]+"))
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
plot_paf <- function(paf_data, xlabel='Target assembly', breakpoints=NULL) {
  start_end_width <- median(paf_data$query_len, na.rm = TRUE)
  
  p <- ggplot(paf_data, aes(x = target_start, y = queryname)) +
    geom_tile(aes(width = queryend - querystart, fill = nm), color = NA) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 1),
      oob = scales::squish,
      breaks = c(0, 1, 2, 3),
      labels = c("0", "1", "2", "3")
    ) +
    geom_tile(data = subset(paf_data, !is.na(highlight_start)), aes(width = start_end_width), fill = "green", size=3, color = NA) +
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

n=6
paf_file = paste0('../res/trio',n,'_paf/trio',n,'.paf')
seqname_x = paste0('Trio',n,'_child')
crossover_detection_window_size <- 100
extremum_detection_threshold <- (crossover_detection_window_size / 2) * 0.8

paf_data <- read_paf(paf_file)
#target='h1tg0000004l'
#query='utg000217l'
#paf_data = paf_data[paf_data$target == target,]
#paf_data = paf_data[paf_data$queryname == query,]
valid_queries <- determine_wellmapping_utigs(paf_data, 0.2,1) #0.5, 5)#0.5, 5)
length(valid_queries)
paf_data <- paf_data %>% filter(queryname %in% valid_queries)
paf_data <- annotate_utig_starts_ends(paf_data)

# sort paf by nm
paf_data = paf_data[order(paf_data$nm, decreasing=T),]
paf_data$threshold <- "no"
breakpoints = c()
plot_paf(paf_data, seqname_x)


all_utigs <- unique(paf_data$queryname)
for (utig_to_analyse in all_utigs){
  # Filter, group by querystart, and sort the dataframe
  filtered_data <- paf_data %>%
    filter(queryname == utig_to_analyse) %>%
    group_by(querystart) %>%
    slice_min(nm, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(querystart)

  # Extract the 'nm' column as a vector
  nm_vector <- filtered_data$nm
  nm_vector[nm_vector > 1] <- 1

  # Find the edges using the convolution function
  conv_result <- find_edges_convolution(nm_vector, crossover_detection_window_size)

  # Detect extrema
  extrema_index <- detect_extrema(conv_result, crossover_detection_window_size, extremum_detection_threshold)

  # if (utig_to_analyse == 'utg000153l'){
  #   # Plotting
  #   plot(1:length(nm_vector), nm_vector, type = "p", col = "blue", main = "NM Vector", xlab = "Index", ylab = "NM")
  #   if (!is.null(extrema_index)) {
  #     abline(v = extrema_index + (crossover_detection_window_size/2), col = "red", lty = 2)
  #   }
  # 
  #   plot(1:length(conv_result), conv_result, type = "l", col = "green", main = "Convolution Result", xlab = "Index", ylab = "Convolution")
  #   if (!is.null(extrema_index)) {
  #     abline(v = extrema_index, col = "red", lty = 2)
  #   }
  #   
  # }
  # Translate to query and target coordinates
  if (!is.null(extrema_index)) {
    breakpoint_index <- extrema_index + floor(crossover_detection_window_size / 2)
    breakpoints <- c(breakpoints, filtered_data$target_start[breakpoint_index])
    
    # Add a 'threshold' column to the original paf_data
    paf_data$threshold[paf_data$queryname == utig_to_analyse & paf_data$querystart == filtered_data$querystart[breakpoint_index]] <- "yes"
  }
}

# Plotting the PAF data with the inferred breakpoint
plot_paf(paf_data, xlabel=seqname_x, breakpoints=breakpoints)

print(paf_data[paf_data$threshold == 'yes',])

get_query_coord <- function(paf, target_seq, target_pos,query, pm_per_side) {
  hits <- paf[startsWith(paf$queryname, query), ]
  hits <- hits[grepl(target_seq, hits$target), ]
  hits <- hits[hits$num_matching_bases > 995,]
  if (nrow(hits) == 0) return(NULL)
  i <- which.min(abs((hits$target_start + hits$target_end) / 2 - target_pos))
  h <- hits[i, ]
  #prop <- (target_pos - h$target_start) / (h$target_end - h$target_start)
  #rel_query <- prop * (h$query_end - h$query_start)
  
  region_to_return=c(h$querystart - pm_per_side, h$queryend + pm_per_side)
  #qc <- if (h$strand == "+") h$query_start + rel_query else h$queryend - rel_query
  return(paste0(query, ':', region_to_return[1], '-', region_to_return[2],'_', h$strand))
}

#breakpoints=c(56375)
target='seq1'
# Breakpoints
#breakpoints=c(106384)
target='h1tg000004l'
extract_windowsize = 50000
pm_per_side = (extract_windowsize - 1000) / 2
# Match and extract
matched <- grep(target, paf_data$target, value = TRUE)[1]
start_num <- as.numeric(sub(paste0(".*", target, ":(\\d+)-.*"), "\\1", matched))
start_num = 0

paste0(target, ':', (breakpoints[1]+start_num)-pm_per_side, '-', (breakpoints[1]+start_num+1000)+pm_per_side, '_+')
get_query_coord(paf_data, target, breakpoints[1], "utg000348l", pm_per_side)
get_query_coord(paf_data, target, breakpoints[1], "utg000354l", pm_per_side)


