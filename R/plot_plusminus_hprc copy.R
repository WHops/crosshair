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
  
  p <- ggplot(paf_data, aes(x = target_start, y = strand)) +
    geom_tile(aes(width = queryend - querystart, fill = nm), color = NA) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 100),
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
paf_file = '../hg38_HPRC_10000.paf'
#paf_file = '../hg38_HG00621.2_10000.paf'
#paf_file = '../HG01978.2_HG00621.2_1000.paf'
seqname_x = 'HPRC_testplot'
chunklen = 10000
crossover_detection_window_size <- 100
limit = 20
extremum_detection_threshold <- (crossover_detection_window_size / 2) * 0.6

paf_data <- read_paf(paf_file)
valid_queries <- determine_wellmapping_utigs(paf_data, 0.0, 1) #0.5, 5)#0.5, 5)
paf_data <- paf_data %>% filter(queryname %in% valid_queries)
paf_data <- annotate_utig_starts_ends(paf_data)

# sort paf by nm
paf_data = paf_data[order(paf_data$nm, decreasing=T),]
paf_data$threshold <- "no"
breakpoints = c()
#plot_paf(paf_data, seqname_x)

# Plotting the PAF data with the inferred breakpoint
#plot_paf(paf_data, xlabel=seqname_x, breakpoints=breakpoints)


df = paf_data

df$target_start = round_custom(df$target_start, chunklen)

df$nm = chunklen - df$num_matching_bases
# Assuming your dataframe is named df
df <- df %>%
  # Group by strand and targetstart
  group_by(strand, target_start) %>%
  # Filter for the minimum 'nm' value for each group
  filter(nm == min(nm)) %>%
  # Ungroup to allow further operations
  ungroup()

# Now, handle missing combinations by completing the dataframe
df <- df %>%
  # Complete all combinations of strand and targetstart
  tidyr::complete(strand, target_start, fill = list(nm = 10000))

df = df %>% group_by(strand, target_start) %>% slice(1)

# View the final dataframe
print(df[,c('strand','target_start','nm')])


#ggplot(df) + geom_point(aes(x=target_start, y=-nm, color=strand), size=0.00000001) +
#  ylim(c(-50,0))


df2 = df[,c('target_start','strand','nm')]

df2[df2$strand == '+','strand'] = 'plus'
df2[df2$strand == '-','strand'] = 'minus'

# Assuming your dataframe is named df
df_merged <- df2 %>%
  tidyr::pivot_wider(names_from = strand, values_from = nm, names_prefix = "nm_")# %>%
df_merged$plus_aln = df_merged$nm_plus < df_merged$nm_minus 

df_merged$diff = df_merged$nm_minus - df_merged$nm_plus
df_merged[df_merged$diff > limit,]$diff = limit
df_merged[df_merged$diff < -limit,]$diff = -limit



genetrack_link = '../../../data/genomes/hg38/15q13_golgas/15q13_GOLGA.bed'
sdtrack_link = '../../../data/segdups/hg38/simple_bed/15q13_ABC.bed'

# Load both
genetrack = read.table(genetrack_link, header=F)
sdtrack = read.table(sdtrack_link, header=F)

colnames(genetrack) = c('chr','start','end','genename')
colnames(sdtrack) = c('chr','start','end','strand', 'SDname')

chrx_startcoord = 28000000

genetrack$start = genetrack$start - chrx_startcoord
genetrack$end = genetrack$end - chrx_startcoord
sdtrack$start = sdtrack$start - chrx_startcoord
sdtrack$end = sdtrack$end - chrx_startcoord


p = ggplot()  +
  # Add the segments (lines) for the lollipop sticks
  geom_segment(data = df_merged, 
               aes(x = target_start, 
                   xend = target_start, 
                   y = 0, 
                   yend = diff,
                   color = plus_aln), 
               size = 1) +  # Adjust the line thickness
  # Add the points at the end of the segments
  geom_point(size = 0.1) +  # Adjust the point size
  # Add titles and labels
  labs(
    x = "Target Start",
    y = "nm_plus - nm_minus",
    color = "Pos-aln") +
  # Enhance the overall look with a minimal theme
  theme_minimal(base_size = 15) +
  # Optional: Rotate x-axis labels for better readability
  theme(axis.text.x = element_text(angle = 45, hjust = 1))# + 
  #ylim(c(-40,40))# +
# Optional: Customize the color palette
#scale_color_brewer(palette = "Set1")



# Plot the gene track on top of p
p = p + geom_rect(data=genetrack, aes(xmin=start, xmax=end, ymin=26, ymax=30), fill="blue", alpha=0.2) +
  geom_rect(data=sdtrack, aes(xmin=start, xmax=end, ymin=22, ymax=25), fill="black", alpha=0.2)

p
