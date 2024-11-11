library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(readr)
library(stringr)
library(GenomicRanges)
library(IRanges)

# Load and process QTL data
df_qtls <- fread("/path/to/map7_all_qtls_cross_control_peaks_pos.csv")

df_qtls <- df_qtls %>%
  arrange(chr, qtl_start) %>%
  group_by(chr) %>%
  mutate(index = row_number())  

# Adjust y-position for QTLs
df_qtls_adjusted <- df_qtls %>%
  mutate(chr = as.numeric(gsub("chr", "", chr))) %>% 
  group_by(chr) %>%
  mutate(index_within_chr = row_number()) %>%
  ungroup() %>%
  mutate(Chromosome = as.numeric(chr))

colnames(df_qtls_adjusted)[1] <- "QTL"

# Load and process combined sweep detection data
directory <- "/path/to/"
csv_files <- list.files(directory, pattern = "combined_.*minwin10000_maxwin100000_grid200kb.csv", full.names = TRUE)
data_list <- lapply(csv_files, read_csv)
combined_data <- bind_rows(data_list, .id = "source")
combined_data$source <- csv_files[as.numeric(combined_data$source)]

combined_data <- combined_data %>%
  mutate(source = basename(source)) %>%
  mutate(minwin = as.numeric(str_extract(source, "(?<=minwin)\\d+")),
         maxwin = as.numeric(str_extract(source, "(?<=maxwin)\\d+")),
         model = ifelse(str_detect(source, "ArabisNemo"), "ArabisNemo", "ArabisSagit"))

combined_data <- combined_data %>%
  mutate(chromosome = as.numeric(gsub("chr", "", chromosome)))

combined_data <- combined_data %>%
  mutate(Likelihood = pmin(Likelihood, 270))

merged_data <- combined_data %>%
  left_join(df_qtls_adjusted, by = c("chromosome" = "Chromosome"))

qtl_count <- merged_data %>%
  group_by(chromosome) %>%
  summarise(max_index_within_chr = max(index_within_chr, na.rm = TRUE))

merged_data <- merged_data %>%
  left_join(qtl_count, by = "chromosome") %>%
  mutate(y_position = 275 + (index_within_chr * 18),
         y_top = y_position + 20,
         max_y_limit = 320 + (max_index_within_chr * 15))

thresholds <- data.frame(
  model = "ArabisSagit",
  threshold = 21.08
)

# Colors for QTL phenotypes
phe_names <- c("Days to Bolting", "Days to Flowering", "Fertility Score", "Inflorescence Height",
               "Lamina Length", "Lamina L/W", "Leaf Length", "Leaf Width",
               "Number of Leaves", "Petal Length", "Petiole Length",
               "Rosette Diameter 1", "Rosette Diameter 2", "Rosette Diameter 3", "Rosette Diameter 4",
               "Side Shoots", "Stem Height", "Stem Leaf Density",
               "Stem Leaf Length", "Stem Leaf Width")

phe_colorp <- c("#b00058", "#ff0000", "#ffc000", "#7f7f7f",
                "#0070c0", "#004272", "#002060", "#87afff",
                "#266f8b", "#ff89c4", "#00b0f0",
                "#502273", "#7030a0", "#b17ed8", "#dac2ec",
                "#7fbd98", "#007b32", "#b9e08c",
                "#533b1e", "#cfa879")

colors2 <- setNames(phe_colorp, phe_names)

colors_sweep <- c("ArabisNemo" = "#B2182B", "ArabisSagit" = "#2166AC")

# Plotting
ggplot() +
  # sweep detection data as line plots
  geom_line(data = combined_data, aes(x = Position, y = Likelihood, color = model), size = 1) +
  # Add a horizontal threshold line
  geom_hline(data = thresholds, aes(yintercept = threshold), color = "purple", linetype = "dashed") +
  # Add QTLs as rectangle bars
  geom_rect(data = merged_data,
            aes(xmin = qtl_start + 1000, xmax = qtl_end - 1000, 
                ymin = y_position, ymax = y_top, fill = QTL),
            alpha = 0.6, inherit.aes = FALSE) +
  # Add black dot for peak positions of QTLs
  geom_point(data = merged_data, 
             aes(x = qtl_peak_pos, y = (y_position + y_top) / 2), 
             color = "black", size = 1.6, inherit.aes = FALSE) +  
  facet_wrap(~chromosome, scales = "free_x") +
  coord_cartesian(ylim = c(0, max(merged_data$max_y_limit))) +
  scale_color_manual(values = colors_sweep) +
  scale_fill_manual(values = colors2) +
  labs(x = "Position (bp)", y = "Likelihood", color = "Model", fill = "Phenotype",
       title = "Sweep detection and QTL overlap across chromosomes") +
  theme_minimal() +
  theme(legend.position = "bottom")

################## overlap QTL peak and sweep 200kb window ################


sweeps_above_threshold <- combined_data %>%
  filter(Likelihood > 21) %>%
  mutate(sweep_window_start = Position - 100000,
         sweep_window_end = Position + 100000)

sweeps_above_threshold <- sweeps_above_threshold %>%
  mutate(sweep_window_start = pmax(sweep_window_start, 0))

# overlaps between sweeps and QTL peaks
overlapping_sweeps <- sweeps_above_threshold %>%
  inner_join(df_qtls_adjusted, by = c("chromosome" = "chr")) %>%
  filter(sweep_window_start <= qtl_peak_pos & sweep_window_end >= qtl_peak_pos)
print(overlapping_sweeps)

# Plotting the overlaps between sweep windows and QTL peaks
ggplot() +
  # Add sweep windows as horizontal lines, colored by model
  geom_segment(data = sweeps_above_threshold, 
               aes(x = sweep_window_start, xend = sweep_window_end, 
                   y = chromosome, yend = chromosome, color = model), 
               size = 2) +
  # Add QTL peaks as points, colored by phenotype (phe)
  geom_point(data = overlapping_sweeps, 
             aes(x = qtl_peak_pos, y = chromosome, color = QTL),
             size = 4) +
  labs(x = "Position (bp)", y = "Chromosome", color = "Legend", 
       title = "Overlaps between sweep windows and QTL peaks") +
  theme_minimal() +
  scale_color_manual(values = c(colors_sweep, colors2)) +
  scale_y_continuous(breaks = seq(min(sweeps_above_threshold$chromosome), max(sweeps_above_threshold$chromosome), by = 1))



################## overlap 10% quantile QTL peak and sweep 200kb window ################


# Adjust QTLs to create the 10% quantile around the peak position
df_qtls_adjusted <- df_qtls %>%
  mutate(quantile_start = pmax(qtl_peak_pos - (qtl_peak_pos - qtl_start) * 0.10, qtl_start),
         quantile_end = pmin(qtl_peak_pos + (qtl_end - qtl_peak_pos) * 0.10, qtl_end))
# Adjust sweeps_above_threshold to include a window of + and - 100,000 bp around each position
sweeps_above_threshold <- combined_data %>%
  filter(Likelihood > 21) %>%
  mutate(sweep_window_start = Position - 100000,
         sweep_window_end = Position + 100000)

# Ensure the sweep window does not go below 0
sweeps_above_threshold <- sweeps_above_threshold %>%
  mutate(sweep_window_start = pmax(sweep_window_start, 0))

# overlaps between sweeps and the 10% quantile region of the QTLs
overlapping_sweeps <- sweeps_above_threshold %>%
  inner_join(df_qtls_adjusted, by = c("chromosome" = "chr")) %>%
  filter(sweep_window_start <= quantile_end & sweep_window_end >= quantile_start)
print(overlapping_sweeps)

# Plotting the overlaps between sweep windows and QTL quantile regions
ggplot() +
  # Add sweep windows as horizontal lines, colored by model
  geom_segment(data = sweeps_above_threshold, 
               aes(x = sweep_window_start, xend = sweep_window_end, 
                   y = chromosome, yend = chromosome, color = model), 
               size = 2) +
  # Add QTL quantile ranges as horizontal bars, colored by phenotype (phe)
  geom_rect(data = overlapping_sweeps, 
            aes(xmin = quantile_start, xmax = quantile_end, 
                ymin = chromosome + 0.1, ymax = chromosome + 0.3, fill = phe), 
            alpha = 0.6) +
  # Add a black vertical bar at the QTL peak position, spanning the height of the QTL bar
  geom_segment(data = overlapping_sweeps, 
               aes(x = qtl_peak_pos, xend = qtl_peak_pos, 
                   y = chromosome + 0.1, yend = chromosome + 0.3), 
               color = "black", size = 1) +
  labs(x = "Position (bp)", y = "Chromosome", color = "Species", fill = "Phenotype",
       title = "Overlaps between sweep windows and QTL 10% quantile regions") +
  theme_minimal() +
  scale_color_manual(values = c(colors_sweep)) +
  scale_fill_manual(values = colors2) +
  scale_y_continuous(breaks = seq(min(sweeps_above_threshold$chromosome), max(sweeps_above_threshold$chromosome), by = 1))


#### plots of 10% quantile but with overlapping sweeps on top of each other as offset ###
sweeps_above_threshold$chromosome <- as.character(sweeps_above_threshold$chromosome)
df_qtls_adjusted$chr <- as.character(df_qtls_adjusted$chr)

sweeps_gr <- GRanges(
  seqnames = sweeps_above_threshold$chromosome,
  ranges = IRanges(
    start = sweeps_above_threshold$sweep_window_start,
    end = sweeps_above_threshold$sweep_window_end
  )
)

qtls_gr <- GRanges(
  seqnames = df_qtls_adjusted$chr,
  ranges = IRanges(
    start = df_qtls_adjusted$quantile_start,
    end = df_qtls_adjusted$quantile_end
  )
)

# Find overlaps between QTLs and sweeps
overlaps <- findOverlaps(qtls_gr, sweeps_gr)
overlapping_qtls <- df_qtls_adjusted[unique(queryHits(overlaps)), ]

# Assign offsets to overlapping QTLs
overlapping_qtls <- overlapping_qtls %>%
  group_by(chr) %>%
  arrange(quantile_start) %>%
  mutate(offset = disjointBins(IRanges(start = quantile_start, end = quantile_end)) - 1) %>%
  ungroup()
overlapping_qtls$offset <- as.numeric(overlapping_qtls$offset)
overlapping_qtls$chr_numeric <- as.numeric(overlapping_qtls$chr)
sweeps_above_threshold$chromosome_numeric <- as.numeric(sweeps_above_threshold$chromosome)

# Plotting the overlaps with adjusted y-positions
ggplot() +
  # Add sweep windows as horizontal lines, colored by model (nemo and allo)
  geom_segment(data = sweeps_above_threshold, 
               aes(x = sweep_window_start, xend = sweep_window_end, 
                   y = chromosome_numeric, yend = chromosome_numeric, color = model), 
               size = 2) +
  
  # Add QTL quantile ranges as horizontal bars, adjusted by offset
  geom_rect(data = overlapping_qtls, 
            aes(xmin = quantile_start, xmax = quantile_end, 
                ymin = chr_numeric + 0.1 + 0.2 * offset, 
                ymax = chr_numeric + 0.3 + 0.2 * offset, fill = phe), 
            alpha = 0.6) +
  
  # Add a black vertical bar at the QTL peak position, adjusted by offset
  geom_segment(data = overlapping_qtls, 
               aes(x = qtl_peak_pos, xend = qtl_peak_pos, 
                   y = chr_numeric + 0.1 + 0.2 * offset, 
                   yend = chr_numeric + 0.3 + 0.2 * offset), 
               color = "black", size = 1) +
  labs(x = "Position (bp)", y = "Chromosome", color = "Species", fill = "Phenotype",
       title = "Overlaps between sweep windows and QTL 10% quantile regions") +
  theme_minimal() +
  scale_color_manual(values = colors_sweep) +
  scale_fill_manual(values = colors2) +
  scale_y_continuous(breaks = seq(min(sweeps_above_threshold$chromosome_numeric), 
                                  max(sweeps_above_threshold$chromosome_numeric), by = 1))

