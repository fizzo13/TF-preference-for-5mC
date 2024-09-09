# meNTT
# Pseudocode
# - Load bismark methylation output
# - Load FIMO motif positions
# - Load motif matrix
# - Check proportion of motifs showing CpG sites
# - Check methylation status on CpG sites
# - Visualization

library(readr)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(ggseqlogo)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(Rsamtools)

setwd("/gpfs/commons/home/fizzo/")
source("/gpfs/commons/home/fizzo/meNTT/Functions/PlotMotifMethylation.R")

PlotMotifMethylation(path_to_bismark_cov = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup.bismark.cov.gz",
                     path_to_fimo_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
                     path_to_meme_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
                     motif_name = "MEME-2",
                     motif_symbol = "GATA1")

PlotMotifMethylation(path_to_bismark_cov = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup.bismark.cov.gz",
                     path_to_fimo_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_1/fimo.tsv",
                     path_to_meme_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
                     motif_name = "MEME-1",
                     motif_symbol = "CTCF",
                     sample_sites = NULL)                   

PlotMotifMethylation(path_to_bismark_cov = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup.bismark.cov.gz",
                     path_to_fimo_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_1/fimo.tsv",
                     path_to_meme_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
                     motif_name = "MEME-1",
                     motif_symbol = "CTCF",
                     sample_sites = 20)

PlotMotifMethylation(path_to_bismark_cov = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/AL02_CTCF_EMseqMERGED_1_S4_L001_R1_001.trimmed_bismark_bt2_pe.sam.bismark.cov.gz",
                     path_to_fimo_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_1/fimo.tsv",
                     path_to_meme_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
                     motif_name = "MEME-2",
                     motif_symbol = "CTCF",
                     sample_sites = 20)

source("/gpfs/commons/home/fizzo/meNTT/Functions/MethylationConcordancePerDistance.R")
MethylationConcordancePerDistance(path_to_bam = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/sorted.bam.gz",
                                  path_to_bam_index = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/sorted.bam.gz.bai",
                                  path_to_fimo_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
                                  motif_symbol = "CTCF",
                                  extended_region = 200,
                                  quantile_number = 99,
                                  ncores = 18)

source("/gpfs/commons/home/fizzo/meNTT/Functions/PatternClassifier.R")
pats = PatternClassifier(path_to_bam = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/sorted.bam.gz",
                         path_to_bam_index = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/sorted.bam.gz.bai",
                         path_to_fimo_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
                         path_to_meme_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
                         motif_symbol = "CTCF",
                         motif_name = "MEME-2",
                         top = 10, 
                         ncores = 18)
for(i in names(pats)){
  plot(pats[[i]])
}

source("/gpfs/commons/home/fizzo/meNTT/Functions/PlotFlankMethylation.R")
PlotFlankingMethylation(path_to_bismark_cov = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/AL02_CTCF_EMseqMERGED_1_S4_L001_R1_001.trimmed_bismark_bt2_pe.sam.bismark.cov.gz",
                        path_to_fimo_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
                        path_to_meme_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
                        motif_symbol = "CTCF",
                        motif_name = "MEME-2",
                        expand_region = 2000,
                        sample_sites = 20
)

source("/gpfs/commons/home/fizzo/meNTT/Functions/PatternHeatmap.R")
PatternHeatmap(path_to_bam = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/sorted.bam",
               path_to_bam_index = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/sorted.bam.bai",
               path_to_fimo_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
               path_to_meme_file = "meNTT/20230523_PBMC/Landau-WC-14411_2023_06_06/Merged_CTCF_EMseq/peak_calling_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
               motif_symbol = "CTCF",
               motif_name = "MEME-2")
