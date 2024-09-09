PlotFlankingMethylation <- function(path_to_bismark_cov = NULL,
                                    path_to_fimo_file = NULL,
                                    path_to_meme_file = NULL,
                                    motif_name = NULL,
                                    motif_symbol = NULL,
                                    expand_region = 250,
                                    sample_sites = 20){
  
  if(is.null(path_to_bismark_cov) | is.null(path_to_fimo_file) | is.null(path_to_meme_file) | is.null(motif_name)){
    stop("Please provide all the paths to the files and the MEME motif name (i.e. MEME-1, MEME-2)")
  }
  if(!is.null(sample_sites)){
    set.seed(10044)
    message("Random seed fixed for reproducible sampling process (seed number = 10044)")
  }
  # - Load bismark methylation output
  bismark.cov <- suppressMessages(read_delim(path_to_bismark_cov, #.cov file from Bismark
                                             delim = "\t", escape_double = FALSE, 
                                             col_names = FALSE, trim_ws = TRUE))
  
  # - Load FIMO motif positions
  motif.sites <-suppressMessages(suppressWarnings(read_delim(path_to_fimo_file, # fimo.tsv file from MEME
                                                             delim = "\t", escape_double = FALSE, 
                                                             trim_ws = TRUE)))
  
  motif.sites = motif.sites[!is.na(motif.sites$sequence_name),]
  
  bismark.cov = GRanges(seqnames = bismark.cov$X1, ranges = IRanges(start = bismark.cov$X2, end = bismark.cov$X3), strand = "*",
                        m.percent = bismark.cov$X4, m.reads = bismark.cov$X5, um.reads = bismark.cov$X6)
  motif.sites = GRanges(seqnames = motif.sites$sequence_name, ranges = IRanges(start = motif.sites$start, end = motif.sites$stop),
                        strand = motif.sites$strand, motif.sequence = motif.sites$matched_sequence, qval = motif.sites$`q-value`)
  
  motif.sites = motif.sites[!is.na(motif.sites$motif.sequence),]
  
  # - Load motif matrix
  motif.pwm  <- as.data.frame(suppressMessages(read_delim(path_to_meme_file, # meme.txt file
                                                          delim = "/t", escape_double = FALSE, 
                                                          trim_ws = TRUE)))
  
  motif.name = paste0(motif_name," position-specific probability matrix")
  motif.length = as.numeric(strsplit(motif.pwm[(grep(motif.pwm[,1],pattern =  motif.name)+2),], split = "w= | nsites")[[1]][2])
  pwm <- motif.pwm[(grep(motif.pwm[,1],pattern = motif.name)+3):((grep(motif.pwm[,1],pattern =  motif.name)+2+(motif.length))),] 
  pwm = do.call(rbind,strsplit(pwm,"  "))
  pwm = t(apply(pwm,2,as.numeric))
  rownames(pwm) = c("A","C","G","T")
  
  # - Check proportion of motifs showing CpG sites
  message("Percentage of motifs with CpG sites: ", round(length(motif.sites$motif.sequence[grep(motif.sites$motif.sequence, pattern = "CG", ignore.case = T)]) / length(motif.sites$motif.sequence) * 100, digits = 2), "%")
  
  # - Check proportion of CpGs contained within motifs
  cpgs.in.motifs = round(length(subsetByOverlaps(bismark.cov,motif.sites)) / length(bismark.cov) * 100, digits = 2)
  message("Percentage of CpGs within motifs: ", cpgs.in.motifs, "%")
  message("Percentage of CpGs outside motifs: ", 100-cpgs.in.motifs,"%") 
  
  # - Check methylation status at CpG sites surrowding the motifs
  message("Expanding motif sites to ",expand_region," bp upstream and downstrem (region size = ",expand_region*2,")")
  expand.region = resize(motif.sites, width = width(motif.sites)+(expand_region*2), fix = "center")
  message("Obtaining CpG sites for measurement across all motif sites...")
  bulk.region.methylation = subsetByOverlaps(bismark.cov,expand.region)
  
  message("Obtaining CpG sites for each individual motif site...")
  per.site.region.methylation = lapply(1:length(expand.region), function(x){
    subsetByOverlaps(bismark.cov,expand.region[x,])
  })
  
  message("Calculating distances from motif center to CpG sites...")
  per.site.region.methylation = lapply(1:length(expand.region), function(x){
    per.site.region.methylation[[x]]$distance = per.site.region.methylation[[x]]@ranges@start - (motif.sites[x,]@ranges@start + ceiling(motif.sites[x,]@ranges@width/2)) 
    return(per.site.region.methylation[[x]])
  })
  
  message("Merging GRanges with CpG distances calculated...")
  bulk = do.call(c, per.site.region.methylation)
  dat = as.data.frame(bulk@elementMetadata)
  
  message("Getting mean methylation values for each distance from motif...")
  mean.dat = as.data.frame(do.call(rbind,lapply(unique(dat$distance), function(x){
    c(mean(dat[dat$distance == x,]$m.percent), x)
  })))
  colnames(mean.dat) = c("m.percent","distance")  
  
  message("Sampling examples for flanking methylation at specific sites...")
  filtered.per.site.region.methylation = lapply(per.site.region.methylation, function(x){
    x = x[x$m.reads+x$um.reads > 4,] # Subset by minimum coverage
    if(length(x) > 1){ # Subset by number of CpGs 
      return(x)  
    }else{
      return(NA)
    }
  })
  
  filtered.per.site.region.methylation = filtered.per.site.region.methylation[!is.na(filtered.per.site.region.methylation)]
  message("Sites remaining: ",length(filtered.per.site.region.methylation))
  
  if(length(filtered.per.site.region.methylation) > sample_sites){
    filtered.per.site.region.methylation = filtered.per.site.region.methylation[sample(1:length(filtered.per.site.region.methylation), size = sample_sites)]
  }
  
  filtered.per.site.region.methylation = lapply(1:length(filtered.per.site.region.methylation), function(x){
    filtered.per.site.region.methylation[[x]]$site = x
    return(filtered.per.site.region.methylation[[x]])
  })
  names(filtered.per.site.region.methylation) = paste0("Site::",1:length(filtered.per.site.region.methylation))
  
  filtered.sites = do.call(rbind,lapply(filtered.per.site.region.methylation, function(x) as.data.frame(x@elementMetadata)))
  
  message("Plotting...")
  # Visualization
  message("Plotting mean methylation around the motif...")
  
  p0 = ggseqlogo(pwm) + 
    labs(title = motif_symbol) + 
    scale_x_continuous(expand = c(0,0), labels = 1:motif.length, breaks = 1:motif.length) + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25, 0, 0, 0), "cm")) 
  
  p0.5 = ggplot(mean.dat) +
    geom_abline(slope = 1, intercept = -20, lty = 2) +
    geom_abline(slope = -1, intercept = -20, lty = 2) +
    xlim(-expand_region,expand_region) +
    ylim(0,expand_region) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))

  p1 = ggplot(mean.dat, aes(x = distance, y = m.percent)) +
    annotate("rect", xmin=-ceiling(motif.sites@ranges@width[1]/2),
             xmax=motif.length - ceiling(motif.sites@ranges@width[1]/2),
             ymin=0,
             ymax=max(mean.dat$m.percent),
             alpha=0.2,
             fill="red") +
    geom_point(size = 1) + 
    geom_smooth() + 
    labs(x = "", y = "CpG methylation (%)") +
    geom_vline(xintercept = c(-ceiling(motif.sites@ranges@width[1]/2),(motif.length - ceiling(motif.sites@ranges@width[1]/2))), lty=2) + 
    theme_classic()+
    xlim(-expand_region,expand_region) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    # Replace the value of sample sites if there are not enough, to keep the scale consistent in plot
    if(length(filtered.per.site.region.methylation) <= sample_sites){
      sample_sites = length(filtered.per.site.region.methylation)
    }
    
  p2 = ggplot(filtered.sites, aes(x = distance, y = site, fill = m.percent)) + 
    geom_hline(yintercept = filtered.sites$site, lwd = 0.1) +
    annotate("rect", xmin=-ceiling(motif.sites@ranges@width[1]/2),
             xmax=motif.length - ceiling(motif.sites@ranges@width[1]/2),
             ymin=1,
             ymax=sample_sites,
             alpha=0.2,
             fill="red") +
    geom_point(data = filtered.sites[filtered.sites$m.percent == 0,], pch=21, size = 2, lwd = 0.1) +
    geom_point(data = filtered.sites[filtered.sites$m.percent > 0,], pch=21, size = 2, lwd = 0.1) +
    scale_fill_gradient(low = "white", high = "black") +
    theme_classic() +
    labs(x = "Distance to motif center (bp)", y = "Motif site example (#)") +
    xlim(-expand_region,expand_region) +
    geom_vline(xintercept = c(-ceiling(motif.sites@ranges@width[1]/2),(motif.length - ceiling(motif.sites@ranges@width[1]/2))), lty=2) + 
    theme(legend.position = "none")
  
  pg = plot_grid(p0,p0.5, p1, p2, ncol = 1, align = 'v', rel_heights = c(7,4,15,7.5))  
  message("Done!")
  pg
}

#### Examples
# PlotFlankingMethylation(path_to_bismark_cov = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup.bismark.cov.gz",
#                         path_to_fimo_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_1/fimo.tsv",
#                         path_to_meme_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
#                         motif_name = "MEME-1",
#                         motif_symbol = "CTCF",
#                         expand_region = 250,
#                         sample_sites = 20)
# 
# PlotFlankingMethylation(path_to_bismark_cov = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup.bismark.cov.gz",
#                         path_to_fimo_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
#                         path_to_meme_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
#                         motif_name = "MEME-2",
#                         motif_symbol = "GATA1",
#                         expand_region = 250,
#                         sample_sites = 20)
# 
