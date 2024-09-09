# Percentage of discordant reads as function of distance to motif
MethylationConcordancePerDistance<- function(path_to_bam = NULL,
                                             path_to_bam_index = NULL,
                                             path_to_fimo_file = NULL,
                                             motif_symbol,
                                             extended_region = 500,
                                             quantile_number = 5,
                                             ncores = 1){
  
  if(is.null(path_to_bam) | is.null(path_to_fimo_file) | is.null(path_to_bam_index)){
    stop("Please provide all the paths to the requested files")
  }
  
  set.seed(10044)
  
  message("Reading in .bam file...")
  bamFile <- BamFile(file = path_to_bam,
                     index = path_to_bam_index)
  bam <- scanBam(bamFile) # Load the bam file into R
  sequence <- bam[[1]]$seq # Extract the read sequence
  methylation <- scanBam(bamFile, param = ScanBamParam(tag = "XM")) # Load the methylation calls per read
  
  # - Load FIMO motif positions
  message("Loading motif site information...")
  motif.sites <-suppressMessages(suppressWarnings(read_delim(path_to_fimo_file, # fimo.tsv file from MEME
                                                             delim = "\t", escape_double = FALSE, 
                                                             trim_ws = TRUE)))
  motif.sites = motif.sites[!is.na(motif.sites$sequence_name),]
  motif.sites = GRanges(seqnames = motif.sites$sequence_name, ranges = IRanges(start = motif.sites$start, end = motif.sites$stop),
                        strand = motif.sites$strand, motif.sequence = motif.sites$matched_sequence, qval = motif.sites$`q-value`)
  motif.sites = motif.sites[!is.na(motif.sites$motif.sequence),]
  
  read.positions = GRanges(seqnames = bam[[1]]$rname,
                           ranges = IRanges(start = bam[[1]]$pos, 
                                            end = bam[[1]]$pos + bam[[1]]$qwidth),
                           strand = bam[[1]]$strand)
  
  # - Get the distance to nearest motif
  message("Calculating distance to nearest motif for each read...")
  distance.to.motif = distanceToNearest(read.positions,motif.sites)
  
  # - Add distance to the GRanges of the read positions
  read.positions$distance.to.motif = NA
  read.positions$distance.to.motif[distance.to.motif@from] = distance.to.motif@elementMetadata$distance
  
  # - Add sequence of nearest motif to the read positions
  message("Adding sequence of nearest motif...")
  read.positions$motif.sequence = NA
  read.positions$motif.sequence[distance.to.motif@from] = motif.sites$motif.sequence[distance.to.motif@to]
  
  # - Add methylation status for each position on the read
  message("Adding methylation status for each CpG in each read...")
  read.positions$methylation = unlist(methylation)
  
  # - Add number of CpGs within each read
  message("Obtaining the number of CpGs for each read...")
  read.positions$number.of.cpgs = unlist(lapply(read.positions$methylation, function(x){
    str_count(x, pattern = "Z|z")
  }))
  
  # - Add fraction of methylated CpGs within each read
  message("Adding methylation percentage for each read...")
  read.positions$m.percent.in.read = unlist(pbmcapply::pbmclapply(read.positions$methylation, function(x){
    unmethylated = str_count(string = x, pattern = "z")
    methylated = str_count(string = x, pattern = "Z")
    proportion = methylated / (unmethylated + methylated)
  }))
  
  # - Add classification on conconrdant/discordant methylation
  message("Classifying concordant / discordant reads...")
  read.positions$concordance = NA
  read.positions$concordance[read.positions$m.percent.in.read != 0 & read.positions$m.percent.in.read != 1] = "Discordant"
  read.positions$concordance[read.positions$m.percent.in.read == 0 | read.positions$m.percent.in.read == 1] = "Concordant"
  
  read.positions$concordance.status = NA
  read.positions$concordance.status[read.positions$m.percent.in.read != 0 & read.positions$m.percent.in.read != 1] = "Discordant"
  read.positions$concordance.status[read.positions$m.percent.in.read == 0] = "Concordant unmethylated"
  read.positions$concordance.status[read.positions$m.percent.in.read == 1] = "Concordant methylated"
  
  read.positions = read.positions[!is.na(read.positions$distance.to.motif),]
  read.positions = read.positions[read.positions$distance.to.motif < extended_region,]
  read.positions$distance.quantile = as.numeric(cut(read.positions$distance.to.motif, breaks = quantile_number))
  
  plot.dat = as.data.frame.matrix(table(read.positions$distance.quantile, read.positions$concordance.status))
  plot.dat$distance.quantile = rownames(plot.dat)
  plot.dat$proportion.of.discordant = plot.dat$`Discordant` / (rowSums(plot.dat[,c("Concordant methylated","Concordant unmethylated","Discordant")]))
  plot.dat$read.support = rowSums(plot.dat[,c("Concordant methylated","Concordant unmethylated","Discordant")])
  plot.dat$m.concordant.percent = plot.dat$`Concordant methylated` / rowSums(plot.dat[,c("Concordant methylated","Concordant unmethylated", "Discordant")]) * 100
  plot.dat$um.concordant.percent = plot.dat$`Concordant unmethylated` / rowSums(plot.dat[,c("Concordant methylated","Concordant unmethylated", "Discordant")]) * 100
  plot.dat$discordant.percent = plot.dat$Discordant / rowSums(plot.dat[,c("Concordant methylated","Concordant unmethylated", "Discordant")]) * 100
  m = melt(plot.dat[,1:6], id.vars = c("distance.quantile","read.support"))
  
  # - Plotting
  message("Plotting...")
  
  s1 = m[m$variable %in% c(c("proportion.of.discordant")),]
  p1 = ggplot(s1, aes(x = as.numeric(distance.quantile), y = value)) +
    geom_bar(stat = "identity", fill = "black", width = 0.75) + 
    #geom_smooth() + 
    theme_classic() + 
    labs(x = "Distance quantile from motif (200 bp)", y = "Proportion of discordant reads", title = motif_symbol) +
    theme(plot.margin = unit(c(0.75, 0.5, 0, 0.5), "cm"), plot.title = element_text(hjust = 0.5))
    
  s2 = m[m$variable %in% c(c("Concordant methylated","Concordant unmethylated","Discordant")),]
  p2  = ggplot(s2, aes(x = as.numeric(distance.quantile), y = value, color = variable)) +
    #geom_point() + 
    geom_smooth() + 
    labs(x = "Distance quantile from motif (200 bp)", y = "Number of reads") + 
    theme_classic() + 
    theme(legend.position = c(0.75, 0.8),
          plot.margin = unit(c(0.5, 2, 0, 2), "cm"),
          legend.key.size = unit(0.25, 'cm'),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank())
  
  
  p3 = ggplot(plot.dat, aes(x = as.numeric(distance.quantile), y = read.support)) + 
    geom_bar(stat = "identity", fill = "black", width = 0.75) + 
    labs(x = "Distance quantile from motif (200 bp)", y = "Total reads in quantile") + 
    theme_classic() +
    theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))
  
  s4 = melt(plot.dat[,4:9], id.vars = c("distance.quantile","read.support"))
  s4 = s4[s4$variable != "proportion.of.discordant",]
  s4$variable = factor(s4$variable, levels = c("discordant.percent","m.concordant.percent","um.concordant.percent"))
  p4 = ggplot(s4, aes(x = as.numeric(distance.quantile), y = value, fill = variable)) +
    geom_bar(stat = "identity", color = "black", width = 0.75) +
    labs(x = "Distance quantile from motif (200 bp)", y = "Percentage of reads (%)") + 
    theme_classic() + 
    scale_fill_manual(values = c("firebrick2","black","white")) +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  pg = plot_grid(p1,p2, p3, ncol = 1, align = 'v', rel_heights = c(3,3,3))  
  pg2 = plot_grid(p3,p4, ncol=1, align = 'v', rel_heights = c(3,3))
  print(pg)
  print(pg2)
  message("Done!")
  
}