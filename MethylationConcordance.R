# - Load read data
# - Define reads that contain the motif site
# - Define the methylation status for each read
# - Calculate per-read metrics:
## - In motif, methylation status
## - In motif, methylation patterns present (possible combinations)
## - Outside motif, proportion of discordant reads

# Nomenclature:
# Z = methylated C in CpG context
# z = unmethylated C in CpG context
# H = methylated C in CHH context
# h = unmethylated C in CHH context
# X = methylated C in CHG context
# x = unmethylated C in CHG context

MethylationConcordance <- function(path_to_bam = NULL,
                                   path_to_bam_index = NULL,
                                   path_to_fimo_file = NULL,
                                   motif_symbol,
                                   sample_reads = 100,
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
  # hist(read.positions$distance.to.motif[read.positions$distance.to.motif <= 2000], breaks = 100)
  
  # - Add sequence of nearest motif to the read positions
  message("Adding sequence of nearest motif...")
  read.positions$motif.sequence = NA
  read.positions$motif.sequence[distance.to.motif@from] = motif.sites$motif.sequence[distance.to.motif@to]
  # hist(table(read.positions$motif.sequence), main = "Number of distinct motif sequences", xlab = "Motif sequence")
  
  # - Add methylation status for each position on the read
  message("Adding methylation status for each CpG in each read...")
  read.positions$methylation = unlist(methylation)
  
  # - Add number of CpGs within each read
  message("Obtaining the number of CpGs for each read...")
  read.positions$number.of.cpgs = unlist(lapply(read.positions$methylation, function(x){
    str_count(x, pattern = "Z|z")
  }))
  # hist(read.positions$number.of.cpgs, main = "CpG sites per read", breaks = 20)
  
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
  
  message("Performing Fisher Exact test for concordant / discordant reads containing or not containing the motif...")
  tab = as.data.frame.matrix(table(read.positions$distance.to.motif, read.positions$concordance))
  tab$Proportion = tab$Discordant / (tab$Discordant+tab$Concordant)
  tab$Read.support = table(read.positions$distance.to.motif)
  tab$Distance = as.numeric(rownames(tab))
  tab$Distance.quantile = as.numeric(cut(tab$Distance, breaks = 100))
  
  motif.tab = read.positions[!is.na(read.positions$distance.to.motif),]
  motif.tab = motif.tab[motif.tab$distance.to.motif == 0,]
  
  fisher.tab = rbind(table(motif.tab$concordance),table(read.positions$concordance[read.positions$distance.to.motif > 0]))
  rownames(fisher.tab) = c("With_motif","Without_motif")
  ft = fisher.test(fisher.tab)
  
  pval = format(ft$p.value, digits = 4)
  or = round(ft$estimate, 2)
  
  m = as.data.frame(apply(fisher.tab,1,function(x) x[2] / (x[2]+x[1]) *100))
  colnames(m) = "Percent.of.discordant.reads"
  m$motif.status = c("Reads with motif", "Reads without motif")
  
  stats = data.frame(y.position = max(m$Percent.of.discordant.reads),
                     group1 = "Reads with motif",
                     group2 = "Reads without motif",
                     n1 = rowSums(fisher.tab)[1],
                     n2 = rowSums(fisher.tab)[2],
                     statistic = NA,
                     df = or,
                     p = pval,
                     p.signif = NA)
  
  # - Extract the positions of methylated and unmethylated CpGs within each read
  message("Extracting the positions for methylated or unmethylated CpGs within each read...")
  in.read.m = pbmclapply(1:length(read.positions$methylation), function(x){
    info = unlist(strsplit(x = read.positions$methylation[[x]],""))
    unmeth.cpg = as.numeric(which(info ==  "z"))
    meth.cpg = as.numeric(which(info ==  "Z"))
    out = data.frame(methylation = c(rep("z", length(unmeth.cpg)),rep("Z", length(meth.cpg))),
                     position = c(unmeth.cpg,meth.cpg))
    if(nrow(out)>0){
      out$tag = names(read.positions$methylation)[x]
      return(out)
    }else{NA}
  }, mc.cores = ncores)
  
  in.read.m = in.read.m[!is.na(in.read.m)]
  in.read.m.out = as.data.frame(do.call(rbind,in.read.m))
  
  # - Sample a number of reads to show as example under the discordant reads plot
  
  ind.1 = sample(names(which(read.positions$distance.to.motif == 0 & read.positions$number.of.cpgs > 1)),sample_reads)
  ind.2 = sample(names(which(read.positions$distance.to.motif %in% 1:2000 & read.positions$number.of.cpgs > 1)),sample_reads)
  
  group1 = in.read.m.out[in.read.m.out$tag %in% ind.1,]
  group2 = in.read.m.out[in.read.m.out$tag %in% ind.2,]
  
  group1$in.motif = "Reads with motif"
  group2$in.motif = "Reads without motif"
  
  concordance1 = apply(table(group1$methylation, group1$tag),2, function(x) x[2]/(x[2]+x[1]))
  concordance2 = apply(table(group2$methylation, group2$tag),2, function(x) x[2]/(x[2]+x[1]))
  
  group1$concordance = concordance1[group1$tag]
  group2$concordance = concordance2[group2$tag]
  
  group1 = group1[order(group1$concordance, decreasing = F),]
  group2 = group2[order(group2$concordance, decreasing = F),]
  
  tag.index1 = 1:length(unique(group1$tag)); names(tag.index1) = unique(group1$tag)
  tag.index2 = 1:length(unique(group2$tag)); names(tag.index2) = unique(group2$tag)
  
  group1$tag.number = tag.index1[group1$tag]
  group2$tag.number = tag.index2[group2$tag]
  
  dat.plot = rbind(group1,group2)
  
  dat.plot$concordance.call = NA
  dat.plot$concordance.call[dat.plot$concordance %in% c(0,1)] = "Concordant"
  dat.plot$concordance.call[!(dat.plot$concordance %in% c(0,1))] = "Disconcordant"
  
  dat.plot$concordance.status = NA
  dat.plot$concordance.status[dat.plot$concordance == 0] = "Concordant unmethylated"
  dat.plot$concordance.status[dat.plot$concordance == 1] = "Concordant methylated"
  dat.plot$concordance.status[!(dat.plot$concordance %in% c(0,1))] = "Disconcordant"
  
  # - Plotting
  message("Plotting...")
  p1 = ggplot(m, aes(x = motif.status, y = Percent.of.discordant.reads)) + 
    geom_bar(stat = "identity", fill = "lightgrey", color ="black") +
    labs(x = "", y = "Discordant reads (%)", title = motif_symbol) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0.5, 0.5, 0, 0.25), "cm")) +
    stat_pvalue_manual(data = stats,label = paste0("P = {p}"),
                       vjust = -1, bracket.nudge.y = 0.5
    ) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  
  p2 = ggplot(dat.plot, aes(x = 1:read.positions[1,]@ranges@width, y = dat.plot$tag.number)) + 
    geom_hline(yintercept = dat.plot$tag.number, lwd = 0.25) +
    geom_point(aes(x = position, fill = methylation),pch=21) + 
    scale_fill_manual(values = c("white", "black")) +
    facet_wrap(~in.motif, scales = "free_y") + 
    labs(x = "Position in read", y = "Read") +
    theme_classic() +
    theme(legend.position = "none",axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), axis.line.y = element_blank(), panel.border = 
            element_blank(),strip.background = element_blank(), strip.text = element_blank(),
          plot.margin = unit(c(0, 1, 0, 1), "cm"))
  
  pg = plot_grid(p1, p2, ncol = 1, align = 'v', rel_heights = c(6,7))  
  message("Done!")
  pg
}






