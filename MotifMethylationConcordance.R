# Motif methylation discordance

# Pseudocode:

# - Extract motif positions
# - Assign motif sequence in read
# - Subset bases that are part of the motif
# - Get methylation status of those bases
# - Calculate concordance or all potential combinations of methylation observed
# - Define metric for methylation preference (i.e. frequency of each methylation pattern pattern)

MotifMethylationConcordance <- function(path_to_bam = NULL,
                                        path_to_bam_index = NULL,
                                        path_to_fimo_file = NULL,
                                        path_to_meme_file = NULL,
                                        motif_symbol,
                                        motif_name,
                                        sample_reads = 10,
                                        ncores = 1){
  
  if(is.null(path_to_bam) | is.null(path_to_fimo_file) | is.null(path_to_bam_index) | is.null(path_to_meme_file)){
    stop("Please provide all the paths to the requested files")
  }
  
  set.seed(10044)
  
  message("Reading in .bam file...")
  bamFile <- BamFile(file = path_to_bam,
                     index = path_to_bam_index)
  bam <- scanBam(bamFile) # Load the bam file into R
  sequence <- bam[[1]]$seq # Extract the read sequence
  methylation <- scanBam(bamFile, param = ScanBamParam(tag = "XM")) # Load the methylation calls per read
  
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
  read.positions$seq = bam[[1]]$seq
  read.positions$tag = names(unlist(methylation))
  
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
  
  read.motif = subsetByOverlaps(read.positions,motif.sites, type = "within")
  idx = findOverlaps(query = motif.sites, subject = read.positions, type = "within")
  positions.start =  motif.sites[idx@from,]@ranges@start - read.positions[idx@to,]@ranges@start
  positions.end =  motif.sites[idx@from,]@ranges@start - read.positions[idx@to,]@ranges@start + motif.sites[idx@from,]@ranges@width-1
  
  # - Restrict the methylation analysis to motif sequence region in each read
  message("")
  read.motif = read.positions[idx@to,]
  motif.meth <- pbmclapply(1:length(read.motif), function(x){
    str_c(unlist(strsplit(read.motif[x,]$methylation,"")[[1]][positions.start[x]:positions.end[x]]),collapse = "")
  })
  read.motif$motif.methylation = unlist(motif.meth)
  
  # - Extract the positions of methylated and unmethylated CpGs within each read
  message("Extracting the positions for methylated or unmethylated CpGs within the motif for each read...")
  in.read.m = pbmclapply(1:length(read.motif$motif.methylation), function(x){
    info = unlist(strsplit(x = read.motif$motif.methylation[[x]],""))
    unmeth.cpg = as.numeric(which(info ==  "z"))
    meth.cpg = as.numeric(which(info ==  "Z"))
    out = data.frame(methylation = c(rep("z", length(unmeth.cpg)),rep("Z", length(meth.cpg))),
                     position = c(unmeth.cpg,meth.cpg))
    if(nrow(out)>0){
      out$tag = read.motif$tag[x]
      return(out)
    }else{NA}
  }, mc.cores = ncores)
  
  in.read.m = in.read.m[!is.na(in.read.m)]
  in.read.m.out = as.data.frame(do.call(rbind,in.read.m))
  
  cpg.coords = split(in.read.m.out, f = in.read.m.out$tag)
  cpg.coords.values = unlist(lapply(cpg.coords, function(x) str_c(as.character(unique(x$position)),collapse = ",")))
  ggplot(as.data.frame(sort(table(cpg.coords.values), decreasing = T)), aes(x = cpg.coords.values, y = Freq)) +
    geom_bar(stat = "identity", fill = "black") +
    labs(x = "CpG position in motif", y = "Frequency") +
    theme_classic() 
  
  
  
  # - Discordant reads plot
  
  concordance = apply(table(in.read.m.out$methylation, in.read.m.out$tag),2, function(x) x[2]/(x[2]+x[1]))
  read.motif$concordance = concordance[read.motif$tag]
  
  read.motif = read.motif[!is.na(read.motif$motif.methylation),]
  
  read.motif$concordance.call = NA
  read.motif$concordance.call[read.motif$concordance %in% c(0,1)] = "Concordant"
  read.motif$concordance.call[!(read.motif$concordance %in% c(0,1))] = "Disconcordant"
  
  read.motif$concordance.status = NA
  read.motif$concordance.status[read.motif$concordance == 0] = "Concordant unmethylated"
  read.motif$concordance.status[read.motif$concordance == 1] = "Concordant methylated"
  read.motif$concordance.status[!(read.motif$concordance %in% c(0,1))] = "Discordant"
  
  dat.plot = in.read.m.out
  dat.plot$concordance = concordance[dat.plot$tag]
  
  concordance.call = read.motif$concordance.call
  names(concordance.call) = read.motif$tag
  
  concordance.status = read.motif$concordance.status
  names(concordance.status) = read.motif$tag
  
  dat.plot$concordance.call = concordance.call[dat.plot$tag]
  dat.plot$concordance.status = concordance.status[dat.plot$tag]
  dat.plot$cpg.number = table(dat.plot$tag)[dat.plot$tag]
  
  
  dat.plot.c.um = dat.plot[dat.plot$concordance.status == "Concordant unmethylated",]
  dat.plot.c.m = dat.plot[dat.plot$concordance.status == "Concordant methylated",]
  dat.plot.dis = dat.plot[dat.plot$concordance.status == "Discordant",]
  
  
  message("Sampling reads for plotting examples...")
  if(length(unique(dat.plot.c.um$tag)) > sample_reads){
    dat.plot.c.um = dat.plot.c.um[dat.plot.c.um$tag %in% c(sample(unique(dat.plot.c.um$tag), size = sample_reads)),]
    dat.plot.c.um = dat.plot.c.um[order(dat.plot$cpg.number, decreasing = T),]
    dat.plot.c.um = dat.plot.c.um[complete.cases(dat.plot.c.um),]
    dat.plot.c.um$tag.number = 1:nrow(dat.plot.c.um)
    tag.number = 1:length(unique(dat.plot.c.um$tag))
    names(tag.number) = unique(dat.plot.c.um$tag)
    dat.plot.c.um$tag.number = tag.number[dat.plot.c.um$tag]
  }else{
    tag.number = 1:length(unique(dat.plot.c.um$tag))
    names(tag.number) = unique(dat.plot.c.um$tag)
    dat.plot.c.um$tag.number = tag.number[dat.plot.c.um$tag]
  }
  if(length(unique(dat.plot.c.m$tag)) > sample_reads){
    dat.plot.c.m = dat.plot.c.m[dat.plot.c.m$tag %in% c(sample(unique(dat.plot.c.m$tag), size = sample_reads)),]
    dat.plot.c.m = dat.plot.c.m[order(dat.plot$cpg.number, decreasing = T),]
    dat.plot.c.m = dat.plot.c.m[complete.cases(dat.plot.c.um),]
    dat.plot.c.m$tag.number = 1:nrow(dat.plot.c.m)
    tag.number = 1:length(unique(dat.plot.c.m$tag))
    names(tag.number) = unique(dat.plot.c.m$tag)
    dat.plot.c.m$tag.number = tag.number[dat.plot.c.m$tag]
  }else{
    tag.number = 1:length(unique(dat.plot.c.m$tag))
    names(tag.number) = unique(dat.plot.c.m$tag)
    dat.plot.c.m$tag.number = tag.number[dat.plot.c.m$tag]
  }
  if(length(unique(dat.plot.dis$tag)) > sample_reads){
    dat.plot.dis = dat.plot.dis[dat.plot.dis$tag %in% c(sample(unique(dat.plot.dis$tag), size = sample_reads)),]
    dat.plot.dis = dat.plot.dis[order(dat.plot$cpg.number, decreasing = T),]
    dat.plot.dis = dat.plot.dis[complete.cases(dat.plot.dis),]
    tag.number = 1:length(unique(dat.plot.dis$tag))
    names(tag.number) = unique(dat.plot.dis$tag)
    dat.plot.dis$tag.number = tag.number[dat.plot.dis$tag]
  }else{
    tag.number = 1:length(unique(dat.plot.dis$tag))
    names(tag.number) = unique(dat.plot.dis$tag)
    dat.plot.dis$tag.number = tag.number[dat.plot.dis$tag]
  }
  
  dat.plot.sub = rbind(dat.plot.c.um, dat.plot.c.m, dat.plot.dis)
  
  # - Plotting
  message("Plotting...")
  p0 = ggseqlogo(pwm) + 
    labs(title = motif_symbol) + 
    scale_x_continuous(expand = c(0,0), labels = 1:motif.length, breaks = 1:motif.length) + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm")) 
  
  
  p1u = ggplot(dat.plot.c.um, aes(x = 1:motif.sites[1,]@ranges@width, y = tag.number)) + 
    geom_hline(yintercept = dat.plot.c.um$tag.number, lwd = 0.25) +
    geom_hline(yintercept = 0, lwd = 0.25, lty = 2) +
    geom_point(aes(x = position, fill = factor(methylation, levels = c("z","Z"))), pch=21) + 
    scale_fill_manual(values = c("white")) +
    labs(x = "Position in read", y = "Concordant unmethylated") +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")) 
  
  p1m = ggplot(dat.plot.c.m, aes(x = 1:motif.sites[1,]@ranges@width, y = dat.plot.c.m$tag.number)) + 
    geom_hline(yintercept = dat.plot.c.m$tag.number, lwd = 0.25) +
    geom_hline(yintercept = 0, lwd = 0.25, lty = 2) +
    geom_point(aes(x = position, fill = factor(dat.plot.c.m$methylation, levels = c("z","Z"))), pch=21) + 
    scale_fill_manual(values = c("black")) +
    labs(x = "Position in read", y = "Concordant methylated") +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")) 
  
  p1d = ggplot(dat.plot.dis, aes(x = 1:motif.sites[1,]@ranges@width, y = dat.plot.dis$tag.number)) + 
    geom_hline(yintercept = dat.plot.dis$tag.number, lwd = 0.25) +
    geom_hline(yintercept = 0, lwd = 0.25, lty = 2) +
    geom_point(aes(x = position, fill = factor(methylation, levels = c("z","Z"))), pch=21) + 
    scale_fill_manual(values = c("white", "black")) +
    labs(x = "Position in read", y = "Discordant") +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")) 
  
  max.read.counts = max(c(max(unique(dat.plot.c.um$tag.number)), max(unique(dat.plot.c.m$tag.number)), max(unique(dat.plot.dis$tag.number))))
  
  rel_values = rescale_max(x = c(max.read.counts/2,max(unique(dat.plot.c.um$tag.number)), max(unique(dat.plot.c.m$tag.number)), max(unique(dat.plot.dis$tag.number))), to = c(0,90))
  rel_values[is.infinite(rel_values)] = 0
  pg = plot_grid(p0, p1u, p1m, p1d, ncol = 1, align = 'v', rel_heights = rel_values)
  message("Done!")
  pg
}
