# Motif pattern classification
PatternHeatmap <- function(path_to_bam = NULL,
                              path_to_bam_index = NULL,
                              path_to_fimo_file = NULL,
                              path_to_meme_file = NULL,
                              motif_symbol,
                              motif_name,
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
  read.positions$tag = names(unlist(methylation)) # Tags are unique <sum(duplicated(read.positions$tag))>
  
  message(sum(duplicated(read.positions$tag))," duplicated read tags...")
  
  # - Add methylation status for each position on the read
  message("Adding methylation status for each CpG in each read...")
  read.positions$methylation = unlist(methylation)
  
  # - Add number of CpGs within each read
  message("Obtaining the number of CpGs for each read...")
  read.positions$number.of.cpgs = unlist(lapply(read.positions$methylation, function(x){
    str_count(x, pattern = "Z|z")
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
  
  # - Restrict methylation data to those reads overlapping motifs
  idx = findOverlaps(query = motif.sites, subject = read.positions, type = "within")
  positions.start =  motif.sites[idx@from,]@ranges@start - read.positions[idx@to,]@ranges@start
  positions.end =  motif.sites[idx@from,]@ranges@start - read.positions[idx@to,]@ranges@start + motif.sites[idx@from,]@ranges@width
  
  motif.abundance.dat = as.data.frame(sort(table(toupper(read.positions$motif.sequence)), decreasing =  T))
  motif.abundance.dat$rank = 1:nrow(motif.abundance.dat)
  
  # - Restrict the methylation analysis to motif sequence region in each read
  read.motif = read.positions[unique(idx@to),]
  motif.meth <- pbmclapply(1:length(read.motif), function(x){
    str_c(unlist(strsplit(read.motif[x,]$methylation,"")[[1]][positions.start[x]:positions.end[x]]),collapse = "")
  })
  
  read.motif$motif.methylation = unlist(motif.meth)
  
  motif.meth = lapply(motif.meth, function(x) gsub(x, pattern = "h|H|x|X", replacement = "."))
  motif.meth = lapply(motif.meth, function(x) gsub(x, pattern = "z|Z", replacement = "o")) # We only care about the position not the methylation status here
  names(motif.meth) = unlist(motif.meth)
  motif.meth = lapply(motif.meth, function(x) unlist(strsplit(x, "")))
  
  mtx = do.call(rbind,lapply(motif.meth, function(x){
    unlist(lapply(motif.meth, function(y){
      sum(x != y)
    }))
  }))
  
  mtx = mtx[!is.na(mtx[,1]),!is.na(mtx[1,])]
  p = pheatmap::pheatmap(mtx, main = "Distance between CpG patterns")
}
