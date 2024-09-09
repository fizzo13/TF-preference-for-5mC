PlotCpGDistribution = function(path_to_fimo_file = NULL,
                               path_to_meme_file = NULL,
                               motif_name = NULL,
                               motif_symbol = NULL,
                               output_file = F){
  
  
  if(is.null(path_to_fimo_file) | is.null(path_to_meme_file) | is.null(motif_name)){
    stop("Please provide all the paths to the files and the MEME motif name (i.e. MEME-1, MEME-2)")
  }
  
  # - Load FIMO motif positions
  motif.sites <-suppressMessages(suppressWarnings(read_delim(path_to_fimo_file, # fimo.tsv file from MEME
                                                             delim = "\t", escape_double = FALSE, 
                                                             trim_ws = TRUE)))
  
  motif.sites = motif.sites[!is.na(motif.sites$sequence_name),]
  
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
  
 counts = as.data.frame(unlist(lapply(motif.sites$motif.sequence, function(x) str_count(string = x, pattern = "CG"))))
 colnames(counts) = "CpG_sites"
 counts$rank = NA
 counts = counts[order(counts$CpG_sites, decreasing = T),]
 counts$rank = 1:nrow(counts)

 p1 = ggseqlogo(pwm) + 
   labs(title = motif_symbol) + 
   scale_x_continuous(expand = c(0,0), labels = 1:motif.length, breaks = 1:motif.length) + 
   theme(plot.title = element_text(hjust = 0.5),
         plot.margin = unit(c(0.25, 0, 0, 0), "cm")) 
 
 p2 = ggplot(counts, aes(x = CpG_sites))  +
   geom_histogram(bins = 20, color = "black", fill = "lightgrey") + 
   labs(x = "Number of CpG sites per motif", y = "Frequency") +
   theme_classic() +
   scale_x_continuous(breaks= seq(0, max(counts$CpG_sites),1))
 
 pg = plot_grid(p1, p2, ncol = 1, align = 'v', rel_heights = c(4,7))  
 message("Done!")
 if(output_file){
   return(counts)
 }
 pg
} 
 