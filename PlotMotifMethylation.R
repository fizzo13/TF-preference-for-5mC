PlotMotifMethylation <- function(path_to_bismark_cov = NULL,
                                 path_to_fimo_file = NULL,
                                 path_to_meme_file = NULL,
                                 motif_name = NULL,
                                 motif_symbol = NULL,
                                 sample_sites = NULL){
  
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
  motif.sites = motif.sites[grep(motif.sites$sequence_name, pattern = "random|alt|Unk",invert = T),] # Remove random or alt chromosomes
  
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
  
  cpgs.in.motifs = round(length(subsetByOverlaps(bismark.cov,motif.sites)) / length(bismark.cov) * 100, digits = 2)
  message("Percentage of CpGs within motifs: ", cpgs.in.motifs, "%")
  message("Percentage of CpGs outside motifs: ", 100-cpgs.in.motifs,"%") 
  
  # - Check methylation status on CpG sites
  motif.sites$cg.position = lapply(motif.sites$motif.sequence, function(x){
    temp = (str_locate_all(pattern = "CG|cg|cG|Cg", x)) # Include all potential CpG sites
    if(length(temp[[1]]) == 0){
      temp = NA # Return NA if no CpG site is found
    }else{
      temp = as.data.frame(temp)$start # Look specifically into the C position of the CpG site
    }
  })
  
  # Fraction of CpG sites per motif position
  fraction.cg = data.frame(table(factor(unlist(motif.sites$cg.position), levels = 1:motif.length)) / sum(table(unlist(motif.sites$cg.position))))
  
  # Methylation status of each CpG across sites per position in the motif
  cg.coords = lapply(1:length(motif.sites@seqnames), function(x){
    coords = list()
    if(!is.na(unlist(motif.sites[x,]$cg.position))){
      
      for(i in unlist(motif.sites[x,]$cg.position)){
        coords[[paste0(x,"_",i)]] = paste0(as.character(motif.sites[x,]@seqnames),"_",motif.sites[x,]@ranges@start+i,"_pos",i) 
      }
    }
    return(coords)
  })
  
  cg.coords = GRanges(seqnames = unlist(lapply(strsplit(unlist(cg.coords),"_"), function(x) x[1])),
                      IRanges(start = as.numeric(unlist(lapply(strsplit(unlist(cg.coords),"_"), function(x) x[2])))),
                      position = unlist(lapply(strsplit(unlist(cg.coords),"_"), function(x) x[3])))
  bismark.cov$motif.position = NA
  ind = findOverlaps(query = bismark.cov, subject = cg.coords)
  bismark.cov$cpg.position = NA
  bismark.cov$cpg.position[ind@from] = cg.coords$position[ind@to]
  
  meth = list()
  average.meth = list()
  for(i in paste0("pos",1:motif.length)){
    meth[[i]] = sort(bismark.cov[bismark.cov$cpg.position %in% i,]$m.percent, decreasing = T)
    average.meth[[i]] = mean(meth[[i]])
  }
  m = melt(meth)
  m$position.numeric = as.numeric(gsub(m$L1,pattern = "pos",replacement = ""))
  m$position.numeric = factor(m$position.numeric, levels = 1:motif.length)
  m$height = NA
  m$height = unlist(lapply(unique(m$position.numeric),function(x){
    values = 1:sum(m$position.numeric == x)
    m$height[m$position.numeric == x] = values
  }))
  
  # For sampling a fixed number of CpG sites per motif position
  if(!is.null(sample_sites)){
    sampled = list()
    for(i in 1:motif.length){
      if(nrow(m[m$position.numeric == i,]) > sample_sites){
        sampled[[i]] = m[sample(rownames(m[m$position.numeric == i,]), size =  sample_sites),]
        sampled[[i]] = sampled[[i]][order(sampled[[i]]$value, decreasing = T),]
        sampled[[i]]$height = 1:sample_sites
      }else{
        sampled[[i]] = m[rownames(m)[m$position.numeric == i],]
      }
    }
    sampled = do.call(rbind,sampled)
  }
  
  average.meth = as.data.frame(unlist(average.meth))
  average.meth$position.numeric = factor(as.numeric(gsub(rownames(average.meth),pattern = "pos", replacement = "")),levels = 1:motif.length)
  colnames(average.meth) = c("average.meth","position.numeric")
  
  # Visualization
  p1 = ggseqlogo(pwm) + 
    labs(title = motif_symbol) + 
    scale_x_continuous(expand = c(0,0), labels = 1:motif.length, breaks = 1:motif.length) + 
    theme(plot.title = element_text(hjust = 0.5)) 
  
  p2 = ggplot(average.meth, aes(x = position.numeric, y = average.meth, fill = average.meth)) +
    geom_bar(stat = "identity", color = "black") + 
    #scale_fill_gradient2(low = "white", high = "black", na.value = "white")+
    scale_fill_gradientn(limits = c(0,100),
                         colours=c("white", "black"), na.value = "white") +
    theme_classic() +
    labs(x = "", y = "Mean methylation (%)")+
    theme(legend.position = "none")  + 
    ylim(0,100)
  
  if(is.null(sample_sites)){
    p3 = ggplot(m, aes(x = position.numeric, y = -height, fill = value)) +
      geom_point(pch=21,size = 3 , color = "black", lwd =0.001) +
      #scale_fill_gradient(low = "white", high = "black", na.value = "grey") +
      labs(x = "", y = "Per-site methylation") +
      theme_classic() +
      scale_x_discrete(drop = FALSE) +
      scale_fill_gradientn(limits = c(0,100),
                           colours=c("white", "black"),
                           breaks=seq(1,100,1)) +
      theme(legend.position = "none",axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 
  }else{
    p3 = ggplot(sampled, aes(x = position.numeric, y = -height, fill = value)) +
      geom_point(pch=21,size = 3 , color = "black", lwd =0.001) +
      #scale_fill_gradient(low = "white", high = "black", na.value = "grey") +
      labs(x = "", y = "Per-site methylation") +
      theme_classic() +
      scale_x_discrete(drop = FALSE) +
      scale_fill_gradientn(limits = c(0,100),
                           colours=c("white", "black"),
                           breaks=seq(1,100,1)) +
      theme(legend.position = "none",axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 
            
  }
  
  p4 = ggplot(fraction.cg, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = "black") +
    labs(x = "Motif position (bp)", y = "Fraction of CpG sites") +
    theme_classic()
  
  if(is.null(sample_sites)){
    pg = plot_grid(p1, p2, p3, p4, ncol = 1, align = 'v', rel_heights = c(5,5,15,5))  
  }else{
    pg = plot_grid(p1, p2, p3, p4, ncol = 1, align = 'v', rel_heights = c(5,5,sample_sites/2.5,5))  
  }
  pg
}

#### Examples ####
# library(readr)
# library(GenomicRanges)
# library(stringr)
# library(ggplot2)
# library(ggseqlogo)
# library(cowplot)
# library(RColorBrewer)
# 
# PlotMotifMethylation(path_to_bismark_cov = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup.bismark.cov.gz",
#                      path_to_fimo_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
#                      path_to_meme_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
#                      motif_name = "MEME-2",
#                      motif_symbol = "GATA1")
# 
# PlotMotifMethylation(path_to_bismark_cov = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup.bismark.cov.gz",
#                      path_to_fimo_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_1/fimo.tsv",
#                      path_to_meme_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
#                      motif_name = "MEME-1",
#                      motif_symbol = "CTCF",
#                      sample_sites = 20)
# 
# PlotMotifMethylation(path_to_bismark_cov = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup.bismark.cov.gz",
#                      path_to_fimo_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_1/fimo.tsv",
#                      path_to_meme_file = "meNTT/20221005_7_CTCF/K562_CTCF_converted_S71_L002_R1_001.trimmed_bismark_bt2_pe_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
#                      motif_name = "MEME-1",
#                      motif_symbol = "CTCF",
#                      sample_sites = NULL)
# 
# PlotMotifMethylation(path_to_bismark_cov = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup.bismark.cov.gz",
#                      path_to_fimo_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/fimo_out_2/fimo.tsv",
#                      path_to_meme_file = "meNTT/20221005_8_GATA1/K562_GATA1_converted_S72_L002_R1_001.trimmed_bismark_bt2_pe_fixmate_sorted_dedup_macs_summits.temp.exp.exp.tab.fa.out_meme/meme_out/meme.txt",
#                      motif_name = "MEME-2",
#                      motif_symbol = "GATA1", 
#                      sample_sites = 5)