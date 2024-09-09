# Motif pattern classification

#Pseudocode:

#1 - Obtain most frequent CpG patterns within the motif
#2 - Select the top-represented CpG patterns in motifs and extract motif methylation status
#3 - For each motif, show the % for each methylation pattern found

PatternClassifier <- function(path_to_bam = NULL,
                              path_to_bam_index = NULL,
                              path_to_fimo_file = NULL,
                              path_to_meme_file = NULL,
                              motif_symbol,
                              motif_name,
                              top =5,
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
  
  # Motif classification based on CpG context positions
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
  
  in.read.m = in.read.m[!is.na(in.read.m)] # Clean up those that don't contain CpG sites
  in.read.m.out = as.data.frame(do.call(rbind,in.read.m))
  # - Obtain label for CpG positions
  cpg.position.in.motif = unlist(lapply(in.read.m, function(x){
    str_c(x[,2], collapse = "_")
  }))
  names(cpg.position.in.motif) = unlist(lapply(in.read.m, function(x) unique(x[,3])))
  
  read.motif$cpg.arrangement = NA
  read.motif$cpg.arrangement =  cpg.position.in.motif[read.motif$tag]
  cpg.in.motif.arrangement = sort(table(cpg.position.in.motif), decreasing = T)
  cpg.in.motif.arrangement.count = as.data.frame(cpg.in.motif.arrangement)
  cpg.in.motif.arrangement.count$rank = 1:nrow(cpg.in.motif.arrangement.count)
  p0 = ggplot(cpg.in.motif.arrangement.count, aes(x = rank, y = Freq)) +
    geom_point() +
    theme_classic() +
    labs(x = "CpG arrangement in motif rank", y = "Frequency")
  
  in.read.m.out$cpg.arrangement = cpg.position.in.motif[in.read.m.out$tag]
  read.motif = read.motif[!is.na(read.motif$cpg.arrangement),]
  
  # - Get proportion of methylation pattterns per motif
  potential.patterns = unique(read.motif$cpg.arrangement)[!is.na(unique(read.motif$cpg.arrangement))]
  patterns = lapply(potential.patterns, function(x){
    temp = in.read.m.out[in.read.m.out$cpg.arrangement == x,]
    bases = c(as.numeric(unlist(strsplit(x,"_"))))
    # print(x)
    if(length(unique(temp$tag))>1){
      out = as.data.frame(do.call(rbind,lapply(unique(temp$tag), function(y){
        strsplit(read.motif[read.motif$tag == y,]$motif.methylation,"")[[1]][bases]
      })))
      colnames(out) = bases
      
      pats = apply(out,1,function(x) str_c(x,collapse = "."))
      pats.freq = table(apply(out,1,function(x) str_c(x,collapse = ".")))/ sum(table(apply(out,1,function(x) str_c(x,collapse = "."))))
      
      out = apply(out, 2, function(z){
        z = gsub(z, pattern = "z", replacement = 0)
        z = gsub(z, pattern = "Z", replacement = 1)
      })
      out = list(m = apply(out,2,as.numeric), pat = pats, freq = pats.freq)
      return(out)
    }else{
      if(length(unique(temp$tag))==1){
        out = unlist(lapply(unique(temp$tag), function(y){
          strsplit(read.motif[read.motif$tag == y,]$motif.methylation,"")[[1]][bases]
        }))
        names(out) = bases
        
        pats = str_c(out,collapse = ".")
        pats.freq = 1
        
        out = gsub(out, pattern = "z", replacement = 0)
        out = gsub(out, pattern = "Z", replacement = 1)
        out = as.data.frame(t(out))
        out = apply(out,2,as.numeric)
        
        out = list(m = out, pat = pats, freq = pats.freq)
        return(out)
      }}
  })
  names(patterns) = potential.patterns
  
  pf = as.data.frame(sort(table(read.motif$cpg.arrangement[!is.na(read.motif$cpg.arrangement)]), decreasing = F))
  patterns = patterns[as.character(pf$Var1)]
  
  top.patterns = patterns[as.character(cpg.in.motif.arrangement.count$cpg.position.in.motif)[1:top]]
  plot.dat = lapply(top.patterns, function(i){
    #print(i)
    if(ncol(as.data.frame(i$m))>1){
      samp = as.data.frame(i$m)
      samp$pat = i$pat
      if(length(unique(samp$pat))>1){
        ind = NULL
        for(j in unique(samp$pat)){
          print(j)
          ind[j] = min(which(samp$pat == j))
        }
        samp = samp[ind,]
        rownames(samp) = samp$pat
        samp = samp[names(sort(i$freq, decreasing = F)),]
        samp$position = 1:nrow(samp)
      }else{
        samp = samp[1,]
        samp$position = 1
      }
      
      counts = cpg.in.motif.arrangement.count$Freq[as.character(cpg.in.motif.arrangement.count$cpg.position.in.motif) == str_c(colnames(i$m), collapse = "_")]
      
      p1 = ggseqlogo(data = pwm)+
        annotate('rect', xmin = as.numeric(colnames(i$m))-0.5,
                 xmax = as.numeric(colnames(i$m))+0.5, ymin = -0.05, ymax = 2, alpha = .1, col='black', fill='yellow')+
        theme_logo() +
        theme(axis.line.y = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(), axis.title.y = element_blank()
        )
      
      m = melt(samp, id.vars = "position")
    }else{
      samp = as.data.frame(i$m)
      samp$pat = i$pat
      
      if(length(unique(samp$pat))>1){
        ind = NULL
        for(j in unique(samp$pat)){
          print(j)
          ind[j] = min(which(samp$pat == j))
        }
        samp = samp[ind,]
        rownames(samp) = samp$pat
        samp = samp[names(sort(i$freq, decreasing = F)),]
        samp$position = 1:nrow(samp)
      }else{
        samp = samp[1,]
        samp$position = 1
      }

      m = melt(samp, id.vars = c("position","pat"))
      
      counts = cpg.in.motif.arrangement.count$Freq[as.character(cpg.in.motif.arrangement.count$cpg.position.in.motif) == str_c(colnames(i$m), collapse = "_")]
      
      p1 = ggseqlogo(data = pwm)+
        annotate('rect', xmin = as.numeric(colnames(i$m))-0.5,
                 xmax = as.numeric(colnames(i$m))+0.5, ymin = -0.05, ymax = 2, alpha = .1, col='black', fill='yellow')+
        theme_logo() +
        theme(axis.line.y = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(), axis.title.y = element_blank()
        )
    }
    
    p2 = ggplot(m, aes(x = 1:motif.sites[1,]@ranges@width, y = m$position)) + 
      geom_hline(yintercept = m$position, lwd = 0.25) +
      geom_point(aes( y = m$position, x = as.numeric(as.character(variable)), fill = factor(value, levels = c(0,1))), pch=21, size = 2) + 
      scale_fill_manual(values = c("white", "black")) +
      labs(x = "Position in motif", y = "Methylation pattern", title = paste0("N = ",counts," reads with CpG pattern")) +
      theme_classic()  + 
      theme(legend.position = "none",
            axis.line.y = element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      scale_x_continuous("",labels =  1:motif.sites[1,]@ranges@width, breaks = 1:motif.sites[1,]@ranges@width,limits = c(1,motif.sites[1,]@ranges@width)) +
      ylim(0, max(m$position)+1)
    
    d = as.data.frame(i$freq)
    if(nrow(d)>1){
      d = d[order(d$Freq, decreasing = F),]
    }
    d$Var1 = factor(d$Var1, levels = d$Var1)  
    p3 = ggplot(d, aes(x = Var1, y = Freq)) +
      geom_bar(stat = "identity", fill = "black") + theme_classic() +
      labs(x = "Methylation pattern", y = "Frequency") +
      ylim(0,1) +
      coord_flip() +
      geom_text(aes(y = Freq-0.025, label=round(Freq,2)),color = "white", position=position_dodge(width=0.9), size =3)
    
    pg = plot_grid(p1, p2,p3, nrow = 1, align = 'h', rel_widths = c(10,10,10))
    message("Done!")
    # print(pg)
    return(as_grob(pg))
  })
  plot.dat[["p0"]] = as_grob(p0)
  return(plot.dat)
}