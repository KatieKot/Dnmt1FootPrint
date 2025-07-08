library(ggthemes, viridis)
get_data <- function(samplas, wdir, outdatadir, dCG, ats=3) {
  if(!file.exists(paste0(outdatadir, samplas, ".RDS"))) {
    rezu <- list()
    reads <- readRDS(paste0(wdir, "output/", samplas, "/reads_CG.RDS"))
    rezu[["dist2best"]] <- table(mcols(reads)$dist2bestCG)

    j <- abs(mcols(reads)$dist2bestCG) < ats
    j[is.na(j)] <- FALSE
    tmp <- mcols(reads[j])$bestCG %>% 
      as.data.table %>% 
      setnames(., c("ID")) %>% 
      .[, ID := gsub("_[+-]_[CG]+", "", ID)] %>% 
      .[, .N, by=ID] %>% 
      .[, chr := as.character(ID) %>% strsplit(., "_") %>% sapply(., `[`, 1)] %>% 
      merge(., (as.data.table(dCG) %>% .[, .(ID)]), by="ID", all=TRUE) %>% 
      .[, chr := as.character(ID) %>% strsplit(., "_") %>% sapply(., `[`, 1)] %>% 
      .[, start := as.character(ID) %>% strsplit(., "_") %>% sapply(., `[`, 2)] %>% 
      .[, end := start] %>%
      .[, chr := paste0("chr", chr)] %>% 
      .[, ID := paste0("chr", ID)] %>%
      setkey(., chr, start) 
    tmp[is.na(tmp)] <- 0 
    rezu[["covTable"]] <- tmp
    rezu[["strandTable"]] <- mcols(reads[j])$bestCG %>% 
      as.data.table %>% 
      setnames(., c("ID")) %>% 
      #.[, ID := gsub("_[+-]_[CG]+", "", ID)] %>% 
      .[, .N, by=ID] %>% 
      .[, chr := strsplit(ID, "_") %>% sapply(., `[`, 1)] %>% 
      .[, start := strsplit(ID, "_") %>% sapply(., `[`, 2)] %>% 
      .[, strand := strsplit(ID, "_") %>% sapply(., `[`, 3)] %>% 
      .[, end := start] %>%
      makeGRangesFromDataFrame(., keep.extra.columns=TRUE) 
    seqlevelsStyle(rezu[["strandTable"]]) <- "UCSC"
      saveRDS(rezu, paste0(outdatadir, samplas, ".RDS"))} else {
        rezu <- readRDS(paste0(outdatadir, samplas, ".RDS"))}   
  return(rezu)
  } 


get_id <- function(x) {return(paste0(seqnames(x), "_", start(x)))}


cols_samples  <- c("WT"="#070D0D", 
                        "Het"="#FF0000", 
                        "Hom"="#008000",
                        "KO"="#606060"
                        )

cols_samplesLong  <- c("WT"="#070D0D", 
                        "hetero"="#FF0000", 
                        "homo"="#008000",
                        "KO"="#606060"
                        )
                        
cols_order  <- c("WT", "hetero", "homo", "KO")


cols_samples_RRBS <- c(
                       "A0D"="#ffc800",
                       "A2D"="#ff9600",
                       "A4D"="#ff6400",
                       "A8D"="#ff0000"
)

cols_samples_4taskai_TOP <- c(
                       "D0azide"="#ffc800",
                       "D2azide"="#ff9600",
                       "D4azide"="#ff6400",
                       "D8azide"="#ff0000",
                       "Ctr2i"="grey"
)

theme_Publication <- function(base_size=10, base_family="Arial") {

      (theme_foundation(base_size=base_size, base_family=base_family) + 
        theme(
              #plot.title = element_text(face = "bold",
              #                           size = rel(1.5), 
              #                           hjust = 0.5),
               text = element_text(family = base_family, size = base_size, color = "black"),

               panel.background = element_rect(fill = "transparent", colour = NA),
               #panel.border = element_rect(colour = "black"),
               #panel.border = element_rect(colour = "black", fill=NA, size=1),
               #panel.border = element_blank(),

               plot.background = element_rect(fill = "transparent", colour = NA),
               
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               
               #axis.title.y = element_text(angle=90, vjust = 2),
               #axis.title.x = element_text(vjust = -0.2),
               axis.line = element_line(colour=NA),
               axis.text = element_text(size=base_size, color="black"),
               #axis.text.x =element_text(size = rel(1), angle = 90, vjust = 0.5, hjust=1),
               #axis.text.y =element_text(size = rel(1)),
               axis.ticks = element_line(linewidth=0.3, color="black"),
               axis.title = element_text(size = rel(1.2), color="black"),
               
               strip.background=element_blank(),
               #strip.text = element_text(face="bold"),
               strip.text.x = element_text(size=base_size, color="black"),
               strip.text.y = element_text(size=base_size, color="black", angle = 0),
               legend.position="bottom",  
               legend.background = element_rect(
                    fill = alpha("white", 0),
                    color = alpha("white", 0)
                  ),
              legend.key = element_rect(color = NA, fill = NA),
              plot.margin=unit(c(2, 2, 2, 2),"mm"),
              legend.key.size = unit(1, "line"),
              legend.text = element_text(size = base_size),

              plot.tag=element_text(size=14, face="bold")

#               legend.key = element_rect(colour = NA),
#               legend.position = "bottom",
#               legend.direction = "horizontal",
#               legend.key.size= unit(0.2, "cm"),
#               legend.margin = margin(t = 0, unit = "cm"),
#               legend.title = element_text(face="italic"),
#               plot.margin=unit(c(5, 5, 5, 5),"mm")
               
          ))      
}

plotCor <- function(x, taitlas) {
  d <- foreach(i=colnames(x[, -1]), .combine="rbind") %do% {
    tmp <- foreach(j=colnames(x[, -1]), .combine="rbind") %do% {
      rezu <- cor.test(x[, get(i)], x[, get(j)])
      return(data.frame(pval=rezu$p.value, estimate=rezu$estimate, s1=i, s2=j))} 
    return(tmp)  
  } %>% 
    as.data.table() %>%  
    .[, s1 := gsub("_met", "", s1)] %>% 
    .[, s2 := gsub("_met", "", s2)] %>% 
    .[, s1 := gsub("_cov", "", s1)] %>% 
    .[, s2 := gsub("_cov", "", s2)] %>% 
    .[, significance := "NotSignificant"] %>% 
    .[pval <= 0.05, significance := "Significant"] 

  p2 <- d %>%  
    ggplot(aes(s1, s2, fill=estimate)) +
      geom_tile() +
      theme_bw() +
      geom_text(data=d, aes(s1, s2, colour=significance, label=round(estimate, 2))) +
      scale_fill_viridis(limits = c(0, 1)) +
      theme(axis.title=element_blank(), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position="bottom") +
      scale_colour_manual(values=c("black", "grey80")) +
      labs(fill="pearson", 
          subtitle=taitlas)    
  return(p2)
}
  
plotPCA <- function(x, taitlas, samp2do) {
  tmp <- x[, paste0(samp2do, "_met"), with=FALSE]
  tmp <- tmp[!is.na(rowSums(tmp)), ]

  pca <- tmp %>% t %>% prcomp()
  summ <- summary(pca)
  p <- pca$x %>% 
    as.data.table(., keep.rownames=TRUE) %>% 
    .[, technology := str_extract(rn, "..BS")] %>% 
    .[, genotype := strsplit(rn, "_") %>% sapply(., `[`, 2)] %>% 
    ggplot(aes(PC1, PC2, shape=technology, colour=genotype)) +
      geom_point(size=3) +
      theme_bw() +
      labs(subtitle=paste0("PCA ", taitlas,  ". All bins")) +
      theme(legend.position="bottom") +
      xlab(paste0("PC 1: ", round(summ$importance[2,1]*100), "%")) +
      ylab(paste0("PC 2: ", round(summ$importance[2,2]*100), "%")) 
  return(p)    
}



### Metilimo grupės 
add_gr_1 <- function(x) {
  tmp <- x %>% 
    as.data.table() %>% 
    .[, group := "KITA" ] %>% 
    .[is.na(x), group := "lowCov"] %>% 
    .[x >= 0, group := "0-20"] %>% 
    .[x >= 20, group := "20-80"] %>% 
    .[x >= 80, group := "80-100"] %>%
    .[, group] %>% 
    table
  return(tmp/length(x))
  }

add_gr_1_noNA <- function(x) {
  tmp <- x %>% 
    as.data.table() %>% 
    .[, group := "KITA" ] %>% 
    .[is.na(x), group := "lowCov"] %>% 
    .[x >= 0, group := "0-20"] %>% 
    .[x >= 20, group := "20-80"] %>% 
    .[x >= 80, group := "80-100"] %>%
    .[group != "lowCov", ] %>% 
    .[, group] %>% 
    table
  return(tmp/length(x[!is.na(x)]))
  }

add_gr_2 <- function(x) {
  tmp <- x %>% 
    as.data.table() %>% 
    .[, group := "KITA" ] %>% 
    .[is.na(x), group := "lowCov"] %>% 
    .[x >= 0, group := "0-10"] %>% 
    .[x >= 10, group := "10-50"] %>% 
    .[x >= 50, group := "50-100"] %>%
    .[, group] %>% 
    table
  return(tmp/length(x))  
  }

add_gr_2_noNA <- function(x) {
  tmp <- x %>% 
    as.data.table() %>% 
    .[, group := "KITA" ] %>% 
    .[is.na(x), group := "lowCov"] %>% 
    .[x >= 0, group := "0-10"] %>% 
    .[x >= 10, group := "10-50"] %>% 
    .[x >= 50, group := "50-100"] %>%
    .[group != "lowCov", ] %>% 
    .[, group] %>% 
    table
  return(tmp/length(x[!is.na(x)]))
  }



get_barPlot_bisulfite <- function(pos2do, bisulfite, namas, samples2do, colSamples, samplesOrder) {
  d <- foreach(i=samples2do, .combine="rbind") %do% {
    tmp_cov <- bisulfite[queryHits(findOverlaps(bisulfite, pos2do))]
    ivertis <- tmp_cov %>% 
      as.data.table() %>% 
      .[, i, with=FALSE]  %>% 
      unlist()

    data.table(mod=ivertis, meginys=i, pos=paste0("pos", 1:length(tmp_cov))) %>% 
      .[, pos := factor(pos, levels=paste0("pos", 1:length(tmp_cov)))] %>% 
      .[] 
    }

  fig <- d %>% 
    as.data.table() %>% 
    .[, spalva := gsub("..BS_", "", meginys) %>% gsub("_met", "", .)] %>% 
    .[, meginys := factor(meginys, levels=samplesOrder)] %>% 
    ggplot(aes(pos, mod)) +
      geom_col(aes(fill=spalva)) +
      theme_bw() +
      ylab("Modified CpGs, methylation") +
      scale_fill_manual(values = colSamples, name="Sample") +
      facet_wrap(~meginys, ncol=1) +   
      theme_Publication() +      
      theme(axis.ticks.x = element_blank(), 
            panel.border = element_blank(), 
            #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            strip.text = element_text(colour = 'black'),
            #strip.background =element_rect(fill=NA),
            strip.background = element_blank(),
            legend.position="bottom",
            axis.title.x = element_blank(),
            axis.text.x = element_blank()
                ) +
      geom_smooth(aes(group="1"), 
            se=F,
            method="glm",
            formula=y~ns(x, 8),
            #family=gaussian(link="log"),
            method.args = list(family = "gaussian"),
            show_guide = FALSE, 
            lwd=0.7) + 
      ggtitle(namas) + 
      coord_cartesian(ylim=c(0, 100))                               
  return(fig)         
} 

do_dif <- function(pos2do, bisulfite, namas, samp1, samp2) {
  #samp1 <- "WGBS_WT_met"
  #samp2 <- "WGBS_hetero_met"
  tmp_cov <- bisulfite[queryHits(findOverlaps(bisulfite, pos2do))]
  ivertis_S1 <- tmp_cov %>% 
      as.data.table() %>% 
      .[, samp1, with=FALSE]  %>% 
      unlist()
  ivertis_S2 <- tmp_cov %>% 
      as.data.table() %>% 
      .[, samp2, with=FALSE]  %>% 
      unlist()
  p <- data.table(S1=ivertis_S1, S2=ivertis_S2) %>% 
    .[, dif := S1-S2] %>% 
    .[, pos := paste0("pos", 1:length(ivertis_S2))] %>% 
    ggplot(aes(pos, dif)) +
      geom_col() +
      theme_bw() +
      ylab("Dif. in methylation") +
      theme_Publication() +      
      theme(axis.ticks.x = element_blank(), 
            panel.border = element_blank(), 
            #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            strip.text = element_text(colour = 'black'),
            #strip.background =element_rect(fill=NA),
            strip.background = element_blank(),
            legend.position="bottom",
            axis.title.x = element_blank(),
            axis.text.x = element_blank()
                ) + 
      geom_smooth(aes(group="1"), 
            se=F,
            method="glm",
            formula=y~ns(x, 8),
            #family=gaussian(link="log"),
            method.args = list(family = "gaussian"),
            show_guide = FALSE, 
            lwd=0.7) +           
      ggtitle(namas) 
  return(p)
}

get_barplot_topSeq <- function(pos2do, coverages_gr_rep, namas, samples2do, colSamples, samplesOrder) {
  samplesInfo <- samples2do %>% 
    as.data.table() %>% 
    setnames(., "rep") %>% 
    .[, samplas := gsub("_R.", "", rep)] %>% 
    .[]

  d <- foreach(i=unique(samplesInfo$samplas), .combine="rbind") %do% {
    tmp_cov <- coverages_gr_rep[queryHits(findOverlaps(coverages_gr_rep, pos2do))]
    vidurkiai <- tmp_cov %>% 
      as.data.table() %>% 
      .[, samplesInfo[samplas == i, rep], with=FALSE] %>% 
      rowMeans
    nuokrypiai <- tmp_cov %>% 
      as.data.table() %>% 
      .[, samplesInfo[samplas == i, rep], with=FALSE] %>% 
      apply(., 1, sd)
    data.table(meanas=vidurkiai, nuokrypiai, meginys=i, pos=paste0("pos", 1:length(vidurkiai))) %>% 
      .[, pos := factor(pos, levels=paste0("pos", 1:length(vidurkiai)))] %>% 
      .[] 
  }

  fig <- d %>% as.data.table() %>% 
      .[, apacia := meanas - (nuokrypiai/2)] %>% 
      .[, virsus := meanas + (nuokrypiai/2)] %>%  
      .[, meginys := factor(meginys, levels=samplesOrder)] %>% 
      ggplot(aes(pos, meanas, fill=meginys)) +
        geom_col() +
        geom_errorbar(aes(ymin = apacia, ymax = virsus), width = 0.2) +
        theme_bw() +
        ylab("Modified CpGs, reads") +
        facet_wrap(~meginys, ncol=1, scales="free_y") +
        theme_Publication() +    
        theme(axis.ticks.x = element_blank(), 
                panel.border = element_blank(), 
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"),
                strip.text = element_text(colour = 'black'),
                #strip.background =element_rect(fill=NA),
                strip.background = element_blank(),
                legend.position="bottom",
                axis.title.x = element_blank(),
                axis.text.x = element_blank()
                ) +
        geom_smooth(aes(group="1"), 
                se=F,
                method="glm",
                formula=y~ns(x, 6),
                #family=gaussian(link="log"),
                method.args = list(family = "gaussian"),
                show_guide = FALSE, 
                lwd=0.7) +     
       scale_fill_manual(values = colSamples, name="Sample") +                       
      ggtitle(namas)                                
  return(fig)                        
} 