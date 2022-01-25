#install.packages("RCircos")
library(RCircos)
# https://cran.r-project.org/web/packages/RCircos/vignettes/Using_RCircos.pdf
# https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html


#BiocManager::install("liftOver")
library(liftOver)
# https://bioconductor.org/packages/devel/workflows/vignettes/liftOver/inst/doc/liftov.html
library(rtracklayer)
library(GenomicRanges)


create_acc_mock_data <- function() {
  
  chromosomes <- c(1:22)
  chromosomes_names <- paste0('chr', chromosomes)

  data(UCSC.HG38.Human.CytoBandIdeogram)
  data(RCircos.Gene.Label.Data)
  data(RCircos.Heatmap.Data)
  
  chr.exclude <- NULL;
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram;
  
  region_len <- 50
  acc_count_by_chr <- 30
  
  result <- data.frame("Chromosome"= character(),
                       "chromStart" = integer(),
                       "chromEnd" = integer())
  
  for (chr in chromosomes_names) {
   
    min_pos <- 0
    mask <- UCSC.HG38.Human.CytoBandIdeogram$Chromosome == chr
    max_pos <- max(UCSC.HG38.Human.CytoBandIdeogram$chromEnd[mask])
    
    acc_starts <- sample(c(min_pos:max_pos), size = acc_count_by_chr, replace = FALSE)
    acc_ends <- acc_starts + region_len - 1
    
    tmp_df <-  data.frame("Chromosome"= rep(chr, acc_count_by_chr),
                          "chromStart" = acc_starts,
                          "chromEnd" = acc_ends)
    
    result <- rbind(result, tmp_df)
    
  }
    
  return(result)

    
}

liftover_acc_data <- function(data, source_specie, target_specie) {
  
  
  # data <- data_3
  # source_specie <- 'hg38' # tolower
  # target_specie <- 'GalGal4' # to upper first
  
  #chain <- import.chain("hg19ToHg18.over.chain")
  #library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  #tx_hg19 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
  #tx_hg18 <- liftOver(tx_hg19, chain)

  chain_file_name <- paste0(source_specie, 'To', target_specie)
    
  chain_file_path <- paste0('/home/paulati/Documents/paula/ingebi/2022/paula/precision/synteny/data/', chain_file_name,'.over.chain')
  chain <- import.chain(chain_file_path)
  
  source_iranges <- IRanges(start = data$chromStart, end = data$chromEnd)
  source_granges <- GRanges(seqnames = data$Chromosome, ranges = source_iranges)
  
  #result <- data
  #colnames(result) <- paste0('source_', colnames(result))
  
  i <- c(1:length(source_granges))
  target_granges_lst <- lapply(i, function(x) {
    target_grange_lst <- liftOver(source_granges[x], chain)
    target_grange <- unlist(target_grange_lst)
    if(length(target_grange) == 0) {
      result <- data.frame('target_Chromosome' = NA,
                           'target_chromStart' = NA,
                           'target_chromEnd' = NA,
                           stringsAsFactors = FALSE)
    } else {
      result <- data.frame('target_Chromosome' = seqnames(target_grange),
                           'target_chromStart' = start(target_grange),
                           'target_chromEnd' = end(target_grange),
                           stringsAsFactors = FALSE)
    }    
    return(result)
  })
  
  #una region de source puede mapear a mas de una region en target con liftover!!

  result <- data.frame('source_Chromosome' = character(),
                       'source_chromStart' = numeric(),
                       'source_chromEnd' = numeric(),
                       'target_Chromosome' = character(),
                       'target_chromStart' = numeric(),
                       'target_chromEnd' = numeric(),
                       stringsAsFactors = FALSE)
  
  source_index <- c(1: length(source_granges))
  
  target_data_lst <- lapply(source_index, function(i) {
    target_data <- target_granges_lst[[i]]
    target_data$source_Chromosome <- data$Chromosome[i]
    target_data$source_chromStart <- data$chromStart[i]
    target_data$source_chromEnd <- data$chromEnd[i]
    return(target_data)
  })
  
  result_raw <- do.call(rbind, target_data_lst)
  
  result <- result_raw[complete.cases(result_raw), ]
  
  
  #hay que checkear que no se pierdan muchos elementos en liftover, en este caso se pierde el 90%
  # pase de 600 a 60 para el chr 1
    
  return(result)
  
}

plot_acc_track <- function(hg38_data, galgal4_data) {

  
  hg38_data$Data <- 100
  hg38_data$Specie <- 'H'
  data_cols <- c("Chromosome", "chromStart", "chromEnd", "Data")
  #hg38_data <- hg38_data[hg38_data$Chromosome == "chr1", data_cols]
  hg38_data <- hg38_data[, data_cols]
  hg38_data$Chromosome <- paste0('H', hg38_data$Chromosome)
  
  colnames(galgal4_data) <- c('Chromosome', 'chromStart', 'chromEnd')
  galgal4_data$Data <- 100
  galgal4_data$Specie <- 'G'
  data_cols <- c("Chromosome", "chromStart", "chromEnd", "Data")
  galgal4_data <- galgal4_data[, data_cols]
  galgal4_data$Chromosome <- paste0('G', galgal4_data$Chromosome)
  
  
  poligon_data <- rbind(hg38_data, galgal4_data)
  #poligon_data <- UCSC.HG38.Human.CytoBandIdeogram[UCSC.HG38.Human.CytoBandIdeogram$Chromosome == 'chr1', ]
  #poligon_data$Data <- 1
  #poligon_data <- poligon_data[, data_cols]
  #RCircos.Polygon.Plot(hg38_data, data.col = 4, genomic.columns = 3, 
  #                    polygon.col="blue",
  #                    track.num, side)
  
  RCircos.Histogram.Plot(hist.data=poligon_data, data.col=4, min.value = 0, max.value = 100,
                         track.num=1, side='in', genomic.columns=3)
  
    
  
}

plot_synteny_track <- function(synteny_data) {
  
  #los bloques de sintenia deberian tener el color del cromosoma de humano
  
  synteny_blocks_count <- nrow(synteny_data)
  
  synteny_blocks_target <- data.frame("Chromosome" = character(synteny_blocks_count),
                                      "chromStart" = numeric(synteny_blocks_count),
                                      "chromEnd" = numeric(synteny_blocks_count))
  
  synteny_blocks_source <- data.frame("Chromosome" = character(synteny_blocks_count),
                                      "chromStart" = numeric(synteny_blocks_count),
                                      "chromEnd" = numeric(synteny_blocks_count))
  
  synteny_blocks_target$Chromosome <- paste0('H', synteny_data$tName)
  synteny_blocks_target$chromStart <- synteny_data$tStart
  synteny_blocks_target$chromEnd <- synteny_data$tEnd
  
  synteny_blocks_source$Chromosome <- paste0('G', synteny_data$qName)
  synteny_blocks_source$chromStart <- synteny_data$qStart
  synteny_blocks_source$chromEnd <- synteny_data$qEnd
  
  synteny_blocks <- rbind(synteny_blocks_target, synteny_blocks_source)
  
  galgal4_chrs <- paste0('chr', c(1:32, 'W', 'Z'))
  valid_query_chrs <- paste0('G', galgal4_chrs)
  valid_target_chrs <- paste0('Hchr', c(1:22))
  valid_chrs <- c(valid_query_chrs, valid_target_chrs)
  valid_mask <- synteny_blocks$Chromosome %in% valid_chrs
  synteny_blocks_clean <- synteny_blocks[valid_mask, ]
  
  synteny_blocks_clean$Data <- 100
  
  RCircos.Histogram.Plot(hist.data=synteny_blocks_clean, data.col=4, min.value = 0, max.value = 100,
                         track.num=2, side='in', genomic.columns=3)
  
  
}

plot_synteny_links <- function(synteny_data) {
  
  synteny_data$qName <- paste0('G', synteny_data$qName)
  synteny_data$tName <- paste0('H', synteny_data$tName)
  
  galgal4_chrs <- paste0('chr', c(1:32, 'W', 'Z'))
  valid_query_chrs <- paste0('G', galgal4_chrs)
  valid_target_chrs <- paste0('Hchr', c(1:22))
  # valid_chrs <- c(valid_query_chrs, valid_target_chrs)
  valid_mask <- (synteny_data$qName %in% valid_query_chrs) & 
                (synteny_data$tName %in% valid_target_chrs)
  
  synteny_data_clean <- synteny_data[valid_mask, ]
  
  
  RCircos.Link.Plot(link.data = synteny_data_clean, track.num = 3, genomic.columns = 3)
  
}

prepare_data <- function(target_chrs) {
  
  # devolver todos los datasets necesarios con las validaciones de chr validos
  
}

init_plot <- function(target_chrs) {

  
  data(UCSC.HG38.Human.CytoBandIdeogram)
  
  col_names <- c('Chromosome', 'chromStart', 'chromEnd', 'Name', 'Stain')
  galgal4.CytoBandIdeogram_data_path <- '/media/paulati/Nuevo vol/paula/ingebi/2022/paula/precision/synteny/data/cytoBandIdeo.txt'
  galgal4.CytoBandIdeogram_all <- read.table(galgal4.CytoBandIdeogram_data_path, sep = "\t", stringsAsFactors = FALSE)
  colnames(galgal4.CytoBandIdeogram_all) <- col_names
  
  galgal4_chrs <- paste0('chr', c(1:32, 'W', 'Z'))
  mask <- galgal4.CytoBandIdeogram_all$Chromosome %in% galgal4_chrs
  galgal4.CytoBandIdeogram <- galgal4.CytoBandIdeogram_all[mask, ]
  galgal4.CytoBandIdeogram$Chromosome <- factor(galgal4.CytoBandIdeogram$Chromosome, levels = galgal4_chrs)
  
  galgal4.CytoBandIdeogram <- galgal4.CytoBandIdeogram[order(galgal4.CytoBandIdeogram$Chromosome), ]
  
  mask <- UCSC.HG38.Human.CytoBandIdeogram$Chromosome %in% target_chrs
  
  cyto.list <- list(UCSC.HG38.Human.CytoBandIdeogram[mask, ],
                    galgal4.CytoBandIdeogram)
  
  species <- c("H", "G")
  RCircos.Multiple.Species.Core.Components(
    cyto.list, species, chr.exclude=NULL,
    tracks.inside=10, tracks.outside=0)
  
  # RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram, NULL, tracks.inside = 10, tracks.outside = 0)
  # RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram[UCSC.HG38.Human.CytoBandIdeogram$Chromosome == 'chr1', ], NULL, tracks.inside = 10, tracks.outside = 0)
  
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.params$radiu.len <- 1.5;
  RCircos.Reset.Plot.Parameters(rcircos.params);
  
  # 6. Open graphic device:
  
  RCircos.Set.Plot.Area();
  
  #or submit your own code. For example: 
  
  par(mai=c(0.25, 0.25, 0.25, 0.25));
  plot.new();
  plot.window(c(-2.5,2.5), c(-2.5, 2.5));
  
  # 7. Call plot function to plot each data track:
  
  RCircos.Chromosome.Ideogram.Plot();
  
}


base_path <- '/media/paulati/Nuevo vol/paula/ingebi/2022/paula/precision/synteny/data/out'

file_name_1 <- 'hg38_galGal4.psl'
file_name_2 <- 'hg38_ornAna1.psl'
file_name_3 <- 'sarc_mammals_chr22_acc_50.bed'


# http://genome.ucsc.edu/FAQ/FAQformat#format2
# matches - Number of bases that match that aren't repeats
# misMatches - Number of bases that don't match
# repMatches - Number of bases that match but are part of repeats
# nCount - Number of "N" bases
# qNumInsert - Number of inserts in query
# qBaseInsert - Number of bases inserted in query
# tNumInsert - Number of inserts in target
# tBaseInsert - Number of bases inserted in target
# strand - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
# qName - Query sequence name
# qSize - Query sequence size.
# qStart - Alignment start position in query
# qEnd - Alignment end position in query
# tName - Target sequence name
# tSize - Target sequence size
# tStart - Alignment start position in target
# tEnd - Alignment end position in target
# blockCount - Number of blocks in the alignment (a block contains no gaps)
# blockSizes - Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
# qStarts - Comma-separated list of starting positions of each block in query
# tStarts - Comma-separated list of starting positions of each block in target
psl_col_names <- c('matches', 'misMatches', 'repMatches', 'nCount', 
                   'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 
                   'strand', 
                   'qName', 'qSize', 'qStart', 'qEnd',
                   'tName', 'tSize', 'tStart', 'tEnd', 
                   'blockCount', 'blockSizes', 'qStarts', 'tStarts')

bed_col_names <- c('Chromosome', 'chromStart', 'chromEnd')

file_path_1 <- file.path(base_path, file_name_1)
file_path_2 <- file.path(base_path, file_name_2)
file_path_3 <- file.path(base_path, file_name_3)

data_1 <- read.table(file_path_1)
colnames(data_1) <- psl_col_names

data_2 <- read.table(file_path_2)
colnames(data_2) <- psl_col_names

#data_3 <- read.table(file_path_3)
data_3 <- create_acc_mock_data()
colnames(data_3) <- bed_col_names

acc_mockup_file_path <- file.path(dirname(base_path), 'acc_mockup_hg38.csv')
write.table(data_3, acc_mockup_file_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

# source_specie <- 'hg38' # tolower
# target_specie <- 'GalGal4' # to upper first
data_4 <- liftover_acc_data(data_3, 'hg38', 'GalGal4')


######################################################

# multiple species
# ----------------

# https://rdrr.io/cran/RCircos/man/RCircos.Multiple.Species.Dataset.html
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3765848/

#target_chrs <- paste0('chr', c(1:22))
target_chrs <- paste0('chr22')

init_plot(target_chrs)


# custom track:
# side <- "in"  # En que lado del ideograma estará: "in" or "out"

#track.num <- 1  # La posición del track 

# mask <- RCircos.Gene.Label.Data$Chromosome != 'chr9'
# sum(mask) / length(mask)
# genes_data <- RCircos.Gene.Label.Data[mask, ]
#RCircos.Gene.Connector.Plot(genes_data, track.num, side)

#name.col <- 4  #indicamos la columna en la que está el nombre
#RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side)

# data(RCircos.Polygon.Data)

hg38_mask <- data_3$Chromosome %in% target_chrs
hg38_data <- data_3[hg38_mask, ]
galgal4_data <- data_4[, c('target_Chromosome', 'target_chromStart', 'target_chromEnd')]
plot_acc_track(hg38_data, galgal4_data) 


# cinta uniendo regiones con color asociado igual al chr de humano

#data_4$source_Chromosome <- paste0('H', data_4$source_Chromosome)
#data_4$target_Chromosome <- paste0('G', data_4$target_Chromosome)
#RCircos.Link.Plot(link.data = data_4, track.num = 2, genomic.columns = 3)

# ver cuales son las limitaciones de liftover poruqe hay regiones que no traduce
# ver la referencia del liftover mas actual que es del paquete hal

#synteny blocks
synteny_data <- data_1[, c('qName', 'qStart', 'qEnd', 'tName', 'tStart', 'tEnd' )]
plot_synteny_track(synteny_data)

plot_synteny_links(synteny_data)


#unique(synteny_blocks$Chromosome)

# synteny?

# creo que solo habia calculado para chr22 de humano (como target)
#unique(data_1$qName)
#unique(data_1$tName)
#synteny_data <- data_1[, c('qName', 'qStart', 'qEnd', 'tName', 'tStart', 'tEnd' )]
#synteny_data$qName <- paste0('H', synteny_data$qName)
#synteny_data$tName <- paste0('G', synteny_data$tName)
#RCircos.Link.Plot(link.data = synteny_data, track.num = 3, genomic.columns = 3)




# dev.off();

######################################################

# https://www.researchgate.net/publication/260196003_OmicCircos_A_Simple-to-Use_R_Package_for_the_Circular_Visualization_of_Multidimensional_Omics_Data/figures?lo=1

