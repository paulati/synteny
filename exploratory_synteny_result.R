#https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomeInfoDb/html/fetchExtendedChromInfoFromUCSC.html
#BiocManager::install(version = "3.14")
#BiocManager::install('GenomeInfoDb')

# Be aware that the coordinates for a negative strand in a dna query PSL line are handled in a special way. 
# In the qStart and qEnd fields, the coordinates indicate the position where the query matches from the point of view of the forward strand, 
# even when the match is on the reverse strand. However, in the qStarts list, the coordinates are reversed.

###########################
# ** First position is 0 **
###########################

library(stringr)
library(GenomeInfoDb)

galGal3_chrominfo <- getChromInfoFromUCSC("galGal3")
galGal3_chrominfo

hg38_chrominfo <- getChromInfoFromUCSC("hg38")
hg38_chrominfo

# no está
#ornAna1_chrominfo <- getChromInfoFromUCSC("ornAna1")
#ornAna1_chrominfo

##########

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


base_path <- '/media/paulati/Nuevo vol/paula/ingebi/2022/paula/precision/synteny/data/out'

file_name_1 <- 'hg38_galGal4.psl'
file_name_2 <- 'hg38_ornAna1.psl'

psl_col_names <- c('matches', 'misMatches', 'repMatches', 'nCount', 
                   'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 
                   'strand', 
                   'qName', 'qSize', 'qStart', 'qEnd',
                   'tName', 'tSize', 'tStart', 'tEnd', 
                   'blockCount', 'blockSizes', 'qStarts', 'tStarts')

file_path_1 <- file.path(base_path, file_name_1)
file_path_2 <- file.path(base_path, file_name_2)


data_1 <- read.table(file_path_1)
colnames(data_1) <- psl_col_names

data_2 <- read.table(file_path_2)
colnames(data_2) <- psl_col_names

# matches - Number of bases that match that aren't repeats
hist(data_1$matches)
hist(data_2$matches)

# misMatches - Number of bases that don't match
unique(data_1$misMatches)
unique(data_2$misMatches)

# repMatches - Number of bases that match but are part of repeats
unique(data_1$repMatches)
unique(data_2$repMatches)

# nCount - Number of "N" bases
unique(data_1$nCount)
unique(data_2$nCount)

# qNumInsert - Number of inserts in query
hist(data_1$qNumInsert)
hist(data_2$qNumInsert)

# qBaseInsert - Number of bases inserted in query
hist(data_1$qBaseInsert)
hist(data_2$qBaseInsert)

plot(data_1$qNumInsert, data_1$qBaseInsert)
plot(data_2$qNumInsert, data_2$qBaseInsert)

# tNumInsert - Number of inserts in target
hist(data_1$tNumInsert)
hist(data_2$tNumInsert)

# tBaseInsert - Number of bases inserted in target
hist(data_1$tBaseInsert)
hist(data_2$tBaseInsert)

plot(data_1$tNumInsert, data_1$tBaseInsert)
plot(data_2$tNumInsert, data_2$tBaseInsert)

# strand - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
table(data_1$strand)
table(data_2$strand)

# qName - Query sequence name
table(data_1$qName)
table(data_2$qName) 
# TODO: que son los contig? los descarto?

# qSize - Query sequence size.
hist(data_1$qSize)
plot(c(1:nrow(data_1), data_1$qSize))
table(data_1$qSize)
nrow(table(data_1$qSize))
hist(data_2$qSize)
plot(c(1:nrow(data_2), data_2$qSize))
table(data_2$qSize)

# cada qSize corresponde a un chromosoma en la query (galGal3 en data_1)
unique_qSize <- unique(data_1$qSize)
qSize_chrom <- data.frame(
  'qName' = character(length(unique_qSize)),
  'qSize' = unique_qSize,
  stringsAsFactors = FALSE
)

for (i in c(1:nrow(qSize_chrom))) {
  mask <- data_1$qSize == qSize_chrom$qSize[i]
  qSize_chrom$qName[i] <- unique(data_1$qName[mask])
}

extended_qSize_chrom <- merge(qSize_chrom, galGal3_chrominfo,
                              by.x = 'qName', by.y = 'chrom', 
                              all.x = TRUE, all.y = FALSE)

extended_qSize_chrom$diff <- extended_qSize_chrom$size - extended_qSize_chrom$qSize
# TODO: no siempre en tamaño de query es menor que el del cromosoma de galGal3
# ver que significa


# qStart - Alignment start position in query
# qEnd - Alignment end position in query

#esto no da -> OK!!!
mask <- data_1$qSize == data_1$qEnd - data_1$qStart
sum(mask) / nrow(data_1)

# qSize no tiene que ver con esto, qSize es parecido al tamaño de cada cromosoma
plot(data_1$qStart, data_1$qEnd)
calc_qSize <- data_1$qEnd - data_1$qStart 
plot(calc_qSize, data_1$qSize)
hist(calc_qSize)

# qSize no tiene que ver con esto, qSize es parecido al tamaño de cada cromosoma
plot(data_2$qStart, data_2$qEnd)
calc_qSize <- data_2$qEnd - data_2$qStart
plot(calc_qSize, data_2$qSize)

# tName - Target sequence name
table(data_1$tName)
table(data_2$tName) 

# tSize - Target sequence size
# tSize es el tamaño del cromosoma de humano en data_1
hist(data_1$tSize)
plot(c(1:nrow(data_1)), data_1$tSize)
table(data_1$tSize)
hist(data_2$tSize)
plot(c(1:nrow(data_2)), data_2$tSize)
table(data_2$tSize)
unique(data_2$tSize) 

#largo del cromosoma 22?
# chr22	50,818,468, coincide con tSize



# tStart - Alignment start position in target
# tEnd - Alignment end position in target


# blockCount - Number of blocks in the alignment (a block contains no gaps)
hist(data_1$blockCount)
hist(data_2$blockCount)


# blockSizes: size of each of the blockCount blocks in the record
calc_blockCount_lst <- sapply(data_1$blockSizes, function(x) {
  x_clean <- str_trim(x, side = 'both')
  parts <- unlist(str_split(x_clean, pattern = ","))
  valid_parts <- unlist(lapply(parts, function(x) nchar(x) > 0))
  result <- sum(valid_parts)
  return(result)
})
calc_blockCount <- as.integer(calc_blockCount_lst)
mask <- calc_blockCount == data_1$blockCount
sum(mask) / length(calc_blockCount)


calc_blockCount_lst <- sapply(data_2$blockSizes, function(x) {
  x_clean <- str_trim(x, side = 'both')
  parts <- unlist(str_split(x_clean, pattern = ","))
  valid_parts <- unlist(lapply(parts, function(x) nchar(x) > 0))
  result <- sum(valid_parts)
  return(result)
})
calc_blockCount <- as.integer(calc_blockCount_lst)
mask <- calc_blockCount == data_2$blockCount
sum(mask) / length(calc_blockCount)



# qStarts - Comma-separated list of starting positions of each block in query

qStarts_count_lst <- sapply(data_1$qStarts, function(x) {
  x_clean <- str_trim(x, side = 'both')
  parts <- unlist(str_split(x_clean, pattern = ","))
  valid_parts <- unlist(lapply(parts, function(x) nchar(x) > 0))
  result <- sum(valid_parts)
  return(result)
})

qStarts_count <- unlist(qStarts_count_lst)
mask <- qStarts_count == data_1$blockCount
sum(mask) / length(data_1$blockCount)

# cada elemento de la lista qStarts esta entre qStart y qEnd:

calc_qStarts_lst <- apply(data_1, MARGIN = 1, function(row) {
  x <- row['qStarts']
  
  x_clean <- str_trim(x, side = 'both')
  parts <- unlist(str_split(x_clean, pattern = ","))
  valid_parts <- unlist(lapply(parts, function(x) nchar(x) > 0))
  positions <- as.integer(parts[valid_parts])

  start <- as.integer(row['qStart'])
  end <- as.integer(row['qEnd'])
  
  
  mask <- sapply(positions, function(e) e <= end & e >= start)
  
  result <- sum(mask) / length(positions)
  return(result)
})

sum(calc_qStarts_lst) / nrow(data_1)

###

calc_qStarts_lst <- apply(data_2, MARGIN = 1, function(row) {
  x <- row['qStarts']
  
  x_clean <- str_trim(x, side = 'both')
  parts <- unlist(str_split(x_clean, pattern = ","))
  valid_parts <- unlist(lapply(parts, function(x) nchar(x) > 0))
  positions <- as.integer(parts[valid_parts])
  
  start <- as.integer(row['qStart'])
  end <- as.integer(row['qEnd'])
  
  
  mask <- sapply(positions, function(e) e <= end & e >= start)
  
  result <- sum(mask) / length(positions)
  return(result)
})

sum(calc_qStarts_lst) / nrow(data_2)



# tStarts - Comma-separated list of starting positions of each block in target

# cada elemento de la lista tStarts esta entre tStart y tEnd: NO!
# hay que transformar las coordenadas de las que están en el strand negativo:
# Be aware that the coordinates for a negative strand in a dna query PSL line are handled in a special way. 
# In the qStart and qEnd fields, the coordinates indicate the position where the query matches from the point of view of the forward strand, 
# even when the match is on the reverse strand. However, in the qStarts list, the coordinates are reversed.


calc_tStarts_lst <- apply(data_1, MARGIN = 1, function(row) {
  x <- row['tStarts']
  
  x_clean <- str_trim(x, side = 'both')
  parts <- unlist(str_split(x_clean, pattern = ","))
  valid_parts <- unlist(lapply(parts, function(x) nchar(x) > 0))
  positions <- as.integer(parts[valid_parts])
  
  strand <- as.character(row['strand'])
  qStrand <- str_sub(strand, start = 1, end = 1)
  tStrand <- str_sub(strand, start = 2, end = 2)
  
  
  tSize <- as.integer(row['tSize'])
  # tStarty tEnd siempre tienen las coordenadas dek strand +:
  start <- as.integer(row['tStart'])
  end <- as.integer(row['tEnd'])
  
  if(tStrand == '-') {
  
    #ahora son tEnds:
    positions <- tSize - positions
    
  } else {
    # do nothing
    
  }
  
  
  mask <- sapply(positions, function(e) e <= end & e >= start)
  
  result <- sum(mask) / length(positions)
  return(result)
})

sum(calc_tStarts_lst) / nrow(data_1)

# con la transformacion de coordenadas da OK!! todos estan comprendidios entre start y end


#Negative-strand-coordinate-qStart = qSize - qEnd   = 61 - 56 =  5
#Negative-strand-coordinate-qEnd   = qSize - qStart = 61 -  4 = 57


###


calc_tStarts_lst <- apply(data_2, MARGIN = 1, function(row) {
  x <- row['tStarts']
  
  x_clean <- str_trim(x, side = 'both')
  parts <- unlist(str_split(x_clean, pattern = ","))
  valid_parts <- unlist(lapply(parts, function(x) nchar(x) > 0))
  positions <- as.integer(parts[valid_parts])
  
  strand <- as.character(row['strand'])
  qStrand <- str_sub(strand, start = 1, end = 1)
  tStrand <- str_sub(strand, start = 2, end = 2)
  
  
  tSize <- as.integer(row['tSize'])
  # tStarty tEnd siempre tienen las coordenadas dek strand +:
  start <- as.integer(row['tStart'])
  end <- as.integer(row['tEnd'])
  
  if(tStrand == '-') {
    
    #ahora son tEnds:
    positions <- tSize - positions
    
  } else {
    # do nothing
    
  }
  
  
  mask <- sapply(positions, function(e) e <= end & e >= start)
  
  result <- sum(mask) / length(positions)
  return(result)
})

sum(calc_tStarts_lst) / nrow(data_2)



###################


# https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomeInfoDb/html/fetchExtendedChromInfoFromUCSC.html


# 32 cromosomas, W y Z


library(BSgenome)
genome <- getBSgenome("GRCh38")  # this loads the BSgenome.Hsapiens.NCBI.GRCh38 package

## A quick look at the GRCh38 seqlevels:
length(seqlevels(genome))

hg38_chrominfo <- getChromInfoFromUCSC("hg38")



