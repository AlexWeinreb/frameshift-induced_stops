# Do frameshifts of 1 bp induce more stop codons than frameshifts of 2 bp?



library(wbData)
tx2g <- wb_load_tx2gene(284)



# Use Wormbase's RESTful API to recover the CDS sequence
get_seq <- function(tx_id){
  cds_id <- stringr::str_remove(tx_id, "\\.[[:digit:]]+$")
  
  xx <- httr::GET(paste0("https://wormbase.org/rest/field/cds/",cds_id,"/cds_sequence"))
  if(xx$status_code != 200) return(NA_character_)
  
  json <- httr::content(xx, as = "text",encoding = "UTF-8") |> jsonlite::fromJSON()
  
  switch(json$cds_sequence$data$strand,
         `+` = json$cds_sequence$data$positive_strand$sequence,
         `-` = json$cds_sequence$data$negative_strand$sequence,
         NA_character_
  )
}

# count nb of stop codons
count_stops <- function(seq, shift){
  seq <- substring(seq, shift+1)
  n <- floor(nchar(seq)/3)
  codons <- substring(seq, 3*(1:n) -2, 3*(1:n))
  sum(codons == "TAA" | codons == "TAG" | codons == "TGA")
}

# Compute nb of stops for different frameshifts
nb_of_stop_codons <- function(seq){
  gene_length <- nchar(seq)
  
  if(count_stops(seq,0) != 1) warning("Incorrect number of codons without shift.")
  
  shift_1 <- count_stops(seq,1)
  
  shift_2 <- count_stops(seq,2)
  
  data.frame(gene_length, shift_1, shift_2)
}




tx_id <- sample(tx2g$transcript_id[tx2g$transcript_biotype == "protein_coding"], 2000)

seqs <- tx_id |>
  lapply(get_seq)

any(is.na(seqs))
# table(is.na(seqs))
# tx_id[is.na(seqs)]
# seqs <- seqs[!is.na(seqs)]

seqs_as_bs <- seqs |>
  lapply(Biostrings::DNAString)

nb_codons <- seqs_as_bs |>
  lapply(nb_of_stop_codons)

shifts <- do.call(what = rbind, args = nb_codons)
shifts$tx_id <- tx_id


library(ggplot2)


ggplot(shifts) + theme_classic() +
  geom_point(aes(shift_1, shift_2)) +
  geom_abline(slope = 1, color = "grey") +
  coord_equal() +
  ggrepel::geom_text_repel(aes(shift_1, shift_2, label = tx_id)) +
  scale_x_log10() + scale_y_log10()

hist(shifts$shift_2/shifts$shift_1, breaks = 100)



# xx2 <- bench::press(n = seq(from = 1, to = 25, by = 5),
#                     bench::mark(a = seqs_as_bs[sample(500, n)] |>
#                                   lapply(nb_of_stop_codons)))
# 
# 
# plot(xx2)
# plot(xx2$n, xx2$median)
# 
