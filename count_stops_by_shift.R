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
  
  shift_0 <- count_stops(seq,0)
  shift_1 <- count_stops(seq,1)
  shift_2 <- count_stops(seq,2)
  
  data.frame(gene_length, shift_0, shift_1, shift_2)
}




tx_id <- sample(tx2g$transcript_id[tx2g$transcript_biotype == "protein_coding"], 2000)

seqs <- tx_id |>
  lapply(get_seq)

any(is.na(seqs))
# table(is.na(seqs))
# tx_id[is.na(seqs)]
# seqs <- seqs[!is.na(seqs)]

shifts <- seqs |>
  lapply(nb_of_stop_codons) |>
  do.call(what = rbind, args = _)
shifts$tx_id <- tx_id
shifts$gene_id <- wb_tx2g(shifts$tx_id, tx2g_tab = tx2g, warn_missing = TRUE)

# Should all have a single stop codon without frameshift
# ctb-1/MTCE.21 weird, there may be others
shifts[shifts$shift_0 != 1,]
shifts <- shifts[shifts$shift_0 == 1,]


library(ggplot2)


ggplot(shifts) + theme_classic() +
  geom_point(aes(shift_1, shift_2)) +
  geom_abline(slope = 1, color = "grey") +
  coord_equal() +
  # ggrepel::geom_text_repel(aes(shift_1, shift_2, label = tx_id)) +
  scale_x_log10() + scale_y_log10()

ggplot(shifts) +theme_classic() +
  geom_histogram(aes(x = shift_2/shift_1), bins = 100, color = "white") +
  scale_x_log10()

hist(log(shifts$shift_2/shifts$shift_1), breaks = 100)

# there is a systematic bias: frameshifts of 2  bp tend to induce more stop codons
# than frameshifts of 1 bp.
# However, not clearly mutually exclusive: in general genes with more of 1 also 
# have more of the other (presumably because longer, see below), but the initial hypothesis
# that 'some genes "prefer" a 1 bp shift whereas others "prefer" 2 bp' is not correct.

ggplot(shifts) + theme_classic() +
  geom_point(aes(gene_length, shift_2)) +
  # ggrepel::geom_text_repel(aes(shift_1, shift_2, label = tx_id)) +
  scale_x_log10() + scale_y_log10()

ggplot(shifts) + theme_classic() +
  geom_point(aes(gene_length, shift_1)) +
  # ggrepel::geom_text_repel(aes(shift_1, shift_2, label = tx_id)) +
  scale_x_log10() + scale_y_log10()


mod <- lm(shift_2 ~ shift_1 + gene_length, data = shifts)

plot(fitted.values(mod), residuals(mod))
plot(log(fitted.values(mod)), residuals(mod))
plot(log(fitted.values(mod)), log(residuals(mod)), xlim=c(-2,6))
mod
summary(mod)

ggplot(shifts) + theme_classic() +
  geom_point(aes(x = shift_1/gene_length, y = shift_2/gene_length),
             alpha = .2) +
  geom_abline(slope = 1, color = "grey")

# Once removing the influence of gene_length, this relationship is even clearer




