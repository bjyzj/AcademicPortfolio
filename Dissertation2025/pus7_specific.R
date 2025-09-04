
library(tidyverse)
library(biomaRt)
library(mclust)
library(readxl)
library(rentrez)
library(biomaRt)
library(rtracklayer)
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38) 
###################
#HeLa DATASETS
###################

# bid pus7 positive site filteration
# match pos form shcontrol w/ shPUS7 to get IDs

#pus7
shpus7_bid <- read_excel("/Users/beyzaerkal/Desktop/datasets_proj/bid/BID_aver_pus7.xlsx", skip = 1, col_names = TRUE)
shpus7_bid <- shpus7_bid %>% rename(pos = position)
colnames(shpus7_bid) <- gsub(" ", "_", colnames(shpus7_bid))

#control 
shcontrol_bid <- read_excel("/Users/beyzaerkal/Desktop/datasets_proj/bid/HeLa_mRNA_shControl_BID-seq.xlsx", skip = 3, col_names = TRUE)
colnames(shcontrol_bid) <- gsub(" ", "_", colnames(shcontrol_bid))

#align
aligned_pos <- merge(shpus7_bid, shcontrol_bid, by = "pos") # 98

# refseq for bid
refseq_id_bid <- unique(aligned_pos$refseq)

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(attributes = c("refseq_mrna", "ensembl_transcript_id_version", "ensembl_gene_id", "external_gene_name"),
  filters = "refseq_mrna", values = refseq_id_bid, mart = mart) # 112

table(mapping$refseq_mrna)
sum(duplicated(mapping$refseq_mrna))

# fasta file of 121nt 
bid_fasta <- readLines("/Users/beyzaerkal/Desktop/datasets_proj/bid/extraction_site_output3rna.fa")

headers <- bid_fasta[grepl("^>", bid_fasta)]
headers <- sub("^>", "", headers)
fasta_ids <- sub(":.*", "", headers)
fasta_ids <- trimws(fasta_ids)

#align and filtering
matched_map <- mapping[mapping$ensembl_transcript_id_version %in% fasta_ids, ]
# keep only match and same order as fasta id
matched_map <- matched_map[match(fasta_ids, matched_map$ensembl_transcript_id_version), ]
#too many 0 so smaple size bigegr here 
# pus7 samples 
matched_map <- matched_map[!is.na(matched_map$ensembl_transcript_id_version), ]
#95

##############################
# set pus7 positive threshold / decision
# symbol with _pct in all column names
colnames(aligned_pos) <- gsub("%", "pct", colnames(aligned_pos))

# diff show the chanegs in pseudouridine
# to tackle the gmm error made a new variable
aligned_pos2 <- aligned_pos
aligned_pos2$diff <- aligned_pos2$Average_Ψ_fraction_in_shControl_pct - aligned_pos2$Average_Ψ_fraction_in_shPUS7_pct

# distribution viz
hist(aligned_pos$diff, breaks=50, main="ΔΨ: shControl - shPUS7", xlab="Ψ Difference", col="skyblue")
abline(v=0, col="red", lty=2)

# GMM application on identifying best threshold
gmm <- Mclust(aligned_pos$diff, G=2)  # Two-component model

# mean and shared variance
mu1 <- gmm$parameters$mean[1]
mu2 <- gmm$parameters$mean[2]
sigma <- sqrt(gmm$parameters$variance$sigmasq[1])  # shared variance

# intersection point
threshold <- uniroot(function(x) {
  dnorm(x, mean = mu1, sd = sigma) - dnorm(x, mean = mu2, sd = sigma)
}, lower = min(aligned_pos$diff), upper = max(aligned_pos$diff))$root
# threshold
cat("Suggested threshold:", round(threshold, 2), "\n") # 17.88


# label
aligned_pos$label <- ifelse(aligned_pos$diff > 17.88, "positive", "background")
positive_ids <- aligned_pos$refseq[aligned_pos$label == "positive"]

matched_map$label <- ifelse(matched_map$refseq_mrna %in% positive_ids, "PUS7_positive", "background")

# save separate variables
positive_sites <- matched_map[matched_map$label == "PUS7_positive", ] #39
background_sites <- matched_map[matched_map$label == "background", ] # 56 # not necessarily negative


# get positives and save as fasta

fasta_headers <- grep("^>", bid_fasta, value = TRUE) # lines start with >
clean_headers <- trimws(sub("^>", "", fasta_headers))# clean >, gaps
clean_headers <- sub(":.*", "", fasta_ids)

positive_ids <- positive_sites$ensembl_transcript_id_version

#clean_positive_ids <- sub("\\..*$", "", positive_ids)
#clean_headers_nover <- sub("\\..*$", "", clean_headers)

pos_header_idx <- which(clean_headers %in% positive_ids) # index no
matched_headers <- clean_headers[pos_header_idx]

length(pos_header_idx) #39
#############

# getting the sequences and other info back function

extract_fasta_seqs <- function(fasta_lines, header_indices) {
  extracted <- character()
  headers <- grep("^>", fasta_lines, value = TRUE)
  for (i in header_indices) {
    start <- which(fasta_lines == headers[i])
    # next header index after start
    next_headers <- which(grepl("^>", fasta_lines) & seq_along(fasta_lines) > start)
    end <- if(length(next_headers) == 0) length(fasta_lines) else min(next_headers) - 1
    extracted <- c(extracted, fasta_lines[start:end])
  }
  return(extracted)
}

positive_fasta_seqs <- extract_fasta_seqs(bid_fasta, pos_header_idx)

# save as fasta
writeLines(positive_fasta_seqs, "/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_positive_sequences.fa")


# backround

background_ids <- background_sites$ensembl_transcript_id_version
# Find header index of ID
background_header_index <- which(clean_headers %in% background_ids)
length(background_header_index) # 56

background_fasta_seqs <- extract_fasta_seqs(bid_fasta, neg_header_index)
writeLines(background_fasta_seqs, "/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_background_sequences.fa")


# checks of U and 121nt

pos_seqs_bid <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_positive_sequences.fa")
pos_lengths <- width(pos_seqs_bid)
table(pos_lengths)# 121 39

pos_mid_base <- substring(as.character(pos_seqs_bid), 61, 61)
table(pos_mid_base)  # 14 is U


bg_seqs <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_background_sequences.fa")

bg_lengths <- width(bg_seqs)
table(bg_lengths) # 121 56

bg_mid_base <- substring(as.character(bg_seqs), 61, 61)
table(bg_mid_base) # 14 U
# create new ones
pos_seqs_bid <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_positive_sequences.fa")
pos_mid_base <- substring(as.character(pos_seqs_bid), 61, 61)
pos_seqs_U <- pos_seqs_bid[pos_mid_base == "U"]
writeXStringSet(pos_seqs_U, "/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_positive_sequences_middleU.fa")

# Load background sequences
bg_seqs <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_background_sequences.fa")
bg_mid_base <- substring(as.character(bg_seqs), 61, 61)
bg_seqs_U <- bg_seqs[bg_mid_base == "U"]
writeXStringSet(bg_seqs_U, "/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_background_sequences_middleU.fa")


pos_seqs_bid <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_positive_sequences_middleU.fa")
pos_lengths <- width(pos_seqs_bid)
table(pos_lengths)# 121 14

pos_mid_base <- substring(as.character(pos_seqs_bid), 61, 61)
table(pos_mid_base)  # 14 is U


bg_seqs <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/bid/PUS7_background_sequences_middleU.fa")
# FINAL DATASET - positive
bg_lengths <- width(bg_seqs)
table(bg_lengths) # 121 14

bg_mid_base <- substring(as.character(bg_seqs), 61, 61)
table(bg_mid_base) # 14 U



# ######## #
#negatiev bid 

shcontrol_bid <- read_excel("/Users/beyzaerkal/Desktop/datasets_proj/bid/HeLa_mRNA_shControl_BID-seq.xlsx", skip = 3, col_names = TRUE)
colnames(shcontrol_bid) <- gsub(" ", "_", colnames(shcontrol_bid))

refseq_id_control <- unique(shcontrol_bid$refseq)
# map enst
control_mapping <- getBM(attributes = c("refseq_mrna", "ensembl_transcript_id_version", 
                                        "ensembl_gene_id", "external_gene_name"),
                         filters = "refseq_mrna", values = refseq_id_control, mart = mart)

# match to fasta id # load FASTA
bid_fasta <- readLines("/Users/beyzaerkal/Desktop/datasets_proj/bid/extraction_site_output3rna.fa")

fasta_headers <- bid_fasta[grepl("^>", bid_fasta)]
fasta_ids <- sub("^>", "", fasta_headers)
fasta_ids <- sub(":.*", "", fasta_ids)
fasta_ids <- trimws(fasta_ids)

matched_control_map <- control_mapping[control_mapping$ensembl_transcript_id_version %in% fasta_ids, ]
matched_control_map <- matched_control_map[match(fasta_ids, matched_control_map$ensembl_transcript_id_version), ]
matched_control_map <- matched_control_map[!is.na(matched_control_map$ensembl_transcript_id_version), ]

# extract seq
control_ids <- matched_control_map$ensembl_transcript_id_version

# Get index of headers in FASTA
control_header_idx <- which(fasta_ids %in% control_ids)

# Use the same function from before:
extract_fasta_seqs <- function(fasta_lines, header_indices) {
  extracted <- character()
  headers <- grep("^>", fasta_lines, value = TRUE)
  for (i in header_indices) {
    start <- which(fasta_lines == headers[i])
    next_headers <- which(grepl("^>", fasta_lines) & seq_along(fasta_lines) > start)
    end <- if(length(next_headers) == 0) length(fasta_lines) else min(next_headers) - 1
    extracted <- c(extracted, fasta_lines[start:end])
  }
  return(extracted)
}

control_fasta_seqs <- extract_fasta_seqs(bid_fasta, control_header_idx)
writeLines(control_fasta_seqs, "/Users/beyzaerkal/Desktop/datasets_proj/bid/shControl_all_sequences.fa")


control_seqs <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/bid/shControl_all_sequences.fa")
control_mid_base <- substring(as.character(control_seqs), 61, 61)

control_seqs_U <- control_seqs[control_mid_base == "U"]
writeXStringSet(control_seqs_U, "/Users/beyzaerkal/Desktop/datasets_proj/bid/shControl_sequences_middleU.fa")

# FINAL DATA

#(bid finish)
#########################################################
########################################################

#library(seqinr) # clashes with biomart
# pseudo-seq the pool is at 69 but some fo them do not have 
# this poses issues as transcript mismatch due to verison/isoform differences 

Pus_data <- read.table("/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/GSE99487_hPUS_Pool_peaks.txt", header=TRUE, sep = "\t") # tab (tabular) separated ones
head(Pus_data)

selected_pus <- Pus_data[Pus_data$PUS_invitro %in% c("PUS7", "TRUB1"), ]
head(selected_pus)

table(Pus_data$PUS_invitro) # PUS7: 23, TRUB1: 61

# only PUS7 # log transformed no.
pus7_selected <- Pus_data[Pus_data$PUS_invitro == "PUS7", ]
pus7_selected <- Pus_data[Pus_data$PUS_invitro == "PUS7", c("pool_name", "pool_pos", "PUS_invitro", "PUS7")]


pus7_selected$transcript_id <- str_extract(pus7_selected$pool_name, "^(NM_[0-9]+|ENST[0-9]+)")

id_parts <- str_match(pus7_selected$pool_name, "^(ENST\\d+|NM_\\d+)_U(\\d+)")
pus7_selected$transcript_id <- id_parts[, 2]
pus7_selected$U_position <- as.numeric(id_parts[, 3])

#transcript_ids <- pus7_selected$pool_name
# Split 
nm_ids <- unique(pus7_selected$transcript_id[grepl("^NM_", pus7_selected$transcript_id)])
enst_ids <- unique(pus7_selected$transcript_id[grepl("^ENST", pus7_selected$transcript_id)])


# get enst
ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

enst_fasta_df <- getSequence(id = enst_ids, type = "ensembl_transcript_id",
                             seqType = "cdna", mart = ensembl_mart)


nm_seqs <- lapply(nm_ids, function(id) {
  Sys.sleep(0.3)
  tryCatch(entrez_fetch(db = "nuccore", id = id, rettype = "fasta", retmode = "text"), error = function(e) NA)
})
names(nm_seqs) <- nm_ids


nm_seqs_clean <- sapply(nm_seqs, function(fa) {
  if (is.na(fa)) return(NA)
  paste0(strsplit(fa, "\n")[[1]][-1], collapse = "")
})

all_seq_map <- c(nm_seqs_clean, setNames(enst_fasta_df$cdna, enst_fasta_df$ensembl_transcript_id))

fasta_entries <- mapply(function(id, pos) {
  if (!(id %in% names(all_seq_map))) return(NULL)  # skip 
  
  seq <- all_seq_map[[id]]
  if (is.na(seq)) return(NULL)  # skip if NA
  
  start <- max(pos - 60, 1)
  end <- min(pos + 60, nchar(seq))
  sub_seq <- substr(seq, start, end)
  
  header <- paste0(">", id, ":", start, "-", end, "(+)")
  seq_wrapped <- paste(strwrap(sub_seq, width = 80), collapse = "\n")
  
  paste(header, seq_wrapped, sep = "\n")
}, pus7_selected$transcript_id, pus7_selected$U_position, SIMPLIFY = TRUE)

####

fetch_enst_window <- function(enst_id, pos, seq_df) {
  start <- max(pos - 60, 1)
  end <- pos + 60
  full_seq <- seq_df$cdna[seq_df$ensembl_transcript_id == enst_id]
  if (length(full_seq) == 0 || is.na(full_seq)) return(NULL)
  if (end > nchar(full_seq)) end <- nchar(full_seq)
  subseq <- substr(full_seq, start, end)
  if (nchar(subseq) < 121) return(NULL) 
  subseq_rna <- chartr("T", "U", subseq)
  header <- paste0(">", enst_id, ":", start, "-", end, "(+)")
  paste(header, subseq_rna, sep = "\n")
}

enst_only <- pus7_selected[grepl("^ENST", pus7_selected$transcript_id), ]

# ENST entries
enst_fasta_list <- mapply(fetch_enst_window,
                          enst_only$transcript_id,
                          enst_only$U_position,
                          MoreArgs = list(seq_df = enst_fasta_df),
                          SIMPLIFY = TRUE)


# get nm
fetch_nm_window <- function(nm_id, pos, nm_seq_map) {
  if (!(nm_id %in% names(nm_seq_map))) return(NULL)
  full_seq <- nm_seq_map[[nm_id]]
  if (is.na(full_seq)) return(NULL)
  
  start <- max(pos - 60, 1)
  end <- min(pos + 60, nchar(full_seq))
  subseq <- substr(full_seq, start, end)
  if (nchar(subseq) < 121) return(NULL) 
  subseq_rna <- chartr("T", "U", subseq)
  header <- paste0(">", nm_id, ":", start, "-", end, "(+)")
  paste(header, subseq_rna, sep = "\n")
}

nm_only <- pus7_selected[grepl("^NM_", pus7_selected$transcript_id), ]
nm_fasta_list <- mapply(fetch_nm_window,
                        nm_only$transcript_id,
                        nm_only$U_position,
                        MoreArgs = list(nm_seq_map = nm_seqs_clean),
                        SIMPLIFY = TRUE)

all_entries <- c(enst_fasta_list, nm_fasta_list)
all_entries <- all_entries[!sapply(all_entries, is.null)]
all_entries <- as.character(all_entries)
writeLines(all_entries, "/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/pus7_pseudo_seq_combined.fasta")


# Load fasta sequences as RNAStringSet
fasta_seqs <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/pus7_pseudo_seq_combined.fasta")

seqs <- as.character(fasta_seqs)
headers <- names(fasta_seqs)

middle_bases <- substring(seqs, 61, 61)
table(middle_bases)  # Check base counts

seq_lengths <- nchar(seqs)
summary(seq_lengths) # Confirm lengths (ideally 121 nt)
idx_u <- which(middle_bases == "U")
headers_u <- headers[idx_u]
seqs_u <- seqs[idx_u]


entries_u <- unlist(mapply(function(h, s) c(h, s), 
                           paste0(">", headers_u), seqs_u, 
                           SIMPLIFY = FALSE))


writeLines(entries_u, "/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/pus7_corrected_pse.fasta")
# use this for training


#############
#############

# splitting verion is not good
# getting nm convetred - 4 smaples left only
# split and map nm and enst:
nm_ids <- pus7_selected$transcript_id[grepl("^NM_", pus7_selected$transcript_id)]

enst_ids <- pus7_selected$transcript_id[grepl("^ENST", pus7_selected$transcript_id)]

# biomart for conversion
ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

nm_to_enst <- getBM(attributes = c("refseq_mrna", "ensembl_transcript_id"), filters = "refseq_mrna",
                    values = nm_ids, mart = ensembl_mart)
# 11
mapped_ensts <- c(nm_to_enst$ensembl_transcript_id, enst_ids)
#mapped_ensts <- unique(mapped_ensts) # for duplicates

# get cdna of these
fasta_df <- getSequence(id = mapped_ensts, type = "ensembl_transcript_id",
  seqType = "cdna", mart = ensembl_mart)

nm_to_enst_map <- setNames(nm_to_enst$ensembl_transcript_id, nm_to_enst$refseq_mrna)
# replace nm to ensts
pus7_selected <- pus7_selected %>% mutate(transcript_id = ifelse(grepl("^NM_", transcript_id), nm_to_enst_map[transcript_id],
      transcript_id)) %>% filter(!is.na(transcript_id)) 



pus7_selected$strand <- "+" 

merged_df <- pus7_selected %>%
  left_join(fasta_df, by = c("transcript_id" = "ensembl_transcript_id"))

merged_df <- merged_df %>%
  mutate(
    seq_length = nchar(cdna),
    start_pos = pmax(pool_pos - 60, 1),
    end_pos = pmin(pool_pos + 60, seq_length),
    subseq_121nt = str_sub(cdna, start = start_pos, end = end_pos),
    subseq_121nt_rna = str_replace_all(subseq_121nt, "T", "U"),
    header = paste0(">", transcript_id, ":", start_pos, "-", end_pos, "(+)")
  )


merged_df <- merged_df %>%
  mutate(middle_base = str_sub(subseq_121nt_rna, 61, 61)) %>%
  filter(middle_base == "U")

fasta_file <- "/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/pus7_pse_seq.fasta"
fasta_entries <- unlist(mapply(function(h, s) c(h, strwrap(s, 80)), merged_df$header, merged_df$subseq_121nt_rna, SIMPLIFY = FALSE))
cat(fasta_entries, file = fasta_file, sep = "\n")

print(dplyr::select(merged_df, transcript_id, pool_pos, middle_base) %>% head())

# take out the ones that do not have U at 69 (61 psoitons since start at 9)

###########
#########

# pseudo-seq negative


pus7_negatives <- Pus_data[Pus_data$PUS_invitro != "PUS7", ]
# extrcat U position
id_parts_neg <- str_match(pus7_negatives$pool_name, "^(ENST\\d+|NM_\\d+)_U(\\d+)")
pus7_negatives$transcript_id <- id_parts_neg[, 2]
pus7_negatives$U_position <- as.numeric(id_parts_neg[, 3])


neg_nm_ids <- unique(pus7_negatives$transcript_id[grepl("^NM_", pus7_negatives$transcript_id)])
neg_enst_ids <- unique(pus7_negatives$transcript_id[grepl("^ENST", pus7_negatives$transcript_id)])

# enst
fasta_df_neg_enst <- getSequence(id = neg_enst_ids,
                                 type = "ensembl_transcript_id",
                                 seqType = "cdna",
                                 mart = ensembl_mart)
# nm
neg_nm_seqs <- lapply(neg_nm_ids, function(id) {
  Sys.sleep(0.3)
  tryCatch(entrez_fetch(db = "nuccore", id = id, rettype = "fasta", retmode = "text"), error = function(e) NA)
})
names(neg_nm_seqs) <- neg_nm_ids

neg_nm_seqs_clean <- sapply(neg_nm_seqs, function(fa) {
  if (is.na(fa)) return(NA)
  paste0(strsplit(fa, "\n")[[1]][-1], collapse = "")
})

# fetch the window
# enst
enst_only_neg <- pus7_negatives[grepl("^ENST", pus7_negatives$transcript_id), ]

enst_fasta_list_neg <- mapply(fetch_enst_window,
                              enst_only_neg$transcript_id,
                              enst_only_neg$U_position,
                              MoreArgs = list(seq_df = fasta_df_neg_enst),
                              SIMPLIFY = TRUE)
#nm
nm_only_neg <- pus7_negatives[grepl("^NM_", pus7_negatives$transcript_id), ]

nm_fasta_list_neg <- mapply(fetch_nm_window,
                            nm_only_neg$transcript_id,
                            nm_only_neg$U_position,
                            MoreArgs = list(nm_seq_map = neg_nm_seqs_clean),
                            SIMPLIFY = TRUE)
#COMBINE
all_entries_neg <- c(enst_fasta_list_neg, nm_fasta_list_neg)
all_entries_neg <- all_entries_neg[!sapply(all_entries_neg, is.null)]
all_entries_neg <- as.character(all_entries_neg)

#raw negatives
writeLines(all_entries_neg, "/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/pus7_neg_raw.fasta")

#  filter for middle U
fasta_seqs_neg <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/pus7_neg_raw.fasta")
seqs_neg <- as.character(fasta_seqs_neg)
headers_neg <- names(fasta_seqs_neg)

middle_bases_neg <- substring(seqs_neg, 61, 61)
idx_u_neg <- which(middle_bases_neg == "U")

entries_neg_u <- unlist(mapply(function(h, s) c(h, s),
                               paste0(">", headers_neg[idx_u_neg]),
                               seqs_neg[idx_u_neg],
                               SIMPLIFY = FALSE))

# Final filtered FASTA (middle U)
writeLines(entries_neg_u, "/Users/beyzaerkal/Desktop/datasets_proj/pseudo_files/pus7_negatives_filtered.fasta")
######### FINAL DATA
# 23!



#get_fasta_sequences(transcript_ids, positions, flank = 60) # 22
# fasta for pseudo-seq (no confidence label)
#(failed_ids) # "ENST00000426131"





#################################
# HEK293 DATASETS
###################################


#####################
# TRUB1 is the predominant pseudouridine synthase acting on mammalian mRNA via a predictable and conserved code

hek_csv <- read.csv("/Users/beyzaerkal/Desktop/datasets_proj/HEK293/supplement_table_1.csv")

hek_csv_pus7 <- hek_csv[hek_csv$group == "pus7", ]
# nrow = 968
# get only mrna
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- unique(hek_csv_pus7$geneSymbol)

gene_biotypes <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                       filters = "external_gene_name", values = gene_symbols, mart = mart)

protein_coding_genes <- gene_biotypes %>% filter(gene_biotype == "protein_coding") %>% pull(external_gene_name)

hek_csv_pus7_mrna <- hek_csv_pus7 %>% filter(geneSymbol %in% protein_coding_genes)


# nrow = 746

# for base-level sequencing - prcise

# criteria

# Ensure correct chromosome format
hek_csv_pus7_mrna <- hek_csv_pus7_mrna %>%mutate(chr = ifelse(grepl("^chr", gcoords), sub(":.*", "", gcoords),
                                                              paste0("chr", sub(":.*", "", gcoords))),
                                                 center = as.integer(sub(".*:", "", gcoords)))


# filter to chromosomes present in BSgenome
valid_chr <- intersect(hek_csv_pus7_mrna$chr, names(seqlengths(genome)))

hek_csv_pus7_mrna <- hek_csv_pus7_mrna %>% filter(chr %in% valid_chr)

chr_lengths <- seqlengths(genome)

# coords to chromosome size
hek_csv_pus7_mrna <- hek_csv_pus7_mrna %>% mutate(chr_len = chr_lengths[chr],
                                                  start = pmax(1, pmin(center - 60, chr_len)),
                                                  end = pmax(1, pmin(center + 60, chr_len)))

# extract sequence
hek_csv_pus7_mrna$seq <- mapply(function(chr, start, end, strand) {
  as.character(getSeq(genome, names = chr, start = start, end = end, strand = strand))
}, hek_csv_pus7_mrna$chr, hek_csv_pus7_mrna$start, hek_csv_pus7_mrna$end, hek_csv_pus7_mrna$strand)

# remove shorter than 121nt
hek_csv_pus7_mrna <- hek_csv_pus7_mrna %>% mutate(seq_length = nchar(seq)) %>% filter(seq_length == 121) 
#convetr t to u
hek_csv_pus7_mrna <- hek_csv_pus7_mrna %>% mutate(seq_rna = chartr("T", "U", seq))

# fasta save
# headers:
fasta_names <- paste0(hek_csv_pus7_mrna$geneSymbol, ":", hek_csv_pus7_mrna$start, "-", hek_csv_pus7_mrna$end, "(", hek_csv_pus7_mrna$strand, ")")
#rnastring object BY GENE SYMBOL
rna_seqs <- RNAStringSet(hek_csv_pus7_mrna$seq_rna)
names(rna_seqs) <- fasta_names

# convert genesymbol to ensembl ID
gene_symbols <- unique(hek_csv_pus7_mrna$geneSymbol)

# Query biomart for mapping
mapping <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),filters = "external_gene_name",
                 values = gene_symbols, mart = mart)

head(mapping) # no version


# keep only one ensembl ID per gene
mapping_unique <- mapping %>% distinct(external_gene_name, .keep_all = TRUE)

hek_csv_pus7_mrna_ensg <- left_join(hek_csv_pus7_mrna, mapping_unique, 
                                    by = c("geneSymbol" = "external_gene_name"))
# 746
sum(duplicated(hek_csv_pus7_mrna$geneSymbol))

# fasta for ENSEMBL ID -ensg

fasta_names_ensg <- paste0(hek_csv_pus7_mrna_ensg$ensembl_gene_id, ":", hek_csv_pus7_mrna_ensg$start, "-", 
                      hek_csv_pus7_mrna_ensg$end, "(", hek_csv_pus7_mrna_ensg$strand, ")")

# RNAStringSet with sequences
rna_seqs <- RNAStringSet(hek_csv_pus7_mrna_ensg$seq_rna)

names(rna_seqs) <- fasta_names_ensg

writeXStringSet(rna_seqs, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/pus7_hek_csv.fasta")
# the sequence is in 2 lines

# convetr ensg to enst
transcript_mapping <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                            filters = "ensembl_gene_id",
                            values = hek_csv_pus7_mrna_ensg$ensembl_gene_id,
                            mart = mart)

transcript_mapping_unique <- transcript_mapping %>% distinct(ensembl_gene_id, .keep_all = TRUE)

hek_csv_pus7_mrna_enst <- left_join(hek_csv_pus7_mrna_ensg, transcript_mapping_unique, 
                                    by = "ensembl_gene_id")
fasta_names_enst <- paste0(hek_csv_pus7_mrna_enst$ensembl_transcript_id, ":", 
                           hek_csv_pus7_mrna_enst$start, "-", 
                           hek_csv_pus7_mrna_enst$end, "(", 
                           hek_csv_pus7_mrna_enst$strand, ")")

rna_seqs <- RNAStringSet(hek_csv_pus7_mrna_enst$seq_rna)
names(rna_seqs) <- fasta_names_enst

writeXStringSet(rna_seqs, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/pus7_hek_csv_ENST.fasta")
# 746

rna_seqs <- RNAStringSet(hek_csv_pus7_mrna_enst$seq_rna)
seqs <- as.character(rna_seqs)

# lengths
seq_lengths <- nchar(seqs)
table(seq_lengths)  # 737 is 121nt and 9 are only 1 nt long(failed ones)
# hek_csv_pus7_mrna <- hek_csv_pus7_mrna %>%mutate(seq_length = nchar(seq)) %>%filter(seq_length == 121)
# after applying above its fixed

# middle base position 61, 1-based
mid_bases <- substring(seqs, 61, 61) # 206 sequences have "U" at position 61
table(mid_bases)  

seqs <- as.character(rna_seqs)
# position 61 (1-based)
mid_bases <- substring(seqs, 61, 61)

idx_U <- which(mid_bases == "U")
rna_seqs_U <- rna_seqs[idx_U]
print(names(rna_seqs_U)[1:5])

cat("Number of sequences with middle base U:", length(rna_seqs_U), "\n")

writeXStringSet(rna_seqs_U, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/pus7_hek_csv_ENST_U_only.fasta")

# get teh U only
seqs <- as.character(hek_csv_pus7_mrna_enst$seq_rna)
#(position 61, 1-based)
mid_bases <- substring(seqs, 61, 61)
rows_with_U <- mid_bases == "U"
hek_csv_pus7_mrna_enst_U <- hek_csv_pus7_mrna_enst[rows_with_U, ]

# fasta
fasta_names_U <- paste0(hek_csv_pus7_mrna_enst_U$ensembl_transcript_id, ":",
                        hek_csv_pus7_mrna_enst_U$start, "-", 
                        hek_csv_pus7_mrna_enst_U$end, "(",
                        hek_csv_pus7_mrna_enst_U$strand, ")")

#  filtered sequences
rna_seqs_U <- RNAStringSet(hek_csv_pus7_mrna_enst_U$seq_rna)
# headers
names(rna_seqs_U) <- fasta_names_U

writeXStringSet(rna_seqs_U, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/pus7_hek_csv_ENST_U_only.fasta")
# use these for U mid only - final
 
# FINAL DATA

##################
# negative of this


hek_csv <- read.csv("/Users/beyzaerkal/Desktop/datasets_proj/HEK293/supplement_table_1.csv")

#  (negatives)
hek_csv_other <- hek_csv[hek_csv$group != "pus7", ]

# Setup biomart
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get unique gene symbols
gene_symbols <- unique(hek_csv_other$geneSymbol)

# Get gene biotypes and filter protein-coding only
gene_biotypes <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                       filters = "external_gene_name", values = gene_symbols, mart = mart)

protein_coding_genes <- gene_biotypes %>% filter(gene_biotype == "protein_coding") %>% pull(external_gene_name)

hek_csv_other_mrna <- hek_csv_other %>% filter(geneSymbol %in% protein_coding_genes)

# Correct chromosome format and center coordinate
hek_csv_other_mrna <- hek_csv_other_mrna %>%
  mutate(chr = ifelse(grepl("^chr", gcoords), sub(":.*", "", gcoords),
                      paste0("chr", sub(":.*", "", gcoords))),
         center = as.integer(sub(".*:", "", gcoords)))

# Validate chromosomes present in genome
genome <- BSgenome.Hsapiens.UCSC.hg38
valid_chr <- intersect(hek_csv_other_mrna$chr, names(seqlengths(genome)))
hek_csv_other_mrna <- hek_csv_other_mrna %>% filter(chr %in% valid_chr)

# Get chromosome lengths from genome
chr_lengths <- seqlengths(genome)

# Calculate start and end coordinates (60 nt up/downstream, clipped to chr length)
hek_csv_other_mrna <- hek_csv_other_mrna %>%
  mutate(chr_len = chr_lengths[chr],
         start = pmax(1, pmin(center - 60, chr_len)),
         end = pmax(1, pmin(center + 60, chr_len)))

# Convert chr and strand to character
hek_csv_other_mrna$chr <- as.character(hek_csv_other_mrna$chr)
hek_csv_other_mrna$strand <- as.character(hek_csv_other_mrna$strand)

# Create GRanges for sequence extraction
gr <- GRanges(
  seqnames = hek_csv_other_mrna$chr,
  ranges = IRanges(start = hek_csv_other_mrna$start, end = hek_csv_other_mrna$end),
  strand = hek_csv_other_mrna$strand
)


seqinfo(gr) <- seqinfo(genome)
# extract genome
hek_csv_other_mrna$seq <- as.character(getSeq(genome, gr))

# length 121nt
hek_csv_other_mrna <- hek_csv_other_mrna %>%
  mutate(seq_length = nchar(seq)) %>% filter(seq_length == 121)

# T to U 
hek_csv_other_mrna <- hek_csv_other_mrna %>% mutate(seq_rna = chartr("T", "U", seq))

# fasta headers by gene symbol
fasta_names <- paste0(hek_csv_other_mrna$geneSymbol, ":", hek_csv_other_mrna$start, "-", hek_csv_other_mrna$end, "(", hek_csv_other_mrna$strand, ")")

rna_seqs <- RNAStringSet(hek_csv_other_mrna$seq_rna)
names(rna_seqs) <- fasta_names

#  initial fasta by gene symbol
writeXStringSet(rna_seqs, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/other_pus_hek_csv.fasta")

# map gene symbol to ensembl gene ID
mapping <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                 filters = "external_gene_name",
                 values = unique(hek_csv_other_mrna$geneSymbol),
                 mart = mart)

mapping_unique <- mapping %>% distinct(external_gene_name, .keep_all = TRUE)

hek_csv_other_mrna_ensg <- left_join(hek_csv_other_mrna, mapping_unique,
                                     by = c("geneSymbol" = "external_gene_name"))

# Fasta headers with ensembl gene IDs
fasta_names_ensg <- paste0(hek_csv_other_mrna_ensg$ensembl_gene_id, ":", hek_csv_other_mrna_ensg$start, "-", hek_csv_other_mrna_ensg$end, "(", hek_csv_other_mrna_ensg$strand, ")")

rna_seqs <- RNAStringSet(hek_csv_other_mrna_ensg$seq_rna)
names(rna_seqs) <- fasta_names_ensg

writeXStringSet(rna_seqs, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/other_pus_hek_csv_ensg.fasta")

# Map gene IDs to transcript IDs
transcript_mapping <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                            filters = "ensembl_gene_id",
                            values = hek_csv_other_mrna_ensg$ensembl_gene_id,
                            mart = mart)

transcript_mapping_unique <- transcript_mapping %>% distinct(ensembl_gene_id, .keep_all = TRUE)

hek_csv_other_mrna_enst <- left_join(hek_csv_other_mrna_ensg, transcript_mapping_unique,
                                     by = "ensembl_gene_id")

# Fasta headers with transcript IDs
fasta_names_enst <- paste0(hek_csv_other_mrna_enst$ensembl_transcript_id, ":", hek_csv_other_mrna_enst$start, "-", hek_csv_other_mrna_enst$end, "(", hek_csv_other_mrna_enst$strand, ")")

rna_seqs <- RNAStringSet(hek_csv_other_mrna_enst$seq_rna)
names(rna_seqs) <- fasta_names_enst

writeXStringSet(rna_seqs, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/other_pus_hek_csv_enst.fasta")

# Filter sequences with 'U' at the middle position (61)
seqs <- as.character(rna_seqs)
mid_bases <- substring(seqs, 61, 61)
rows_with_U <- mid_bases == "U"
hek_csv_other_mrna_enst_U <- hek_csv_other_mrna_enst[rows_with_U, ]

fasta_names_U <- paste0(hek_csv_other_mrna_enst_U$ensembl_transcript_id, ":", hek_csv_other_mrna_enst_U$start, "-", hek_csv_other_mrna_enst_U$end, "(", hek_csv_other_mrna_enst_U$strand, ")")

rna_seqs_U <- RNAStringSet(hek_csv_other_mrna_enst_U$seq_rna)
names(rna_seqs_U) <- fasta_names_U

writeXStringSet(rna_seqs_U, filepath = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/other_pus_hek_csv_enst_U_only.fasta")
# FINAL DATA

### cut the sample size for balancing
neg_fasta <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/HEK293/other_pus_hek_csv_enst_U_only.fasta")
set.seed(42) # reproducibility
subset_neg <- neg_fasta[sample(length(neg_fasta), 412)]
writeXStringSet(subset_neg, "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/other_pus_hek_csv_enst_downsampled.fasta")



#######################
# Quantitative profiling of pseudouridylation landscape in the human transcriptome

hek_qp <- read_excel("/Users/beyzaerkal/Desktop/datasets_proj/HEK293/Q41589_2023_1304_MOESM3_ESM.xlsx", sheet = 3, skip = 2, col_names = TRUE)
# only get mRNA - also contain ncRNA
# only M: mRNA and curated: N -> NM data only : since all experiemntally done all of these are high confidence confirmed sites

curated_mrna <- hek_qp[grepl("^NM_", hek_qp$chr_name), ]
# pus7 only 
qp_pus7 <- curated_mrna[curated_mrna$enzyme_dependency == "PUS7", ]

# label by confidence - Deletion ratio low or zero, p-value is < 0.05

# get diff
qp_pus7$diff_rep1 <- qp_pus7$wt_treated_deletion_ratio - qp_pus7$ko_rep1_deletion_ratio
qp_pus7$diff_rep2 <- qp_pus7$wt_treated_deletion_ratio - qp_pus7$ko_rep2_deletion_ratio

qp_pus7$decrease_rate_rep1 <- ifelse(qp_pus7$wt_treated_deletion_ratio == 0, 0,
                                     qp_pus7$diff_rep1 / qp_pus7$wt_treated_deletion_ratio)
qp_pus7$decrease_rate_rep2 <- ifelse(qp_pus7$wt_treated_deletion_ratio == 0, 0,
                                     qp_pus7$diff_rep2 / qp_pus7$wt_treated_deletion_ratio)


qp_pus7 <- qp_pus7 %>% mutate(confidence = case_when(
    ko_rep1_total_counts >= 10 &
      ko_rep2_total_counts >= 10 &
      diff_rep1 > 0.05 &
      decrease_rate_rep1 > 0.05 &
      diff_rep2 > 0.05 &
      decrease_rate_rep2 > 0.05 ~ "high",
    
    ko_rep1_total_counts >= 10 &
      ko_rep2_total_counts >= 10 &
      diff_rep1 < 0.01 &
      decrease_rate_rep1 < 0.01 &
      diff_rep2 < 0.01 &
      decrease_rate_rep2 < 0.01 ~ "background",
    
    TRUE ~ NA_character_
  ))

head(qp_pus7[, c("chr_site", "confidence")]) #105


# getting fasta
# rentrez use loop as it only takes only one sample per run - fetch_nm_sequence

hek_qp_nm_ids <- unique(qp_pus7$chr_name)

fetch_sequence <- function(nm_id) {
  tryCatch({
    fasta <- entrez_fetch(db = "nuccore", id = nm_id, rettype = "fasta", retmode = "text")
    seq_lines <- unlist(strsplit(fasta, "\n"))
    seq <- paste(seq_lines[-1], collapse = "")
    seq <- chartr("T", "U", seq)
    return(seq)
  }, error = function(e) {
    message(paste("Error:", nm_id))
    return(NA)
  })
}

# for multiple ones: "5913-5914" "3024-3026" "1707-1710" "663-666"   "2448-2450"
# take mid-range
extract_mid_pos <- function(site) {
  if (grepl("-", site)) {
    parts <- strsplit(site, "-")[[1]]
    mean(as.numeric(parts))
  } else {
    as.numeric(site)
  }
}

positions <- sapply(qp_pus7$site, extract_mid_pos)


get_fasta_sequences <- function(transcript_ids, positions, flank = 60, output_file = "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/hek_qb.fasta") {
  fasta_list <- c()
  
  for (i in seq_along(transcript_ids)) {
    id <- transcript_ids[i]
    mod_pos <- positions[i]
    
    if (is.na(id) || is.na(mod_pos)) next # for invalids skip
    
    seq <- fetch_sequence(id)
    
    if (is.na(seq)) next  # skip when fetch failed
    
    seq_len <- nchar(seq)
    start_pos <- max(mod_pos - flank, 1)
    end_pos <- min(mod_pos + flank, seq_len)  # all window
    
  
    mod_seq <- substr(seq, start_pos, end_pos)
    if (nchar(mod_seq) != 121) next
    header <- paste0(">", id, ":", start_pos, "-", end_pos, "(+)")
    
    fasta_list <- c(fasta_list, header, mod_seq)
  }
  
  writeLines(fasta_list, con = output_file)
}

# run
get_fasta_sequences(qp_pus7$chr_name, positions, flank = 60)


fasta_seqs <- readRNAStringSet("/Users/beyzaerkal/Desktop/datasets_proj/HEK293/hek_qb.fasta")
seqs <- as.character(fasta_seqs)
# FINAL DATA

seq_lengths <- nchar(seqs)
table(seq_lengths)
# first check
# 104 sequences that are 121 nt long
# 1 sequence that is 110 nt long
# second check
# all 121nt - 104 in total
mid_bases <- substring(seqs, 61, 61)
table(mid_bases)
# all 105 sequences have a "U"
#################################

##################

# Filter: Not PUS7
hek_qp <- read_excel("/Users/beyzaerkal/Desktop/datasets_proj/HEK293/Q41589_2023_1304_MOESM3_ESM.xlsx", sheet = 3, skip = 2, col_names = TRUE)
# only get mRNA - also contain ncRNA
# only M: mRNA and curated: N -> NM data only : since all experiemntally done all of these are high confidence confirmed sites

curated_mrna <- hek_qp[grepl("^NM_", hek_qp$chr_name), ]

qp_negatives <- curated_mrna[curated_mrna$enzyme_dependency != "PUS7", ]

qp_negatives$chr_name <- as.character(qp_negatives$chr_name)
qp_negatives$site <- as.character(qp_negatives$site)

# Extract midpoint of 'site' column 
extract_mid_pos <- function(site) {
  if (grepl("-", site)) {
    parts <- strsplit(site, "-")[[1]]
    return(mean(as.numeric(parts)))
  } else {
    return(as.numeric(site))
  }
}

positions <- sapply(qp_negatives$site, extract_mid_pos)

# Sequence fetcher
fetch_sequence <- function(nm_id) {
  tryCatch({
    fasta <- entrez_fetch(db = "nuccore", id = nm_id, rettype = "fasta", retmode = "text")
    seq_lines <- unlist(strsplit(fasta, "\n"))
    seq <- paste(seq_lines[-1], collapse = "")
    seq <- chartr("T", "U", seq)
    return(seq)
  }, error = function(e) {
    message(paste("Error fetching", nm_id, ":", e$message))
    return(NA)
  })
}

# FASTA writer
get_fasta_sequences <- function(transcript_ids, positions, flank = 60, output_file) {
  fasta_list <- character()
  
  n <- length(transcript_ids)
  for (i in seq_len(n)) {
    id <- transcript_ids[i]
    mod_pos <- positions[i]
    
    cat("Processing", i, "of", n, "-", id, "\r")  #progress
    
    if (is.na(id) || is.na(mod_pos)) next
    
    seq <- fetch_sequence(id)
    if (is.na(seq)) next
    
    seq_len <- nchar(seq)
    start_pos <- max(mod_pos - flank, 1)
    end_pos <- min(mod_pos + flank, seq_len)
    
    mod_seq <- substr(seq, start_pos, end_pos)
    
    if (nchar(mod_seq) != 121) next
    if (substr(mod_seq, 61, 61) != "U") next  # central base must be U
    
    header <- paste0(">", id, ":", start_pos, "-", end_pos, "(+)")
    fasta_list <- c(fasta_list, header, mod_seq)
  }
  
  if (length(fasta_list) == 0) {
    message("no valid sequences found")
    return(NULL)
  }
  
  writeLines(fasta_list, con = output_file)
  cat("\nFinished writing to:", output_file, "\n")
}



output_file <- "/Users/beyzaerkal/Desktop/datasets_proj/HEK293/hek_qb_negatives.fasta"
get_fasta_sequences(qp_negatives$chr_name, positions, output_file = output_file)
# FINAL DATA
















