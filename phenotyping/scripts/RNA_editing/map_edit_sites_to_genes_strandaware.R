#!/usr/bin/env Rscript
#
# map_edit_sites_to_genes_strandaware.R
#
# 1) if a site overlaps exactly one gene body → map to that gene  
# 2) if it overlaps >1 gene → pick the gene whose TSS is closest to the site  
#
# Usage:
#   Rscript map_edit_sites_to_genes_strandaware.R <sites.bed> <genes.gtf> <out.tsv>

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(glue)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: map_edit_sites_to_genes_strandaware.R <sites.bed> <genes.gtf> <out.tsv>")
}
sites_bed <- args[1]
gtf_file  <- args[2]
out_tsv   <- args[3]

message("1) Reading sites BED…")
sites_df <- read_table(
  sites_bed,
  col_names = c("seqname","start0","end0","dot","score","strand"),
  col_types = "ciicic"
) %>%
  mutate(
    start1  = start0 + 1,
    site_id = str_glue("{seqname}_{end0}")
  )

sites_gr <- GRanges(
  seqnames = sites_df$seqname,
  ranges   = IRanges(start = sites_df$start1, width = 1),
  strand   = sites_df$strand,
  site_id  = sites_df$site_id
)

# keep only autosomes 1–22
keep_autos <- as.character(1:22)
sites_gr   <- keepSeqlevels(sites_gr, keep_autos, pruning.mode="coarse")

message("2) Importing GTF and extracting gene bodies…")
genes   <- import(gtf_file, format="gtf", features="gene")
gene_body_gr <- keepSeqlevels(genes, keep_autos, pruning.mode="coarse")
mcols(gene_body_gr)$gene_id <- mcols(genes)$gene_id

# build TSS GRanges
tss_pos <- ifelse(strand(gene_body_gr) == "+",
                  start(gene_body_gr),
                  end(gene_body_gr))
tss_gr <- GRanges(
  seqnames = seqnames(gene_body_gr),
  ranges   = IRanges(start = tss_pos, width = 1),
  strand   = strand(gene_body_gr),
  gene_id  = mcols(gene_body_gr)$gene_id
)

message("3) Finding overlaps…")
hits <- findOverlaps(sites_gr, gene_body_gr, ignore.strand = FALSE)
if (length(hits) == 0) {
  stop("No overlaps found: check BED vs. GTF chromosome names")
}

message("4) Building small hits table…")
hits_df <- tibble(
  site_idx = queryHits(hits),
  gene_idx = subjectHits(hits)
) %>%
  mutate(
    site_id = mcols(sites_gr)$site_id[site_idx],
    gene_id = mcols(gene_body_gr)$gene_id[gene_idx],
    tss_pos = start(tss_gr)[gene_idx],
    site_pos = start(sites_gr)[site_idx]
  ) %>%
  select(site_id, gene_id, tss_pos, site_pos)

# split single- vs multi-overlaps
counts <- hits_df %>% count(site_id, name="n_hits")
singles    <- counts %>% filter(n_hits == 1) %>% pull(site_id)
multiples  <- counts %>% filter(n_hits >  1) %>% pull(site_id)

message("  • single‑hit sites: ", length(singles))
message("  • multi‑hit sites:  ", length(multiples))

# 5a) single‑gene sites
mapped_single <- hits_df %>%
  filter(site_id %in% singles) %>%
  select(site_id, gene_id)

# 5b) multi‑gene sites → nearest TSS
mapped_multi <- hits_df %>%
  filter(site_id %in% multiples) %>%
  group_by(site_id) %>%
  slice_min(abs(tss_pos - site_pos), with_ties = FALSE) %>%
  ungroup() %>%
  select(site_id, gene_id)

# 6) combine & write
mapped <- bind_rows(mapped_single, mapped_multi)
message("7) Writing ", nrow(mapped), " total mappings to ", out_tsv)
write_tsv(mapped, out_tsv, col_names = FALSE)
