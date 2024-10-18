library("MiscMetabar")
library("targets")

here::i_am("_targets.R")
if (!dir.exists("output")) {
  dir.create("output")
}

list(
  tar_target(
    name = file_sam_data_csv,
    command = "data/samples_data/sam_data_english.csv",
    format = "file"
  ),
  tar_target(
    name = file_maarjam_dada2,
    command = "data/maarjam/maarjam_dada2.fasta.gz",
    format = "file"
  ),
  tar_target(
    name = file_pr2_dada2,
    command = "data/PR2/pr2_version_5.0.0_SSU_dada2.fasta",
    format = "file"
  ),
  tar_target(
    name = file_maarjam_dada2_species,
    command = "data/maarjam/maarjam_dada2_species.fasta.gz",
    format = "file"
  ),
  tar_target(
    name = file_cutadapt_sh,
    command = "code/cutadapt.sh",
    format = "file"
  ),
  tar_target(
    name = folder_wo_primers,
    command = "data/data_raw/wo_primers",
    format = "file"
  ),
  tar_target(
    name = file_maarjam_blast,
    command = "data/maarjam/maarjam_blast.fasta",
    format = "file"
  ),
  # The folder wo_primers is obtain using the script cutadapt.sh
  tar_target(cutadapt, system(paste0(
    "bash ", file_cutadapt_sh, " --paired_end true"
  ))),
  tar_target(data_raw, {
    cutadapt
    list_fastq_files(path = here::here(folder_wo_primers),
                     paired_end = TRUE)
  }),
  tar_target(data_fnfs, data_raw$fnfs),
  tar_target(data_fnrs, data_raw$fnrs),
  ### Pre-filtered data with low stringency
  tar_target(
    filtered,
    filter_trim(
      output_fw = paste(getwd(), "/output/filterAndTrim_fwd", sep = ""),
      output_rev = paste(getwd(), "/output/filterAndTrim_rev", sep = ""),
      fw = data_fnfs,
      rev = data_fnrs,
      multithread = 4,
      compress = TRUE
    )
  ),

  ### Dereplicate fastq files
  tar_target(derep_fs, derepFastq(filtered[[1]]), format = "qs"),
  tar_target(derep_rs, derepFastq(filtered[[2]]), format = "qs"),
  ### Learns the error rates
  tar_target(err_fs, learnErrors(derep_fs, multithread = 4), format = "qs"),
  tar_target(err_rs, learnErrors(derep_rs, multithread = 4), format = "qs"),
  ### Make amplicon sequence variants
  tar_target(ddF, dada(derep_fs, err_fs, multithread = 4), format = "qs"),
  tar_target(ddR, dada(derep_rs, err_rs, multithread = 4), format = "qs"),
  ### Merge paired sequences
  tar_target(
    merged_seq,
    mergePairs(
      dadaF = ddF,
      dadaR = ddR,
      derepF = derep_fs,
      derepR = derep_rs,
      minOverlap = 8,
      maxMismatch = 1
    ),
    format = "qs"
  ),
  ### Build a a table of ASV x Samples
  tar_target(seq_tab_Pairs, makeSequenceTable(merged_seq)),

  ### Remove sequences < 200 bp length.
  tar_target(seqtab_wo_chimera, chimera_removal_vs(seq_tab_Pairs)),
  tar_target(seqtab,
             seqtab_wo_chimera[, nchar(colnames(seqtab_wo_chimera)) >= 200]),
  tar_target(
    tax_tab_maarjam,
    assignTaxonomy(
      seqtab,
      refFasta = paste0(here::here(), "/", file_maarjam_dada2),
      taxLevels = c(
        "Kingdom",
        "Phyla",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "VT"
      ),
      multithread = 3,
      minBoot = 50
    )
  ),
  tar_target(tax_tab_maarjam_species,
             setNames(
               data.frame(cbind(
                 tax_tab_maarjam,
                 assignSpecies(
                   seqtab,
                   refFasta = paste0(here::here(), "/", file_maarjam_dada2_species),
                   allowMultiple = TRUE
                 )
               )),
               c(
                 "Kingdom",
                 "Phyla",
                 "Class",
                 "Order",
                 "Family",
                 "Genus",
                 "Species",
                 "VT",
                 "Genus_100",
                 "Sp_100"
               )
             )),
  tar_target(id_sam,
             read.delim(
               paste0(here::here(), "/", file_sam_data_csv)
             )$Ref_for_bioinformatic),
  tar_target(id_asv_table,
             unlist(regmatches(
               basename(data_fnfs), regexec(pattern = "(.*?)_S", basename(data_fnfs))
             ))
             [seq(2, length(data_fnfs) * 2, by = 2)]),
  tar_target(
    sam_tab,
    sample_data_with_new_names(
      paste0(here::here(), "/", file_sam_data_csv),
      names_of_samples = id_sam,
      samples_order = na.omit(match(id_asv_table, id_sam)),
      dec = ","
    )
  ),
  tar_target(asv_tab,
             otu_table(
               rename_samples(
                 otu_table(seqtab, taxa_are_rows = FALSE),
                 names_of_samples = id_asv_table
               ),
               taxa_are_rows = FALSE
             )),
  ### Create the phyloseq object 'data_phyloseq' with
  ###   (i) table of asv,
  ###   ii) taxonomic table,
  ###   (iii) sample data and
  ###   (iv) references sequences
  tar_target(
    data_phyloseq,
   add_dna_to_phyloseq(phyloseq(asv_tab,
                     sample_data(sam_tab),
                     tax_table(
                       as.matrix(tax_tab_maarjam_species,
                                 dimnames = rownames(tax_tab_maarjam_species))
                     )))
  ),
  tar_target(
    data_phyloseq_PR2,
   add_new_taxonomy_pq(
      data_phyloseq,
      file_pr2_dada2,
      suffix = "PR2",
      taxLevels = c(
        "Kingdom",
        "Supergroup",
        "Division",
        "Subdivision",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      )
    )
  ),
  tar_target(
    d_vs,
    asv2otu(data_phyloseq_PR2, method = "vsearch", tax_adjust = 0)
  ),
  tar_target(
    d_blast,
    filter_asv_blast(
      data_phyloseq_PR2,
      fasta_for_db = paste0(here::here(), "/", file_maarjam_blast),
      nproc = 4
    )
  ),
  tar_target(
    d_vs_blast,
    filter_asv_blast(
      d_vs,
      fasta_for_db = paste0(here::here(), "/", file_maarjam_blast),
      nproc = 4
    )
  ),
  tar_target(d_vs_blast_rarefy,
             rarefy_even_depth(d_vs_blast, rngseed = 22)),
  tar_target(track_sequences_samples_clusters,
             track_wkflow(
               list(
                 "Paired sequences" = seq_tab_Pairs,
                 "Paired sequences without chimera" = seqtab_wo_chimera,
                 "Paired sequences without chimera and longer than 200bp" = seqtab,
                 "ASV denoising" = data_phyloseq_PR2,
                 "OTU after vsearch reclustering at 97%" = d_vs,
                 "OTU after blast filter without reclustering" = d_blast,
                 "OTU after blast filter with reclustering" = d_vs_blast
               )
             )),
  tar_target(track_by_samples,
             track_wkflow_samples(
               list(
                 "ASV denoising" = data_phyloseq_PR2,
                 "OTU after vsearch reclustering at 97%" = d_vs,
                 "OTU after blast filter without reclustering" = d_blast,
                 "OTU after blast filter with reclustering" = d_vs_blast
               )
             ))
)
