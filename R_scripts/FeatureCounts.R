library(Rsubread)

# List your BAM files
bam_files <- list.files(path = "/fh/fast/mian_o/grp/Mianlab/NEBCProject/StarWDL/", pattern = "\\.bam$", full.names = TRUE)

# Run featureCounts
fc_out <- featureCounts(files = bam_files,
                        annot.ext = "/fh/fast/mian_o/grp/Mianlab/Bioinformatics/INDEXS/GCF_000001635.27_GRCm39_genomic.gtf",  # use your GTF file
                        isGTFAnnotationFile = TRUE,
                        isPairedEnd = TRUE,  # set FALSE if single-end
                        nthreads = 8)

# This gives you fc_out$counts and fc_out$annotation
write.csv(fc_out$counts, "WDLSTAR_counts_matrix.csv")


