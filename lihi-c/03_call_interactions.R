#!/usr/bin/env Rscript
# One replicate (or a separately supplied pooled disomy map), all three resolutions.
suppressPackageStartupMessages({library(HiCDCPlus);library(data.table);library(ggplot2)})
# Read the contact map, output directory, donor label and optional resources.
a <- commandArgs(TRUE)
if(length(a)<3L) stop("Usage: Rscript 03_call_interactions.R input.hic outdir sample [cores=4] [Dmax_bp=package_default]")
hic <- normalizePath(a[1],mustWork=TRUE);out <- a[2];sample <- a[3]
if(!sample %in% c("16216","17759","16216_17759")) stop("Disomy sample labels only")
cores <- if(length(a)>=4L) as.integer(a[4]) else 4L
if(is.na(cores)||cores<1L) stop("cores must be positive")
# Keep historical package defaults unless an explicit distance limit is supplied.
dmax <- if(length(a)>=5L) as.numeric(a[5]) else NULL
if(!is.null(dmax)&&(!is.finite(dmax)||dmax<=0)) stop("Invalid Dmax")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Save the package version and its distance defaults alongside the results.
capture.output({print(packageDescription("HiCDCPlus"));print(sessionInfo());print(formals(HiCDCPlus_parallel));print(formals(generate_bintolen_gi_list))},file=file.path(out,paste0(sample,"_software_and_defaults.txt")))
# Analyze the same map independently at each of the three bin sizes.
for(bin in c(5000L,10000L,25000L)) {
 prefix <- file.path(out,paste0(sample,"_hg38_",bin,"_AAGCTT"))
 result <- paste0(prefix,"_combinedresult.txt.gz")
 if(file.exists(result)) stop("Output exists: ",result)
 # Locate HindIII sites in hg38 and calculate GC content and effective length.
 # No mappability track is supplied.
 construct_features(output_path=prefix,gen="Hsapiens",gen_ver="hg38",sig="AAGCTT",bin_type="Bins-uniform",binsize=bin)
 gi_args <- list(bintolen_path=paste0(prefix,"_bintolen.txt.gz"),gen_ver="hg38")
 if(!is.null(dmax)) gi_args$Dthreshold <- dmax
 # Create candidate cis bin pairs and attach the observed counts from the map.
 gi <- do.call(generate_bintolen_gi_list,gi_args)
 gi <- add_hic_counts(gi,hic_path=hic);gi <- expand_1D_features(gi)
 # Combine the two anchors' features and fit the negative-binomial contact model.
 # A fixed random seed makes the model's sampling step reproducible.
 set.seed(1010)
 model_args <- list(gi_list=gi,ncore=cores)
 if(!is.null(dmax)) model_args$Dmax <- dmax
 gi <- do.call(HiCDCPlus_parallel,model_args)
 # Write chromosome-wise to avoid gi_list_write's large-table append failure.
 temp <- paste0(result,".partial");con <- gzfile(temp,"wt");header <- TRUE
 summary <- list();distances <- list()
 tryCatch({
  for(chr in names(gi)) {
   x <- gi[[chr]]
   i <- as.data.frame(InteractionSet::anchors(x,type="first"))[,c("seqnames","start","end")]
   j <- as.data.frame(InteractionSet::anchors(x,type="second"))[,c("seqnames","start","end")]
   names(i)<-c("chrI","startI","endI");names(j)<-c("chrJ","startJ","endJ")
   z <- as.data.frame(S4Vectors::mcols(x))[,c("D","counts","pvalue","qvalue","mu","sdev")]
   write.table(cbind(i,j,z),con,sep="\t",quote=FALSE,row.names=FALSE,col.names=header);header<-FALSE
   # HiC-DC+ supplies within-chromosome BH-adjusted q-values.
   sig <- !is.na(z$qvalue)&z$qvalue<=0.05
   summary[[chr]] <- data.table(sample=sample,bin_size=bin,chromosome=chr,tested_with_qvalue=sum(!is.na(z$qvalue)),significant_q_le_0.05=sum(sig))
   distances[[chr]] <- z$D[sig & !is.na(z$D) & z$D>0]
  }
 },finally=close(con))
 if(!file.rename(temp,result)) stop("Failed to finalize ",result)
 # Save per-chromosome call counts, contact distances, and a distance histogram.
 fwrite(rbindlist(summary),paste0(prefix,"_summary.tsv"),sep="\t")
 d <- data.table(distance_bp=unlist(distances,use.names=FALSE))
 fwrite(d,paste0(prefix,"_significant_distances.tsv"),sep="\t")
 if(nrow(d)) ggsave(paste0(prefix,"_significant_distance_histogram.pdf"),ggplot(d,aes(distance_bp/1000))+geom_histogram(bins=60)+scale_x_log10()+labs(x="Contact distance (kb)",y="Significant bin pairs")+theme_classic(),width=6,height=4)
 rm(gi);gc()
}
