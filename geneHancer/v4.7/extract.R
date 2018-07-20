tbl.elite <- read.table("enhancer_elite_ids.txt", sep="\t", as.is=TRUE, header=TRUE)            #  243281      6
tbl.geneScoresAll <- read.table("enhancer_gene_scores.txt", sep="\t", as.is=TRUE, header=TRUE)  #  934287      9
tbl.tfs <- read.table("enhancer_tfs.txt", sep="\t", as.is=TRUE, header=TRUE)
tbl.tissues <- read.table("enhancer_tissues.txt", sep="\t", as.is=TRUE, header=TRUE)

clusters.elite <- sort(unique(tbl.elite$cluster_id))
clusters.geneScores <- sort(unique(tbl.geneScoresAll$cluster_id))
clusters.tfs <- sort(unique(tbl.tfs$enhancer_cluster_id))
clusters.tissues <- sort(unique(tbl.tissues$cluster_id))

printf("clusters.elite: %d", length(clusters.elite))
printf("clusters.geneScores: %d", length(clusters.geneScores))
printf("clusters.tfs: %d", length(clusters.tfs))
printf("clusters.tissues: %d", length(clusters.tissues))

all.clusterIDs <- sort(unique(c(clusters.elite, clusters.geneScores, clusters.tfs, clusters.tissues)))
printf("all clusters: %d", length(all.clusterIDs))

dim(tbl.elite); dim(tbl.geneScoresAll)
   # [1]  250733      7
   # [1] 1086552      9


targetGene <- "TREM2"
tbl.scores <- subset(tbl.geneScoresAll, symbol==targetGene)
clusters <- tbl.scores$cluster_id
tbl.enhancers <- subset(tbl.elite, cluster_id %in% clusters)[, c("chr", "enhancer_start", "enhancer_end")]
colnames(tbl.enhancers) <- c("chrom", "start", "end")
tbl.enhancers$chrom <- paste("chr", tbl.enhancers$chrom, sep="")
rownames(tbl.enhancers) <- NULL
filename <- sprintf("enhancers.v47.%s.RData", targetGene)
save(tbl.enhancers, file=filename)


# explore these four tables using the pre-term birth snp, rs7594852, provided by alison
genome <- "hg38"
targetGene <- "CKAP2L"
#sgm <- trenaSGM(genome, targetGene, quiet=FALSE)
chromosome <- "chr2"
tss <- 112764686
snp <- "rs7594852"
snp.bp <- 112764177
snp.loc <- sprintf("%s:%d-%d", chromosome, snp.bp, snp.bp)

snp.clusterID <- subset(tbl.elite, chr==2 & enhancer_start <= snp.bp & enhancer_end >= snp.bp)$cluster_id
#       chr enhancer_start enhancer_end cluster_id        GHid is_elite regulatory_element_type
# 25397   2      112762481    112765802      78883 GH02I112762        1       Promoter/Enhancer

head(subset(tbl.tfs, enhancer_cluster_id == snp.clusterID))
#         enhancer_cluster_id     TF              tissues
# 1803477               78883   HDGF              GM12878
# 1803478               78883 PKNOX1 K562;HEK293T;GM12878
# 1803479               78883  CLOCK                MCF-7
# 1803480               78883  SMAD1              GM12878
# 1803481               78883 ARID4B                HepG2
# 1803482               78883  SIN3A   MCF-7;H1-hESC;A549

subset(tbl.geneScoresAll, cluster_id == snp.clusterID)
#         cluster_id          symbol  eqtl_score erna_score chic_score expression_score distance_score combined_score is_elite
# 6932         78883          CKAP2L         NaN        NaN        NaN              NaN     0.77140794     550.771408        0
# 351416       78883     GC02P113658         NaN        NaN        NaN              NaN     0.33995102       0.339951        0
# 363923       78883            IL1A 6.27948e-06        NaN        NaN              NaN     0.33995102       5.542027        0
# 529643       78883          NT5DC4 1.81319e-05        NaN        NaN              NaN     0.26274150       5.004298        0
# 782088       78883        FLJ42351 9.21176e-05        NaN        NaN              NaN     0.12729262       4.162950        0
# 909257       78883          POLR1B         NaN        NaN        NaN     2.810298e-26     0.07446023      25.625708        0
# 966304       78883           RGPD8         NaN        NaN        NaN     5.735886e-13     0.05292404      12.294324        0
# 989773       78883            PSD4 9.77664e-07        NaN        NaN              NaN     0.04976234       6.059573        0
# 1072209      78883          ANAPC1         NaN        NaN        NaN     1.849851e-33     0.02382727      32.756691        0
# 1079290      78883 ENSG00000227359         NaN        NaN        NaN     4.257254e-11     0.02478953      10.395660        0

dim(subset(tbl.tissues, cluster_id == snp.clusterID))  # [1] 104   5
length(unique(subset(tbl.tissues, cluster_id == snp.clusterID)$tissue)) [1] 97

