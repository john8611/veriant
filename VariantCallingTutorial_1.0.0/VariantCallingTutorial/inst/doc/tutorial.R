
library(VariantCallingTutorial)

data(tallies)

calling.filters <- VariantCallingFilters()

post.filters <- VariantPostFilters()

variants <- callVariants(tallies,
                         calling.filters,
                         post.filters)

snvs <- variants[isSNV(variants)]

indels <- variants[isIndel(variants)]

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")
color.scheme <- cbPalette
my.theme <- lattice::trellis.par.get()
my.theme <- within(my.theme, {
  box.rectangle$col <- "black"
  box.rectangle$fill <- "gray90"
  box.umbrella$col <- "black"
  box.umbrella$lty <- 1
  strip.background$col <- rep("gray90", 7)
  plot.symbol$col <- "black"
  plot.symbol$pch <- 19
  box.dot$cex <- 1
  box.dot$pch <- "|"
  add.line$col <- "red"
  plot.polygon$col <- "gray90"
  superpose.symbol$fill <- color.scheme
  superpose.symbol$col <- color.scheme
  superpose.line$fill <- color.scheme
  superpose.line$col <- color.scheme
  superpose.polygon$fill <- color.scheme
  superpose.polygon$col <- color.scheme
})
lattice::lattice.options(default.theme = my.theme)

pdf(file="fig/cov.pdf",height=4,width=4)
densityplot(~ totalDepth(variants),
            xlim=c(0, 2*median(totalDepth(variants))),
            plot.points=FALSE, n=200)
dev.off()

pdf(file="fig/freq-cov-bin.pdf",height=4,width=8)
variants$coverage.bin <- cut(totalDepth(variants), c(0, 20, 80, Inf))
densityplot(~ altDepth/totalDepth | coverage.bin,
            as.data.frame(variants),
            plot.points=FALSE, layout=c(3, 1),
            xlab="variant frequency by coverage bin")
dev.off()

indel.windows <- indels + 10

snvs$near.indel <- ifelse(snvs %over% indel.windows,
                          "over indel", "off indel")

pdf(file="fig/freq-near-indel.pdf",height=4,width=8)
densityplot(~ altDepth/totalDepth | near.indel,
            as.data.frame(snvs),
            plot.points=FALSE, layout=c(2,1),
            xlab="variant frequency by coverage bin")
dev.off()

chr20.sequence <- getSeq(Hsapiens, "chr20")
chr20.hp <- ranges(Rle(as.raw(chr20.sequence)))

chr20.hp <- chr20.hp[width(chr20.hp) > 6L]

snvs$over.hp <- ifelse(ranges(snvs) %over% chr20.hp,
                       "over homopolymer",
                       "off homopolymer")

pdf(file="fig/freq-near-hp.pdf",height=4,width=8)
densityplot(~ altDepth/totalDepth | over.hp,
            as.data.frame(snvs),
            plot.points=FALSE, layout=c(2,1),
            xlab="variant frequency by coverage bin")
dev.off()

data(selfChains, package="VariantCallingTutorial")

unchained.filter <-
  FilterRules(list(unchained = function(x) {
    x %outside% selfChains
  }))
variants <- softFilter(variants, unchained.filter)
summary(softFilterMatrix(variants))

SRAdb::startIGV("lm")
sock <- SRAdb::IGVsocket()

mcols(variants) <- NULL
sampleNames(variants) <- "NA12878"
vcf <- writeVcf(sort(variants),
                "variants.vcf")
vcf <- tools::file_path_as_absolute(vcf)
vcf_gz <- paste0(tools::file_path_sans_ext(vcf), ".gz")
indexTabix(bgzip(vcf, vcf_gz, overwrite=TRUE),
           format="vcf4")

bam <- tools::file_path_as_absolute(bam)
session <- SRAdb::IGVsession(c(bam, vcf_gz),
                             "session.xml",
                             "hg19")

SRAdb::IGVload(sock, session)

rtracklayer::export(variants, "roi.bed")

vcf.file <- NA12878_pg.chr20.vcf.bgz
header <- scanVcfHeader(vcf.file)
header

geno(header)

info(header)["END",]

vcf <- readVcf(vcf.file, genome="hg19")

print(object.size(vcf), unit="auto")

ranges.chr20 <- as(seqinfo(Hsapiens)["chr20"],
                   "GRanges")

param <- ScanVcfParam(which=ranges.chr20)
vcf.chr20 <- readVcf(vcf.file, genome="hg19",
                     param=param)

colnames(mcols(rowData(vcf)))

head(colnames(info(vcf)))

names(geno(vcf))

illumina_variants <- vcf[geno(vcf)$GT[,1] != "0/0",]

head(alt(illumina_variants))

selectSNVs <- isSNV(illumina_variants)
illumina_snvs <- subset(illumina_variants, selectSNVs)

class(alt(illumina_snvs))

illumina_snvs <- expand(illumina_snvs)
class(alt(illumina_snvs))

print(object.size(illumina_snvs), unit="auto")

prefilters <-
  FilterRules(list(onlyVariants=function(text) {
    !grepl("0/0", text, fixed=TRUE)
  }))

filters <- FilterRules(list(onlySNVs=isSNV))

filterVcf(vcf.file, genome="hg19", "snvs.vcf", index=TRUE,
          prefilters=prefilters, filters=filters,
          param=ScanVcfParam(info=NA))

fixed(header(illumina_snvs))$FILTER

table(unlist(strsplit(filt(illumina_snvs), ";")))

passed <- grep("PASS", filt(illumina_snvs), fixed=TRUE)
illumina_snvs <- illumina_snvs[passed,]

rowData(illumina_snvs)$coverage.bin <-
  cut(geno(illumina_snvs)$DP, c(0, 20, 80, Inf))
table(rowData(illumina_snvs)$coverage.bin)

illumina_vr <- as(illumina_snvs, "VRanges")

seqlevelsStyle(illumina_vr) <- "NCBI"
illumina_vr <- dropSeqlevels(illumina_vr, "MT")
genome(illumina_vr) <- "GRCh37"

illumina_vr$in.vt <- illumina_vr %in% snvs
mean(illumina_vr$in.vt)

snvs$in.illumina <- snvs %in% illumina_vr
mean(snvs$in.illumina)

illumina_vr$vt.freq <-
  altFraction(snvs)[match(illumina_vr, snvs)]

pdf(file="fig/scatter.pdf",height=4,width=4)
xyplot(vt.freq ~ altFraction(illumina_vr),
       as.data.frame(illumina_vr),
       panel=panel.smoothScatter,
       xlab="Illumina frequency",
       ylab="VariantTools frequency")
dev.off()

runs <- vcf[!is.na(info(vcf)$END),]
end(rowData(runs)) <- info(runs)$END

calling.filters <- hardFilters(snvs)[3:5]
tallies <- resetFilter(tallies)
tallies <- softFilter(tallies, calling.filters,
                      serial=TRUE)
fn <- tallies[tallies %in% subset(illumina_vr, !in.vt)]
t(summary(softFilterMatrix(fn)))

gene.models <- TxDb.Hsapiens.UCSC.hg19.knownGene
snvs <- keepSeqlevels(snvs, "20")
seqlevelsStyle(snvs) <- "UCSC"
genome(snvs) <- "hg19"
locations <- locateVariants(snvs, gene.models,
                            CodingVariants())

colnames(mcols(locations))

snvs$coding.tx <- NA_integer_
snvs$coding.tx[locations$QUERYID] <- locations$TXID

gene_ids <- sub("GeneID:", "",
                locations$GENEID[!is.na(locations$GENEID)])
syms <- unlist(mget(gene_ids,
                    org.Hs.egSYMBOL,
                    ifnotfound=NA))
locations$SYMBOL[!is.na(locations$GENEID)] <- syms

coding <- predictCoding(snvs, gene.models, Hsapiens,
                        varAllele = DNAStringSet(alt(snvs)))

setdiff(colnames(mcols(coding)), colnames(mcols(snvs)))

table(coding$CONSEQUENCE)
