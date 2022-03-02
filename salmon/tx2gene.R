#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("rhdf5")
library('tximport')
library("argparse")

# arguments parsing
parser <- ArgumentParser()
parser$add_argument("--input",default=TRUE)
parser$add_argument("--tx2gene",default=TRUE,help="")
parser$add_argument("--tool", type="character",default=TRUE,help="Specify expression counting tool")
parser$add_argument("--output",default=TRUE)
args <- parser$parse_args()

# transcript 2 gene annotation
tx2gene <- read.table(args$tx2gene,skip=5,sep='\t',header=FALSE)
names(tx2gene) <- c('TXNAME','GENEID')

# readCounts
txi <- tximport(args$input, type=args$tool, tx2gene=tx2gene, ignoreAfterBar=TRUE)

txi$counts <- cbind(gene = rownames(txi$counts), txi$counts)
rownames(txi$counts) <- 1:nrow(txi$counts)
colnames(txi$counts) <- c('gene','readCounts')

txi$length <- cbind(gene=rownames(txi$length), txi$length)
rownames(txi$length) <- 1:nrow(txi$length)
colnames(txi$length) <- c('gene','length')

# output 
txi_out <- merge(txi$length, txi$counts,key='gene')
write.table(txi_out,file=args$output,row.names=FALSE,sep='\t',quote=FALSE)
