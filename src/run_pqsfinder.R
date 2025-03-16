library(seqinr)
library(pqsfinder)
library(rtracklayer)
library(Biostrings)

# Process command line arguments
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen != 4) {
  stop("Specify the fasta file, output file, minimum score, and overlapping flag (1 or 0, for True or False)", call.=FALSE)
}

FastaFile = paste0(args[1])
OutputFile = paste0(args[2])
MinScore = as.numeric(paste0(args[3]))
Overlapping = as.numeric(paste0(args[4]))

sprintf("FastaFile: %s", FastaFile)
sprintf("OutputFile: %s", OutputFile)
sprintf("MinScore: %s", MinScore)
sprintf("Overlapping: %s", Overlapping)

# Read the sequence in the Fasat file
print("Reading the fasta file")
seq_list<-read.fasta(FastaFile, as.string=TRUE)
seq = paste(seq_list[[1]], collapse="")

# Convert seq to DNAString
dnaString = DNAString(seq)

print("Running pqsfinder")
# Run pqsfinder
#pqs <- pqsfinder(dnaString, min_score=MinScore, overlapping=(Overlapping == 1)) # Prints the status of the search
pqs <- suppressMessages(suppressWarnings(pqsfinder(dnaString, min_score=MinScore, overlapping=(Overlapping == 1))))

# Export all PQS into a GFF3-formatted file
print("Exporting pqsfinder output")
writeXStringSet(as(pqs, "DNAStringSet"), file = OutputFile, format = "fasta")
