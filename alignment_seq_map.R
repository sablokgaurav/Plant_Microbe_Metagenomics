# a streamline analysis workflow for generating all
# the alignments and the alignment maps and the sequence
# distribution
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#SBATCH -J self.name 
#SBATCH -p constraint="snb|hsw 
#SBATCH -p self.queue 
#SBATCH -n self.threads 
#SBATCH -c self.core 
#SBATCH --mem=self.memory 
#SBATCH --workdir = self.change 
#SBATCH --mail = self.mail 
#SBATCH --mail-type=END'

for f in *.fasta; do echo muscle -in \
       $f -out $f.muscle.fasta; done
for f in *.fasta; do echo mafft --$thread --threadtb $thread \
    --threadit 0 --reorder --auto $f > $f.mafft.fasta; done
for f in *.fasta; do echo prank -d=$f \
    -o=$f.prank.fasta; done 

library(ggmsa)
getwd()
dir <- setwd("set_your_working_directory")
alignment_mafft <-list.file(path=dir, pattern="*.mafft.fasta", full.names=TRUE)
for (i in alignment_mafft){
   print(i)
}
for (f in alignment_mafft){
seq <- system.file("extdata", "$f", package = "ggmsa") \
ggmsa(seq, start = 221, end = 280, \
   char_width = 0.5, seq_name = TRUE) + geom_seqlogo() + geom_msaBar() \
}
alignment_muscle <-list.file(path=dir, pattern="*.muscle.fasta", full.names=TRUE)
for (i in alignment_muscle){
   print(i)
}
for (f in alignment_mafft){
seq <- system.file("extdata", "$f", package = "ggmsa") \
ggmsa(seq, start = 221, end = 280, \
   char_width = 0.5, seq_name = TRUE) + geom_seqlogo() + geom_msaBar() \
}
alignment_prank <-list.file(path=dir, pattern="*.prank.fasta", full.names=TRUE)
for (i in alignment_mafft){
   print(i)
}
for (f in alignment_mafft){
seq <- system.file("extdata", "$f", package = "ggmsa") \
ggmsa(seq, start = 221, end = 280, \
   char_width = 0.5, seq_name = TRUE) + geom_seqlogo() + geom_msaBar() \
}