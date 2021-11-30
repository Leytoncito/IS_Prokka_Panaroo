# Suporting code for manuscript: "Beyond to the stable: role of the insertion sequences as epidemiological descriptors in *Corynebacterium striatum*"

## 1. Genome recovery *C. striatum*.

There are several ways to retrieve genomes in NCBI, one way is through the command line. NCBI provides [command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/quickstarts/command-line-tools/) that allow you to download genomes from databases, including GenBank.
For this study we used the following:

`./datasets download genome taxon "Corynebacterium striatum" --exclude-gff3 --exclude-protein --exclude-rna --assembly-source genbank`

Then to sort the files:This is especially useful for using the access ID, such as a fasta filename.

`unzip *.zip`
```
for f in ./ncbi_dataset/data/*; do
   if [ -d "$f" ]; then
      echo "$f"
      for i in ./${f%}/*.fna; do
          mv $i ${f%}.fna
      done
  fi
 done
 ```

The metadata was retrieved from [NCBI pathogens](https://www.ncbi.nlm.nih.gov/pathogens//isolates/#taxgroup_name:%22Corynebacterium%20striatum%22), then with a bit of R code, we input the information to each previously downloaded genome.
 
```
for f in *.fna; do echo ${f%%.fna}; done > genomas.txt
R
```
```
t1<-read.table("genomas.txt", header = FALSE, sep = "\t")
t2<-read.csv("isolates.csv", sep = ",")
t2$Assembly=as.factor(t2$Assembly)
t1$V1=as.factor(t1$V1)
t3<-merge(t1, t2, by, by.x = "V1", by.y = "Assembly", all.x=TRUE, sort=FALSE)
write.csv(t3, "metadatos.csv")
```
After this, we already have genomes and their respective neat and clean information.

## 2. Study of insertion sequences using Prokka and Panaroo.

```
for file in *.fna; do
    ../../benjamin_leyton/prokka/bin/prokka $file --cpus 0 --kingdom Bacteria --proteins --prefix ${file%%.fna} --locustag ${file%%.fna} --centre Bioren --compliant --addgenes --outdir ${file%%.fna}
done
```

`panaroo-qc -t 8 --graph_type all -i *.gff --ref_db ../../benjamin_leyton/refseq.genomes.k21s1000.msh -o panaroo-qc`

```
library(xlsx)
panaroo_striatum <- read.csv(file = "gene_presence_absence_roary.csv")
filas <- grep("transposase", panaroo_striatum$Annotation, ignore.case = TRUE)
IS_striatum <- panaroo_diphtheriae[c(filas_diph),]
write.xlsx(IS_striatum, "IS_diphtheriae.xlsx")
```
=SI(ESBLANCO(H93);"";LARGO(H93)-LARGO(SUSTITUIR(H93;" ";""))+1)

## 3. Coregenome alignment using parsnp  phylogenetic tree using iqtree

parsnp -d ./genomas/fna/ -r ./referencia/GCA_002803965.1.fna -c -p 8 -o ./parsnp/
harvesttools -i parsnp.ggr -M parsnp.aln
iqtree -s parsnp.aln -m GTR -pre arbol -bb 1000 -nt auto

## 4. Recombination Analysis Using ClonalFrame and Gubbins

ClonalFrameML newick_file aln_file output_file 
run_gubbins aln_file 

11. Classification of lineages using Rhierbaps.

library(rhierbaps)
library(ggtree)
library(phytools)
library(ape)
library(xlsx)

set.seed(1234)
fasta.file.name <- "parsnp.aln" #Aqui tambien se puede usar un aln enmascarado libre de recombinacion producido por gubbins o clonal_frame
snp.matrix <- load_fasta(fasta.file.name)
hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 20, quiet = TRUE)
#hb.results  <- hierBAPS ( snp.matrix , max.depth  =  2 , n.pops  =  20 , n.extra.rounds  =  Inf , 
     #quiet  =  TRUE )
head(hb.results$partition.df)
newick.file.name  <- arlbol.tree
iqtree  <-  phytools :: read.newick ( newick.file.name )
pdf("rhierbaps_results_circular.pdf")
gg  <- ggtree ( iqtree , layout  =  " circular " )
gg  <-  gg % < + % hb.results $ partición.df 
gg  <-  gg  + geom_tippoint (aes ( color  =  factor ( `level 1` )))
gg
dev.off()
write.xlsx(hb.results$partition.df, "levels.xlsx")


13. Distance matrix

15. Constrained PCoA analysis

17. Random Forest

