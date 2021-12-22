# Supporting code for manuscript: "Beyond to the stable: role of the insertion sequences as epidemiological descriptors in *Corynebacterium striatum*"

## 1. Assemblies recovery of *C. striatum* genomes.

NCBI provides [command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/quickstarts/command-line-tools/) that allow you to download genomes from databases, including GenBank.
For this study we used the following:

`./datasets download genome taxon "Corynebacterium striatum" --exclude-gff3 --exclude-protein --exclude-rna --assembly-source genbank`

Then, to sort the files.

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

The metadata was retrieved from [NCBI pathogens](https://www.ncbi.nlm.nih.gov/pathogens//isolates/#taxgroup_name:%22Corynebacterium%20striatum%22), then with R code, we enter the information for each previously downloaded genome.
 
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


## 2. Study of insertion sequences using Prokka and Panaroo.

Our aproach cosidered to use [Prokka](https://github.com/tseemann/prokka) for the annotation. 

```
for file in *.fna; do
    ../../benjamin_leyton/prokka/bin/prokka $file --cpus 0 --kingdom Bacteria --proteins --prefix ${file%%.fna} --locustag ${file%%.fna} --centre Bioren --compliant --addgenes --outdir ${file%%.fna};
done
```

To rebuild the pangenome we use [Panaro](https://github.com/gtonkinhill/panaroo). 

`panaroo-qc -t 8 --graph_type all -i *.gff --ref_db ../../benjamin_leyton/refseq.genomes.k21s1000.msh -o panaroo-qc` # for QC

`nohup panaroo -i ./gff/*.gff -o ./panaroo/ --clean-mode moderate -a core -- core_threshold 0.95 -t 8 -f 0.6 --refind_prop_match 0.5 --search_radius 6000 > panaroo.log &`

We then got a subset of all transposases related with IS of Panroo's results with a R code:

```
library(xlsx)
panaroo_striatum <- read.csv(file = "gene_presence_absence_roary.csv")
filas <- grep("transposase", panaroo_striatum$Annotation, ignore.case = TRUE)
IS_striatum <- panaroo_diphtheriae[c(filas_diph),]
write.xlsx(IS_striatum, "IS_diphtheriae.xlsx")
```
Finally, to create the matrix of presence and abundance of IS (transposases), it is enough to "count" the number of loci in the Exel's cells. A Excel formula like this will work:

=SI(ESBLANCO(H93);"";LARGO(H93)-LARGO(SUSTITUIR(H93;" ";""))+1)

=IF(ISBLANK(H93);"";LARGE(H93)-LARGE(SUBSTITUTE(H93;""; ""))+1)

## 3. Classify the studied genomes into lineages

SNP-based alignment was performed using Parsnp. [Pasnp](https://harvest.readthedocs.io/en/latest/content/parsnp.html) makes a global alignment of the core genome and identifies the SNPs

`parsnp -d ./genomas/fna/ -r ./referencia/GCA_002803965.1.fna -c -p 8 -o ./parsnp/`

then we transform the output of parsnp to a conventional alignment output.

`harvesttools -i parsnp.ggr -M parsnp.aln`

to build a tree based on the core alignment: 

`iqtree -s parsnp.aln -m GTR -pre arbol -bb 1000 -nt auto`

To lineages detection [RhierBAPS](https://github.com/gtonkinhill/rhierbaps) ( hierarchical Bayesian Analysis of Population Structure) was used.

```
library(ggtree)
library(phytools)
library(ape)
library(xlsx)

set.seed(1234)
fasta.file.name <- "parsnp.aln"
snp.matrix <- load_fasta(fasta.file.name)
hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 20, quiet = TRUE)
#hb.results  <- hierBAPS ( snp.matrix , max.depth  =  2 , n.pops  =  20 , n.extra.rounds  =  Inf , 
     #quiet  =  TRUE )
head(hb.results$partition.df)
newick.file.name  <- arbol.tree
iqtree  <-  phytools :: read.newick ( newick.file.name )
pdf("rhierbaps_results_circular.pdf")
gg  <- ggtree ( iqtree , layout  =  " circular " )
gg  <-  gg % < + % hb.results $ partición.df 
gg  <-  gg  + geom_tippoint (aes ( color  =  factor ( `level 1` )))
gg
dev.off()
write.xlsx(hb.results$partition.df, "levels.xlsx")
```

This generates two outputs, 1) a colorful tree with hierarchical levels and 2) an excel table with the hierarchical level for each genome studied.
 

## 4. Recombinacion/Mutacion.

Recombination / mutation analyzes were carried out using [Gubbins](https://github.com/nickjcroucher/gubbins) and [ClonalFrame](https://github.com/xavierdidelot/ClonalFrameML).

`run_gubbins parsnp.aln`

`ClonalFrameML arbol.tree parsnp.aln output_CFML`

## 5. Constrained PCoA analysis

To find and explore the relationships of the characteristics of *C. striatum* with ISs, CPCoA was used.

We define the lineages predicted by RhierBAPS as a restrictive variable. Our restricted PCoA analysis was validated using the permanova test. Permanova (Permutacional multivariable analysis of variance) is a non-parametric test based on dissimilarities. Permanova has the null hypothesis that groups do not differ in spread or position in multivariable space. In this study we did a permanova test with 10,000 permutations.

 You can download the Rmd and metadata [here](https://github.com/Leytoncito/IS_Prokka_Panaroo/tree/main/Constrained_PCoA) and run it in Rstudio.

## 6. Random Forest

As we mentioned in the manuscript Random Forest has the advantage of incorporating a vast and diverse number of data to predict characteristics. We define the response variable as the number of ARGs, lineages and country. Then we defend insertion sequences as predictor variables. We rank the variables according to their importance in the classification and then we graph the 5% of the most important variables for the classification of characteristics.

We use the kappa concordance test to measure the performance of the classification, and Cross validation to evaluate the Random Forest model. You can see the code [here](https://github.com/Leytoncito/IS_Prokka_Panaroo/blob/main/RandomForest)

