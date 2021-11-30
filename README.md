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

3. Annotation using Prokka.
4. Construction of pangenome using Panaroo.
5. Coregenome alignment using parsnp
6. phylogenetic tree using iqtree
7. Recombination Analysis Using ClonalFrame and Gubbins
8. Classification of lineages using Rhierbaps.
9. Distance matrix
10. Constrained PCoA analysis
11. Random Forest

