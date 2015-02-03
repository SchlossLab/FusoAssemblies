---
layout: post
title:  "What I've done so far"
date:   2014-11-05
comments: true
output:
  html_document:
    keep_md: yes
---



Since Summer 2014 I've been working on creating a pipeline to assemble metagenomic shotgun reads from the Human Microbiome Project. So far I've been working with synthetic metagenomic data to try different methods. Eventually I will have a pipeline that I can use to look at variation in Fusobacterium nucleatum genomes between individuals and body sites. 

###Practice with test data

CONCOCT was recently published in *Nature methods* by [Alneburg et al.](http://www-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pubmed/?term=binning+metagenomic+contigs+by+coverage+and+composition) They introduce a new metagenomic assembly method that bins contigs based on read coverage and kmer frequency information. 

In the paper they use two synthetic mock metagenomic data sets. The first called "species" has 101 species that were identified from HMP data using 16S. Once the organisms were chosen, they compiled all the chromosome and plasmid information from NCBI for those 101 genomes. They also created a second data set that included fewer organisms (20 organisms abundant in the gut), but these include multiple strains of E.coli, Bacteroides, and Clostridium. In both datasets the frequencies of the organisms were similar to frequencies in the HMP.

From these large genome compilations, Alneburg et al. similulated reads that would be obtained from an Illumina HiSeq. The species mock has 96 samples with 7.75 million paired-end reads each. The strain mock has 64 samples with 11.75 million paired-end reads each (~2GB per sample). This is small compared to the HMP samples which can reach 69 million reads (need to look up exact range). In the simulated reads they included error profiles from a real data set. 

```
working directory: /mnt/EXT/Schloss-data/amanda/Fuso/concoct/testdata/
Data obtained from: https://export.uppmax.uu.se/b2010008/projects-public/concoct-paper-data/
```

I used the velvet script shuffleSequences_fasta to interleave all of the samples combined into one big file. This ended up taking a really long time (2 days on axiom).

```
cd $CONCOCT_SPECIES
zcat Sample*R1*.fasta.gz > run1/All_R1.fa #unzips all files and cats to single file
zcat Sample*R2*.fasta.gz > run1/All_R2.fa 
quicksubmit "perl /mnt/EXT/Schloss-data/amanda/velvet/velvet_1.2.10/contrib/shuffleSequences_fasta/shuffleSequences_fasta.pl All_R1.fa All_R2.fa All.fa" --pm nodes=1:ppn=8,mem=48GB --jobname shuffle5 --walltime 999:00:00 --cput 999:00:00
```

###Digital normalization

Digital normalization is a method created by Titas Brown which will remove highly redundant sequences to decrease the file size, remove errors, and make it easier to assemble using a single-genome assembler like velvet. There is pretty good documentation on the [read the docs](http://khmer.readthedocs.org/en/v0.5/scripts.html). I had trouble getting khmer to work with the permissions that I have on axiom, so I followed the directions to create a khmer environment. I made the alias in my bash_profile so that khmerEnv will activate the environment. This works when running jobs in the queue as well.

```
quicksubmit "khmerEnv; python2.7 /mnt/EXT/Schloss-data/amanda/Fuso/khmer/khmerEnv/bin/normalize-by-median.py -k 41 All.fa -o normalized" --pm nodes=1:ppn=8,mem=48GB --jobname DN_k41 --walltime 999:00:00 --cput 999:00:00
```

Output:

```
DONE with All.fa; kept 190653 of 1488000000 or  0%
output in All.fa.keep
fp rate estimated to be 1.000
```

There was a problem with the -p command. It said that the paired-end reads weren't in the correct format. Went ahead and ran DN without the paired option. Will need to later do DN again and include a cutoff of 10 or 20 to speed up the normalization. It cut out way too many of the reads which will affect the velvet assembly.

###Velvet assembly

Velvet is the assembler that Alneburg et al. used in their online example, so I figured I would try it the same way that they did. In the online [example](https://concoct.readthedocs.org/en/latest/complete_example.html) that includes only 4 species, they suggest to just assemble the paired end reads with a kmer of 71. This wouldn't work with my 250GB file, so I am forced to do DN first. The test data in the paper was assembled using k=41. I will assemble the normalized dataset. I would also like to be able to assemble the raw reads if I can figure out how.


```
velveth velveth2_k32 71 All.fa.keep
velvetg velveth2_k32 -cov_cutoff auto

cd $CONCOCT_SPECIES/run1/velveth_k41

awk '{/>/&&++a||b+=length()}END{print b/a}' contigs.fa #find average sequence length
292.795 bp
grep -c '>' contigs.fa #counts number of contigs
2556

awk '!/^>/ {next} {getline s} length(s) >= 500 { print $0 "\n" s }' contigs.fa > contigs.500.fa #pulls out all the contigs greater than 500bp

grep -c '>' contigs.500.fa
0 contigs greater than 500bp

Log:
Median coverage depth = 2.023810
Final graph has 2573 nodes and n50 of 331, max 1837, total 646309, using 31194/1
90653 reads
```



###Megahit assembly

Titas brown says that megahit is [pretty good](http://ivory.idyll.org/blog/2014-how-good-is-megahit.html) and it is faster and more memory efficient. Digital normalization is not needed. I cloned [megahit](https://github.com/voutcn/megahit) into /mnt/EXT/Schloss-data/amanda/Fuso/megahit/megahit. It's pretty simple to use. It will take fasta, fastq, or zipped files. I don't think it can take paired read information.

```
quicksubmit "python ./megahit -m 45e9 --input-cmd 'cat /mnt/EXT/Schloss-data/amanda/Fuso/concoct/testdata/species/run1/All_*' --cpu-only -l 100 -o /mnt/EXT/Schloss-data/amanda/Fuso/concoct/testdata/species/run1/megahit" --pm nodes=1:ppn=8,mem=48GB --jobname megahit2 --walltime 999:00:00 --cput 999:00:00

awk '{/>/&&++a||b+=length()}END{print b/a}' final.contigs.fa #find average sequence length
3650.39 bp
grep -c '>' final.contigs.fa #counts number of contigs
100419

awk '!/^>/ {next} {getline s} length(s) >= 500 { print $0 "\n" s }' final.contigs.fa > final.contigs.500.fa #pulls out all the contigs greater than 500bp

grep -c '>' final.contigs.500.fa
19663 contigs greater than 500bp

awk '!/^>/ {next} {getline s} length(s) >= 1000 { print $0 "\n" s }' final.contigs.fa > final.contigs.1000.fa
grep -c '>' final.contigs.1000.fa
12741 contigs greater than 1kb
```

There are 12741 contigs greater than 1kb. In the CONCOCT paper they processed the contigs to filter out any less than 1kb because the composition vectors with the frequencies of each kmer isn't accurate when contigs are shorter than 1kb. This may throw away some of the rare bugs that couldn't assemble to greater than 1kb, but it will also get rid of errors from too short contigs. 

##Ray assembly

[Ray](http://denovoassembler.sourceforge.net/) is the assembler used in the CONCOCT paper, so I also tried to assemble the reads using this method. The method requires an MPI to run, which is what makes this method superior as it can run in parallel on multiple nodes. Since we don't have the ability to run MPI on multiple nodes on axiom, I assembled using just one node.

```
mpiexec -n 1 ./Ray -k 32 -p /mnt/EXT/Schloss-data/amanda/Fuso/concoct/testdata/species/run1/All_R1.fasta /mnt/EXT/Schloss-data/amanda/Fuso/concoct/testdata/species/run1/All_R2.fasta -o /mnt/EXT/Schloss-data/amanda/Fuso/concoct/testdata/species/run1/ray
```

##CONCOCT processing

Once the prelimnary assembly is done, the coverage table can be generated. First, cut up the contigs into chunks less than 10kb so that there isn't a bias towards mapping onto long contigs. Then, map the reads onto this "reference" of the contigs.

To run concoct on axiom, I have to enter the concoct environment. I made an alias so just type "concoctenv".

```
cd $CONCOCT_SPECIES/run1
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m megahit/final.contigs.fa > megahit/megahit.congits_c10K.fa
bowtie2-build megahit/megahit.congits_c10K.fa megahit/megahit.congits._c10K.fa
quicksubmit "bowtie2-build megahit/megahit.congits_c10K.fa megahit/megahit.congits._c10K.fa" --pm nodes=1:ppn=8,mem=48GB --jobname bowtie_megahit --walltime 999:00:00 --cput 999:00:00 
```

The next step is to map all of the reads to the assembled contigs. I was having a hard time trouble shooting this huge for loop, so I put the whole thing into a bash script (/Users/Amanda/Documents/Schloss/Fuso/concoct/paperdata/bash script/markdup_bash4.sh) so I could run quicksubmit "sh markdup_bash4.sh" to submit the job. Here is the final script:

```
#!/bin/bash
for f in $CONCOCT_SPECIES/Sample*R1*.fasta.gz; do 
  gunzip -c $f > $(echo $f | sed s/".gz"/""/)
	gunzip -c $(echo $f | sed s/"R1"/"R2"/) > $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/)
	mkdir -p $CONCOCT_SPECIES/map/$(basename $f .gz)
	cd $CONCOCT_SPECIES/map/$(basename $f .gz)
	$CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $(echo $f | sed s/".gz"/""/) $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/) pair $CONCOCT_SPECIES/run1/megahit/megahit.congits_c10K.fa asm_megahit bowtie2_megahit
	rm $(echo $f | sed s/".gz"/""/)
	rm $(echo $(echo $f | sed s/".gz"/""/) | sed s/"R1"/"R2"/)
	cd ../..
done
```

Once the bam tables have been made, I used the gen_input_table.py script to make a single table the includes all the coverage information from each sample per contig.

```
cd $CONCOCT_SPECIES/map
python $CONCOCT/scripts/gen_input_table.py --isbedfiles --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) $CONCOCT_SPECIES/run1/megahit/megahit.congits_c10K.fa */bowtie2_megahit/asm_megahit_pair-smds.coverage> concoct_inputtable.tsv
mkdir $CONCOCT_SPECIES/run1/concoct-input
mv concoct_inputtable.tsv $CONCOCT_SPECIES/run1/concoct-input/
```

Then cut the file to input into concoct.

'''
cut -f1,3-26 concoct-input/concoct_inputtable.tsv > concoct-input/concoct_inputtableR.tsv
'''

Finally, run concoct. I started by using these parameters

'''
concoct -c 40 --coverage_file concoct-input/concoct_inputtableR.tsv --composition_file $CONCOCT_SPECIES/run1/megahit/megahit.congits_c10K.fa -b concoct-output/
'''

It finished! And pretty quickly, too. Here is the logfile:

'''
2014-11-12 16:52:40,484:INFO:root:Results created at /mnt/EXT/Schloss-data/amand
a/Fuso/concoct/testdata/species/run1/concoct-output
2014-11-12 16:58:44,634:INFO:root:Successfully loaded composition data.
2014-11-12 16:58:46,519:INFO:root:Successfully loaded coverage data.
2014-11-12 16:58:48,360:INFO:root:Performed PCA, resulted in 0.9 dimentions
2014-11-12 16:58:55,050:INFO:root:Wrote original filtered data file.
2014-11-12 16:58:56,083:INFO:root:Wrote PCA transformed file.
2014-11-12 16:58:56,087:INFO:root:Wrote PCA components file.
2014-11-12 16:58:56,087:INFO:root:PCA transformed data.
2014-11-12 16:58:56,087:INFO:root:Will call vbgmm with parameters: concoct-outpu
t/, 40, 1000
2014-11-12 17:02:33,220:INFO:root:CONCOCT Finished
'''

















