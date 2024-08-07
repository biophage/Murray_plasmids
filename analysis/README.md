# Analysis of PAE plasmids and their modern relatives

This document describes the software and commands used to analyse the plasmid sequences reported in our study. Further details on the rational and different steps of the analysis can be found in the methods section of our manuscript.  

> ___Software list___:
> * enaBrowserTools: https://github.com/enasequence/enaBrowserTools
> * unicycler: https://github.com/rrwick/Unicycler
> * quast: https://github.com/ablab/quast
> * kraken2: https://github.com/DerrickWood/kraken2
> * bracken: https://github.com/jenniferlu717/Bracken
> * plasmidverify: https://github.com/ablab/plasmidVerify
> * circlator: https://github.com/sanger-pathogens/circlator
> * abricate: https://github.com/tseemann/abricate
> * mob_suite: https://github.com/phac-nml/mob-suite
> * dedupe: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/dedupe-guide/
> * minimap2: https://github.com/lh3/minimap2
> * CD-HIT: https://sites.google.com/view/cd-hit
> * prokka: https://github.com/tseemann/prokka
> * Integron finder: https://github.com/gem-pasteur/Integron_Finder
> * mash: https://github.com/marbl/Mash
> * BLAST: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
> * Panaroo: https://github.com/gtonkinhill/panaroo
> * Clinker: https://github.com/gamcil/clinker

## Murray plasmids

### Retrieving Murray sequencing reads

_Software_: enaBrowserTools v1.6  
_Command_:  
```
enaGroupGet -f fastq -d OUTPUT_DIR/ ACCESSION
```

### Genome Assembly

_Software_: unicycler v0.4.7  
_Command_:  
```
unicycler -1 READS_1.fastq.gz -2 READS_2.fastq.gz -o OUTPUT_DIR
```

### Assembly metrics

_Software_: quast v4.5  
_Command_:  
```
quast -o OUTPUT_DIR/ unicycler/ASSEMBLY.fasta
```

### Taxonomic classification

_Software_: kraken2 v2.1.2 & bracken v2.5  
_Commands_:  
___kraken___
```
kraken2 --use-names --threads 8 --db k2_standard_8gb_20210517 --fastq-input --report OUTPUT.REPORT-reads.txt --gzip-compressed --paired READS_1.fastq.gz  READS_2.fastq.gz --output OUTPUT-reads.kraken
```
___bracken___
```
bracken -d k2_standard_8gb_20210517 -i OUTPUT.REPORT-reads.txt -o OUTPUT-reads.bracken -l S
```

### Chromosomal contigs depletion
_Software_: plasmidverify.py  
_Command_:  
```
python2.7 ~/plasmidVerify/plasmidverify.py -f unicycler/ASSEMBLY.fasta -o PVERIFY/OUTPUT_DIR --hmm ~/Pfam-HMM/Pfam-A.hmm.gz -t 8
```

### Circularisation of plasmid sequences
_Software_: circlator v1.5.5  
_Command_:  
```
circlator all --threads 8 --b2r_min_read_length 50 --merge_min_length 100 --merge_min_length_merge 200 --assemble_spades_k 107,97,87,77,67,57,51 --b2r_discard_unmapped --clean_min_contig_length 100 /data/${d}-unicycler_noncircularised.fasta interleaved_reads.fastq circlator/OUTPUT_DIR
```

### Typing
_Software_: abricate v1.0.1, mob_suite v3.0.1  
_Commands_:  
___abricate - plasmidfinder___
```
abricate -db plasmidfinder --minid 80 --mincov 60 --quiet --nopath --threads 8 Post-circlator/PLASMID_SEQUENCES.fna > PLASMIDFINDER/Post-circlator_seqs-plasmidfinder_80id_60cov.txt
```
___mob_typer___
```
mob_typer --infile Post-circlator/PLASMID_SEQUENCES.fna --out_file mob_typer/OUTPUT_report.txt --num_threads 8 --multi
```

### Plasmid identification 
_Software_: dedupe (BBTools package v38.90), minimap2 v2.17-r941, mob_suite v3.0.1  
_Commands_:  
___mapping against known plasmids___  
- Deduplicate reference database
```
bash ~/bbmap/dedupe.sh in=iPlasmidDB-core_all.fna threads=8 sort=id out=dedupe_out.fasta ac=f outd=duplicates.fasta
```
- Sequences mapping  
```
minimap2 -x asm5 -t 8 -c dedupe/dedupe_out.fasta PLASMID_SEQUENCES.fna -o minimap2/OUTPUT.paf
```
___mob_recon___  
```
mob_recon --infile Post-circlator/PLASMID_SEQUENCES.fna --outdir mob_recon/OUTPUT_DIR -t -u -c -n 8
```


## Plasmids sequence characterisation

### AMR/Virulence databases clustering
_Software_: CD-HIT v4.8.1  
_Commands_:  
___Concatenate databases___
```
cat card/sequences ncbi/sequences argannot/sequences resfinder/sequences > cd-hit_res/sequences-all
```
___cd-hit-est___  
```
cd-hit-est -i sequences-all -o sequences -c 1.00 -n 10 -T 8 -s 0.9 -sc 1
```
_The above are examples of concatenating & clustering the AMR databases; the same was done for the virulence databases VFDB and Ecoli_VF._

### AMR/Virulence genes identification
_Software_: abricate v1.0.1  
_Command_:  
```
abricate -db cd-hit_res --minid 80 --mincov 80 --quiet --nopath --threads 8 --fofn sequences/input_fna_files.txt > abricate/amr/OUTPUT_REPORT.80cov-cd-hit_res.tsv
```
_The above is an example of searching for AMR genes. A similar command was used to search for virulence factors using the relevant database._

### Plasmid annotation
_Software_: prokka v1.14.6  
_Command_:
```
prokka --outdir prokka/OUTPUT_DIR --prefix PLASMID_NAME --locustag PLASMID_NAME --plasmid PLASMID_NAME --kingdom Bacteria --gcode 11 --cpus 8 PLASMID_SEQ.fna
```

### Integrons detection
_Software_: Integron finder v2.0.1  
_Command_:  
```
integron_finder --local-max --func-annot --promoter-attI --circ --cpu 8 --outdir integron_finder/OUTPUT_DIR PLASMID_SEQ.fna
```


## Identification of Murray plasmid relatives

### Genome distance estimation
_Software_: mash v2.3  
_Commands_:  
___Make sequence sketches___  
```
mash sketch -k 14 -p 8 -o mash/sketches/OUTPUT_FILE PLASMID_SEQS.fna
```
_Mash sketches were created for both the Murray and iPlasmidDB plasmid nucleotide sequences._  
  
___Estimate distance between sketches___  
```
mash dist -p 8 -d 0.05 mash/sketches/iPlasmidDB.msh MURRAY_PLASMID_SEQ.fna > mash/dist/OUTPUT_FILE.txt
```

### Nucleotide sequence similarity - BLAST
_Software_: BLAST v2.10.1+  
_Commands_:  
___Make a nucleotide blast database___
```
makeblastdb -dbtype nucl -in iPlasmidDB.fna -out blastndb/iPlasmidDB
```
___Search for homologous sequences in the BLAST database___  
```
blastn -query MURRAY_PLASMID_SEQ.fna -db blastndb/iPlasmidDB -evalue 1e-05 -num_threads 8 -subject_besthit -outfmt "6 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore qcovhsp qcovs" | awk -F"\t" '{if ($NF >= 70) {print $0}}' > blastn/OUTPUT_FILE.txt
```
_To retrieve all matches, as we did to identify all plasmids from iPlasmidDB related to Murray plasmid sequences, remove the part of the code filtering out matches with query coverage below 70% (`| awk -F"\t" '{if ($NF >= 70) {print $0}}'`)._


## Network clustering and comparative analysis

### All-vs-all mash distance estimation
_Software_: mash v2.3  
_Commands_:  
___Combine Murray and iPlasmidDB plasmid sketches into a single file___
```
find mash/sketches/ -name *msh  > input_msh_files.txt
```
```
mash paste mash/sketches/Murray-iPlasmidDB_plasmids -l input_msh_files.txt
```
___Calculate all-vs-all mash distance___
```
mash dist -p 8 -d 0.05 mash/sketches/Murray-iPlasmidDB_plasmids.msh -l input_msh_files.txt > mash_network/dist/Murray-iPlasmidDB_plasmids-dist_210shashes.tsv
```
_The above output file was used as input to create the plasmids network in Cytoscape_

### Pangenome analysis
_Software_: Panaroo v1.2.9  
_Command_:  
```
panaroo -i input_files/*_prokka.gff -o panaroo/OUTPUT_DIR --threshold 0.7 --len_dif_percent 0.7 -t 8 --clean-mode sensitive --remove-invalid-genes --no_clean_edges
```

### Pairwise comparisons
_Software_: Clinker v0.0.28  
_Command_:  
```
clinker input_files/*_prokka.gbk -p clinker/OUTPUT.html
```
