# Analysis of PAE plasmids and their modern relatives

This document describes the software and commands used to analyse the plasmid sequences reported in our study.  


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
The above are examples of concatenating & clustering the AMR databases; the same was done for the virulence databases VFDB and Ecoli_VF.

### AMR/Virulence genes identification
_Software_: abricate v1.0.1  
_Command_:  
```
abricate -db cd-hit_res --minid 80 --mincov 80 --quiet --nopath --threads 8 --fofn sequences/input_fna_files.txt > abricate/amr/OUTPUT_REPORT.80cov-cd-hit_res.tsv
```

### Plasmid annotation
_Software_: prokka v1.14.6
_Command_:
```
prokka --outdir prokka/OUTPUT_DIR --prefix ${bname} --locustag ${bname} --plasmid ${bname} --kingdom Bacteria --gcode 11 --cpus 8 sequences/FNA/individual_files/${l}
```


#===== Integrons detection =====#
# Software: Integron finder v2.0.1
# Command:
integron_finder --local-max --func-annot --promoter-attI --circ --cpu 8 --outdir ~/Plasmids/Murray/Murray_collection/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/annotation/INTEGRON_FINDER/output_files ~/Plasmids/Murray/Murray_collection/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/sequences/FNAs/${l}.fna


#===== Genome distance estimation =====#
# Software: mash v2.3
# Command:
# Make sequence sketches
mash sketch -k 14 -p 8 -o MASH/sketches/${fname} ${l}
# Estimate distance between sketches
mash dist -p 8 -d 0.05 ~/Plasmids/Murray/iPlasmidDB/v1.1/mash/sketches/iPlasmidDB-genbank.msh ${l} > MASH/vs_iPlasmidDBv1.1/output_files/${fname}.out.txt


#===== BLAST comparison =====#
# Software: BLAST v2.10.1+
# Command:
# Make a nucleotide blast database
makeblastdb -dbtype nucl -in tmp -out blast/blastndb/iPlasmidDB-genbank
# Search for homologous sequences in the BLAST database
blastn -query FNAs/${l} -db ~/Plasmids/Murray/iPlasmidDB/v1.1/blast/blastndb/iPlasmidDB-genbank -evalue 1e-05 -num_threads 8 -subject_besthit > BLAST/vs_iPlasmidDBv1.1/BLASTN/output_files/${fname}.out.txt;


cd /home/linuxbrew/.linuxbrew/Cellar/abricate/1.0.1/libexec/db/
mkdir cd-hit_res
cat card/sequences ncbi/sequences argannot/sequences resfinder/sequences > cd-hit_res/sequences-all
cd cd-hit_res
cd-hit-est -i sequences-all -o sequences -c 1.00 -n 10 -T 8 -s 0.9 -sc 1
rm sequences-all
mv sequences.clstr ../abricate/sequences-amr.clstr
makeblastdb -in sequences -title cd-hit_res -dbtype nucl -hash_index





#===== Update Murray-iPlasmidDB plasmid (n=9977) master table =====#
# From: ~/Plasmids/Murray/Murray_collection/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1

## Clustering
perl -pe 's{",}{\t}g; s {,"}{\t}g; s{"}{}g; s{\t\t}{\t0\t}' clustering/MASH-CLUSTERMAKER/mcl_clusters-400shashes_threshold.csv | awk -F"\t" '{ if (NR != 1) {print $11,$1,$2}}' OFS="\t" | sort -k1,1 > clustering/MASH-CLUSTERMAKER/plasmids_cluster_and_lineage.tsv
sed -i '1s/^/ID\tCluster\tLineage\n/' clustering/MASH-CLUSTERMAKER/plasmids_cluster_and_lineage.tsv
# Check all lines have the same number of fields in the joined file
join -t $'\t' --header clustering/MASH-CLUSTERMAKER/plasmids_cluster_and_lineage.tsv Murray-iPlasmidDB_plasmids_master_table.txt | awk -F"\t" '{print NF}' | sort -u
join -t $'\t' --header clustering/MASH-CLUSTERMAKER/plasmids_cluster_and_lineage.tsv Murray-iPlasmidDB_plasmids_master_table.txt > tmp_joined; mv tmp_joined Murray-iPlasmidDB_plasmids_master_table.txt

## Number of CDS (prokka)
sort -k1,1 metadata/Murray-iPlasmidDB_plasmids-cds_count.txt > tmp_join
sed -i '1s/^/ID\tCDS_prokka\n/' tmp_join
# Check all lines have the same number of fields in the joined file
join -t $'\t' --header Murray-iPlasmidDB_plasmids_master_table.txt tmp_join | awk -F"\t" '{print NF}' | sort -u
join -t $'\t' --header Murray-iPlasmidDB_plasmids_master_table.txt tmp_join > tmp_joined; mv tmp_joined Murray-iPlasmidDB_plasmids_master_table.txt
rm tmp_join

## Transposases and integrases
sort -k1,1 annotation/PROKKA_products/transposases_and_integrases_per_plasmid.txt > tmp_join
sed -i '1s/^/ID\tTransposases\tIntegrases\n/' tmp_join
# Check all lines have the same number of fields in the joined file
join -t $'\t' --header Murray-iPlasmidDB_plasmids_master_table.txt tmp_join | awk -F"\t" '{print NF}' | sort -u
join -t $'\t' --header Murray-iPlasmidDB_plasmids_master_table.txt tmp_join > tmp_joined; mv tmp_joined Murray-iPlasmidDB_plasmids_master_table.txt
rm tmp_join

## BioSample accession
sort -k1,1 metadata/Murray-iPlasmidDB_plasmids-biosample_acc.txt > tmp_join
sed -i '1s/^/ID\tBioSample\n/' tmp_join
# Check all lines have the same number of fields in the joined file
join -t $'\t' --header tmp_join Murray-iPlasmidDB_plasmids_master_table.txt | awk -F"\t" '{print NF}' | sort -u
join -t $'\t' --header tmp_join Murray-iPlasmidDB_plasmids_master_table.txt > tmp_joined; mv tmp_joined Murray-iPlasmidDB_plasmids_master_table.txt

## Host and date of isolation
sort -k1,1 metadata/Murray-iPlasmidDB_plasmids-host_year.txt > tmp_join
sed -i '1s/^/ID\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tCollection_date\n/' tmp_join
# Check all lines have the same number of fields in the joined file
join -t $'\t' --header tmp_join Murray-iPlasmidDB_plasmids_master_table.txt | awk -F"\t" '{print NF}' | sort -u
join -t $'\t' --header tmp_join Murray-iPlasmidDB_plasmids_master_table.txt > tmp_joined; mv tmp_joined Murray-iPlasmidDB_plasmids_master_table.txt
rm tmp_join

## Integrons (Integron_finder)
sort -k1,1 metadata/Murray-iPlasmidDB_plasmids-integrons.txt > tmp_join
sed -i '1s/^/ID\tCALIN\tcomplete\tIn0\ttopology\n/' tmp_join
# Check all lines have the same number of fields in the joined file
join -t $'\t' --header Murray-iPlasmidDB_plasmids_master_table.txt tmp_join | awk -F"\t" '{print NF}' | sort -u
join -t $'\t' --header Murray-iPlasmidDB_plasmids_master_table.txt tmp_join > tmp_joined; mv tmp_joined Murray-iPlasmidDB_plasmids_master_table.txt
rm tmp_join

## Master table fields list
head -n1 Murray-iPlasmidDB_plasmids_master_table.txt | perl -pe 's{\t}{\n}g' | nl -w 1 > Murray-iPlasmidDB_plasmids_master_table-fields.txt

#___________________________________#



#===== Network manipulation in codon =====#

# Transfer files to codon
# From: ~/Dropbox/ESPOD_fell/03_Project/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network
scp Murray-iPlasmidDB_plasmids-dist_210shashes_figs.cys codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network

# Launching cytosacpe
# Read for details: https://github.com/leoisl/cytoscape_on_cluster/blob/master/README.md#usage
# Add the following to ~/.ssh/config
#Host *
#    ForwardX11Trusted yes
#    ForwardX11 yes
# From: /homes/acaza
module load openjdk-11.0.1-gcc-9.3.0-unymjzh
MEM_IN_GB=20
bsub -q gui -XF -I -R "select[mem>$((MEM_IN_GB*1024))] rusage[mem=$((MEM_IN_GB*1024))]" -M$((MEM_IN_GB*1024)) -o cytoscape_gui_.o -e cytoscape_gui_.e -J cytoscape_gui /hps/nobackup/iqbal/leandro/adrian/cytoscape/cytoscape_installation/cytoscape.sh
bjobs -w | grep cytoscape_gui | awk '{print $1}' | xargs bkill
rm cytoscape_gui_.*

scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_mobility.png .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_mobility.svg .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_AMR.png .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_AMR.svg .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_lin.svg .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_lin.png .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_linlarge.png .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_lin.png .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_PAE.png .
scp codon:/nfs/research/zi/acaza/projects/Murray_collection/Plasmid_identification/Illumina_sequencing/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/mash_network/Murray-iPlasmidDB_plasmids-network_PAE.svg .

#___________________________________#
