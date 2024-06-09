# Analysis of PAE plasmids and their modern relatives

This document describes the software and commands used to analyse the plasmid sequences reported in our study.


## Murray plasmids


### Retrieving Murray sequencing reads

__Software__: enaBrowserTools v1.6
___Command___:
enaGroupGet -f fastq -d Illumina/ ${acc}


### Genome Assembly

# Software: unicycler v0.4.7
# Command:
unicycler -1 ${d}/*_1.fastq.gz -2 ${d}/*_2.fastq.gz -o /home/ubuntu/Plasmids/Murray/Murray_collection/UNICYCLER/${d}


#===== Assembly metrics =====#

# Software: quast v4.5
# Command:
quast -o QUAST/ UNICYCLER/*/assembly.fasta &


#===== Taxonomic classification =====#

# Software: kraken2 v2.1.2 & bracken v2.5
# Commands:
# kraken
kraken2 --use-names --threads 8 --db k2_standard_8gb_20210517 --fastq-input --report M266_test.report-reads.txt --gzip-compressed --paired ../Murray_genomes/Illumina/Reads/M266-ERS222761/ERR316957_1.fastq.gz  ../Murray_genomes/Illumina/Reads/M266-ERS222761/ERR316957_2.fastq.gz --output ../Murray_genomes/Illumina/sandbox/KrakenBracken_Reads/M266_test-reads.kraken
# bracken
bracken -d k2_standard_8gb_20210517 -i M266_test.report-reads.txt -o M266_reads_bracken -l S


#===== Chromosomal contigs depletion =====#
# Software: plasmidverify.py
# Command:
python2.7 ~/plasmidVerify/plasmidverify.py -f ${d}/assembly.fasta -o ~/Plasmids/Murray/Murray_collection/PVERIFY/${d} --hmm ~/Pfam-HMM/Pfam-A.hmm.gz -t 8


#===== Contigs circularisation =====#
# Software: circlator v1.5.5
# Command:
circlator all --threads 8 --b2r_min_read_length 50 --merge_min_length 100 --merge_min_length_merge 200 --assemble_spades_k 107,97,87,77,67,57,51 --b2r_discard_unmapped --clean_min_contig_length 100 /data/${d}-noncir_uni.fasta /data/${d}_reads.fastq /data/circlator_out


#===== Contigs mapping =====#
# Software: minimap2 v2.17-r941
# Command:
minimap2 -x asm5 -t 8 -c ~/Plasmids/Murray/iPlasmidDB/v1/sequences/FNA/individual_files-nr/DEDUPE/dedupe_out.fasta ${d}/*.fna -o ~/Plasmids/Murray/Murray_collection/Predicted_plasmids/PVERIFY/MINIMAP2/Post-circlator_contigs-iPlasmidDB_all/${d}.paf


#===== Plasmids typing =====#
# Software: abricate v1.0.1, mob_suite v3.0.1,
# Command:
# abricate - plasmidfinder
abricate -db plasmidfinder --minid 80 --mincov 60 --quiet --nopath --threads 8 Post-circlator_contigs/*/*fna > PLASMIDFINDER/Post-circlator_contigs-plasmidfinder_80id_60cov.txt &
# mob_recon
mob_recon --infile /mnt/${d}/${d}_post-*.fna --outdir /mnt/${d}/mob_recon_output -t -u -c -n 8
# mob_typer
mob_typer --infile ${f}.fna --out_file ~/Plasmids/Murray/Murray_collection/Identified_plasmids/Post-PVERIFY/comparison/vs_iPlasmidDBv1.1/annotation/MOB_TYPER/tmp_reports/${f}.txt --num_threads 8 --multi


#===== Sequence deduplication =====#
# Software: Dedupe (BBTools package v38.90)
# Command:
bash ~/bbmap/dedupe.sh in=iPlasmidDB-core_all.fna threads=8 sort=id out=dedupe_out.fasta ac=f outd=duplicates.fasta &


#===== Plasmids NCBI metadata =====#
# Software: plasmid_querier
# Command:
Ask LEANDRO


#===== AMR/Virulence databases clustering =====#
# Software: CD-HIT v4.8.1
# Command:
# Concatenate databases
cat card/sequences ncbi/sequences argannot/sequences resfinder/sequences > cd-hit_res/sequences-all
# cd-hit-est
cd-hit-est -i sequences-all -o sequences -c 1.00 -n 10 -T 8 -s 0.9 -sc 1


#===== AMR/Virulence genes identification =====#
# Software: abricate v1.0.1
# Command:
abricate -db cd-hit_res --minid 80 --mincov 80 --quiet --nopath --threads 8 --fofn sequences/input_fna_files.txt > annotation/ABRICATE/amr/Murray-iPlasmidDB_plasmids.80cov-cd-hit_res.tsv &


#===== Plasmid annotation =====#
# Software: prokka v1.14.6
# Command:
prokka --outdir annotation/PROKKA/${bname} --prefix ${bname} --locustag ${bname} --plasmid ${bname} --kingdom Bacteria --gcode 11 --cpus 8 sequences/FNA/individual_files/${l}


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
