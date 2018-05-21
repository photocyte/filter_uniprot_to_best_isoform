##Drosophila melanogaster Uniprot reference proteome
##http://www.uniprot.org/proteomes/UP000000803
echo "Manually download the Drosophila melanogaster Uniprot reference proteome from Uniprot"
echo "Find it here: http://www.uniprot.org/proteomes/UP00000080"
echo ""
##Tribolium Uniprot reference proteome
##lynx http://www.uniprot.org/proteomes/UP000007266
echo "Manually download the Tribolium castaneum Uniprot reference proteome from Uniprot"
echo "Find it here: http://www.uniprot.org/proteomes/UP000007266"
echo ""
echo "I've manually downloaded both these FASTA files and added them to the library for convinience"
echo "Copyright etc. reserved by Uniprot"
echo ""

##Can run this to see the # of Swissprot proteins for a given species
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#zcat uniprot_sprot.fasta.gz | grep ">" | grep -oP "OS=.+OX=[0-9]+" | sort | uniq -c | sort -nr | grep -P "Drosophila melanogaster|Caenorhabditis elegans|Apis melifera|Bombyx mori|Anopheles gambiae|Tribolium castaneum"
echo ""
echo "Now downloading the NCBI proteome FASTA file for Drosophila melanogaster"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz
echo "Now downloading the NCBI proteome FASTA file for Tribolium castaneum"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_protein.faa.gz
echo "Now downloading the NCBI proteome FASTA file for Apis melifera"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/195/GCF_000002195.4_Amel_4.5/GCF_000002195.4_Amel_4.5_protein.faa.gz
echo "Now downloading the NCBI proteome FASTA file for Bombyx mori"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/625/GCF_000151625.1_ASM15162v1/GCF_000151625.1_ASM15162v1_protein.faa.gz
echo "Now downloading the NCBI proteome FASTA file for Caenorhabditis elegans"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.faa.gz
echo "Now downloading the NCBI proteome FASTA file for Anopheles gambiae"
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/575/GCF_000005575.2_AgamP3/GCF_000005575.2_AgamP3_protein.faa.gz

echo "Producing BLASTP database for testing Drosophila isoforms."
zcat GCF_000002335.3_Tcas5.2_protein.faa.gz\
 GCF_000002195.4_Amel_4.5_protein.faa.gz\
 GCF_000151625.1_ASM15162v1_protein.faa.gz\
 GCF_000002985.6_WBcel235_protein.faa.gz\
 GCF_000005575.2_AgamP3_protein.faa.gz > for_drosophila.faa
makeblastdb -dbtype prot -in for_drosophila.faa

echo "Producing BLASTP database for testing Tribolium isoforms."
zcat GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz\
 GCF_000002195.4_Amel_4.5_protein.faa.gz\
 GCF_000151625.1_ASM15162v1_protein.faa.gz\
 GCF_000002985.6_WBcel235_protein.faa.gz\
 GCF_000005575.2_AgamP3_protein.faa.gz > for_tribolium.faa
makeblastdb -dbtype prot -in for_tribolium.faa

echo "Now running the filtering script on Drosophila, but exiting before the BLAST step"
python3 ../filter_uniprot_to_best_isoform.py -nb <(zcat 2018_05_uniprot-proteome_Dmelanogaster_GCA_000001215.4.fasta.gz)
mv tmp.query.fa Drosophila.tmp.query.fa
echo "Running blastp in parallel"
/lab/solexa_weng/testtube/parallel_blast/parallel_blastp_map.sh Drosophila.tmp.query.fa for_drosophila.faa Dmel_blastp.results.cache.tsv
echo "Running the final filtering step... Find the results in isoform-filtered-*"
python3 ../filter_uniprot_to_best_isoform.py -cr Dmel_blastp.results.cache.tsv -sb <(zcat 2018_05_uniprot-proteome_Dmelanogaster_GCA_000001215.4.fasta.gz) > isoform-filtered_2018_05_uniprot-proteome_Dmelanogaster_GCA_000001215.4.fasta
echo "Done."
echo ""

echo "Now running the filtering script on Tribolium, but exiting before the BLAST step"
python3 ../filter_uniprot_to_best_isoform.py -nb <(zcat 2018_05_uniprot-proteome_Tcastaneum_GCA_000002335.3.fasta.gz)
mv tmp.query.fa Tribolium.tmp.query.fa
echo "Running blastp in parallel"
/lab/solexa_weng/testtube/parallel_blast/parallel_blastp_map.sh Tribolium.tmp.query.fa for_tribolium.faa Tcas_blastp.results.cache.tsv
echo "Running the final filtering step... Find the results in isoform-filtered-*"
python3 ../filter_uniprot_to_best_isoform.py -cr Tcas_blastp.results.cache.tsv -sb <(zcat 2018_05_uniprot-proteome_Tcastaneum_GCA_000002335.3.fasta.gz) > isoform-filtered_2018_05_uniprot-proteome_Tcastaneum_GCA_000002335.3.fasta
echo "Done."
echo ""

seqkit stat ./*.fasta ./*.fasta.gz
