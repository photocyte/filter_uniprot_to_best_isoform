This is a Python3 script which is designed to take a [Uniprot reference proteome](https://www.uniprot.org/help/reference_proteome) in FASTA format (which are non-redundant in peptide sequence, but may have multiple isoforms per gene), and attempt to filter down the protemome to a single "canonical isoform" per gene.

If your organism is already well-manually reviewed on Uniprot (e.g. the protein is included in the SwissProt database), it has already undergone such reduction to a single canonical isoform per gene (barring some [exceptions](https://www.uniprot.org/help/canonical_and_isoforms) where the definition of "gene" gets a bit tricky). 

This "well-manually reviewed" criteria really only holds true for human, budding-yeast, *Arabidopsis thaliana*, and *E. coli*.  Everything else (e.g. *D. melanogaster*) is "barely" manually reviewed, and this script can help reduce redundancy.

The script filters on three criteria:

 1. Is there a single manually-reviewed SwissProt protein per gene-isoform cluster? If so, pass only that protein, and do not pass the rest.  If there is >1, then it is reported as a "special-case", and all proteins in that gene-cluster are passed.  If it is 0, then that gene-cluster goes into the next steps for further evaluation/
 2. Filter based on the [Uniprot protein existence annotation](https://www.uniprot.org/help/protein_existence). If there is 1 isoform that has a better protein existence than all the others in the gene-isoform cluster, then it is passed, and the rest are not passed. If there is >1 proteins that are equal in terms of their existence annotation, then put those proteins into the 3rd step.
 3. Pass the isoform which has the greatest bitscore, when BLASTPed to the proteomes (non-redundant, non filtered) of closely related species. If the isoform BLASTP bitscores are tied, pass all isoforms.

See the [examples](examples) folder for more details.

Representative results are below:

* *Drosophila melanogaster* **(unfiltered)**
  * 21,957 proteins (for 14,213 protein coding genes)
  * BUSCO score: C:99.8%[S:67.0%,**D:32.8%**],F:0.2%,M:0.0%,n:2442
  
* *Drosophila melanogaster* **(script filtered)**
  * 17,684 proteins (for 14,213 protein coding genes)
  * BUSCO score: C:99.6%[S:84.4%,**D:15.2%**],F:0.3%,M:0.1%,n:2442
  * Checksum: 889f251c848b103e1d892c53dc79ae26

* *Tribolium castaneum* **(unfiltered)**
  * 18,505 proteins (for 16,578 protein coding genes)
  * BUSCO score: C:98.1%[S:89.4%,**D:8.7%**],F:1.5%,M:0.4%,n:2442

* *Tribolium castaneum* **(script filtered)**
  * 18,436 proteins (for 16,578 protein coding genes)
  * BUSCO score: C:98.1%[S:89.4%,**D:8.7%**],F:1.5%,M:0.4%,n:2442
  * Checksum: 05d499b2f68714cef040267d45c532c2

 
