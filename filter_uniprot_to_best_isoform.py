#! /usr/bin/python 
import Bio
import Bio.SeqIO
import argparse
import re
import sys
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument('fasta',type=str,help='The path to the fasta file to clean')
parser.add_argument('-nb',default=False,action='store_true',help='Exit before blastp-ing but after writing the query fasta file')
parser.add_argument('-sb',default=False,action='store_true',help='Skip the blast step')
parser.add_argument('-cr',type=str,required=False,default="blastp.results.cache.tsv",help='Path to the cached BLASTP results file (output fmt 6)')
args = parser.parse_args()

gene_isoform_clusters = dict()
record_dict = dict()

uniprot_description_match = re.compile("(tr|sp)\|([0-9a-zA-Z_]+)\|([0-9a-zA-Z_]+) (.+) OS=([0-9a-zA-Z_]+) ([0-9a-zA-Z_]+) OX=([0-9]+) GN=(.+) PE=([1-5]) SV=([1-9])")

sys.stderr.write("Loading FASTA records into memory...")
handle = open(args.fasta)
for record in Bio.SeqIO.parse(handle, "fasta"):
        result = uniprot_description_match.match(record.description)
        if result == None:
            sys.stderr.write("This shouldn't happen.")
            sys.stderr.write(record.description)
            exit()
        record.db = result.group(1)
        record.uni_id = result.group(2)
        record.entry_name = result.group(3)
        record.protein_name = result.group(4)
        record.genus = result.group(5)
        record.species = result.group(6)
        record.taxonid = result.group(7)
        record.geneid = result.group(8)
        record.pe = int(result.group(9))
        record.sv = int(result.group(10))
        record_dict[record.uni_id] = record
        record.similarity_score = 0.0
        if record.geneid not in gene_isoform_clusters.keys():
            gene_isoform_clusters[record.geneid] = [record.uni_id]
        else:
            gene_isoform_clusters[record.geneid].append(record.uni_id)
handle.close()

sys.stderr.write("Loaded "+str(len(record_dict))+" FASTA records from:"+args.fasta+"\n")
sys.stderr.write("Produced "+str(len(gene_isoform_clusters))+" gene clusters\n")
#for key in record_dict.keys():
#    r = record_dict[key] 

sys.stderr.write("Evalutating based on swissprot/reviewed anotation...\n")
IDs_to_pass = []
special_cases = []
further_evaluation = []
for gene in gene_isoform_clusters.keys():
    #sys.stderr.write(gene,gene_isoform_clusters[gene])
    if len(gene_isoform_clusters[gene]) == 1:
        p = gene_isoform_clusters[gene][0]
        IDs_to_pass.append(p)
    elif len(gene_isoform_clusters[gene]) > 1:
        swissprot_ids = []
        for p in gene_isoform_clusters[gene]:
            if record_dict[p].db == "sp":
                swissprot_ids.append(p)
        if len(swissprot_ids) == 0:
            further_evaluation.append(gene)
        elif len(swissprot_ids) == 1:
            IDs_to_pass.append(swissprot_ids[0])
        elif len(swissprot_ids) > 1:
            special_cases.append(gene)
            #sys.stderr.write("Special case",len(swissprot_ids))
            #sys.stderr.write(gene_isoform_clusters[gene])

sys.stderr.write("IDs to pass:"+str(len(IDs_to_pass))+"\n")
sys.stderr.write("Gene-isoform clusters that are special cases:"+str(len(special_cases))+"\n")
sys.stderr.write("Gene-isoform clusters that are special cases:"+str(special_cases)+"\n")
sys.stderr.write("Gene-isoform clusters that require further evaluation:"+str(len(further_evaluation))+"\n")

#sys.stderr.write(further_evaluation)
###Filter down
###Here are the rules.

###Warning: the Gene Name used in the fasta description may not be unique. E.g. fly "Rpb4" which is
##two genes that share the same promoter, and terminator, but generate different polypeptides via
##alternative splicing which are considered different enough that people consider them 2 different genes

##Ideally, one could lookup the FlyBase or Tribolium gene name from the Uniprot ID mapping service, but that seems too complicated
##Also really only Drosophila is the one with this problem, as only Drosophila has had enough "fine-toothed" annotation fixing
##to produce such weird edge cases

##Requires further evaluation, means none of them are swissprot, so they could all be equally bad.
##Lets use a couple different metrics:
##The protein with the best evidence (Uniprot evidence score)
##If that is a tie, amongst the tied members
##The protein with the best BLASTP score amongst the different isoforms.

##Protein evidence scores
sys.stderr.write("Evaluating based on protein evidence scores...\n")
sys.stderr.write("Only keeping those proteins with the best (lowest) scores in their cluster.\n")
overall_pass = []
for g in further_evaluation:
    can_pass = []
    best_score = 6 ##To start, set it worse than the worst of the scale
    for p in gene_isoform_clusters[g]:
        if record_dict[p].pe == best_score:
            can_pass.append(p)
        elif record_dict[p].pe < best_score:
            can_pass = [p]
            best_score = record_dict[p].pe
        elif record_dict[p].pe > best_score:
            #sys.stderr.write(str(record_dict[p].pe)+","+str(best_score)+"\n")
            pass
    overall_pass+=can_pass

i=0
for g in further_evaluation:
    for p in list(gene_isoform_clusters[g]):
        if p not in overall_pass:
            gene_isoform_clusters[g].remove(p)
            i+=1
        else:
            pass
sys.stderr.write("Removed "+str(i)+" proteins\n")

class Blast_Cache:
    def __init__(self,path):
        try:
            handle = open(path,"rU")
        except:
            self.cache = dict()
            return None
        self.cache = dict()
        for line in handle.readlines():
            splitline = line.rstrip().split("\t")
            if len(splitline) != 12:
                ##Corrupted line, skip
                continue
            uni_id = splitline[0].split("|")[1]
            bitscore = float(splitline[11])
            if uni_id not in self.cache.keys():
                self.cache[uni_id] = bitscore
            else:
                self.cache[uni_id] += bitscore
        handle.close()
    def check(self,uni_ID):
        if uni_ID in self.cache.keys():
            return True
        if uni_ID not in self.cache.keys():
            return False   

    def get_score(self,uni_ID):
        return self.cache[uni_ID]         

cache = Blast_Cache(args.cr)

##Write the query file for blasting
query_fasta_path="tmp.query.fa"
qf_handle = open(query_fasta_path,"wt")
j=0
k=0
for g in further_evaluation:
    for p in gene_isoform_clusters[g]:
        if cache.check(p) == False:
            qf_handle.write(record_dict[p].format("fasta"))
        elif cache.check(p) == True:
            record_dict[p].similarity_score = cache.get_score(p)
            k+=1
        j+=1
qf_handle.close()
sys.stderr.write("Had "+str(j)+" total proteins to lookup, and "+str(k)+" cache hits.\n")
sys.stderr.write("Wrote "+str(j-k)+" to tmp.query.fa.\n")
if args.nb == True:
    exit()

if args.sb == False:
    ##Or, just blastp it to a related proteome and pick the isoform with the best scores
    handle = open(args.cr,"at")
    sys.stderr.write("Now blastp-ing all the different isoforms...\n")
    sys.stderr.write(str(j-k)+" different isoforms to evaluate\n")
    sys.stderr.write("This could take awhile... Recommend you use the -nb/-sb parameter and parallelize the BLASTP search on your own hardware\n")
    target_fasta_path="./GCF_000002335.3_Tcas5.2_genomic.fna"
    cmd_string = "/lab/solexa_weng/testtube/ncbi-blast-2.6.0+/bin/blastp -query "+query_fasta_path+" -subject "+target_fasta_path+" -culling_limit 1 -evalue 1e-5 -outfmt 6 -gapextend 1" 
    blast_process = subprocess.Popen(cmd_string.split(" "),stdout=subprocess.PIPE)
    i=0
    for line in iter(blast_process.stdout.readline,''):
        handle.write(line)
        splitline = line.rstrip().split("\t")
        uni_id = splitline[0].split("|")[1]
        bitscore = float(splitline[11])
        record_dict[uni_id].similarity_score += bitscore
        i +=1
        if i % 100 == 0:
            approx_percent_complete = round(((i*1.0)/(j*5.0))*100,1)
            sys.stderr.write(str(i)+" done. ~"+str(approx_percent_complete)+"%\n")
    handle.close()

##Evaluate which is better:

i=0
for g in list(further_evaluation):
    max_proteins = [] 
    max_score = 0.0
    for p in gene_isoform_clusters[g]:
       if record_dict[p].similarity_score > max_score:
           max_proteins= [p]
           max_score = record_dict[p].similarity_score
       elif record_dict[p].similarity_score == max_score:
           max_proteins.append(p)
       elif record_dict[p].similarity_score < max_score:
           pass
    if len(max_proteins) == 1:
        IDs_to_pass.append(max_proteins[0])
        further_evaluation.remove(g)
        i+=1
sys.stderr.write("Resolved "+str(i)+" gene-isoform clusters.\n")
sys.stderr.write("Gene-isoform clusters that are still unclear:"+str(len(further_evaluation))+"\n")

for p in IDs_to_pass:
    sys.stdout.write(record_dict[p].format("fasta"))

for g in further_evaluation:
    for p in gene_isoform_clusters[g]:
        sys.stdout.write(record_dict[p].format("fasta"))

for g in special_cases:
    for p in gene_isoform_clusters[g]:
        sys.stdout.write(record_dict[p].format("fasta"))
