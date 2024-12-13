from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import os
from argparse import ArgumentParser
import time
import datetime 

date = datetime.datetime.today().strftime('%d/%m/%Y')

date_dict ={"01":"JAN","02":"FEB","03":"MAR","04":"APR","05":"MAY","06":"JUN","07":"JUL","08":"AUG","09":"SEP","10":"OCT","11":"NOV","12":"DEC"}

start_time = time.time()

parser = ArgumentParser(prog='Meta GBK generator',
                    description='Creates a GBK file from an assembly file, a cds file and an annotation file')


parser.usage = 'Use this program as follow:\n\npython meta_gbk_generator.py --asembly ASSEMBLY.fasta --faa ASSEMBLY_PROT_GENES.faa --annot ANNOTATION.tsv --source "SOURCE" --organism "ORGANISM" --taxonomy "TAXONOMY" --isolate "ISOLATE" --isolation_source "ISOLATION SOURCE" --taxon_id "TAXON ID" --country "COUNTRY" --type "TYPE" --gbk "GBK_OUTPUT"'

parser.add_argument('--assembly',help='Assembly file that contains the contigs',required=True)
parser.add_argument('--faa',help='Fasta file containing all CDS in amionacidic format',required=True)
parser.add_argument('--annot',help='Annotation file in EggNogMapper format',required=True)
parser.add_argument('--source',help='Genome Source (For Metagenomes use : Metagenome, for Genome/MAGs use GTDB or NCBI names\nExample: unclassified_Gammaprotobacteria | "Gammaprotobacteria Bacterium"',required=True)
parser.add_argument('--organism',help='Organism, could be the same as Source',required=True)
parser.add_argument('--taxonomy',help='Organism taxonomy. For Metagenomes use "Bacteria", for Genomes use the Format:\n"Phylum:Class:Order:Family:Genus:Species" ',required=True)
parser.add_argument('--isolate',help='MAG / Assembly isolate',required=True)
parser.add_argument('--isolation_source',help='Description of your sample extraction site: Might be soil, desert_soil, marine, environmental etc.',required=True)
parser.add_argument('--taxon_id',help='Taxon ID of the Organism / Isolate',required=True)
parser.add_argument('--country',help='Country of Origin',required=True)
parser.add_argument("--type",help='Create a GBK for usage in NCBI or Metage2Metabo. Default [ncbi]',choices=['m2m','ncbi'],required=True,default='ncbi')

parser.add_argument('--gbk',help='Output file name (GBK)',required=True)

args = parser.parse_args()

# Define all the args used in this script, I may add support to read all these from a tsv file in teh future


if (args):   
    fasta_assembly= args.assembly
    fasta_genes = args.faa
    annotation_file = args.annot
    genome_source = args.source.replace("_"," ")
    genome_organism = args.organism.replace("_"," ")
    taxonomy=args.taxonomy
    genome_isolate = args.isolate
    isolation_source = args.isolation_source
    taxon_id = args.taxon_id
    gbk_output_name = args.gbk
    country = args.country
    gbk_type= args.type

# Reads the CDS / Gene file containing all genes from an assembly. 
# The format of the fasta file must be identical to the ouput of Prodigal v2.6.3


fna_gene_data_dict = dict()   

print("Opening ",fasta_genes )

for record_fna in SeqIO.parse(fasta_genes,"fasta"):
    fna_data = []
    aux_for_each_fna = []
    aux_fna_features = (record_fna.description).split(" ")  
    fna_identifier = record_fna.id
#   print (fna_identifier)
    aux_name_parent = (record_fna.id).split("_")
    name_parent = "_".join(aux_name_parent[:-1])
    #print (name_parent)
    fna_start = aux_fna_features[2]
    fna_stop  = aux_fna_features[4]
    fna_strand= aux_fna_features[6]

    fna_data.append(fna_identifier)     # 2 FASTA CDN SEQ ID
    fna_data.append(fna_start)          # 3 START
    fna_data.append(fna_stop)           # 4 STOP
    fna_data.append(fna_strand)         # 5 STRAND (1 , -1)
    fna_data.append(record_fna.seq.rstrip("\*"))

    if(name_parent not in fna_gene_data_dict.keys()):
        fna_gene_data_dict.setdefault(name_parent,[])
        fna_gene_data_dict[name_parent].append(fna_data)
    else: 
        fna_gene_data_dict[name_parent].append(fna_data)
        
    #print(fna_data)

print ("FNA file CDS", len(fna_gene_data_dict))   
print("---Time for create gene cds dict--- %s seconds ---" % (time.time() - start_time))

    
# EMAPPER V 2.x.x FIELDS 

# 0 query   
# 1 seed_ortholog   
# 2 evalue  
# 3 score   
# 4 eggNOG_OGs  
# 5 max_annot_lvl   
# 6 COG_category    
# 7 Description 
# 8 Preferred_name  
# 9 GOs 
# 10 EC 
# 11 KEGG_ko    
# 12 KEGG_Pathway   
# 13 KEGG_Module    
# 14 KEGG_Reaction  
# 15 KEGG_rclass    
# 16 BRITE  
# 17 KEGG_TC    
# 18 CAZy   
# 19 BiGG_Reaction  
# 20 PFAMs
    
    

all_anotation_features = open(annotation_file,"r")
all_data_anotation_cluster_dict = dict()  


for lines in all_anotation_features:

    if ("#" not in lines):
        line_data = lines.split("\t")
        all_data_anotation = []
        aux_parent_anot = line_data[0].split("_")
        name_parent = "_".join(aux_parent_anot[:-1])
        all_data_anotation.append(line_data[0]) # 1 ID
        all_data_anotation.append(line_data[8]) # 2 GENE
        all_data_anotation.append(line_data[7]) # 2 ANNOT
        all_data_anotation.append(line_data[10])# 6 EC
        all_data_anotation.append(line_data[9]) # 7 GO:TERMS
        if(name_parent not in all_data_anotation_cluster_dict.keys()):
            all_data_anotation_cluster_dict.setdefault(name_parent,[])
            all_data_anotation_cluster_dict[name_parent].append(all_data_anotation)
        else:
            all_data_anotation_cluster_dict[name_parent].append(all_data_anotation)
        #print(all_data_anotation)

print ("ALL data from anotation", len(all_data_anotation_cluster_dict))
print("---Time for creating annotation dict--- %s seconds ---" % (time.time() - start_time))




aux_array_cluster = []

cds_position_by_parent=dict()


for key, values in fna_gene_data_dict.items():

    if(key in all_data_anotation_cluster_dict.keys()):
        if(len(fna_gene_data_dict[key])>=0):
            for entries in fna_gene_data_dict[key]:
                empty_annot = ["-","-","-","-"]
                sub_entries = entries
                sub_entries = sub_entries + empty_annot
                gene_index = fna_gene_data_dict[key].index(entries)
                fna_gene_data_dict[key][gene_index]=sub_entries
                for anot in all_data_anotation_cluster_dict[key]:
                    if(sub_entries[0] == anot[0]):
                        sub_entries[5] = anot[1]
                        sub_entries[6] = anot[2]
                        sub_entries[7] = anot[3]
                        sub_entries[8] = anot[4]
                        fna_gene_data_dict[key][gene_index]=sub_entries

    else:
        for entries in fna_gene_data_dict[key]:
            empty_annot = ["-","-","-","-"]
            sub_entries = entries
            sub_entries = sub_entries + empty_annot
            gene_index = fna_gene_data_dict[key].index(entries)
            fna_gene_data_dict[key][gene_index]=sub_entries
            #print(key ," >>> ",fna_gene_data_dict[key])
            #print("No annotation")
    #break

print ("All data combined", len(all_data_anotation_cluster_dict))
print("---Time for creating merged info --- %s seconds ---" % (time.time() - start_time))
genome_data = []


output_file = open (gbk_output_name+".gbk",'w')


for record in SeqIO.parse(fasta_assembly,"fasta"):

    
    aux_for_each = []
    aux_identifier = 'unknown'
    aux_for_each.append(record.id)
    aux_for_each.append(record.seq)
    genome_data.append(aux_for_each)
    
    sequence_string = str(record.seq)
   
    taxonomy_tree=args.taxonomy
    taxonomy_tree=taxonomy_tree.split(":")
    
    aux_seq_all = Seq(sequence_string)

    name_record = (record.id).split("_")

    all_annotations=dict()  
    all_annotations["molecule_type"]="DNA"
    all_annotations["source"]=genome_source.replace("_"," ")
    all_annotations["organism"]=genome_organism.replace("_"," ")
    all_annotations["topology"]="linear"
    all_annotations["data_file_division"]="PLN"
    
    date_time_parsed=date.split("/")
    all_annotations["date"]=date_time_parsed[0]+"-"+date_dict[date_time_parsed[1]]+"-"+date_time_parsed[2]
    all_annotations["taxonomy"]=taxonomy_tree
    
    new_custom_name = "_".join(name_record[:-2])
    record_gbk = SeqRecord(aux_seq_all,
                id=new_custom_name,
                name=new_custom_name,
                annotations=all_annotations)
   
    qualifier_source = dict()
    qualifier_source['organism']=args.organism.replace("_"," ") 
    qualifier_source['db_xref']="taxon:"+taxon_id
    qualifier_source['mol_type']="genome DNA" 
    qualifier_source['isolate']=genome_isolate 
    qualifier_source['isolation_source']=isolation_source.replace("_"," ") 
    qualifier_source['country']=country 
    feature_source = SeqFeature(FeatureLocation(start=int(0),end=int(len(record_gbk.seq)-1), strand=int(1)), type="source", qualifiers=qualifier_source)
    record_gbk.features.append(feature_source)
    
    if(record.id in fna_gene_data_dict.keys()):    

        for each_fna in fna_gene_data_dict[record.id]:
      
            qualifier_dict = dict()
            name_gbk_locus = (record.id).split("_")
          
            new_custom_locus = str(record.id)
            qualifier_dict['locus_tag'] = each_fna[0]

            if(each_fna[5]!="-" and each_fna[5]!=""):

                qualifier_dict['gene'] = str(each_fna[5]).upper()

            if (each_fna[8]!="" and each_fna[8]!="-"): 
                
                if(gbk_type=="m2m"):
                    qualifier_dict['go_component']=[]
                    if("," in each_fna[8]):
                        all_go_terms=each_fna[8].split(",")
                        qualifier_dict['go_component']=list(all_go_terms)
                    else:
                        qualifier_dict['go_component'] = each_fna[8]

                elif(gbk_type=="ncbi"):
                    qualifier_dict['db_xref']=[]
                    if("," in each_fna[8]):
                        all_go_terms=each_fna[8].split(",")
                        qualifier_dict['db_xref']=list(all_go_terms)
                    else:
                        qualifier_dict['db_xref'] = each_fna[8]


            if (each_fna[5]!="" and each_fna[5]!="-"):
                
                qualifier_dict['product'] = each_fna[5].capitalize()

            if (each_fna[6]!="" and each_fna[6]!="-"):
                
                qualifier_dict['function'] = each_fna[6]
            
            if (each_fna[7]!="" and each_fna[7]!="-"): 
                if("," in each_fna[7]): 
                    all_ec_terms=each_fna[7].split(",") 
                    qualifier_dict['EC_number']=list(all_ec_terms)
                else:
                    qualifier_dict["EC_number"] = each_fna[7]
            
            qualifier_dict['translation'] = each_fna[4]
            feature = SeqFeature(FeatureLocation(start=int(each_fna[1])-1,end=int(each_fna[2]), strand=int(each_fna[3])), type="CDS", qualifiers=qualifier_dict) ##

            qualifier_2_dict = {}
            if(each_fna[5]!="-" and each_fna[5]!=""):
                qualifier_2_dict['gene'] = str(each_fna[5]).upper()
            qualifier_2_dict['locus_tag'] = each_fna[0]

            feature_2 = SeqFeature(FeatureLocation(start=int(each_fna[1])-1,end=int(each_fna[2]), strand=int(each_fna[3])), type="gene", qualifiers=qualifier_2_dict)
            record_gbk.features.append(feature_2)
            record_gbk.features.append(feature)

    else:
        pass
        # In case we Add contigs without genes / features non automatically.


    SeqIO.write(record_gbk,output_file,'genbank')                                          

print ("###########")
print("Completed GBK file in --- %s seconds ---" % (time.time() - start_time))
