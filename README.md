# Meta gbk generator

This script was done to work in junction with outputs coming from : MEGAHit, Prodigal and EggnogMapper, as : Assembly, Gene Prediction And Gene Annotation pipeline
If you have any other input/output file, you are encouraged to convert them in a suitable format.


## Getting Started

### Dependencies

This program uses the following python libraries:

* Biopython
* datetime

### Input Files:

The script requires 3 main inputs. 
1) An Fasta formated Assembly file with contigs from a Metagenome or MAG. This script has been specifically made for MEGAHit contigs names, although you can use the following format :
  * [Sample_name]\_[extra-description]\_[contig-id]
 
 #### Example:
 ```
 >[MAG_01]_[k141]_[3199166]
 GACTAGTTGACGATGTACGTAGCTAGACGTATCGATCTGA
 ATGCTGACGTAGTCATCTGTAGTCATGACTGATGCATGAA
 ...
 ```
2) A Fasta formated Gene File with all the CDS|ORF in Amionacidic / Protein format. This script has been specifically made for Prodigal 2.6.3 names, although you can use the following format:
  * [Sample_name]\_[extra-description]\_[contig-id]\_[gene-id] # [CDS-START] # [CDS-STOP] # [STRAND(1 OR -1)] #

#### Example:
 ```
 >MAG_01_k141_3199166_gene4 # 1695 # 2837 # -1 # 
 MDTFALDNLFENHAIDAVIHFAALKAVGESALQPLQYYQTNVHGSLCLLEAMAKAGVNNFV
 YSSSATVYGESNPSPYCETMALGSPSSPYGASKVMVERILQDKAKANSEFRAVSLRYFNP
 ...
 ```
3) A Tab-sepparated formated File with all the CDS|ORF annotations and or descriptions. This script has been specifically made for EggNogMapper v 2.x.x fields, although you can use the following format:
```
#query  seed_ortholog   evalue  score   eggNOG_OGs      max_annot_lvl   COG_category    Description     Preferred_name  GOs     EC      KEGG_ko KEGG_Pathway      KEGG_Module     KEGG_Reaction   KEGG_rclass     BRITE   KEGG_TC CAZy    BiGG_Reaction   PFAMs
```

  
  * #query Must be the same format as point 2).

#### Example:
 ```
MAG_01_k141_3199166_gene4   1535422.ND16A_1633      3.4e-150        437.0   COG0153@1|root,COG0153@2|Bacteria,1MVQD@1224|Proteobacteria,1RQ0C@1236|Gammaproteobacteria,2Q5XE@267889|Colwelliaceae   1236|Gammaproteobacteria        G       Galactokinase galactose-binding signature       galK    -       2.7.1.6 ko:K00849       ko00052,ko00520,ko01100,map00052,map00520,map01100       M00554,M00632   R01092  RC00002,RC00078 ko00000,ko00001,ko00002,ko01000,ko04147 -       -       -       GHMP_kinases_C,GHMP_kinases_N,GalKase_gal_bdg
 ```

## Example Scripts

* Once you have cloned/downloaded the repository, you can run 2 example scripts. You can also use them as a guide for your own run.

* For MAGs / BINs use:
```
bash bash_example_MAG.sh
```

* For entire Metagenomes:
```
bash_example_Metagenome.sh
```



## Usage

* Once you have all 3 core input files you can run this program adding the inputs and the listed extra info:

```
python meta_gbk_generator.py --help
options:
  -h, --help            show this help message and exit
  --assembly ASSEMBLY   Assembly file that contains the contigs
  --faa FAA             Fasta file containing all CDS in amionacidic format
  --annot ANNOT         Annotation file in EggNogMapper format
  --source SOURCE       Genome Source (For Metagenomes use : Metagenome, for Genome/MAGs use GTDB or NCBI names Example: unclassified_Gammaprotobacteria
                        | "Gammaprotobacteria Bacterium"
  --organism ORGANISM   Organism, could be the same as Source
  --taxonomy TAXONOMY   Organism taxonomy. For Metagenomes use "Bacteria", for Genomes use the Format: "Phylum:Class:Order:Family:Genus:Species"
  --isolate ISOLATE     MAG / Assembly isolate
  --isolation_source ISOLATION_SOURCE
                        Description of your sample extraction site: Might be soil, desert_soil, marine, environmental etc.
  --taxon_id TAXON_ID   Taxon ID of the Organism / Isolate
  --country COUNTRY     Country of Origin
  --type {m2m,ncbi}     Create a GBK for usage in NCBI or Metage2Metabo. Default [ncbi]
  --gbk GBK             Output file name (GBK)
```

## Filling the fields

If you are having trouble filling all the field of the script, here are some roundabout you can take.

* --source : This field should contain the Taxonomic Clasification of your Genomic Data plus an extra description if you wish, ranging from a simple "Soil Bacteria" or "Metagenomic Bacteria" to Specie level descriptors, such as "Isolated E. Coli 1E14"

* --organism : This field should contain only Taxonomic Clasifications accepted by NCBI/Taxonomy or tools such as GTDB, MetaPhlan, etc. Using the same example we can use either, "Bacteria", "Archaea" or specific taxas such as "Eschericia coli 1E14"
  
* --taxonomy : This field should contain the whole naming tree of your "--organism". For metagenomes use only "Bacteria", "Archaea", for isolates, use the maximum clasification available: "Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia"

* --taxon_id : This field contains the TAXONOMIC ID of your "--organism". It can be as simple as "2" for Bacteria, or "1255577" for our "Eschericia coli 1E14". All these IDs can be found on the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/taxonomy)

* --isolate : This field should contain and ID representative to you or your team. Usually names such as "Metagenome_01", "BIN-XX-water", are preferred. You can also use your [Sample_name] from point 1). 

## Authors

Contributors names and contact info

ex. Ricardo Palma  
ex. rpalmavejares@gmail.com

## License

This project is licensed under the GPL-3.0 License - see the LICENSE.md file for details
