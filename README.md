# Meta_gbk_generator

This script was done to work in junction with outputs coming from : *MEGAHit, Prodigal and EggnogMapper, as : Assembly, Gene Prediction And Gene Annotation pipeline
If you have any other input/output file, you are encouraged to convert them in a suitable format.


## Getting Started

### Dependencies

This program uses the next python libraries:

* Biopython
* argparse
* datetime

## How to Use:

The script requires 3 main inputs. 
* 1) An Fasta formated Assembly file with contigs from a Metagenome or MAG. This script has been specifically made for MEGAHit contigs names, although you can use the following format :
  * [Sample_name]_[extra-description]_[contig-id]
 ```
 >[MAG_01]_[k141]_[3199166]
 GACTAGTTGACGATGTACGTAGCTAGACGTATCGATCTGA
 ATGCTGACGTAGTCATCTGTAGTCATGACTGATGCATGAA
 ...
 ```
* 2) A Fasta formated Gene File with all the CDS|ORF in Amionacidic / Protein format. This script has been specifically made for Prodigal 2.6.3 names, although you can use the following format:
  * [Sample_name]_[extra-description]_[contig-id]_[gene-id] # [CDS-START] # [CDS-STOP] # [STRAND(1 OR -1)] #
 ```
 >MAG_01_k141_3199166_gene4 # 1695 # 2837 # -1 # 
 MDTFALDNLFENHAIDAVIHFAALKAVGESALQPLQYYQTNVHGSLCLLEAMAKAGVNNFV
 YSSSATVYGESNPSPYCETMALGSPSSPYGASKVMVERILQDKAKANSEFRAVSLRYFNP
 ...
 ```

## Authors

Contributors names and contact info

ex. Ricardo Palma  
ex. rpalmavejares@gmail.com

## License

This project is licensed under the GNU V3 License - see the LICENSE.md file for details
