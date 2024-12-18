tar -xvzf example_Metagenome.tar.gz
python meta_gbk_generator.py --assembly example_dataset/Metagenome_assembly.fasta --faa example_dataset/Metagenome_genes.faa --annot example_dataset/Metagenome_annotation.tsv --source Bacteria --organism Bacteria --taxonomy Bacteria --isolate Meta-01 --isolation_source "Marine Megagenome" --taxon_id 2 --country Chile --gbk Meta-01_genbank --type m2m
