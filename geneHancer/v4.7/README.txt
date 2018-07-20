README for GeneHancer dump

The GeneHancer dump consists of 4 data files and this README.

Data files:
enhancer_elite_ids.txt
enhancer_gene_scores.txt
enhancer_tfs.txt
enhancer_tissues.txt

File descriptions

enhancer_elite_ids.txt has 7 columns:

chr--chromosome where enhancer is located
enhancer_start--enhancer start position (hg38 coordinates)
enhancer_end--enhancer end position (hg38 coordinates)
cluster_id--enhancer cluster id, which can be used to connect enhancers between files
GHid--GeneHancer identifier
is_elite--elite status of enhancer (boolean)
regulatory_element_type--element type, such as enhancer or promoter

enhancer_gene_scores.txt has 9 columns:

cluster_id--enhancer cluster id, which can be used to connect enhancers between files
symbol--gene symbol for gene associated with enhancer
eqtl_score--GTEx eQTL p-value
erna_score--eRNA co-expression p-value
chic_score--CHiC score
expression_score--TF co-expression p-value
distance_score--score based on distance between enhancer midpoint and gene TSS
combined_score--enhancer-gene association score based on combination of above scores
is_elite--elite status of enhancer-gene association (boolean)

enhancer_tfs.txt has 3 columns:

enhancer_cluster_id--enhancer cluster id, which can be used to connect enhancers between files
TF--transcription factor
tissues--tissues in which TFBS was identified, delimited by a ';'

enhancer_tissues.txt has 5 columns:

cluster_id--enhancer cluster id, which can be used to connect enhancers between files
gh_id--GeneHancer identifier
source--source for the tissue annotation
tissue--tissue in which enhancer was identified
category--tissue biosample type (for ENCODE, FANTOM5)
