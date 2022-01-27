
s1_mapped1.txt: Domain 1 of Simulation 1, branching tree (originally from Liu et al 2019, MMD-MA paper)
s1_mapped2.txt: Domain 2 of Simulation 1, branching tree (originally from Liu et al 2019, MMD-MA paper)
s1_label1.txt: Group labels for Domain 1 of Simulation 1
s1_label2.txt: Group labels for Domain 2 of Simulation 1

s2_mapped1.txt: Domain 1 of Simulation 2, Swiss roll (originally from Liu et al 2019, MMD-MA paper)
s2_mapped2.txt: Domain 2 of Simulation 2, Swiss roll (originally from Liu et al 2019, MMD-MA paper)
s2_label1.txt: Group labels for Domain 1 of Simulation 2
s2_label2.txt: Group labels for Domain 2 of Simulation 2

s3_mapped1.txt: Domain 1 of Simulation 3, circular frustum (originally from Liu et al 2019, MMD-MA paper)
s3_mapped2.txt: Domain 2 of Simulation 3, circular frustum (originally from Liu et al 2019, MMD-MA paper)
s3_label1.txt: Group labels for Domain 1 of Simulation 3
s3_label2.txt: Group labels for Domain 2 of Simulation 3

scGEM_expression.csv:	Gene expression (scRNA-seq) domain of the scGEM dataset
scGEM_typeExpression.txt:	Cell-type labels for the gene expression (scRNA-seq) domain of the scGEM dataset
scGEM_methylation.csv:	DNA methylation (scMethyl-seq) domain of the scGEM dataset
scGEM_typeMethylation.txt:	Cell-type labels for the DNA methylation (scMethyl-seq) domain of the scGEM dataset
Label to cell-type mapping:	1: BJ, 2:d8, 3: d16T+, 4: d24T+, 5:iPSc


scatac_feat.npy:  Chromatin accessibility (snATAC-seq) domain of the SNARE-seq dataset (in numpy array format)
SNAREseq_atac_types.txt:  Cell-type labels for the gene expression domain of the SNARE-seq dataset
scrna_feat.npy:	  Gene expression (scRNA-seq) domain of the SNARE-seq dataset (in numpy array format)
SNAREseq_rna_types.txt:	  Cell-type labels for the chromatin accessibility domain of the SNARE-seq dataset (same as the gene expression domain because gene expression domain proved more reliable for cell-type annotation)
SNAREseq_metadata.txt: 	Some metadata on the SNARE-seq datasets, such as the label number -- cell-type mapping and the barcode names for cells (in order)
