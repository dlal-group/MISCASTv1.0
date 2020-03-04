# MISCASTv1.0

MISCAST (MIssense variant to protein StruCture Analysis web SuiTe; http://miscast.broadinstitute.org/) is a web server to interactively visualize and analyze missense variants in protein sequence and structure space. 

The dataset in the current version of MISCAST spans 1,330 human genes, with 406,449 population variants from gnomAD (release 2.1.1, https://gnomad.broadinstitute.org/) and 54,137 variants from the ClinVar (February, 2019 release, pathogenic and likely-pathogenic variants, https://www.ncbi.nlm.nih.gov/clinvar/) and HGMD (Professional release 2018.4, disease mutations, http://www.hgmd.cf.ac.uk/ac/index.php) databases.

The variants are mapped from Ensembl transcript to 17,953 human protein 3D structures from Protein Data Bank (https://www.rcsb.org/).

Further, a comprehensive set of per annotations of protein structural, physicochemical, and functional features per residues were collected from multiple resources spanning DSSP, PDBsum, PhosphoSitePlus, PANTHER, UniProt. For the details about feature set annotation and mining, we refer the user to read the documentation page of MISCAST web server.

All annotation tracks (pathogenic and population missense variants and protein features) were subsequently mapped and displayed in protein sequence and structure in the web server. 

# MISCASTv1.0 GitHub Repo Content

1. Genes1330_crossrefrences.txt
-- The list of 1,330 genes included in the web server and associated refrences

2. Protein_class_to_gene.txt
-- The list of twenty-four protei functional classe. The genes are grouped into these classes based on the molecular function of the encoded proteins.

3. Gene wise annotation tracks (2 directories)
-- Annotation tracks (pathogenic and population missense variants and forty different protein features) for 1,330 genes

4. app-source-code
-- The R codes to implement the web server
