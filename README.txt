This is the documentation for the annovar/annotation pipeline.

11/07/14

8/30/12
# Added Cosmic sift polyphen NHLBI database annotations 
# Upgraded 1000g version to Apr 2012 

#####################################
How to run the prorgam:

/data/khanlab/apps/annovar/annovar.pl -infile <FILE.txt> -cs 4


-cs swith is added to incorporate the OPTION 1 from below. 
-cs default value is 5 meands your input file have Start and End columns.


### Output file will be FILE.txt.annotations ###


for this version of annovar the input file can have N number of clumns but there should be folluwing columns at the given column number.

#### OPTION 1 #### 
Column  Value 
1 	Chr
2 	Pos
3	Ref  
4	Alt
### OPTION 2 ###
Column  Value
1       Chr
2       Start
3       End
4       Ref
5 	Alt

Chr	Start	End	Ref	Obs	All the other columns 
chr1	161003087	161003087	C	T	a SNP in Illumina SNP arrays
chr1	13133880	13133881	TC	-	a 2-bp deletion
chr1	11326183	11326183	-	AT	a 2-bp insertion
chr1	105293754	105293754	A	ATAAA	a block substitution
chr16	49303427	49303427	C	T	a non-synonymous SNP 
chr16	49314041	49314041	G	C	a non-synonymous SNP 
chr16	49321279	49321279	-	C	a frameshift SNP 
chr13	19661686	19661686	G	-	a frameshift mutation 
chr13	19695176	20003944	0	-	a 342kb deletion

Output of the Module:

The output file name will be file.txt_annovar.
the file will be in the directory of input file.

Output file will have all the columns present in the input file and at the end of the RHS there will be 20 columns added.

Description of the annovar.pl output.


################
gene.Region
        Region in which the position contained on hg18 build (refSeq).
Gene
        Gene name.
        in cases of splice, intergenic ncRNA regions it contain the 2 gene names.
gene.Variant
	this is refSeq Gene Annotation and the column values can be 
	frameshift insertion		an insertion of one or more nucleotides that cause frameshift changes in protein coding sequence
	frameshift deletion		a deletion of one or more nucleotides that cause frameshift changes in protein coding sequence
	frameshift block substitution	a block substitution of one or more nucleotides that cause frameshift changes in protein coding sequence
	stopgain			a nonsynonymous SNV, frameshift insertion/deletion, nonframeshift insertion/deletion or block substitution that lead to the immediate creation of stop codon at the variant site. 
					For frameshift mutations, the creation of stop codon downstream of the variant will not be counted as "stopgain"!
	stoploss			a nonsynonymous SNV, frameshift insertion/deletion, nonframeshift insertion/deletion or block substitution that lead to the immediate elimination of stop codon at the variant site
	nonframeshift insertion		an insertion of 3 or multiples of 3 nucleotides that do not cause frameshift changes in protein coding sequence
	nonframeshift deletion		a deletion of 3 or mutliples of 3 nucleotides that do not cause frameshift changes in protein coding sequence
	nonframeshift block substitution a block substitution of one or more nucleotides that do not cause frameshift changes in protein coding sequence
	nonsynonymous 			a single nucleotide change that cause an amino acid change
	synonymous 			a single nucleotide change that cause an amino acid change
	unknown				unknown function (due to various errors in the gene structure definition in the database file)

gene.AminoAcidChange
	This is the amino acid change due the the change in the base. (refseq).
#################
CytoBand.annovar
        Cytogenic band.
        When a variant spans multiple bands, they will be connected by a dash (for example, 1q21.1-q23.3)
#################
dbsnp138
	dbsnp138 annotation. if present then rs id will be present else 0
#################
1000g08.2014
1000g08.2014_eur
1000g08.2014_afr
1000g08.2014_amr
1000g08.2014_eas
1000g08.2014_sas
	Frequencies from 1000 genome new release (Please mail me if you find a newer release)
	http://www.1000genomes.org/category/frequently-asked-questions/population
#################
NHLBI6500_ALL	
NHLBI6500_EA	
NHLBI6500_AA
	Frequencies from Exome Sequencing Project new release (Please mail me if you find a newer release)
	http://evs.gs.washington.edu/EVS/
#################
cg69	
	Frequency in 69 samples sequenced at Complete Genomics 
	http://www.completegenomics.com/public-data/69-Genomes/
#################	
nci60
	Frequency in 60 Cell lines sequenced at NCI
	http://discover.nci.nih.gov/cellminer/analysis.do
#################
caddRAW	
caddPHRED
	Combined Annotation Dependent Depletion (CADD) Score
#################
clinvar_20140702	
cosmic70
	Membership to Clinvar and cosmic databases 
	http://www.ncbi.nlm.nih.gov/clinvar/
	http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/
################
SIFT Prediction	
	Prediction from SIFT (DAMAGING, TOLERATED, Not scored, Damaging due to stop, N/A, DAMAGING *Warning! Low confidence.)
SIFT Score
	0 = DAMAGING, 1 = TOLERATED, N/A = Damaging due to stop, Not scored, N/A 
	http://sift.jcvi.org/
################
PPH2 Prediction	
	(benign, possibly damaging, probably damaging)
PPH2 Class
	(neutral, deleterious)
PPH2 Probability
	Score
	http://genetics.bwh.harvard.edu/pph2/
################
Acc_No.hgmd2014.3	
Gene.hgmd2014.3	
GeneName.hgmd2014.3	
Disease.hgmd2014.3	
Category.hgmd2014.3	
Reference_PMID.hgmd2014.3
	if HGMD member, accession number, gene name, defination, disease related to in HGMD, Category as defined by HGMD, and Reference Pubmed ID.
	http://www.hgmd.cf.ac.uk/ac/hahaha.php
	https://portal.biobase-international.com/hgmd/pro/global.php#other
	Access is provided at https://portal.biobase-international.com/cgi-bin/portal/login.cgi>. Individual registration and login required. Enter the code 1881-6975-97565225 in the license field during the account registration process.
################
CancerGeneCensus
	Yes if gene is a member of cancer gene census 
	http://cancer.sanger.ac.uk/cancergenome/projects/census/
################
ACMG_reportableGenes
	Yes if gene is a member of ACMG Reportable gene list.
	https://www.acmg.net/docs/IF_Statement_Final_7.24.13.pdf
################
MATCH_v1_08_2014Purpose	
MATCH_v1_08_2014aMOI	
MATCH_v1_08_2014Drug	
MATCH_v1_08_2014Level	
MATCH_v1_08_2014References	
MATCH_2 Purpose
	Molecular Analysis for Therapy Choice Program
	http://www.cancer.gov/clinicaltrials/noteworthy-trials/match#match
################
MCG.Diagnosis	
MCG.Targeted Therapy
MCG.Other Implication
	http://www.mycancergenome.org/
################
Total_ICGC
	https://icgc.org/ (~6500 Patients, 50 different tumor types)
Count_TCGA_ALL 
	(~16K Samples Downloded from http://www.cbioportal.org/public-portal/index.do)
2008.Parsons.Glio.Multiforme
	http://www.sciencemag.org/content/321/5897/1807
2010Barretina.SoftTissueSarcoma
	http://www.ncbi.nlm.nih.gov/pubmed/20601955
2011Heravi.DICER1
	http://www.ncbi.nlm.nih.gov/pubmed/22187960
2011Rausch.Medulloblastoma
	http://www.ncbi.nlm.nih.gov/pubmed/22265402
2011.Zhang.Retinoblastoma
	http://www.nature.com/nature/journal/v481/n7381/full/nature10733.html
2012Cheung.NB
	http://www.ncbi.nlm.nih.gov/pubmed/22416102
2012Gruber.AML
	http://www.ncbi.nlm.nih.gov/pubmed/23153540
2012Harrison.PoorRiskLeukemia
	http://www.sciencedirect.com/science/article/pii/S1535610812003042
2012.Jones Medulloblastoma
	http://www.nature.com/nature/journal/v488/n7409/full/nature11284.html
2012Kannar.LowerGradeGlioma
	http://www.ncbi.nlm.nih.gov/pubmed/23104868
2012.Lee.Rhabdoid
	http://www.ncbi.nlm.nih.gov/pubmed/22797305
2012Molenaar.NB
	http://www.nature.com/nature/journal/v483/n7391/full/nature10910.html
2012Roberts.HighRiskALL
	http://www.ncbi.nlm.nih.gov/pubmed/22897847
2012Robinson.Medulloblastoma.Germline
	http://www.ncbi.nlm.nih.gov/pubmed/22722829
2012Robinson.Medulloblastoma Somatic
	http://www.ncbi.nlm.nih.gov/pubmed/22722829
2012Wu.Glioblastoma
	http://www.ncbi.nlm.nih.gov/pubmed/22286216
2012Zhang.ALL
	http://www.nature.com/nature/journal/v481/n7380/full/nature10725.html
2013.Chen Rhabdo 
	http://www.ncbi.nlm.nih.gov/pubmed/24332040
2013Holmfeldt.HyperdiploidALL
	http://www.ncbi.nlm.nih.gov/pubmed/23334668
2013Loh.ALL
	http://www.ncbi.nlm.nih.gov/pubmed/23212523
2013.Pugh.NB
	http://www.ncbi.nlm.nih.gov/pubmed/23334666
2013Sausen
	http://www.nature.com/ng/journal/v45/n1/full/ng.2493.html
2013Zhang.Glioma
	http://www.nature.com/ng/journal/v45/n6/full/ng.2611.html	
2014Chen.Osteo	
	http://www.cell.com/cell-reports/abstract/S2211-1247(14)00165-X
2014.Huether.1000PediatricCancerGenomes
	http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4119022/
2014.Shern Rhabdo (44)
	http://cancerdiscovery.aacrjournals.org/content/4/2/216.figures-only?cited-by=yes&legid=candisc;4/2/216
2014.Wu High Grade Glioma(116)
	http://www.nature.com/ng/journal/v46/n5/full/ng.2938.html
2014.Brohl.EWS
	http://www.ncbi.nlm.nih.gov/pubmed/25010205
2013.Dorschner.Actionable	
	http://www.sciencedirect.com/science/article/pii/S0002929713003819
2013.Wei.DW
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0077731
2014.Shern.RMSExpressed
        http://cancerdiscovery.aacrjournals.org/content/4/2/216.figures-only?cited-by=yes&legid=candisc;4/2/216
2014Tirode.EWS
	http://cancerdiscovery.aacrjournals.org/content/early/2014/09/13/2159-8290.CD-14-0622.abstract
2014Crompton.EWS
	http://cancerdiscovery.aacrjournals.org/content/early/2014/09/03/2159-8290.CD-13-1037.abstract
PCG_Total 
	(Total of all the above studies except TCGA)
Grand Total
	(Total of PCG and TCGA)
###############
Clinseqc_genotypes
Clinseqc_homref
Clinseqc_het
Clinseqc_homvar
Clinseqc_hetnonref
Clinseqc_other
Clinseqc_refallele
Clinseqc_varallele
Clinseqfreq_homref
Clinseqfreq_het
Clinseqfreq_homvar
Clinseqfreq_hetnonref
Clinseqfreq_refallele
Clinseqfreq_varallele
Clinseqref_is_minor
Clinseqc_major
Clinseqc_minor
Clinseqmaf
Clinseqchisquare
	Column defination on http://trek.nhgri.nih.gov/wiki/index.php/VarSifter
###############
ExAC_Freq
ExAC_Total
ExAC_AFR
ExAC_AMR
ExAC_EAS
ExAC_FIN
ExAC_NFE
ExAC_OTH
ExAC_SAS
ExHet_Total
ExHet_AFR
ExHet_AMR
ExHet_EAS
ExHet_FIN
ExHet_NFE
ExHet_OTH
ExHet_SAS
ExHom_Total
ExHom_AFR
ExHom_AMR
ExHom_EAS
ExHom_FIN
ExHom_NFE
ExHom_OTH
ExHom_SAS
	Exome Aggregation Consortium data.
	http://exac.broadinstitute.org/
	ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.2/
##############
cBIO Link	
	Link to cBIo main site for the gene
COSMIC Link	
	Link to COSMIC site for the gene


Exome Outputs.
Chr	Start	End	Ref	Alt	gene.Region	Gene	gene.VariantType	gene.AminoAcidChange	CytoBand	dbsnp138	1000g08.2014	1000g08.2014_eur	1000g08.2014_afr	1000g08.2014_amr	1000g08.2014_eas	1000g08.2014_sas	NHLBI6500_ALL	NHLBI6500_EA	NHLBI6500_AA	cg69	nci60	caddRAW	caddPHRED	clinvar_20140702	cosmic70	SIFT Prediction	SIFT Score	PPH2 Prediction	PPH2 Class	PPH2 Probability	Acc_No.hgmd2014.3	Gene.hgmd2014.3	GeneName.hgmd2014.3	Disease.hgmd2014.3	Category.hgmd2014.3	Reference_PMID.hgmd2014.3	CancerGeneCensus	ACMG_reportableGenes	MATCH_v1_08_2014Purpose	MATCH_v1_08_2014aMOI	MATCH_v1_08_2014Drug	MATCH_v1_08_2014Level	MATCH_v1_08_2014References	MATCH_2 Purpose	MCG.Diagnosis	MCG.Targeted Therapy	MCG.Other Implication	Count_TCGA_ALL	2008.Parsons.Glio.Multiforme	2010Barretina.SoftTissueSarcoma	2011Heravi.DICER1	2011Rausch.Medulloblastoma	2011.Zhang Retinoblastoma	2012Cheung.NB	2012Gruber.AML	2012Harrison.PoorRiskLeukemia	2012.Jones Medulloblastoma	2012Kannar.LowerGradeGlioma	2012.Lee.Rhabdoid	2012Molenaar.NB	2012Roberts.HighRiskALL	2012Robinson.Medulloblastoma.Germline	2012Robinson.Medulloblastoma Somatic	2012Wu.Glioblastoma	2012Zhang.ALL	2013.Chen Rhabdo 	2013Holmfeldt.HyperdiploidALL	2013Loh.ALL	2013.Pugh.NB	2013Sausen	2013Zhang.Glioma	2014Chen.Osteo	2014.Huether.1000PediatricCancerGenomes	2014.Shern Rhabdo (44)	2014.Wu High Grade Glioma(116)	2014.Brohl.EWS	Actionable	2013.Wei.DW	2014.Shern.RMSExpressed	2014Tiroda.EWS	2014Crompton.EWS	PCG_Total	Total	Clinseqc_genotypes	Clinseqc_homref	Clinseqc_het	Clinseqc_homvar	Clinseqc_hetnonref	Clinseqc_other	Clinseqc_refallele	Clinseqc_varallele	Clinseqfreq_homref	Clinseqfreq_het	Clinseqfreq_homvar	Clinseqfreq_hetnonref	Clinseqfreq_refallele	Clinseqfreq_varallele	Clinseqref_is_minor	Clinseqc_major	Clinseqc_minor	Clinseqmaf	Clinseqchisquare	ExAC_Freq	ExAC_Total	ExAC_AFR	ExAC_AMR	ExAC_EAS	ExAC_FIN	ExAC_NFE	ExAC_OTH	ExAC_SAS	ExHet_Total	ExHet_AFR	ExHet_AMR	ExHet_EAS	ExHet_FIN	ExHet_NFE	ExHet_OTH	ExHet_SAS	ExHom_Total	ExHom_AFR	ExHom_AMR	ExHom_EAS	ExHom_FIN	ExHom_NFE	ExHom_OTH	ExHom_SAS	cBIO Link	COSMIC Link
chr2	29432664	29432664	C	T	exonic	ALK	nonsynonymous SNV	ALK:NM_004304:exon25:c.G3824A:p.R1275Q,	2p23.2	rs113994087	0	0	0	0	0	0	0	0	0	0	0	5.724913	36	CLINSIG=other;CLNDBN=Neuroblastoma_3;CLNACC=RCV000019709.2	ID=COSM28056;OCCURENCE=1(breast),1(large_intestine),72(autonomic_ganglia)	DAMAGING	0	probably damaging	deleterious	1	CM085230	ALK	Anaplastic lymphoma receptor tyrosine kinase	Neuroblastoma	Diseases Causing Mutation	18923525	yes	-	Hotspot	0	0	0	0	Hotspot	Neuroblastoma	Response to crizotinib,other ALK inhibitors	-	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	8	7	0	0	0	0	0	0	0	0	0	0	0	17	17	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	http://www.cbioportal.org/public-portal/cross_cancer.do?cancer_study_id=all&data_priority=0&case_ids=&gene_set_choice=user-defined-list&gene_list=ALK%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit 	http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=ALK#histo 
chr2	48018182	48018182	C	A	exonic	MSH6	stopgain	MSH6:NM_000179:exon2:c.C377A:p.S126X,	2p16.3	0	0	0	0	0	0	0	0	0	0	0	0	5.534283	35	0		Damaging due to stop	N/A	-	-	-	-	-	-	-	-	-	yes	yes	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	http://www.cbioportal.org/public-portal/cross_cancer.do?cancer_study_id=all&data_priority=0&case_ids=&gene_set_choice=user-defined-list&gene_list=MSH6%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit 	http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=MSH6#histo 
chr5	156184679	156184679	C	A	exonic	SGCD	stopgain	SGCD:NM_000337:exon8:c.C663A:p.C221X,SGCD:NM_001128209:exon7:c.C660A:p.C220X,SGCD:NM_172244:exon8:c.C663A:p.C221X,	5q33.3	0	0	0	0	0	0	0	0	0	0	0	0	7.045342	38	0		Damaging due to stop	N/A	-	-	-	-	-	-	-	-	-	-	-	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	http://www.cbioportal.org/public-portal/cross_cancer.do?cancer_study_id=all&data_priority=0&case_ids=&gene_set_choice=user-defined-list&gene_list=SGCD%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit 	http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=SGCD#histo 
chr1	1372806	1372806	C	T	exonic	VWA1	nonsynonymous SNV	VWA1:NM_199121:exon2:c.C178T:p.L60F,	1p36.33	0	0	0	0	0	0	0	0.000077	0	0.000227	0	0	0.540596	6.926	0		N/A	N/A	-	-	-	-	-	-	-	-	-	-	-	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	http://www.cbioportal.org/public-portal/cross_cancer.do?cancer_study_id=all&data_priority=0&case_ids=&gene_set_choice=user-defined-list&gene_list=VWA1%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit 	http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=VWA1#histo 
chr1	2493172	2493172	G	A	exonic	TNFRSF14	stopgain	TNFRSF14:NM_003820:exon6:c.G612A:p.W204X,	1p36.32	0	0	0	0	0	0	0	0	0	0	0	0	2.852355	15.5	0		Damaging due to stop	N/A	-	-	-	-	-	-	-	-	-	yes	-	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	http://www.cbioportal.org/public-portal/cross_cancer.do?cancer_study_id=all&data_priority=0&case_ids=&gene_set_choice=user-defined-list&gene_list=TNFRSF14%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit 	http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=TNFRSF14#histo 
chr1	6638978	6638978	T	G	exonic	TAS1R1	stopgain	TAS1R1:NM_177540:exon5:c.T1098G:p.Y366X,TAS1R1:NM_138697:exon6:c.T1860G:p.Y620X,	1p36.31	0	0	0	0	0	0	0	0	0	0	0	0	4.24462	22.1	0		Damaging due to stop	N/A	-	-	-	-	-	-	-	-	-	-	-	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	http://www.cbioportal.org/public-portal/cross_cancer.do?cancer_study_id=all&data_priority=0&case_ids=&gene_set_choice=user-defined-list&gene_list=TAS1R1%0D%0A&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit 	http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=TAS1R1#histo 

