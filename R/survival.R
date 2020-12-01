library(TCGAbiolinks)
library(SummarizedExperiment)

library(TCGAbiolinks)
library(SummarizedExperiment)

clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
save(clinical, file="/home/minzhang/TCGA/clinical.RData")
write.table(clinical,file = "/home/minzhang/TCGA/clinical.txt", sep='\t')

query <- GDCquery(project = "TCGA-BRCA",
				  data.category = "Transcriptome Profiling",
				  data.type = "Gene Expression Quantification",
				  workflow.type = "HTSeq - Counts")
GDCdownload(query, files.per.chunk = 20, directory="/home/minzhang/TCGA")
expr <- GDCprepare(query = query, 
                   save = TRUE, 
                   directory =  "/home/minzhang/TCGA",
                   save.filename = "/home/minzhang/TCGA/BRCA_rnaseq_counts.RData")
write.table(expr,file = "/home/minzhang/TCGA/BRCA_rnaseq_counts.txt", sep='\t')
save.image(file = "/home/minzhang/TCGA/BRCA_rnaseq_counts_all.RData")

query <- GDCquery(project = "TCGA-BRCA",
				  data.category = "Transcriptome Profiling",
				  data.type = "Gene Expression Quantification",
				  workflow.type = "HTSeq - FPKM")
GDCdownload(query, files.per.chunk = 20, directory="/home/minzhang/TCGA")
expr <- GDCprepare(query = query, 
                   save = TRUE, 
                   directory =  "/home/minzhang/TCGA",
                   save.filename = "/home/minzhang/TCGA/BRCA_rnaseq_fpkm.RData")
write.table(expr,file = "/home/minzhang/TCGA/BRCA_rnaseq_fpkm.txt", sep='\t')
save.image(file = "/home/minzhang/TCGA/BRCA_rnaseq_fpkm_all.RData")

query <- GDCquery(project = "TCGA-BRCA",
				  data.category = "Transcriptome Profiling",
				  data.type = "Gene Expression Quantification",
				  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query, files.per.chunk = 20, directory="/home/minzhang/TCGA")
expr <- GDCprepare(query = query, 
                   save = TRUE, 
                   directory =  "/home/minzhang/TCGA",
                   save.filename = "/home/minzhang/TCGA/BRCA_rnaseq_fpkm_uq.RData")
write.table(expr,file = "/home/minzhang/TCGA/BRCA_rnaseq_fpkm_uq.txt", sep='\t')
save.image(file = "/home/minzhang/TCGA/BRCA_rnaseq_fpkm_uq_all.RData")

query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query, files.per.chunk = 20, directory="/home/minzhang/TCGA")
clinical <- GDCprepare_clinic(query, clinical.info = "patient", directory="/home/minzhang/TCGA")
save(clinical, file="/home/minzhang/TCGA/clinical_patient.RData")
write.table(clinical,file = "/home/minzhang/TCGA/clinical_patient.txt", sep='\t')
clinical <- GDCprepare_clinic(query, clinical.info = "follow_up", directory="/home/minzhang/TCGA")
save(clinical, file="/home/minzhang/TCGA/clinical_follow_up.RData")
write.table(clinical,file = "/home/minzhang/TCGA/clinical_follow_up.txt", sep='\t')
save.image(file = "/home/minzhang/TCGA/clinical_xml_all.RData")


library(TCGAbiolinks)
load("/home/minzhang/TCGA/BRCA_rnaseq_counts.RData")
expr <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6,
									  datatype = "HTSeq - Counts")
rownames(expr) <- data@rowRanges@elementMetadata@listData$external_gene_name
write.table(expr, file='/home/minzhang/TCGA/expr_htseq_counts.xls', sep='\t')

load("/home/minzhang/TCGA/clinical_patient.RData")
write.table(clinical, file='/home/minzhang/TCGA/clinical_patient.xls', sep='\t')


library(TCGAbiolinks)
load("/home/minzhang/TCGA/BRCA_rnaseq_counts.RData")
expr <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6,
									  datatype = "HTSeq - Counts")
rownames(expr) <- data@rowRanges@elementMetadata@listData$external_gene_name
write.table(expr, file='/home/minzhang/TCGA/expr_htseq_counts.xls', sep='\t')

load("/home/minzhang/TCGA/clinical_patient.RData")
write.table(clinical, file='/home/minzhang/TCGA/clinical_patient.xls', sep='\t')

library(survminer)
library(survival)
library("doFuture")
genes <- c("CRISP3", "AFF3", "ESR1", "CLIC6", "ST8SIA6.AS1", "AGR3", "NPY1R", "PIP", "STC2", "PDZK1", "FDCSP", "SOX11", "TFF1", "SUSD3", "PGR", "SCUBE2", "S100A7", "CLCA2", "TPRG1", "CHPT1", "SERPINA3", "S100A8", "LRP2", "MYBPC1", "MMP1", "FGG", "NAT1", "ABAT", "SLC39A6", "GREB1", "RLN2", "RBM24", "CXCL13", "TCN1", "HOTAIR", "BCL2", "IL6ST", "ADAMTS15", "PPAPDC1A", "ADH1B", "FGB", "PDK4", "TMPRSS3", "SLC27A2", "S100A9", "KRT80", "C7", "GREM1", "COL11A1", "C1orf168", "ERBB4", "GFRA1", "CLSTN2", "ANLN", "SERPINA1", "SEC14L2", "FSIP1", "MUC19", "CPB1", "BRINP3", "NANOS1", "C1orf106", "CEL", "CP", "KCND3", "COMP", "LINC00472", "MS4A1", "AGTR1", "MRPS30", "PARD6B", "CHAD", "FMO5", "COL10A1", "RTN1", "CEMIP", "CDCA7", "FAM83D", "ELP2", "ARHGAP36", "RP11.53O19.3", "HRASLS", "NEK10", "CD1E", "C5orf46", "P2RY12", "KIF14", "ANKRD30A", "HEPACAM2", "MUCL1", "PI15", "ELOVL2", "SLC26A3", "LTF", "KCNE4", "LINC01087", "NBPF4", "MMP11", "SLC7A2", "KRT6A", "PRAME", "SCN7A", "BIRC5", "MCM10", "EVL", "KIF20A", "ATP1A4", "TPX2", "PLAU", "TMEM65", "GLYATL2", "ABCA8", "SLC16A6", "NTRK2", "MYB", "TBC1D9", "GRB7", "CBX2", "IBSP", "DYNLRB2", "LAPTM4B", "CXADR", "REPS2", "FGFR2", "FCGBP", "FLJ13744", "SPOCK1", "LAMA3", "ASPM", "LIFR", "CYBRD1", "NPR3", "INHBA", "CDC20", "RGS4", "EFNA2", "TXNIP", "GLUL", "P2RY13", "TFAP2B", "PSAT1", "DLX2", "CFB", "HSPA2", "NME5", "MAOA", "PCSK6", "IL20RA", "ACADSB", "POSTN", "CORIN", "SLC7A5", "C1orf21", "DLGAP5", "E2F8", "PGM2L1", "FAM129A", "THBS2", "CALML5", "TFF3", "DNAJC12", "CAPN8", "RERG", "PKIB", "ACKR1", "KLHDC7A", "CX3CR1", "FKBP5", "CELF2", "LYPD6", "CCNB2", "AEBP1", "ADAMTS12", "ADIPOQ", "OGN", "DACH1", "CXCL8", "DEFB1", "IL33", "DUSP4", "FCER1A", "MT1M", "TMEM45A", "GSTM3", "PITX1", "EPPK1", "CCDC170", "SDC1", "PCDH7", "RGS13", "RHOH", "EIF4EBP1", "ITPR1", "TMEM40", "EPHX2", "ENC1", "LINC00993", "LOC101926959", "CFD", "C6orf141", "HOXC10", "BIRC3", "SPAG6", "CEP55", "NCAPG", "KIF4A", "NOX4", "CASC1", "FGL2", "CENPF", "LBP", "CCL19", "PRR11", "ADM", "ITM2A", "FN1", "F2RL1", "TREM1", "CCDC176", "RP11.28F1.2", "MMP13", "AKR1B10", "FBN2", "NOVA1", "THSD4", "IL17RB", "FNDC1", "CD69", "ASS1", "C2orf54", "SLITRK5", "BUB1B", "MELK", "COL5A1", "KIF2C", "DEPDC1B", "PITPNM3", "KRT81", "CHRDL1", "SDPR", "KANK4", "IGF1R", "LOC145837", "LY6D", "ARNT2", "EGLN3", "SCEL", "RRM2", "ARSG", "NEK2", "ENTPD5", "LOXL2", "COL13A1", "NXPE3", "SRPX2", "N4BP2L1", "CA12", "KLK5", "GAD1", "SCARA5", "ANKRD29", "BMP7", "FAM83B", "LY75", "CCNB1", "AURKA", "CDKN3", "PDZD2", "FGD3", "H2BFXP", "WDR96", "IQGAP3", "SIAH2", "ZNF83", "GLA", "GAS2L3", "WFDC1", "NDP", "DNALI1", "SLC1A2", "S100A7A", "PMAIP1", "GGH", "COL1A2", "SLC6A16", "PLOD2", "DACT1", "ZNF655", "LRIG1", "MAOB", "ZNF385B", "LRRN1", "ENPP5", "KCNK1", "CST6", "DCDC2", "CCNE2", "AK5", "PLCH1", "TSPAN7", "UBE2C", "BUB1", "GTSE1", "PFKP", "MGP", "ZBTB16", "CD24", "MAPT", "VAV3", "C1orf64", "TRDV3", "TGFB2", "TTK", "CELSR1", "FA2H", "HMMR", "FOXM1", "KIF18B", "ABCD3", "GSTM2", "PXDN", "DNMT3B", "CST3", "SMARCA2", "LYZ", "LONRF2", "RAI2", "PCAT18", "PTPRT", "DTL", "PREX1", "LOC400768", "GINS1", "KCTD15", "SLC1A1", "GJB2", "LRRC15", "MMP3", "SYNPO2", "NMU", "PTPRD", "FAM107A", "SIX2", "CDK1", "NUSAP1", "PRC1", "AURKB", "MAP3K1", "UBXN10", "CENPI", "LOC286052", "THBS1", "BTBD19", "C2CD4A", "ASPN", "KIAA1324", "COL5A2", "SPC25", "PHYHD1", "RACGAP1", "LOC100506119", "PTP4A3", "DBN1", "GZMK", "CLU", "GIMAP7", "MTFR2", "KIF23", "CASC5", "TANC2", "DIAPH3", "KLHL5", "PCDH17", "JAK2", "COL18A1", "IGJ", "ADIRF", "PPP4R4", "GPR126", "TTC39A", "C1orf116", "PTGER3", "AIF1L", "MND1", "SH3BGRL", "TTC36", "CALD1", "PCM1", "CENPN", "CENPE", "LRRC48", "FAM105A", "YWHAZ", "HMGCS2", "AGR2", "FRAS1", "C15orf48", "MYCN", "FAM198A", "TPH1", "C4orf32", "ITPR2", "NUCB2", "SYCP3", "PNRC2", "C2orf40", "DCLK1", "FAM196A", "IGFBP5", "BAIAP2L1", "GSTM1", "RABEP1", "BGN", "CLEC2D", "PTK7", "PEG10", "GATA3.AS1", "ERBB2", "LOC101930067", "CACNA1D", "TGFBR3", "ZNF652", "ATP6V1C2", "GP2", "NKX3.1", "SLC12A1", "FAM134B", "COL12A1", "LYPD6B", "DEPDC1", "C6orf211", "SULF1", "RSPH1", "CDON", "MBOAT1", "SQLE", "LINC01016", "ZNF689", "PTTG1", "TNFRSF21", "CD74", "CCNE1", "ASAH1", "TMEM101", "CD55", "KLHL7", "ELOVL5", "RAB8B", "HUNK", "KIF13B", "SERPINH1", "FAM161B", "FABP4", "MKX", "LINC01105", "TTC18", "WDR78", "MSI2", "FHOD1", "L1CAM", "IRAK3", "COL3A1", "PGLYRP2", "KIF15", "SUGCT", "CREBL2", "WNK4", "ABI3BP", "EDIL3", "PBK", "HSD17B6", "CDH2", "CDH11", "VGLL3", "TK1", "BTG2", "PGBD5", "ESM1", "PRICKLE1", "KDM4B", "POLQ", "AGBL2", "SAPCD2", "WFDC2", "MAPT.AS1", "LOX", "SCARA3", "MYBL2", "EN1", "GUCY1A2", "KIAA1467", "COL1A1", "GABRB3", "MEX3A", "ANTXR1", "CRISP2", "C1orf226", "SLC44A4", "CALB2", "PAPPA", "UBE2T", "MKI67", "GSTM4", "TTC39C", "GPR171", "STARD3", "CDO1", "LINC00173", "FAM198B", "ULBP2", "CPEB3", "HPSE", "ABCA12", "SOWAHA", "FUT9", "UGCG", "GPR56", "KCNN2", "PSD3", "HLF", "EYA2", "CTTN", "THY1", "DENND2D", "KLK8", "PTGDS", "SLC15A2", "ADAMTS2", "XBP1", "ZNF396", "VGLL1", "TMEM150C", "INPP4B", "BCAT1", "HOXC9", "IL18R1", "ARHGDIB", "DENND5B", "ITK", "PLAC8", "EFCAB4A", "SGIP1", "ZWINT", "RAD54B", "PSTPIP2", "SMYD2", "ARHGEF6", "ARHGEF3", "GATA3", "TYRP1", "CDH3", "EDN3", "RP11.111M22.3", "SMIM14", "ABCA6", "UBE2S", "CD44", "PCOLCE2", "BANK1", "BC032415", "RAB23", "ZBED5.AS1", "PGAP3", "SLC26A7", "ADAM12", "FBP1", "SYNPO2L", "MAP1B", "PCSK1N", "HSD11B1", "LOC101929122", "MAD2L1", "GIMAP2", "LOC102723927", "CCR6", "CDCA3", "ITGA11", "GIMAP6", "FGD6", "PPFIBP1", "CHSY3", "NREP", "STK38L", "MRVI1", "RP11.803D5.4", "SNAP23", "IER5L", "RUNDC1", "CSNK1E", "CYP4X1", "GZMA", "PVR", "HLA.DQB1", "ZIC2", "FAM155A", "UHRF1", "MPP7", "ANGPTL2", "RFX3", "TLDC1", "PDLIM7", "GRB14", "ADCY1", "ST6GAL2", "TMEM139", "GAL", "FOS", "JPH1", "LDLRAD3", "OGFRL1", "TAGLN", "CLIC3", "KRT16", "COL8A1", "LRG1", "MALAT1", "LHX8", "CCNG2", "TRIP13", "FAM214A", "PTK2", "PCAT6", "ZNF540", "PROM1", "FGFR4", "AHNAK2", "SPAG5", "CDC14A", "MED13L", "PLEKHA8", "EDN2", "LZTFL1", "GABRP", "ATP13A5", "CITED1", "CITED2", "MIR210HG", "SPG20", "KIAA1551", "E2F7", "HOPX", "GALNT14", "GULP1", "ZNF44", "TRIM45", "HJURP", "PRKAA2", "GLRB", "MIEN1", "C16orf54", "RBM20", "MEG3", "PLA2G16", "LRRC6", "SYNDIG1", "SCG5", "NT5DC2", "RPL7", "DCD", "PALM2", "ENPEP", "SLC22A4", "RNF125", "MAATS1", "PTPRF", "PLGRKT", "UGT2B28", "SLC28A3", "DNAH5", "PTPLAD2", "CENPK", "ALDH6A1", "MYO1D", "RAB6B", "APOD", "KIF5C", "ONECUT2", "SESN3", "ACKR4", "TNFAIP6", "PLAUR", "LCN2", "SYBU", "RASGRP1", "FRZB", "CELSR3", "HRASLS5", "FAM83A", "LTBP1", "OSBPL1A", "PTPRC", "CHDH", "SCN3A", "NID2", "TMEM26", "CDK12", "EGOT", "LMCD1", "RP11.111M22.4", "CLIC2", "HOXC13", "FLNB", "P4HA3", "SEMA6A", "SPRR1A", "GPC6", "GALNT7", "SPTSSA", "RMI2", "ZBTB18", "ALDH1A1", "CYP2J2", "CDT1", "SKAP1", "DHRS2", "KCNMA1", "RAC2", "LRP12", "ABHD3", "AVEN", "SLC4A7", "KCNK3", "PRKG1", "PTP4A2", "GBP3", "TNNT1", "SLC6A4", "CAPS", "PLD1", "IL18", "ABCB1", "VCAN", "ADAMTS6", "MDH1B", "TGFB1I1", "GRIA2", "PRKCB", "RDH10", "FGF10", "FUCA1", "CSTA", "SLC4A11", "CRTAM", "LRMP", "FBN1", "PDE9A", "TMEM144", "HIPK2", "ADRBK2", "STON1", "TAPBPL", "CLDN11", "FGFR1", "PNPLA4", "ANKRD6", "TSC22D2", "TFPI2", "DKK1", "CXCL17", "CDC6", "SIM2", "CDKN2B", "MS4A14", "KIF11", "RP5.1092A3.4", "SLC19A2", "C11orf96", "ACSF2", "MEGF6", "DDX17", "NXN", "NR2F1.AS1", "GK", "JADE2", "SIK3", "MYO1B", "HAUS1", "PIEZO1", "ALG13", "EDN1", "XPOT", "SHMT2", "DNER", "RGS22", "CTHRC1", "GSDMB", "RGS16", "ACTA2", "SLC6A1", "VWDE", "FAM120AOS", "ZNF711", "PLCB4", "KIF26B", "C11orf80", "GJA5", "MSRB3", "COL4A1", "DONSON", "KLC3", "RECQL4", "RASAL2", "KIT", "RPS16P5", "MFAP2", "SCRN1", "RP11.274H2.5", "GARS", "BEX1", "RIMS2", "HOMER1", "GPX8", "PERP", "DSCC1", "KIAA1161", "SMTN", "STXBP3")
genes2 <- intersect(genes, colnames(expr))
for (gene in genes2) {

	res.cut <- surv_cutpoint(expr, time="time", event="event", variables=c(gene))
	summary(res.cut)
	res.cat <- surv_categorize(res.cut)
	cn <- colnames(res.cat)
	cn[3] <- 'group'
	colnames(res.cat) <- cn
	fit <- survfit(Surv(time, event) ~ group, data = res.cat)
		
	# ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE, pval=TRUE)
	log_p <- surv_pvalue(fit)
	log_p <- as.numeric(log_p$pval)
	res.cox <- coxph(Surv(time, event) ~ group, data = res.cat)
	s <- summary(res.cox)
	hr <- as.numeric(s$conf.int[1])
	hr_951 <- as.numeric(s$conf.int[3])
	hr_952 <- as.numeric(s$conf.int[4])
	wald_p <- as.numeric(s$waldtest[3])
	sc_p <- as.numeric(s$sctest[3])
	print <- c(gene, hr, hr_951, hr_952, wald_p, sc_p, log_p)
	df <- rbind(print, df)
}