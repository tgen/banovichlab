# Annotation for the four cell populations:
# PTPRC+ (immune cells), EPCAM+ (epithelial cells), PECAM1+/PTPRC− 
# (endothelial cells), and PTPRC−/EPCAM−/PECAM1− (mesenchymal cells)
cellpop_markers <- c("PTPRC", "EPCAM", "PECAM")

## Epithelial markers

# AT1 markers
at1_markers <- c("AGER", "PDPN", "CAV1", "EMP2")

# AT2 markers
at2_markers <- c("SFTPC", "ABCA3", "LAMP3", "AGER")

# Secretory markers
secretory_markers <- c("SCGB1A1", "SCGB3A2", "MGP", "MUC5B", "MUC5AC")

# Basal markers
basal_markers <- c("KRT5", "KRT17", "CHGA", "CALCA", "FOXI1", "SERPINB3")

# Ciliated markers
ciliated_markers <- c("FOXJ1", "SFTPB", "OMG", "PPIL6", "NME5", "NWD1", "SLAIN2", "SPAG16")

# Proliferation markers
proliferating_markers <- c("MKI67", "CDK1")

# PNEC: SYP, CHGA
pnec_markers <- c("EPCAM", "SYP", "CHGA")

## Immune markers

# B cell markers
bcell_markers <- c("MS4A1", "CD19", "CD79A")

# Plasma cell markers
plasma_markers <- c("JCHAIN", "IGHG1", "IGLL5")

# Mast cell markers
mast_markers <- c("CPA3", "KIT")

# Monocyte markers
# Inflammatory monocytes: IL-6, IL-8, CCL2, CCL3, and CCL5. CCR2, GR1
monocyte_markers <- c("CD14", "CD16", "S100A12", "FCN1", "S100A9", "LYZ", "CD14", "CCL2", "CCL3", "CCL5", "IL6", "IL8", "CCR2", "GR1")

# cDC markers
# cDC1: CD11C (ITGAX)
# cDC2: IRF4
cdc_markers <- c("FCER1A", "CD1C", "CLEC9A", "IRF4", "CD11C")

# pDC markers
pdc_markers <- c("LILRA4", "CLEC4C", "JCHAIN")

# Macrophage markers
macrophage_markers <- c("LYZ", "MARCO", "FCGR1A", "C1QA", "APOC1", "SPP1")

# Proliferating Macrophages
prolif_macrophage_markers <- c("MKI67", "CD1", "LYZ")

# NK cell markers. NKG7 (high), CD8A-
nk_markers <- c("NCR1", "KLRB1", "NKG7", "CD8A", "GNLY")

# T cell markers
tcell_markers <- c("CD3E")

# Proliferating T Cells
prolif_tcell_markers <- c("MKI67", "CD1", "CD3E")

# TReg markers
treg_markers <- c("FOXP3")

# CD8 T cells
cd8_tcell_markers <- c("CD8A")

# CD4 T cells, CD8A-
cd4_tcell_markers <- c("CD4", "IL7R", "CD8A")

## Mesenchymal markers

# Mesothelial
mesothelial_markers <- c("MSLN", "UPK3B", "HP", "WT1")

# Smooth Muscle Cells (SMCs)
smc_markers <- c("ACTA2", "PDGFRB", "MYH11", "TAGLN", "DES", "ACTG2")

# Pericytes
pericyte_markers <- c("ACTA2", "PDGFRB", "RGS5", "HIGD1B", "GJA4")

# Fibroblasts
fibroblast_markers <- c("LUM", "DCN", "PDGFRA")

# Myofibroblasts
myofb_markers <- c("MYLK", "ACTA2", "COL8A1", "COL1A1", "WNT2")

# Activated MyoFB
activated_myofb_markers <- c("COL1A1", "POSTN", "CTHRC1")

# WNT2+ Fibroblasts
WNT2_fibro_markers <- c("MYLK", "WNT2", "A2M", "GPC3", "MACF1", "CES1", "LIMCH1")

# Matrix Fibroblasts
matrixfb_markers <- c("SFRP2", "CLU", "APOD", "FBLN1", "CST3", "IGFBP6", "SCARA5", "CD34")

# PLIN2+ Fibroblasts
PLIN2_fibro_markers <- c("PLIN2", "HAS1")

# HAS1 High Fibroblasts
HAS1_fibro_markers <- c("HAS1", "TWIST1", "PLIN2")

# Activated HAS1 High FB
act_HAS1_fibro_markers <- c("LIF", "ICAM1", "CXCL2", "HAS1", "PLIN2")

## Endothelial markers

# Vascular
vascular_markers <- c("VWF")

# Capillary
capillary_markers <- c("CA4", "VIPR1", "RGCC", "CYB5A", "ADGRL2")

# Capillary G1
g1_markers <- c("VWF", "HPGD", "EDNRB", "EMCN")

# Capillary G2
g2_markers <- c("VWF", "FCN3", "IL7R", "SLC6A4")

# Arterial
arterial_markers <- c("DKK2", "GJA5", "BMX", "HEY1")

# Venous
venous_markers <- c("PLA1A", "CPE", "PTGDS", "ACKR1", "SPRY1")

# Bronchial
bronchial_markers <- c("SPRY1", "PLVAP", "COL15A1", "VWA1", "MYC")

# Bronchial G1
bronchial_g1_markers <- c("POSTN", "ACKR1", "SPRY1")

# Bronchial G2
bronchial_g2_markers <- c("HBEGF", "ACKR1", "SPRY1")

# Lymphatic
lymphatic_markers <- c("CCL21", "PROX1", "PDPN")




