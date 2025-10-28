
## ------- Male & Female
# from https://www.sciencedirect.com/science/article/pii/S009286742300908X
M_F_markers <- list("Female" = c("XIST", "PSPHP1"),
                    "Male" = c("TTTY15", "KDM5D", "TTTY14",
                               "ZFY", "PRKY", "GYG2P1", "PSMA6P1",
                               "USP9Y", "TXLNGY", "UTY", "RPS4Y1",
                               "DDX3Y", "EIF1AY"))

## add more markers
blood_markers_more = c(
  "CD34", "SPINK2", "MLLT3", "HLF", "MECOM", "RUNX1", "HOXA9",
  "CDK6", "SELL", "CD52", "PROM1", "MEIS1", "MYB", "ITGA6", # HSC/MPP
  "MPO", "AZU1", "SPI1", "LYZ", # Granulocyte
  "GATA1", "GATA2", "TESPA1", "KLF1", "CTNNB1", # MEMPs (megakaryocyte-erythroid-mast cell progenitor)
  "FAM178B", "BLVRB", "TFRC", "AHSP", "ALAS2", "HBA1", "HBB", "GYPA", "BPGM", # Erythroid
  "ITGA2B", "GP9", "PLEK", 'MPL', 'PECAM1', 'CXCR4', "PPBP", "PF4", # Megakaryocytes
  "HDC", "CPA3", "LMO4", "CD63", 'ENPP3', "TPSAB1", "TPSB2", # Mast cells
  "CD14", "FCGR3A", "FCN1", "S100A9", "CD68", "MNDA", # Monocytes
  "CD163", "MS4A7", "C1QA", "MRC1", "CTSB", "MARCO", "CD5L", "VCAM1", # Kupffer cells
  "CLEC9A", "THBD", "XCR1", "BATF3", # cDC1
  "CD1C", "CLEC4A", "CLEC10A", # cDC2
  "FLT3", "VCAN", # cDC3
  "ACY3", "IRF8", "CLEC4C", "IL3RA", "MPEG1", # pDCs
  "AXL", "SIGLEC6", # ASDC
  "IL7R", "JCHAIN", "LTB", "CD7", # LP
  "EBF1", "PAX5", "CD79A", "CD79B", "MME", "IGLL1", "IGHM", "IGHD",
  "CD19", "MS4A1", "IRF4", "DNTT", "RAG1", "RAG2", "CD24", "CD38", # B cells
  "CD5", "CD27", "SPN", "CCR10", # B1
  "IL2RB", "KLRD1", "KLRF1", "NCR1", "NCAM1", "PRF1", "GZMA", "GNLY", "NKG7", # NK
  "TCF7", "RORC", "AHR", "ID2", "NCR2", # ILC3
  "CD3D", "CD3E", "CD3G", "TRAC", "FOXP3", "TIGIT", "CD4", "CD8A", "CD8B", # T
  "KIT", "GATA3", "IL1A", "IL1B",
  "PTPRC", # CD45
  "ALB", "AFP", # Hepatocytes
  "CDH5", "KDR", # endothelial cells
  "STAB1", "STAB2", "LYVE1", "DCN", # LSECs
  "COL1A1", "COL3A1", "RBP1", # stellate cells
  "KRT19",
  'MKI67', "TOP2A" # cycling
)

# less markers to show for dotplot
blood_markers_less = c(
  "CD34", "SPINK2", "HLF",
  "MPO", "AZU1",
  "GATA2", "TESPA1", "KLF1",
  "FAM178B", "BLVRB", "ALAS2", # "SLC25A21", # also a good one (https://pmc.ncbi.nlm.nih.gov/articles/PMC7793295)
  "ITGA2B", 'PECAM1',
  "HDC", "CPA3", "ENPP3",
  "SPI1", "MRC1",
  "FCN1", "VCAN",
  "CD163", "C1QA",
  "CLEC9A", "BATF3",
  "CD1C", "CLEC10A",
  "CLEC4C", "IL3RA",
  "AXL", "SIGLEC6",
  "ACY3", "IRF8",
  "IL7R", "JCHAIN",
  "EBF1", "PAX5", "IGHM", "IGHD",
  "NCR1", "NCAM1",
  "RORC",
  "CD3D",
  "AFP", "CDH5",
  "MKI67"
)

# HSC markers
HSC_sanity_check = c(
  # HSC
  "CD34", "SPINK2", "MLLT3", "HLF", "MECOM", "RUNX1", "HOXA9",
  "CDK6", "SELL", "CD52", "PROM1", "MEIS1", "MYB", "ITGA6",
  # GP
  "MPO", "AZU1", "SPI1", "LYZ",
  # MEMP
  "GATA1", "GATA2", "TESPA1", "KLF1", "CTNNB1",
  # LP
  "IL7R", "JCHAIN", "LTB", "CD7", "IL2RG",
  # B-lin
  "EBF1", "PAX5", "CD19",
  # ILC
  "TCF7", "RORC",
  # NK/T-lin
  "IL2RB", "KLRD1", "CD3D", "BCL11B",
  # Cycling
  "MKI67", "TOP2A"
)

# Mapping human haematopoietic stem cells from haemogenic endothelium to birth
Hanna_HSC_signature <- c("RUNX1", "HLF", "HOXA9", "MLLT3", "MECOM", "SPINK2")
# Fig. 2f
Hanna_HSC_maturation <- list(
  Up = rev(c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M",
             "HLA-DMA", "HLA-DPB1", "HLA-DRA", "HLA-DQA1", "HLA-DPA1", "HLA-DQB1",
             "MLLT3", "HLF", "MALAT1", "MSI2", "EVI2B", "SOCS2", "HEMGN", "HOPX", "SPINK2", "CD52", "SELL", "PROM1")),
  Dn = rev(c("CDH5", "MEIS2", "RUNX1T1", "HOXB9", "ESAM", "PLVAP", "SELP", "ITGA2B", "GP9",
             "RAB27B", "GMPR", "GFI1", "GBP4", "MECOM", "IL3RA", "IL6R", "IL11RA", "IFNAR1",
             "IFNAR2", "IFNGR1", "CSF1R", "CSF2RA", "IGFBP2", "HMGA2", "LIN28B", "MKI67", "TOP2A", "AURKB"))
)
# Supplementary Table 4
Hanna_HSC_maturation_more <- list(
  Dn = c("RASGRP3", "UBL5", "MTHFD1L", "ID1", "SRSF10", "SAAL1", "CCND1", "NINJ2", "GCSH", "FYB", "TUSC3", "MMP2", "ENY2", "HNRNPA2B1", "UQCR10", "ALDOA", "RRP7A", "CKAP4", "EPB41L2", "MYEOV2", "CTPS1", "SNCG", "DEXI", "BOD1", "NDUFB9", "WDR1", "ITGB1", "EIF3A", "LAPTM4B", "S100A16", "UQCR11", "KIAA0101", "AURKAIP1", "PGK1", "PFDN2", "ESAM", "SLC9A3R2", "FASN", "NDUFB2", "AHCY", "MRPL20", "CALM3", "HSP90B1", "XRCC5", "EIF3J", "YWHAE", "VDAC1", "CTSC", "TMEM97", "DDX18", "COX6C", "HIST1H4C", "PSMB3", "H3F3A", "ACOT7", "TCEB2", "SEC11A", "SCD", "COPRS", "PPP1CC", "PSMD14", "GNB2", "PIN1", "CHCHD2", "S100A13", "ATP6V0B", "FUS", "HSPA8", "CMTM3", "MTA1", "MAZ", "CALU", "CFP", "FUCA2", "TRAPPC2L", "GNAS", "SAE1", "NT5DC2", "SON", "KTN1", "HSPD1", "MRPL51", "EIF5", "FRAT2", "PROCR", "ISYNA1", "ID2", "SSBP4", "PNP", "HMGB3", "TIMM8B", "SLC1A5", "IDH2", "PSMA4", "CDV3", "CLTA", "MDH2", "TMEM98", "C2orf88", "DDX21", "SNRPE", "SFPQ", "TRIM28", "KHDRBS1", "TUBB4B", "LAMTOR5", "SSB", "NDUFV3", "RDX", "MRPS12", "BAZ1B", "ATP5F1", "COX7A2", "NDUFA4", "H2AFX", "TMEM258", "CNRIP1", "RAMP2", "FABP5", "SYNCRIP", "NDUFAB1", "MEIS2", "NASP", "IGF2BP1", "C17orf89", "HOXB8", "HDGF", "CCT6A", "RPL41", "LOXL1", "PHPT1", "DYNLL1", "EIF3B", "SAC3D1", "PRDX2", "ACTR2", "PSMC4", "ECI1", "ZWINT", "SMARCC1", "C4orf48", "NETO2", "GOLM1", "RNF145", "SRSF3", "RPS26", "GTPBP4", "ATP5E", "TIMM13", "QPRT", "KRT18", "ERH", "NDUFS6", "PSIP1", "ATP6V1B2", "S100A11", "HSPA5", "RUNX1T1", "PPIB", "HSP90AA1", "RPL6", "ATP5J2", "ATP5B", "YWHAH", "MINOS1", "PPIF", "COX5A", "CENPV", "ANP32A", "PTP4A2", "HSBP1", "CDH5", "KIF5B", "COX7B", "LSM4", "TMSB10", "KPNA2", "RPL37", "CYC1", "USMG5", "RHOC", "COX8A", "TUFM", "RNF187", "PPP1CA", "PDLIM7", "TYMS", "TRAP1", "SPI1", "NDUFA2", "LMNB1", "SNRPD1", "SRSF7", "PPM1G", "BRD7", "EIF4EBP1", "ARPC1A", "CPVL", "EIF4A3", "KDELR2", "UBE2S", "NDUFB1", "DDX3Y", "TUBA1B", "CALM1", "ACTB", "HSPE1", "SNRNP25", "SLC25A5", "SNRPB", "PGP", "MGST1", "CEBPB", "PSMA7", "SMARCA4", "CALCRL", "NME4", "CST3", "RAC1", "PRELID1", "NUCKS1", "FAM101B", "PKM", "GNAI2", "NCL", "PDLIM1", "CHD7", "CHCHD10", "PA2G4", "RCN1", "CDT1", "SLC25A6", "HOXB9", "ATP5G1", "H1FX", "SRRM2", "MAP2K2", "NME1", "HN1", "CD63", "HBZ", "PRDX4", "SNRPF", "PPP1R14B", "ATOX1", "TMA7", "CCT5", "SET", "PSMB5", "NDUFS5", "MARCKS", "SLC25A3", "TPI1", "FBLN1", "SRSF9", "C14orf2", "C1QBP", "EIF5A", "SNRPG", "YBX1", "RPL23", "PFN1", "H2AFZ", "FKBP10", "TOMM40", "PGAM1", "POLR2L", "SRM", "ODC1", "HBE1", "IGFBP4", "ATP5I", "NAA38", "RANBP1", "CCDC85B", "HNRNPD", "EIF4A1", "GLRX5", "SERBP1", "RAN", "RPS21", "ECSCR.1", "EIF1AY", "HNRNPAB", "TXNDC17", "CFL1", "TPM2", "PLVAP", "CSF1R", "GAPDH", "SLIRP", "CD24", "RPL37A", "MDFI", "MARCKSL1", "IGFBP2", "COTL1", "HMGA1", "RPS4Y1"),
  Up = c("RPS4X", "CD74", "HLA-DRA", "CD37", "HLA-DPB1", "HLA-DPA1", "SELL", "RPL13", "HLA-A", "B2M", "HLA-DRB1", "PLAC8", "HLA-B", "HINT1", "PRSS57", "GNB2L1", "RPS3A", "RPL10", "XIST", "RPS24", "ANKRD28", "FXYD5", "TPT1", "RPL11", "RPL31", "HLA-DQB1", "HLA-DMA", "MLLT3", "EEF1A1", "CD52", "HOPX", "SPINK2", "EVI2B", "RPL9", "AVP", "RPL32", "RPS12", "MALAT1", "SOCS2", "EIF3E", "CRHBP", "CLEC2B", "ICAM3", "HLA-DQA1", "EEF2", "CASP4", "COMMD6", "HEMGN", "KIAA0125", "TKT", "ITM2C", "LSP1", "S1PR4", "RPL3", "RPLP0", "GMFG", "MZB1", "C6orf48", "MEG3", "RPS18", "RPS9", "MSRB3", "H2AFY", "SERPINB1", "RARRES3", "GIMAP2", "GNA15", "RPS15A", "HLF", "RPL18", "SEPT6", "SRGN", "RPL13A", "SARAF", "RPL30", "HLA-E", "PROM1", "LRBA", "RP11-620J15.3", "RPS6", "GABPB1-AS1", "PNISR", "HSD17B11", "PBXIP1", "GBP2", "RPS5", "SORL1", "HLA-C", "RAB37", "LITAF", "RPL10A", "HTR1F", "ARHGAP15", "CD48", "TRBC2", "PSMB9", "CRYGD", "CARD16", "SMDT1", "ACAP1", "TSTD1", "SMIM3", "RPS14", "IGBP1", "CLEC9A", "MSI2", "LGALS9", "TNFRSF14", "NEAT1", "MAP7", "MYO1G", "CTSW", "FHL1", "PCBP2", "AJ006998.2")
)
# Fig. 2e
Hanna_HSC_TFs <- c("RUNX1", "SPI1", "TAL1", "LYL1", "LMO2", "MYC", "ETV6", "BCL11A", "GFI1", "GFI1B",
                   "GATA2", "GATA3", "MYB", "MLLT3", "HLF", "HOXA7", "HOXA9", "PBX1", "PRDM16", "MECOM")

# table 1 of "The genesis of human hematopoietic stem cells"
Hanna_HSC_genesis <- list(
  EmbryonicLiver_wk6_wk7 = c("RUNX1", "MLLT3", "HOXA9", "SPINK2", "HLF", "MECOM", "IGFBP2", "LIN28B", "HOXB9"),
  FetalLiver_wk8_wk12 = c("RUNX1", "MLLT3", "HOXA9", "SPINK2", "HLF", "MECOM", "HEMGN"),
  FetalLiver_wk13_wk20 = c("RUNX1", "MLLT3", "HOXA9", "SPINK2", "HLF", "MECOM", "HEMGN", "MSI2")
)

# Mapping the developing human immune system across organs
# progenitors
suo_science_prog <- c(
  "CD34", "SPINK2", "MLLT3", "HLF", "CLEC9A", # HSC/MPP
  "BCL11A", "IL7R", "IL2RG", # MLP
  "CD19", "VPREB1", "PAX5", # B lineage
  "CD3D", "BCL11B", # T lineage
  "GATA2", # MEMP
  "KLF1", "ITGA2B", # MEP
  "HBB", "HBG1", # Ery lineage
  "MPO", "CSF1R", "CEBPA", # Myeloid
  "MKI67", "TOP2A"
)

# committed
suo_science_megak_ery <- c(
  "CD34", "SPINK2",
  "TESPA1", "GATA2", "FCER1A", # MegaK/Ery precursors
  "KLF1", "APOE", "FAM178B", # early ery
  "BLVRB", "CD36", "OAT", # mid ery
  "GYPA", "GYPB", "SLC4A1", "HBZ", "HBE1", # late ery
  "HBD", "PF4", # early megaK
  "ITGA2B", "ITGB3", "CLK1", # late megaK
  "MKI67", "TOP2A"
)
suo_science_b <- c(
  "CD34", "SPINK2", "IL7R", "KIT", # lymphoid prog.
  "FLT3", "CD19", "VPREB1", # pre pro B
  "MME", "CDC45", "DHFR", # pro B
  "CD27", "RAG1", "DNTT", "VPREB3", # late pro B
  "CD24", "TNFRSF17", # pro -> pre
  "IDH2", "SPIB", "IL4R", "IGHM", # pre -> immature
  "IGHD", "MS4A1", "CD40", "FCER2",  # mature
  "CD5", "SPN", "CCR10", # B1
  "MKI67", "TOP2A"
)
suo_science_mye <- c(
  "CD34", "SPINK2", "MLLT3", # early prog.
  "PRSS57", "PRTN3", "AZU1", # mye prog.
  "ELANE", "DEFA4", "LCN2", "LTF", "ORM1", # neutrophil
  "CD52", "S100A8", "MS4A6A", "CD14", "CXCR4", "CCR2", "IL1B", "CD300E", # mono
  "ACY3", "TIFAB", "KIF17", # DC prog.
  "CLEC4C", "JCHAIN", "IRF7", # pDC
  "SIGLEC6", "AXL", # ASDC
  "CLEC10A", "CD1C", # DC2
  "CLEC9A", "BATF3", # DC1
  "CCR7", "LAMP3", # migratory DC
  "IDO1", "CD207", "CD1A", # Langerhans cells
  "CLC", "KIT", "TPSAB1", # EO_BASO_Mast
  "F13A1", "LYVE1", "SPP1", # mac LYVE1 high
  "CD5L", "APOE", "VCAM1", # mac iron recycling
  "HLA-DRA", "HLA-DPA1", "CLEC7A", # mac MHCII high
  "ENG", "KDR", "CAV1", # mac kupffer like
  "TREM2", "P2RY12", # mac trem2
  "TIMD4", "FOLR2", # mac TLF+
  "MKI67", "TOP2A"
)


# MHC-I/II
# from: https://www.genenames.org/data/genegroup/#!/group/588
# rm pseudo-gene
MHC_genes = list(
  MHC1 = c("HLA-E", "HLA-F", "HLA-G"),
  MHC2 = c("HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DQB3", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB"),
  MHC2_related <- c("CD74", "CIITA")
)

# IFN-stimulated gene (ISG)
# from: https://www.nature.com/articles/s41591-018-0302-5/figures/5
isg_genes <- c(
  "ADAR", "DDX60", "HERC6", "IRF7", "OASL", "PSME2", "STAT2", "TRIM25",
  "BST2", "DHX58", "IFI35", "ISG15", "OGFR", "RSAD2", "TDRD7", "UBE2L6",
  "CASP1", "EIF2AK2", "IFIH1", "ISG20", "PARP12", "RTP4", "TRAFD1", "USP18",
  "CMPK2", "EPSTI1", "IFIT2", "MX1", "PARP14", "SAMD9L", "TRIM14",
  "CXCL10", "GBP4", "IFIT3", "NMI", "PNPT1", "SP110", "TRIM21"
)
