#!/usr/bin/env Rscript

##---------------------------------------------##
##------------------Colors---------------------##
##---------------------------------------------##

annot2col <- c(
  "HSC"="#E41A1C",
  "GP"="#E0FFFF",
  "Granulocyte"="#B3CDE3",
  "MEMP-t"="#E6AB02",
  "MEMP"="#FF7F00",
  "MEP"="#CD661D",
  "MEMP-Mast-Ery"="#FDCDAC",
  "MEMP-Ery"="#E9967A",
  "Early-Ery"="#CD5555",
  "Late-Ery"="#8B0000",
  "MEMP-MK"="#663C1F",
  "MK"="#40E0D0",
  "MastP-t"="#1E90FF",
  "MastP"="#1F78B4",
  "Mast"="#253494",
  "MDP"="#E6F5C9",
  "Monocyte"="#005A32",
  "Kupffer"="#00EE00",
  "cDC1"="#ADFF2F",
  "cDC2"="#B3DE69",
  "pDC"="#4DAF4A",
  "ASDC"="#CDC673",
  "LMPP"="#FFF2AE",
  "LP"="#FFD92F",
  "Cycling-LP"="#FFFF33",
  "PreProB"="#FFF0F5",
  "ProB-1"="#FFB5C5",
  "ProB-2"="#E78AC3",
  "Large-PreB"="#CD1076",
  "Small-PreB"="#FF3E96",
  "IM-B"="#FF00FF",
  "NK"="#A020F0",
  "ILCP"="#49006A",
  "T"="#984EA3",
  "Hepatocyte"="#666666",
  "Endothelia"="#000000",
  "LowQ"="#CCCCCC"
)

pcw2col <- c(
  "5" = "#ffbf00",
  "6" = "#ff986d",
  "7" = "#ff8782",
  "8" = "#ff7b9b",
  "9" = "#f176b4",
  "10" = "#d97fce",
  "11" = "#b989e2",
  "12" = "#9094ed",
  "13" = "#60a4f2",
  "14" = "#32b1eb",
  "15" = "#29badb",
  "16" = "#4bc0c8",
  "17" = "#8dd7db",
  "18" = "#c0dedf",
  "Missing" = "#CCCCCC",
  "Mixed" = "#CCCCCC"
)

sex2col <- c(
  "M" = "#718DBF",
  "F" = "#E84D60"
)

tissue2col <- c(
  "FL" = "#FED766",
  "CB" = "#FE4A49",
  "BM" = "#009FB7"
)

cycle2col <- c(
  "G1" = "#66C2A5",
  "S" = "#FC8D62",
  "G2M" = "#8DA0CB"
)

