## @knitr exploreGene

palette(c("black", "tan", "red"))
geneID <- "Ndufa12"

eQTLB6[geneID, ]
eQTLBXD[geneID, ]

variantID <- eQTLB6[geneID, "IDvariant"]

# plot variant-gene expression pair
plotGeneExpressionVariant(ref="B6", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")
plotGeneExpressionVariant(ref="BXD", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")

# plot variant-gene count pair
plotGeneCountsVariant(ref="B6", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")
plotGeneCountsVariant(ref="BXD", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")

if(eQTLB6[geneID, "IDvariant"]!=eQTLBXD[geneID, "IDvariant"]){
  
  variantID <- eQTLBXD[geneID, "IDvariant"]
  
  # plot variant-gene expression pair
  plotGeneExpressionVariant(ref="B6", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")
  plotGeneExpressionVariant(ref="BXD", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")
  
  # plot variant-gene count pair
  plotGeneCountsVariant(ref="B6", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")
  plotGeneCountsVariant(ref="BXD", gene=geneID, variant=variantID, tissue="Cortex", condition="NSD")
}
