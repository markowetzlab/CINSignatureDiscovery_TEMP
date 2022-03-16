#### Summarise multiple runs of BayesNMF

# You can run this script and it will simply chose the best solution according to the 
# Kullback-Leibler divergence (if option "div" chosen). Or you can manually select the most 
# representative solution based on visual inspection of network graphs of all solutions with
# optimal number of K's (rank of NMF, number of signatures). For the manual inspection, run this
# script until line 240 (ish, have a look) and then run this function to see the network graph
# with the currently chosen solution highlighted:
# testSolution(<<solution ID>>) 
# If you want this solution, simply transfer this solution to this variable and continue script:
# bestSolution = lData[[ <<solution ID>> ]]

rm(list=ls(all=TRUE))

## Options
# BASE="/Users/drews01/phd/prjcts/cnsigs2/cnsigs2_revisions"
# BASE="/home/drews01/cnsigs2_revisions"
BASE="/Users/drews01/cnsigs2_revisions"
OUT=file.path(BASE, "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PATHTOFILES=file.path(BASE, "out/5_NMF_CxS_K0-43_TCGA_WES_Two_callers_CNVkit_Standard")
SIGMATPATTERN="W.txt"
EXPMATPATTERN="H.txt"
LOGPATTERN="log.txt"
# See Alexandrov et al., 2018 (bioRxiv)
corThreshold = 0.85
# If null, then mode of K is chosen, otherwise the x-most frequent K is chosen:
WHICHK=NULL
METADATA=file.path(BASE, "data/ASCAT_metadata_TCGA_PCAWG_with_dCIN.txt")
CANCERCOLS=file.path(BASE, "data/TCGA_colour_scheme_Lydia.txt")
OUTPUTDIR=PATHTOFILES
dir.create(OUTPUTDIR, showWarnings = FALSE, recursive = TRUE)
TRANSPOSED=FALSE
# Can be either "cost" or "div"
DECISION="div"


library(data.table)
library(ComplexHeatmap)
library(ggplot2)
library(igraph)
library(qgraph)


# From: https://stats.stackexchange.com/questions/97051/building-the-connection-between-cosine-similarity-and-correlation-in-r
# or just use:
# library(lsa)
# cosine(matrix)
cosineSimForMatrices <- function(X, corr=FALSE){
    if(corr){ X = apply(X, 2, function(x){ x-mean(x) }) }
    denom = solve(diag(sqrt(diag(t(X)%*%X))))
    return( denom%*%(t(X)%*%X)%*%denom )
} 

plotCosineHeat = function(lData, allSamples, thisK, corThreshold = 0.9) {
    
    lFilt = list()
    for(thisSample in allSamples) {
        thisMat = lData[[ thisSample ]]
        if(nrow(thisMat) == thisK) {
            rownames(thisMat) = paste0( thisSample, "_", rownames(thisMat) )
            lFilt[[ thisSample ]] = thisMat
        }
    }
    
    
    matFilt = do.call(rbind, lFilt)
    corMat = cosineSimForMatrices(t(matFilt), corr = FALSE)
    corMat[ corMat < corThreshold ] = 0
    return( Heatmap(corMat, show_row_names = FALSE, show_column_names = FALSE, name = paste("K =", thisK)) )
    
}

makeGraphLayout = function(corMat, bestSample, thisK) {
    net = graph_from_adjacency_matrix(corMat, weighted = TRUE, mode = "undirected", diag = FALSE)
    
    # Get samples names
    theseSamples = sapply(rownames(corMat), function(thisSig) { strsplit(thisSig, "_")[[1]][1] })
    V(net)$Samples = theseSamples
    
    V(net)$color <- "#FFFFFF"
    V(net)$color[ V(net)$Samples == bestSample ] = "#228B22"
    
    # make TCGA larger
    V(net)$size = 3
    V(net)$size[ V(net)$Samples == bestSample ] = 10
    
    # label only TCGA
    V(net)$label = NA
    V(net)$label[ V(net)$Samples == bestSample ] = seq(1:thisK)
    
    # layout for spacing out vertices
    e = get.edgelist(net, names=FALSE)
    l = qgraph.layout.fruchtermanreingold(e, vcount = vcount(net),
                                          area=8*(vcount(net)^2),
                                          repulse.rad=(vcount(net)^3.1))
    
    lOut = list(net = net, l = l)
    return(lOut)
}

makeGraph = function(net, l) {
    plot(net, edge.width = 1, vertex.label.font = 2, vertex.label.color="black", layout=l)    
}

## For manual solution chosing
testSolution = function(bestSample) {
    
    bestSolution = lData[[ bestSample ]]
    
    # For later
    dtScores=data.table(do.call(rbind, lScores))
    dtScores$V2=as.numeric(dtScores$V2)
    dtScores=dtScores[order(dtScores$V2, decreasing = FALSE),]
    colnames(dtScores) = c("Sample", as.character(DECISION))
    
    # Hairball plot with highlighted best solution
    
    lFilt = list()
    for(thisSample in allSamples) {
        thisMat = lData[[ thisSample ]]
        if(nrow(thisMat) == thisK) {
            rownames(thisMat) = paste0( thisSample, "_", rownames(thisMat) )
            lFilt[[ thisSample ]] = thisMat
        }
    }
    
    matFilt = do.call(rbind, lFilt)
    corMat = cosineSimForMatrices(t(matFilt), corr = FALSE)
    corMat[ corMat < corThreshold ] = 0
    rownames(corMat) = rownames(matFilt)
    colnames(corMat) = rownames(matFilt)
    
    # Convert matrix to graph
    lGraph = makeGraphLayout(corMat, bestSample, thisK)
    
    # Plotting heatmap of best solution
    plotHeatBest = Heatmap(bestSolution, name = bestSample, cluster_rows = FALSE, cluster_columns = FALSE)
    makeGraph(lGraph$net, lGraph$l)
}


# Use faster version of hclust for heatmaps (not needed as performance increase was small)
# library(fastcluster)
# ht_global_opt(fast_hclust = TRUE)

allSigFiles = list.files( path = PATHTOFILES, pattern = SIGMATPATTERN, recursive = TRUE, full.names = TRUE )
allExpFiles = list.files( path = PATHTOFILES, pattern = EXPMATPATTERN, recursive = TRUE, full.names = TRUE )
allLogs = list.files( path = PATHTOFILES, pattern = LOGPATTERN, recursive = TRUE, full.names = TRUE )


## Reading in all heatmaps
allSamples = sapply(allSigFiles, function(x) strsplit( basename(x), "_")[[1]][1] )
lData = lapply(allSamples, function(thisSample) {
    
    print(thisSample)
    thisW = fread( allSigFiles[ grep(thisSample, allSigFiles) ] )
    theseRownames = thisW$V1
    thisW$V1 = NULL
    thisSigs = as.matrix(thisW)
    rownames(thisSigs) = theseRownames
    
    # For easier downstream analysis, transpose the matrix, so the signatures are in the rows
    if( TRANSPOSED ) {
        thisSigs = t( apply(thisSigs, 1, function(x) x/sum(x)) )
        return( thisSigs )
    } else {
        return( t(thisSigs) )    
    }
    
        
} )

names(lData) = allSamples


## Histogram of Ks
dfKs = data.frame(sapply(lData, function(x) nrow(x)))
colnames(dfKs) = "K"
plotHist = ggplot(dfKs, aes(x=K)) + geom_histogram() + theme_bw() + ylab("Frequency")


## How close are the signatures of the same K?
## Actual code:
thisK = as.numeric( names(table(dfKs))[ which.max( table(dfKs) ) ] )
## Hacked for WP1 (setting K=10 to compare to pancancer sigs)
# thisK = 10
plotHeat1 = plotCosineHeat(lData, allSamples, thisK, corThreshold = corThreshold)
# Determine second K
allK = names(sort(table(dfKs), decreasing=T))
# 
# if(length(allK) >= 2) {
#     thisSecondK = as.numeric( allK[2] )
#     plotHeat2 = plotCosineHeat(lData, allSamples, thisK = thisSecondK, corThreshold = corThreshold) 
# }
# 
# 
# # Determine third K
# if(length(allK) >= 3) {
#     thisThirdK = as.numeric( allK[3] )
#     plotHeat3 = plotCosineHeat(lData, allSamples, thisK = thisThirdK, corThreshold = corThreshold) 
# }
# 
# 
# ## Add manual step for deciding which K to choose
# if (! is.null(WHICHK)) {
#     print(paste("Go for this K:", WHICHK))
#     # thisK = as.numeric(WHICHK)
#     thisK = as.numeric( allK[as.numeric(WHICHK)] )
# }


## Loop over samples from K=9 and chose solution with best beta_divergence.
lScores=list()
bestCostFun = Inf
for(thisSample in allSamples) {

    thisMat = lData[[ thisSample ]]
    if(nrow(thisMat) == thisK) {
        thisLog = fread( allLogs[ grep(pattern = thisSample, x = allLogs) ] )
        
        if( DECISION == "div" ) {
            thisCost = as.numeric( strsplit( thisLog[ nrow(thisLog) ]$V3, "=")[[1]][2] )
            lScores[[length(lScores)+1]] = c(thisSample, thisCost)
            #print(thisCost)
            if(thisCost < bestCostFun) {
                print(paste("bestSample:", thisSample,"| bestDiv:", thisCost))
                bestSample = thisSample
                bestCostFun = thisCost
                }    
        } else {
            thisCost = as.numeric( strsplit( thisLog[ nrow(thisLog) ]$V2, "=")[[1]][2] )
            #print(thisCost)
            if(thisCost < bestCostFun) {
                print(paste("bestSample:", thisSample,"| bestCost:", thisCost))
                bestSample = thisSample
                bestCostFun = thisCost
            }    
        }
    }
}

# For later
dtScores=data.table(do.call(rbind, lScores))
dtScores$V2=as.numeric(dtScores$V2)
dtScores=dtScores[order(dtScores$V2, decreasing = FALSE),]
colnames(dtScores) = c("Sample", as.character(DECISION))

# Save all suitable solutions and their value
write.table(dtScores, file.path(OUTPUTDIR, paste0("6_summary_BayesNMF_runs_all_Solutions_", bestSample, "_K", thisK, ".txt")), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


# Extract best solution
bestSolution = lData[[ bestSample ]]



###########################################
#### At this point manual inspection of solutions possible by running:
# testSolution(<<solution ID>>) 
# testSolution(bestSample)
# testSolution(dtScores$Sample[1])
#### If settled on a best solution, transfer the best solution to bestSolution and continue script.
###########################################


## Failsafe
bestSolution = lData[[ bestSample ]]

# Hairball plot with highlighted best solution
lFilt = list()
for(thisSample in allSamples) {
    thisMat = lData[[ thisSample ]]
    if(nrow(thisMat) == thisK) {
        rownames(thisMat) = paste0( thisSample, "_", rownames(thisMat) )
        lFilt[[ thisSample ]] = thisMat
    }
}

matFilt = do.call(rbind, lFilt)
corMat = cosineSimForMatrices(t(matFilt), corr = FALSE)
corMat[ corMat < corThreshold ] = 0
rownames(corMat) = rownames(matFilt)
colnames(corMat) = rownames(matFilt)

# Convert matrix to graph
lGraph = makeGraphLayout(corMat, bestSample, thisK)

# Plotting heatmap of best solution
plotHeatBest = Heatmap(bestSolution, name = bestSample, cluster_rows = FALSE, cluster_columns = FALSE)
# makeGraph(lGraph$net, lGraph$l)


## Exposures of bestSample
bestExpRaw = fread( allExpFiles[ grep( bestSample, allExpFiles ) ])
expRownames = bestExpRaw$V1
bestExpRaw$V1 = NULL
matBestExpRaw = as.matrix(bestExpRaw)
rownames(matBestExpRaw) = expRownames
if( TRANSPOSED ) {
    matBest = t( apply(matBestExpRaw, 1, function(x) x/sum(x)) )
} else {
    matBest = t( apply(matBestExpRaw, 2, function(x) x/sum(x)) )
}
plotBestExp = Heatmap(t(matBest), show_row_names = TRUE, show_column_names = FALSE)


## Exposures by cancer type
metadata = fread(METADATA)
cancerCols = read.table(CANCERCOLS, header = FALSE, comment.char = "|")
# Filter and sort like the exposure matrix
metadata = metadata[ metadata$name %in% rownames(matBest), ]
metadata = metadata[ match(rownames(matBest), metadata$name), ]
rownames(metadata) = metadata$name

md2 = data.table( Cancer = factor(metadata$cancer_type) )
vCancerCols = as.character( cancerCols$V2[ match(sort(unique(md2$Cancer)), cancerCols$V1) ] )
names(vCancerCols) = cancerCols$V1[ match(sort(unique(md2$Cancer)), cancerCols$V1) ]

numCancers = length(unique(md2$Cancer))
if(numCancers > 1) {
    plotBestExpCancer = Heatmap(matBest, split = md2[,1], cluster_columns = FALSE, show_row_names = FALSE, show_column_names = TRUE) + 
        rowAnnotation(df = md2, col = list( Cancer = vCancerCols ))
}



## Boxplots of exposures
dtBestExp = as.data.table(matBest)
library(reshape2)
mBestExp = melt(dtBestExp)
plotBoxBestExp = ggplot(mBestExp, aes(x = variable, y = value)) + geom_boxplot()

## Combine all previous plots
# Histogram of Ks
# Heatmap of cosine similarity of K=12 and K=13 solutions
# Hairball of cosine similarity of all solutions with the best one highlighted
# Heatmap of exposures of best solution
# Boxplot of signature exposures


pdf(file = file.path(OUTPUTDIR, paste0("6_summary_BayesNMF_runs_", bestSample, "_K", thisK, ".pdf")), width = 12, height = 8)
plotHist
plotHeat1
if(length(allK) >= 2) plotHeat2
if(length(allK) >= 3) plotHeat3
makeGraph(lGraph$net, lGraph$l)
plotHeatBest
plotBestExp
if(numCancers > 1) plotBestExpCancer
plotBoxBestExp
dev.off()

## Move best solution to output directory
file.copy( allSigFiles[ grep(bestSample, allSigFiles)], file.path(OUTPUTDIR, paste0("6_Signatures_", bestSample, ".txt")))
file.copy( allExpFiles[ grep(bestSample, allExpFiles)], file.path(OUTPUTDIR, paste0("6_Exposures_", bestSample, ".txt")))

# Save signature matrix (normalised by default in BayesNMF algorithm)
saveRDS(object = bestSolution, file = file.path(OUTPUTDIR, paste0("6_Signatures_", bestSample, "_normalised.rds")))

# Save exposure matrix
saveRDS(object = matBestExpRaw, file = file.path(OUTPUTDIR, paste0("6_Exposures_", bestSample, ".rds")))
saveRDS(object = matBest, file = file.path(OUTPUTDIR, paste0("6_Exposures_", bestSample, "_normalised.rds")))
