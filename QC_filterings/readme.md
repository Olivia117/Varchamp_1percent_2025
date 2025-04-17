Our main pipeline for Step-0 (basic QC):

1. Load data files (CSV) per gene

2. Filter data

- Remove bulk-level data (Cell Type = "aggregated" or "cell")
- Remove small clusters (No of cells < 1000) = Remove low-expression clusters (%Cells Expressing Gene < 10% OR below 25th percentile per tissue) => moved this to _S1 (Hypothesis 1 test)
- Remove genes with extremely low absolute expression (Expression < 0.1 in >90% of clusters)
- Compute the Z-score per tissue (Saved as a new column). However, normalized expression is the most useful for us.

3. Compute scRNA-score = 𝑊1×robustness + 𝑊2×expression specificity+ 𝑊3×cluster diversity
- W1 (Robustness) = 0.4 (slightly increased, as it indicates how widely a gene is expressed)
- W2 (Expression Specificity) = 0.4 (remains high, as specificity is important for biological function)
- W3 (Cluster Diversity) = 0.2 (reduced weight, as broad cluster distribution, might not always be informative) scRNA-score=0.4×Robustness + 0.4× Expression Specificity + 0.2×Cluster Diversity where, Robustness = %Cells Expressing Gene Expression specificity = Variance-normalized expression level per tissue Cluster diversity = Number of unique Cell Types (cell types/states) per gene per tissue (normalized by max across dataset)

Here, we reproduced all the findings from previous analysis tested on a subset of genes (n=43; edgeotype responsive genes from Dr. Florent Laval). Now we run this on all the gene associated PPIs. 

Please note that POIs (e.g., AC012254.2) which are not expressed in the single-cell datasets of Cell*Gene & hence excluded from this pipeline.
