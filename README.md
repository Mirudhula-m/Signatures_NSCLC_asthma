# Key Signatures between Non-small Cell Lung Cancer and Severe Asthma


Common perception indicates smoking to be the likely reason for lung cancer but that does not explain for the 25% of the cases attributed to non-smoking related lung cancer. We identified several mate-analysis linking lung cancer to severe cases of asthma. So, we sought to find common molecular signatures that might be existing between the two most common yet deadly diseases. 

We performed a Differential Gene Analysis (DGE) for datasets corresponding to Non-small cell lung cancer (NSCLC) and Severe Asthma independently against a cohort of normal patients. The data was extracted and filtered from Gene Expression Omnibus (GEO). The result was used to perform a Gene Set Enrichment Analysis (GSEA) which will help us determine the statistical significance of one gene set being a part of the DGE results of another. We found that over-expressed genes corresponding to NSCLC were being overrepresented in the asthma dataset.

![enplot_LC_GENE_SET_13.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/de422c31-e7db-4f21-82d4-148ef35a0456/enplot_LC_GENE_SET_13.png)

We further identified the top genes that were responsible for this enrichment and found that these were being differentially expressed in both severe asthma and lung cancer. We used Decision Classification tree for validating this result with another pair of datasets from GEO, and found the best 8 genes that was able to classify these diseases.

![Loss_curve_Dtrees.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/2f7ea4ce-997e-4f72-8a36-186288ef0eaa/Loss_curve_Dtrees.png)

Out of these 8 genes, all of which seemed to have a relationship with either disease in literature, PPARD had a unique aspect. This gene is under-expressed in asthma and over-expressed in NSCLC. But, it was observed that patients with mixed cases of asthma and lung cancer had an increased expression. We hypothesize that this gene, and possibly other genes, could be used as a biomarker for early detection of lung cancer in severe asthma patients.
