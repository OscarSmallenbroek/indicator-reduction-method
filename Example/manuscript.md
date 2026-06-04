Proposed paper structure

1\.	Introduction: explaining the need

When creating indicator dashboards and/or composite indices, indicator selection often involves a gradual process, narrowing down an initially larger set into smaller subsets. Besides the usual data quality, relevance, potential policy priority considerations, statistical considerations are also employed. In case of composites, for example, a key consideration is whether the indicators within a well-defined component/dimension are sufficiently aligned/homogenous.

Especially in the case of dashboards, an additional statistical issue emerges: how much of the information content (“variation”) of the initial subset is preserved in a potential subset? While the predictive power (R2 of a regression) of the subset for the first (or first n) principal component of the full set is indicative, it does not necessarily reflect the multivariate explanatory power.

A specific example was the process of the creation of the two sustainable and inclusive wellbeing dashboards by the European Commission/JRC. From a 140-variable, comprehensive and balanced but far from small subset, a dominantly prioritisation and data quality assessment process has led to a final list of 50 indicators. It was not possible, however, to assess whether and by how much the small dashboard has remained representative for the large one dimension by dimension.

In the compilation of their beyond GDP dashboard, the UN HLEG will be facing a similar issue: it may be relatively easy to agree on a set of 50 or 100 indicators for all dimensions, but for a subset of at most 20 indicators, their degree of comprehensiveness comes up.

The paper/project aims at exploring a handful of potential options to address this issue. The first group of approaches works in the “observation domain” (country perspective): how similar are the observations (countries) if one uses all of their attributes (the full set of indicators) versus the subset only? The second group, instead, takes an “indicator” perspective: how much of the internal variation is captured (a generalisation of the principal component approach)? How much loss of predictive power one faces when using the subset instead of the full set in regressions (not for a specific LHS variable, but for a broader, even if not fully general set of LHS variables)? 

Both methods can be applied to give a number of the relative information content, so that any two potential subsets can be compared. It can also be possible (depending on the number of indicators in the full set and the subset) to select the “best n indicators”. Or, even if the best n may not meet some other criteria (policy priority, establishedness of an indicator etc), the comparison of the metric of the best and the proposed subset can indicate the goodness of the proposal.

2\.	Proposed methodological approaches 

We approach the problem of identifying structure in the dataset from two different perspectives.  First, we consider what we term the “country perspective” where we assess how the measured variables compare and contrast countries between each other.  Second, we consider what we term the “variable perspective” in which we assess what the variables have to say about each other and how well a subset of variables captures the information in the whole group.



2.1 Country Perspective Methods 



Our primary tool for assessing the differences and similarities between countries with regard to the measure variable will be a country-by-country distance matrix using the Euclidean distance between the vectors of scaled variables by country.  Countries with lower distance with regard to the full set of variables will be more similar with regard to the latent states the variables are intended to capture.  Agglomerative clustering combined with a statistical test (SIMPROF) for significant groupings by country will help determine which countries are indistinguishable based on the measure variables and, by association, the latent state.  Such methods are well-suited to datasets where the number of measured variables is high versus the number of samples.



To assess how well a reduced subset of variables preserves the information in the original difference matrix, we will create a new distance matrix just based on the subset and use rank correlation to compare the two matrices.  Rank correlation is used based on the assessment that we are more interested in determining which countries are close versus the absolute measure.  Using rank correlation as a measure of distance between distance matrices, we will also be able to compare different variable subsets to each other as well assess how the performance of subset varies over time.



Further, we will also cluster countries using the second matrix and compare significant clusters between the two measures.  Clusters will be compared by using the clusters to define an equivalence relationship and measuring the percentage of equivalency relationships in the clusters derived from the full variable set that are preserved by the relationship relative to the subset of variables, or what is known as the Rand Index (need citation).



To determine how meaningful the performance of a subset of variables is, we will compare it to the performance of a large number of random subsets of equal size.  This will give us a measure of significance.



2.2 Variable Perspective Methods 



Our approach for the variable perspective is rooted in Principle Component Analysis (PCA).  A PCA ordination yields a set of orthogonal vectors which maximize the variation between samples.  However, each vector is a linear combination of the complete set of variables, and thus no reduction in the number of variables is achieved.  However, following need citation here, we can use the PCAs to assess the amount of information captured by a single variable by looking at the magnitude of its projection onto each of the PCs weighted by the magnitude of each PC’s eigenvalue.  Thus, our measure of performance will be given by:



Equation1:  r\_m=\\sqrt{\\left(\\sum\_{i=1}^{p}{\\lambda\_i\\left(r\_m\\right)\_i^2}/(\\sum\_{j=1}^{p}\\lambda\_j)\\right)\\ } 



where \\left(r\_m\\right)\_i^2 is the squared multiple correlation of the i^{th} eigenvector with the subset of variables and \\lambda\_i is the i^{th} eigenvalue.  r\_m^2 can be interpreted as the percentage of total variation explained by the subset of variables.  When optimized for a given subset size, it should produce the subset of variables that best matches the principle component ordination.  





2.3 Applications and specific use case



The two methods outlined above can be used both diagnostically and generatively.  Diagnostically, they can be used to assess the performance of a proposed subset of variables from each perspective.  Generatively, they can be used to optimize the selection of a subset of variables of a given size with the potential addition of other constraints such as required particular variables to be included in the subset.  In particular, we will apply these methods to the sustainable and inclusive wellbeing indicator set proposed by the Joint Research Center of the European Union.  The full dataset is comprised of 141 variables for the 27 EU countries from 2011 to 2022.  Say more about the data here.



The European Commission has proposed a reduced set of 44 indicators, the result of an iterative, expert-selection process, with the goal of having the dashboard be more interpretable and data tracking demands less intensive.  It is a natural question, addressed by this work, to ask how this reduced set compares to the full set in terms of statistical representativeness and loss of information.  It is also a natural question to wonder whether an even smaller subset might do just as well or better given that 45 is a still a significant number.    



The data for the twelve years were handled two different ways.  First, for each country and each variable, we averaged the values over the twelve years to produce a data matrix with only one row per country.  We also, used a second matrix with 384 rows, one for each country for each year.  In the second case, we are not only assessing differences between countries but also across time.





3\.	Results

3.1.	 The country perspective

3.1.1 Time-averaged dataset, diagnostic application

The distance matrix for the full time-averaged dataset (FTD) and the distance matrix for the reduced 44-variable set (RTD) have a rank correlation of s = 0.908 implying a high level of correlation between the ranked distances.   However, this is not very meaningful for a subset of 44 variables.  Figure 1 shows a histogram of correlation values for 1000 random subsets of size 44.  0.908 is in the 64th percentile.  Figure 2 shows the maximum correlation achieved for 1000 random subsets of size k, showing that a subset of size 44 could have a correlation as high as 0.96 and that a correlation of 0.91 can be achieved with as few as 19 variables.





Figure 1 – Histogram of correlation values for the difference matrices derived from 1000 random variable subsets of size 44.





Figure 2 - Maximum observed correlation with the full similarity matrix for 1000 matrices derived from random subsets of variables of size k as a function of k.



Using agglomerative clustering with the SIMPROF test (α=0.05) (need citation), FTD results in 14 significant clusters for the 27 countries while RTD yielded only 11.  Figure 3 shows the dendrogram for FTD while below each sample is a symbol showing its RTD cluster.  With one exception, the clusters are generally consistent meaning that the RTD clustering is a coarser version of the FTD clusters. The Rand Index between the two sets of clustering was 95.7%.  Slovenia (SI) is an exception having been clustered separately by RTD when FTD includes it with four other countries.  The most significant loss of specificity is seen between Luxembourg, Finland, Sweden, Denmark and the Netherlands which form three separate clusters in FTD but only one in RTD.  



Figure 3 – Dendrogram from agglomerative clustering for the full time-averaged dataset.  Significant divisions based on the SIMPROF test (α=0.05) are shown with solid lines.  Countries joined by dashed red lines did not differ significantly based on the full dataset.  The symbols below show which clusters countries belong to based on the reduced dataset.



3.1.2 Time-averaged dataset, generative applications

Two exercises were conducted to try and produce a reduced dataset that preserves the distances between countries based on the whole dataset.  The indicator variables are divided into six domains which are further divided in to 19 subdomains (please check this, Peter, and maybe say some more).  It is desirable to have a subset include variables representative of each of these, and to this end, we searched for an optimal subset of variables that included either 3 variables per domain (exercise 1) or 1 variable from each subdomain (exercise 2).

Exercise 1 produced a reduced set of 18 variables which had a rank correlation of s = 0.904 which was highly significant (p \~ 0.001 when compared to random subsets of size 18).  The identified variables are reported in Table 1 along with the rank correlation scores for each subset of three variables when compared to all of the variables within the component.  In other words, this is a measure of how well the difference matrices compare just within in the given component.  

Figure 4, similar to Figure 3, shows the dendrogram derived from the clustering for the full dataset with the clusters for the 18-variable subset plotted at the bottom.  The Rand Index for this set of clusters was 82.3% relative to the full set.  The most significant difference with the clusters produced by FTD is the production of a large cluster of nine countries (blue triangles on the left) that was divided into four significant clusters by FTD.



Domain	Variables	Number of Variables in Domain	Within Domain Correlation

Wellbeing Today	"AROPE for children, 0-17, %"             	35	0.828

&#x09;  "Inability to make ends meet"             		

&#x09;  "UHC service coverage index"              		

Resources for the future	  "Financial net worth of the total economy"	18	0.689

&#x09;  "Life expectancy"                         		

&#x09;  "Spending on prevention"                  		

Societal resilience	  "Digital public services for citizens"    	36	0.805

&#x09;  "Net International Investment Position"   		

&#x09;  "Population change"                       		

Nature	"Share of forest area"                    	35	0.723

&#x09;"Net greenhouse gas emissions"            		

&#x09;"Consumption footprint per capita"        		

Inclusiveness	"Feeling discriminated"                   	8	0.885

&#x09;"Gender employment gap"                   		

&#x09;"Income quintile share ratio (S80/S20)"   		

Institutions	"Control of corruption"                   	7	0.959

&#x09;"Political stability"                     		

&#x09;"Rule of Law"         		

Table 1 – Three best variables within each domain





Figure 4 - Dendrogram from agglomerative clustering for the full time-averaged dataset.  Significant divisions based on the SIMPROF test (α=0.05) are shown with solid lines.  Countries joined by dashed red lines did not differ significantly based on the full dataset.  The symbols below show which clusters countries belong to based on the 18-variable subset.

For exercise two, the following variables were selected:

&#x20;\[1] "Severe material and social deprivation rate"                 

&#x20;\[2] "Employment rate"                                             

&#x20;\[3] "UHC service coverage index"                                  

&#x20;\[4] "Average trust in national government and national parliament"

&#x20;\[5] "Pollution, grime or other environmental problems"            

&#x20;\[6] "Financial net worth of the total economy"                    

&#x20;\[7] "Life expectancy"                                             

&#x20;\[8] "Average rating of trust"                                     

&#x20;\[9] "Government spending on health, education, social protection" 

\[10] "Digital public services for citizens"                        

\[11] "Trade openness"                                              

\[12] "Projected old-age dependency ratio"                          

\[13] "Natural and semi-natural vegetated land"                     

\[14] "Net greenhouse gas emissions"                                

\[15] "Resource productivity"                                       

\[16] "Employment in the environmental goods and services sector"   

\[17] "Income quintile share ratio (S80/S20)"                       

\[18] "Worldwide Governance Index"           

This needs to be a table, but I want to get everything else done first.

This subset had a rank correlation 0.920 which was highly significant.  The Rand Index was 88.0% which was the 13th percentile based on a random assessment of subsets of 18 variables.



3.1.3 Year x Country dataset, diagnostic application

The rank correlation between the difference matrix for the full year x country dataset (FYC) with 141 variables and the reduced dataset (RYC) with 44 variables was 0.91.  Note that each of these is a 384 x 384 matrix (27 countries x 12 years).  This was in the 89th percentile versus random subsets of size 44.  Given the large sample size, we did not apply agglomerative clustering to either matrix.



3.1.4 Year x Country dataset, generative applications

We again did an exhaustive search for the best three variables within each component, and the identified variables are in Table 2 (dataset B18).  This set had a rank correlation of s = 0.877 which was in the 99.9th percentile.





Domain	Variables	Number of Variables in Domain	Within Domain Correlation

Wellbeing Today	"AROPE for children, 0-17, %"             	35	0.832

&#x09;  "Inability to make ends meet"             		

&#x09;  "UHC service coverage index"              		

Resources for the future	  "Financial net worth of the total economy"	18	0.681

&#x09;  "Life expectancy"                         		

&#x09;  “Gross domestic expenditure on R\&D”		

Societal resilience	"SMEs with at least a basic level of digital intensity"	36	0.767

&#x09;  "Net International Investment Position"   		

&#x09;  "Trade openness"                       		

Nature	"Share of forest area"                    	35	0.697

&#x09;"Energy productivity”            		

&#x09;"Carbon footprint"        		

Inclusiveness	"Feeling discriminated"                   	8	0.894

&#x09;"Gender employment gap"                   		

&#x09;"Income quintile share ratio (S80/S20)"   		

Institutions	"Control of corruption"                   	7	0.953

&#x09;"Worldwide Governance Index"                     		

&#x09;"Rule of Law"         		

Table 2 – Variables in the reduced subset B18 representing the best variables within each domain

We also repeated the search for the best variable in each subdimension.  This yielded the following variables and had a rank correlation of 0.883.

\[1] "Severe material and social deprivation rate"                                                                         

\[2] "Young people neither in employment nor in education and training (NEETS)"                                            

&#x20;\[3] "UHC service coverage index"                                                                                          

&#x20;\[4] "Average trust in national government and national parliament"                                                        

&#x20;\[5] "Pollution, grime or other environmental problems"                                                                    

&#x20;\[6] "Total fixed assets"                                                                                                  

&#x20;\[7] "Life expectancy"                                                                                                     

&#x20;\[8] "Average rating of trust"                                                                                             

&#x20;\[9] "Difference in GINI coefficient before and after taxes and social transfers 

&#x09;(pensions excluded from social transfers)"

\[10] "Digital public services for citizens"                                                                                

\[11] "Trade openness"                                                                                                      

\[12] "Projected old-age dependency ratio"                                                                                  

\[13] "Natural and semi-natural vegetated land"                                                                             

\[14] "Net greenhouse gas emissions"                                                                                        

\[15] "Resource productivity"                                                                                               

\[16] "Employment in the environmental goods and services sector"                                                           

\[17] "Income quintile share ratio (S80/S20)"                                                                               

\[18] "Worldwide Governance Index" 



3.1.4 Change in country difference matrix across time

For three of the year x country datasets, FYC, RYC, and the best three variables by domain (B3D), we tested for the degree to which year and country were factors in determining the distance between samples.  We used non-metric multidimensional scaling to visually assess trends in distance between samples based on country and year for the dataset FYC.  The results are shown in Figure 5 where clear clustering can be seen by country as well as, for most countries, a fairly consistent trend over time (Ireland and Montenegro being two outliers in this regard).  The overall trend is consistent with the other two datasets (results not shown).



Figure 5 – Non-metric multidimensional scaling ordination of the 384 year x country samples.  Samples are coded by their country, and arrows connect consecutive years within each country.  Most countries cluster heavily and have a clear time trend.

To compare how strong these trends were in each dataset, we ran ANOSIM with a two-way crossed setup with Year as an ordered factor.  The results are in Table 3.

DATASET	ANOSIM R FOR YEAR	ANOSIM R FOR COUNTRY

FYC	0.890	0.921

RYC	0.910	0.947

B18	0.871	0.929

Table 3 – ANOSIM results for the factors Year and Country for the full FYC dataset and the two reduced datasets, RYC and B18.

For all three datasets, both Country and Year were very strong determinants of distance after controlling for the other.  There is clear consistency in the strength of these factors amongst all three datasets.

To further assess the trend over time, we subsetted each dataset by Year and created separate distance matrices for each year.  These matrices were then compared pairwise based on the rank correlation of distances between countries and the resulting correlations were used to create a distance matrix for the 36 distance matrices (12 years x 3 datasets).  We were especially interested in comparing the first (2011) and last (2022) distance matrices for each dataset to determine how much stability there was in the distance matrix representations.

For FYC, the correlation between 2011 and 2022 was s = 0.846 suggesting some change over time.  For B18, the coefficient was s = 0.815, and for RYC it was s = 0.897.  Overall, the worst correlation was between the 2011 matrix for RYC and the 2021 matrix for FYC (s = 0.726).  The nMDS plot derived from the distance matrix of matrices (Figure 6) shows clearly that the three datasets are distinct, but it also shows that the reduction in correlation seen when going from 2011 to 2022 within one dataset is less than the loss in going between datasets.  Further, the overall time trend seen in FYC is generally mirrored by the trends in B18 and RYC.

Figure 6 – nMDS representation of the difference matrix comparing 36 distance matrices derived from the three datasets, FYC, RYC, and B18, subsetted by each of the 12 years from 2011 – 2022.  A minimum spanning tree has been overlaid based on the distance matrix which clearly shows that the apparent separation between datasets and the gradient across years are consistent with the distance matrix.



This consistency is also seen in the dendrogram in Figure 7.  (Note that because the distance matrix is not derived from an underlying dataset, we cannot use SIMPROF to determine cluster significance.)



Figure 7 – Dendrogram based on the matrix distance matrix.  Note that it clearly shows that the matrices derived from the NARROW dataset are generally closer to the FULL dataset than are the BEST matrices, and the difference in average similarity is approximately 0.04.  The data for 2021 in the FULL dataset are a bit of an outlier.

Again, this suggests a high degree of consistency between the three datasets both within countries and across years.  And given the significant variable reduction achieved by BEST, it suggests this is an efficient set of variables for comparing countries based on the larger datasets.



&#x09; The variable perspective 

We used our variable-perspective criterion both to generate potential subsets of variables that capture much of the structure for both the time-averaged dataset and the country by year dataset.  We also used it diagnostically to assess the reduced country by year (RYC) dataset.  It is not a meaningful criterion for assessing the RTD because any variable subset with more variables than samples with have r\_m=1 unless the set is linearly dependent which is highly unlikely.  This is because the space spanned by the set will be identical to the space spanned by the PCs.



3.2.1 Diagnostic application of the variable perspective 

The RYC dataset has r\_m=0.958 which means these 44 variables account for 91.9% of the variability.  The first 44 PCs account for 97.7% of the variability implying that RYC achieves approximately 94% of the optimum.  This is approximately the 50th percentile compared to random subsets of size 44 meaning that, in general, subsets of that size are usually sufficient for capturing most of the variance. 



3.2.2 Generative application of the variable perspective 

For both the time-averaged dataset and the country by year dataset, we identified the best three-variable subset within each domain for capturing the variability within that domain (not the entire dataset).  The criterion for evaluating a subset was identical to Equation 1 except the PCs are for the reduced dataset composed only of variables within the domain.  The results are in Table 3.  For all 18 variables, r\_m=0.940.  The first 18 PCs capture 94.2% of the variability suggesting that this is a highly effective subset for representing the variability within the full dataset.  No random subset of size 18 was seen to have a greater r\_m.



Domain	Variables	Number of Variables in Domain	Within Domain Correlation

Wellbeing Today	 "Severe material and social deprivation rate"                	35	0.746

&#x09; "Long-term unemployment rate"                                		

&#x09; "Standardised death rate due to homicide"                    		

Resources for the future	 "Gross domestic expenditure on R\&D"                          	18	0.730

&#x09; "Lower-secondary completion only"                            		

&#x09; "Tertiary education attainment (25-34)"                      		

Societal resilience	 "Standardised preventable and treatable mortality (low rate)"	36	0.706

&#x09; "Digital public services for citizens"                       		

&#x09; "Concentration of value chain partners"                      		

Nature	 "GHG emissions intensity of the economy"                     	35	0.638

&#x09; "Net greenhouse gas emissions"                               		

&#x09; "Share of renewable energy in gross final energy consumption"		

Inclusiveness	 "Disability employment gap"                                  	8	0.902

&#x09; "Gender employment gap"                                      		

&#x09; "Income quintile share ratio (S80/S20)"                      		

Institutions	 "Political stability"                                        	7	0.991

&#x09; "Regulatory quality"                                         		

&#x09; "Worldwide Governance Index"  		

Table 3 – Best three variables in each domain for capturing variability within the domain.

We also searched for the best singular variable within each domain, producing the following set of 18 variables which had r\_m=0.900.  Again, no random subset of size 18 had a greater r\_m (maximum out of 1000 was 0.871).

&#x20;\[1] "Severe material and social deprivation rate"                 

&#x20;\[2] "Employment rate"                                             

&#x20;\[3] "UHC service coverage index"                                  

&#x20;\[4] "Average trust in national government and national parliament"

&#x20;\[5] "Pollution, noise, grime or other environmental problems"     

&#x20;\[6] "Gross domestic expenditure on R\&D"                           

&#x20;\[7] "Children aged less then 3 years in formal childcare"         

&#x20;\[8] "Participation in voluntary activities"                       

&#x20;\[9] "Government spending on health, education, social protection" 

\[10] "Digital public services for citizens"                        

\[11] "Trade openness"                                              

\[12] "Projected old-age dependency ratio"                          

\[13] "Share of forest area"                                        

\[14] "Net greenhouse gas emissions"                                

\[15] "Resource productivity"                                       

\[16] "Gross added value of environmental goods and services sector"

\[17] "Income quintile share ratio (S80/S20)"                       

\[18] "Worldwide Governance Index"      



We repeated the search within each component for the year by country dataset.  The results are in Table 4.  For this subset, r\_m=0.850 versus a total variance explained by the first 18 PCs of 0.857, again suggesting a performance very close to the optimum.



Domain	Variables	Number of Variables in Domain	Within Domain Correlation

Wellbeing Today	 "Severe material and social deprivation rate"                	35	0.721

&#x09; "Long-term unemployment rate"                                		

&#x09; "Standardised death rate due to homicide"                    		

Resources for the future	 "Gross domestic expenditure on R\&D"                          	18	0.686

&#x09; "Lower-secondary completion only"                            		

&#x09; "Tertiary education attainment (25-34)"                      		

Societal resilience	 "Adult participation in learning"                            	36	0.650

&#x09; "Standardised preventable and treatable mortality (low rate)"		

&#x09; "Trade openness"                                             		

Nature	 "Primary energy consumption"                                 	35	0.603

&#x09; "Share of renewable energy in gross final energy consumption"		

&#x09; "Resource productivity"                                      		

Inclusiveness	 "Disability employment gap"                                  	8	0.896

&#x09; "Feeling discriminated"                                      		

&#x09; "Income quintile share ratio (S80/S20)"                      		

Institutions	 "Political stability"                                        	7	0.986

&#x09; "Regulatory quality"                                         		

&#x09; "Worldwide Governance Index"  		

Table 4 – Best three variables in each domain for capturing variability within the domain for the year by country dataset.

Finally, we broke the data out by year and repeated the search for the best three variables by year.  The resulting subsets of variables were very consistent in explaining most of the total variation with the mean r\_m=0.932 and a range from 0.924 to 0.936.  However, this level was only significant for 4 out of 11 years compared to random subsets of 18 variables.

Peter:  How to describe the variation in variables?



&#x09;Conclusions and open issues



