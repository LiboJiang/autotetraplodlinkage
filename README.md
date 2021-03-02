# Autotetraploidlinkage
A General Framework for Statistical Linkage Analysis in Autotetraplod

1.	Introduction
At present, this package is used to linkage analysis in the full-sib of autotetraploid based on the two-point and three-point model. This guide gives some brief instructions on how to perform the tasks of linkage analysis by this package. The outline of this guide is as follows: 
2.	Data format

ID                        P1   P2  IND_1  IND_10  IND_100  IND_101
solcap_snp_c2_23780       2    0     9       0       1        1
solcap_snp_c2_23781       0    2     2       1       1        1
solcap_snp_c2_23803       1    1     2       1       2        1
solcap_snp_c2_23804       0    2     1       1       1        1
solcap_snp_c2_238045      0    2     1       1       1        1
solcap_snp_c2_7632        0    2     9       1       1        1
solcap_snp_c2_23643       1    0     0       1       0        1
solcap_snp_c2_23669       3    2     9       2       2        3
solcap_snp_c2_23678       1    2     2       2       2        1
solcap_snp_c2_23717       1    1     2       1       2        1
Five genotypes (aaaa=0, Aaaa=1,AAaa=2, AAAa=3, AAAA=4) and missing data (coded as 9 ) are valid marker values. 



