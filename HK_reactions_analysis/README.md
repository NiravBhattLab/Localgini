# Details of the code files in the folder

1) getCoreRxns.m: To get core reactions for the cancer cell line and HPA data using three different thresholds and three different GEMs (Recon 2.2, Recon3D and Human1)
    - The files generated using this code is stored in core_rxns_cancer_data.mat and core_rxns_tissue_data.mat
2) getHKRxns.m: Maps the housekeeping genes from the hkg_sym.txt file to the reactions in the GEMS - Recon2.2,Recon3D and Human1
    - The files generated using this code is stored in HK_rxns_Recon2_2.mat, HK_rxns_Recon3D.mat and HK_rxns_iHuman.mat
3) getHKRxnsInCoreCount.m: To get the count of housekeeping reactions rectified in each of the context across three distinct GEMs and three distinct thresholding methods
    - The files generated using this code is stored in HK_in_core_cancer.mat and HK_in_core_tissue.mat. These files are further used to generate the Figure 2B
4) getHKspecificCount.m: To get the count of contexts on which a particular housekeeping reaction is captured by the thresholding methods
    - The files generated using this code is stored in HK_tissue_Countdata.xlsx. This file is further used to generate the Figures S2, S3 and S4
