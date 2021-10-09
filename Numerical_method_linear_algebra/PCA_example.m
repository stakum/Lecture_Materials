clear all; close all; clc;
%--------------------------------------------------------------------------
% Example: Principal Component Analysis (PCA)
%--------------------------------------------------------------------------

DATA=[86 79 67 68; 
   71 75 78 84;
   42 43 39 44; 
   62 58 98 95;
   96 97 61 63; 
   39 33 45 50;
   50 53 64 72; 
   78 66 52 47;
   51 44 76 72; 
   89 92 93 91];

R=corrcoef(DATA);

[coeff,score,latent] = pca(zscore(DATA));

fl=sqrt(latent).*coeff';