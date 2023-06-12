clear;

%% Load Data
load data\SimVesselData;
%% Define fdct Parameters
fdctPara.fdct_is_real = 0;
fdctPara.ifdct_is_real = 0;
fdctPara.fdct_finest = 1;
fdctPara.fdct_nbscales = ceil(log2(min(size(dataLoss))) - 3);
fdctPara.fdct_nbangles_coarse = 16;
fdctPara.M = size(dataLoss,1);
fdctPara.N = size(dataLoss,2);

%% Display Curvelet domain
figure;imshow(dataLoss);colormap(gray)
fdctDisp(dataLoss, fdctPara)

%% Iteration 

iterPara.fitL1 = fitL1;
iterPara.fitL2 = fitL2;
iterPara.outerloops = 20;
iterPara.innerloops = 1;
iterPara.mu = 0.5;
reconImage = iterateFunc(dataLoss, fdctPara, iterPara);

%% Display Recon Image
reconImagedB = 20*log10(reconImage./max(reconImage(:)));
figure; imshow(reconImagedB,[-100, 0]); colormap(gray)
fdctDisp(reconImage, fdctPara)