clear;


load ('20210211Exp\SRresults_200fs','xi','zi');

% load ('CTSPCamData\mb_imageEnFov50framesAverCM');
load ('CTSPCamData\mb_imageEnFov50framesAverEnb');

load CTSPPaper\CAM\Movie\dataloss_50; % regularization parameter


MBall = mb_image_tot{32};
opticalImg = imread('CTSPPaper\CAM\CAM4_Region01_crop.png');
opticalImg = imresize(opticalImg,size(MBall),'bilinear');
opticalImg = opticalImg(600+25:1200+25,350+22:700+22,:);
opticalHSV = rgb2hsv(opticalImg);
test = opticalHSV(:,:,1);
OpBin = zeros(size(test));
OpBin(test<0.1) = 1;
OpBin(test>0.9) = 1;


for imageidx = 1:size(mb_image_tot,2)
    
disp(['imageidx ',num2str(imageidx)]);    
dataloss = mb_image_tot{imageidx};
dataloss = dataloss(600:1200,350:700);

M=zeros(size(dataloss));
M(dataloss~=0)=1;

%% dataloss stat
    l2errorloss = mean((dataloss-OpBin).^2,'all');
    TPloss = zeros(size(dataloss));
    TPloss((OpBin>0)&(dataloss>0))=1;
    FPloss = zeros(size(dataloss));
    FPloss((OpBin==0)&(dataloss>0))=1;
    TNloss = zeros(size(dataloss));
    TNloss((OpBin==0)&(dataloss==0))=1;
    FNloss = zeros(size(dataloss));
    FNloss((OpBin>0)&(dataloss==0))=1;

    accuracyloss = ( sum(TPloss(:))+sum(TNloss(:)) )/(size(dataloss,1)*size(dataloss,2));
    precisionloss = sum(TPloss(:))/( sum(TPloss(:))+sum(FPloss(:)) );
    sensloss = sum(TPloss(:))/( sum(TPloss(:))+sum(FNloss(:)) );
    specloss = sum(TNloss(:))/( sum(TNloss(:))+sum(FPloss(:)) );

recovered_data(imageidx).accuracyloss = accuracyloss;
recovered_data(imageidx).precisionloss = precisionloss;
recovered_data(imageidx).sensloss = sensloss;
recovered_data(imageidx).specloss = specloss;
recovered_data(imageidx).l2errorloss = l2errorloss;





fdct_is_real = 0;
ifdct_is_real = 0;
fdct_finest = 1;
fdct_nbscales = ceil(log2(min(size(dataloss))) - 3);
fdct_nbangles_coarse = 16;


%%
C = fdct_wrapping(dataloss,fdct_is_real,fdct_finest,fdct_nbscales,fdct_nbangles_coarse); % A^T(y) forward curvelet transform
coeffdisp_original = fdct_wrapping_dispcoef(C);
nscales = length(C);

%% Find the first several largest curvelet coefficients
coeff = 0;
for i = 1:nscales
    for j = 1:length(C{i})
        coeff = max(coeff,max(C{i}{j}(:)));
    end
end

m = mean(dataloss,'all');
L1 = fitL1.a*m^fitL1.b+fitL1.c;
L2 = fitL2.a*m^fitL2.b+fitL2.c;
%% DownSampled
x = C;
u = dataloss;
counter = 1;
outerloops = 10;
innerloops = 1;

CurDecR = L1*(1:-1/outerloops:1/outerloops);
SpatDecR = L2*(1/outerloops:1/outerloops:1);

lambda = coeff(1)*CurDecR(1);
lambda2 = max(dataloss(:))*SpatDecR(1);
mu = 0.5;
tic;
while counter <= outerloops
    disp(['loop ',num2str(counter)]);
    
    coeffSum(counter) = 0;
    for j = 1:nscales
        for k = 1:length(C{j}) 
            coeffSum(counter) = coeffSum(counter) + sum(abs(x{j}{k}),'all');         
        end
    end
    l1costCur(counter) = lambda*coeffSum(counter);
    l1costSpat(counter) = lambda2*sum(abs(u),'all');
    l2cost(counter) = norm(dataloss-u);
    cost(counter) = l1costCur(counter)+l1costSpat(counter)+ l2cost(counter);
    
    for i = 1:innerloops 
        
        coeff = 0;
        temp1 = ifdct_wrapping(x,ifdct_is_real,size(dataloss,1),size(dataloss,2)).*M;
        temp2 = dataloss-temp1;
        temp3 = fdct_wrapping((temp2.*M),fdct_is_real,fdct_finest,fdct_nbscales,fdct_nbangles_coarse);
        
        for j = 1:nscales
            for k = 1:length(C{j}) 
                dummy = x{j}{k} + temp3{j}{k};
                x{j}{k} = sign(dummy).*max(0,abs(dummy) - abs(lambda));
                coeff = max(coeff,max(x{j}{k}(:)));
            end
        end
        
        u = ifdct_wrapping(x,ifdct_is_real,size(dataloss,1),size(dataloss,2));
        temp4 = u-M.*u;
        u = M.*u+sign(temp4).*max(0,abs(temp4) - abs(lambda2));
        
        [ux,uy] = gradient(u);
        tot = sqrt(sum(ux.^2+uy.^2,'all'));
        [uxx,~] = gradient(ux./tot);
        [~,uyy] = gradient(uy./tot);
        TV = uxx+uyy;
        u = u-mu*TV;
            
        x = fdct_wrapping(u,fdct_is_real,fdct_finest,fdct_nbscales,fdct_nbangles_coarse);
        
    end
    coeff = sort(coeff,'descend');
    lambda = coeff(1)*CurDecR(counter);
    lambda2 = max(u(:))*SpatDecR(counter);
    counter = counter + 1;

end
toc;

I_re = ifdct_wrapping(x,ifdct_is_real,size(dataloss,1),size(dataloss,2));
I_redb = 20*log10(abs(I_re)./max(abs(I_re(:))));


  TPRe = zeros(size(I_redb));
    TPRe((OpBin>0)&(I_redb>-60))=1;
    FPRe = zeros(size(I_redb));
    FPRe((OpBin==0)&(I_redb>-60))=1;
    TNRe = zeros(size(I_redb));
    TNRe((OpBin==0)&(I_redb<-60))=1;
    FNRe = zeros(size(I_redb));
    FNRe((OpBin>0)&(I_redb<-60))=1;

    accuracyRe = ( sum(TPRe(:))+sum(TNRe(:)) )/(size(I_redb,1)*size(I_redb,2));
    precisionRe = sum(TPRe(:))/( sum(TPRe(:))+sum(FPRe(:)) );
    sensRe = sum(TPRe(:))/( sum(TPRe(:))+sum(FNRe(:)) );
    specRe = sum(TNRe(:))/( sum(TNRe(:))+sum(FPRe(:)) );
    
    bin_recovered = zeros(size(I_redb));

    bin_recovered((I_redb>-60))=1;
    l2error = mean((bin_recovered-OpBin).^2,'all');

recovered_data(imageidx).data_recovered = I_re;
recovered_data(imageidx).accuracyRe = accuracyRe;
recovered_data(imageidx).precisionRe = precisionRe;
recovered_data(imageidx).sensRe = sensRe;
recovered_data(imageidx).specRe = specRe;
recovered_data(imageidx).l2error = l2error;
recovered_data(imageidx).l1costCur = l1costCur;
recovered_data(imageidx).l1costSpat = l1costSpat;
recovered_data(imageidx).l2cost = l2cost;
recovered_data(imageidx).cost = cost;
recovered_data(imageidx).dataloss = dataloss;        
    
end

% save ('SimulationLatest210310\MovieGenerate\mbImageRecoveredEnFov50framesAccuMu05','data_recovered_tot')
% generateMovie;
data_recovered = recovered_data(1).data_recovered;
data_recovereddb = 20*log10(abs(data_recovered)./max(abs(data_recovered(:))));
data_recovereddb(data_recovereddb<-55) = -1000;
data_recovereddb(data_recovereddb>-55) = data_recovereddb(data_recovereddb>-55)+40;
figure;imagesc(data_recovereddb,[-60,0]);colormap(gray);xticks(zeros(1,0));yticks(zeros(1,0));axis image


dataloss = recovered_data(15).dataloss;
datalossdb = 20*log10(abs(dataloss)./max(abs(dataloss(:))));
datalossdb(datalossdb<-55)=-1000;
datalossdb(datalossdb>-55) = datalossdb(datalossdb>-55)+40;
figure;imagesc(datalossdb,[-55,0]);colormap(gray);xticks(zeros(1,0));yticks(zeros(1,0));axis image
% axis image;xlabel('x (mm)');ylabel('z (mm)');

save CTSPCamData\recovered_data_Enb_50_256 recovered_data;