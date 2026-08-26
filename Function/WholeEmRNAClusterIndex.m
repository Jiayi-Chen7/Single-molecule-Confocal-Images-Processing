clc
clear all
close all



%initial parameters
%folder
% in_folder = 'E:/CJY/20220314/EM1/';
% mid_sur_in_folder = 'E:/CJY/20220314/EM1/MID SUR/';
% out_folder =  'E:/CJY/20220314/EM1/output/';
%the image names of zstack
dapi_channel = '20250424def 9632 hb SunS1 Em9.lif_Em9Large_z*_ch00.tif';
gfp_channel = '20250424def 9632 hb SunS1 Em9.lif_Em9Large_z*_ch02.tif';
mRNA_channel = '20250424def 9632 hb SunS1 Em9.lif_Em9Large_z*_ch01.tif';
% image list
imlist_dapi = dir(dapi_channel);
imlist_gfp = dir(gfp_channel);
imlist_mRNA = dir(mRNA_channel);
% nuclei segmentation parameters
cir_th = 0.75;
low_th = 1500;
ExVal = 3;%3-nc14
IxVal = 40;%30-nc14
Nuc_th =15000;
%the size of the image in one z step
size_in_merge_image = [2048,2048];
size_in_plane_image = [1024,1024];
% low_th = 1000;    %%% lower area limit for circular mask algorithm
% low_th = 100;    %%% lower area limit for circular mask algorithm
high_th = 50000;   %%% higher area limit for circular mask algorithm
CountScale = 1;
%register parametes
HighMag = 63;
LowMag = 20;
ZoomFactor = 1.3;
RegScaleFactor = LowMag./(HighMag.*ZoomFactor);
TrueNucVol = 20000;% 3D volume threshold to remove the yolk nuclei


%% read in raw data

dapi_raw_3D = zeros(size_in_merge_image(1),size_in_merge_image(2),numel(imlist_dapi));
gfp_raw_3D = zeros(size_in_merge_image(1),size_in_merge_image(2),numel(imlist_gfp));
mRNA_raw_3D =zeros(size_in_merge_image(1),size_in_merge_image(2),numel(imlist_mRNA));

parfor image_I = 1:length(imlist_dapi)
    dapi_raw_3D(:,:,image_I) = imread([imlist_dapi(image_I).name]);
    gfp_raw_3D(:,:,image_I) = imread([imlist_gfp(image_I).name]);
    mRNA_raw_3D(:,:,image_I) = imread([imlist_mRNA(image_I).name]);
end


%%image mask of nuclei in high magnitude image
RefineZ = (numel(imlist_dapi));
Zstart = 28;
Zend = 38;

DapiRawMask = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
GFPRawMask = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
mRNARawMask =zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
mRNAFociRawMask = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
mRNAGauFil = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
GFPGauFil = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
thresholdvalue = zeros(1,RefineZ);
%estimate the nuclei cycle of the embryo
NCEst = 14;
grayth = 1;
graythAdj = 0;

tic 
fig = 0;
parfor ImageMaskIndex = 1:RefineZ
EmSingleOri = dapi_raw_3D(:,:,ImageMaskIndex);
% [Mask T] = Nuclei_seg_cjy(EmSingleOri, ExVal, IxVal,low_th,cir_th,ImageMaskIndex,NCEst);
% [Mask] = Nuclei_seg_cjy(EmSingleOri,5, 25,[],0.6,ImageMaskIndex,NCEst,500,5,fig,5);
[Mask] = NucleiSeg(EmSingleOri,2,[],200,0.3,ImageMaskIndex,NCEst,[],10,fig,10,grayth,graythAdj,6500);
DapiRawMask(:,:,ImageMaskIndex) = Mask;
end
toc
% % % 
% % % 
% % % 
% for j = 1:RefineZ
% figure();
% [B,L] = bwboundaries(DapiRawMask(:,:,j));
% imshow(dapi_raw_3D(:,:,j),[]);
% hold on
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
%    hold on
% end
%      hold off     
% end
% % %      saveas(gcf, sprintf('RawNucleiWithMask%d.tif',j));
% % % end

%remove yolk nuclei and refine
BWDapiRawMask = bwlabeln(DapiRawMask);
StsDapiRawMask = regionprops3(BWDapiRawMask,dapi_raw_3D(:,:,1:RefineZ),'Volume','MeanIntensity','Centroid');
NucVol = [StsDapiRawMask.Volume];
NucMeanInt = [StsDapiRawMask.MeanIntensity];
NucCenRaw = [StsDapiRawMask.Centroid];
figure();
histogram(NucMeanInt);
%the threshold is estimate from the prior knowledge of the size of a nuclei
%the diameter of ~ 5um: 5 ./ 0.069 um/pixel ~ 72
%should be observed in at least 5 layers of the z stacks
TrueNucVolLB = 30000;%20000
TrueNucVolHB = 150000;%20000
NucIntenThLB = 12000;
NucIntenThHB = 17000;
TrueNucVolLB = 0;%20000
TrueNucVolHB = 900000;%20000
NucIntenThLB = 12000;
NucIntenThHB = 20000;
TrueNucVol2 = 0;%20000
NucIntenTh2 = 0;
TrueNucInd = find(NucVol > TrueNucVol2& NucMeanInt>NucIntenTh2);
Nuc_true = ismember(BWDapiRawMask,TrueNucInd);
LLNuc_true = bwlabeln(Nuc_true);
stsNuc_true = regionprops3(LLNuc_true,dapi_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid');
NuccVol = [stsNuc_true.Volume];
NucWeiCen = [stsNuc_true.WeightedCentroid];
% length(NuccVol);
%shrink the muclei mask for more accurate nuclei assignment of the mRNA and
%GFP spot
ShrinkNucPara = 5;
se = strel('cube',ShrinkNucPara);
ErodeNuc_true = imerode(Nuc_true,se);
LLErodeNuc_true = bwlabeln(ErodeNuc_true);
stsLLErodeNuc_true = regionprops3(LLErodeNuc_true,dapi_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid');
ErodeNucVol = [stsLLErodeNuc_true.Volume];
ErodeNucWeiCen = [stsLLErodeNuc_true.WeightedCentroid];


load('RNAMat_EM9_20250424_hb_nc14_GCN4.mat');
NorNucMeanRNAClusterIndex = zeros(length(ErodeNucWeiCen),1);
NucMeanRNADistance = zeros(length(ErodeNucWeiCen),1);
load('OnRNAMat_EM9_20250424_hb_nc14_GCN4.mat');
NorNucMean_On_RNAClusterIndex = zeros(length(ErodeNucWeiCen),1);
NucMeanRNA_On_Distance = zeros(length(ErodeNucWeiCen),1);
load('OFFRNAMat_EM9_20250424_hb_nc14_GCN4.mat');
NorNucMean_Off_RNAClusterIndex = zeros(length(ErodeNucWeiCen),1);
NucMeanRNA_Off_Distance = zeros(length(ErodeNucWeiCen),1);
TrueDistanceScale = 69;

parfor i = 1:length(ErodeNucWeiCen)
 
    % RNA2NucIndSelZ = find(MatrixOfRNASpotSelZ(:,6) == i);
    % RNANumForNucSelZ(i,:) = sum(MatrixOfRNASpotSelZ(RNA2NucIndSelZ,5));
    RNA2NucInd = find(MatrixOfRNASpot(:,6) == i);
    [min_d] =  NearestNeighborDistance(MatrixOfRNASpot(RNA2NucInd,1:3));
     MeanRNAAdjDis = mean(min_d).*TrueDistanceScale;
     NucMeanRNADistance(i) = MeanRNAAdjDis;
NorNucMeanRNAClusterIndex(i) = MeanRNAAdjDis./length(RNA2NucInd);

%on
 On_RNA2NucInd = find(OnMatrixOfRNASpot(:,6) == i);
    [min_d_on] =  NearestNeighborDistance(OnMatrixOfRNASpot(On_RNA2NucInd,1:3));
     Mean_On_RNAAdjDis = mean(min_d_on).*TrueDistanceScale;
     NucMeanRNA_On_Distance(i) = Mean_On_RNAAdjDis;
NorNucMean_On_RNAClusterIndex(i) = Mean_On_RNAAdjDis./length(On_RNA2NucInd);

%off
 Off_RNA2NucInd = find(OFFMatrixOfRNASpot(:,6) == i);
    [min_d_off] =  NearestNeighborDistance(OFFMatrixOfRNASpot(Off_RNA2NucInd,1:3));
     Mean_Off_RNAAdjDis = mean(min_d_off).*TrueDistanceScale;
     NucMeanRNA_Off_Distance(i) = Mean_Off_RNAAdjDis;
NorNucMean_Off_RNAClusterIndex(i) = Mean_Off_RNAAdjDis./length(Off_RNA2NucInd);

end


%%
parfor i = 1:length(ErodeNucWeiCen)
 
    % RNA2NucIndSelZ = find(MatrixOfRNASpotSelZ(:,6) == i);
    % RNANumForNucSelZ(i,:) = sum(MatrixOfRNASpotSelZ(RNA2NucIndSelZ,5));
    RNA2NucInd = find(MatrixOfRNASpot(:,6) == i);
    [min_d] =  pdist(MatrixOfRNASpot(RNA2NucInd,1:3));
     MeanRNAAdjDis = var(min_d);
     NucMeanRNADistance(i) = MeanRNAAdjDis;
NorNucMeanRNAClusterIndex(i) = MeanRNAAdjDis./length(RNA2NucInd);

%on
 On_RNA2NucInd = find(OnMatrixOfRNASpot(:,6) == i);
    [min_d_on] =  pdist(OnMatrixOfRNASpot(On_RNA2NucInd,1:3));
     Mean_On_RNAAdjDis = var(min_d_on);
     NucMeanRNA_On_Distance(i) = Mean_On_RNAAdjDis;
NorNucMean_On_RNAClusterIndex(i) = Mean_On_RNAAdjDis./length(On_RNA2NucInd);

%off
 Off_RNA2NucInd = find(OFFMatrixOfRNASpot(:,6) == i);
    [min_d_off] =  pdist(OFFMatrixOfRNASpot(Off_RNA2NucInd,1:3));
     Mean_Off_RNAAdjDis = var(min_d_off);
     NucMeanRNA_Off_Distance(i) = Mean_Off_RNAAdjDis;
NorNucMean_Off_RNAClusterIndex(i) = Mean_Off_RNAAdjDis./length(Off_RNA2NucInd);

end

%%
MidEm = imread(['20250423def 9623 hb SunTag S1.lif_Em9 nc14 hb Mid_ch00.tif']);
SurEm = imread(['20250423def 9623 hb SunTag S1.lif_Em9 nc14 hb Sur_ch00.tif']);

RotateMidEm = rot90(MidEm,-1);
RotateSurEm = rot90(SurEm,-1);
ScaleFactor = 20./(63*1.3);
%shrink the nuclei image for accurate feature extraction
MIPNuc = max(dapi_raw_3D(:,:,1:RefineZ),[],3);
% % % MIPErodeNuc_true = max(ErodeNuc_true,[],3);
% % % AA = MIPNuc.*MIPErodeNuc_true;
AA = MIPNuc;
ScaleMIPNuc = imresize(AA,ScaleFactor);
% figure();
% imshow(ScaleMIPNuc,[]);
%shoule transform the ScaleMIPNuc to uint16(ScaleMIPNuc)
% MIPNuc = max(dapi_raw_3D(:,:,1:RefineZ),[],3);
% % NucPlane = dapi_raw_3D(:,:,54);
% MIPErodeNuc_true = max(Nuc_true,[],3);
% % AA = MIPNuc.*MIPErodeNuc_true;
% ScaleMIPNuc = imresize(NucPlane,ScaleFactor);
% ScaleMIPNuc = (ScaleMIPNuc - min(min(ScaleMIPNuc)))./(max(max(ScaleMIPNuc)) - min(min(ScaleMIPNuc))).*65535;
% RotateSurEm = (RotateSurEm - min(min(RotateSurEm)))./(max(max(RotateSurEm)) - min(min(RotateSurEm))).*255;
% figure();
% imshow(ScaleMIPNuc,[]);
[MOVINGREG,ErrorMteric] = RegisterImages(uint16(ScaleMIPNuc),RotateSurEm);
% not valid: [MOVINGREG,ErrorMteric] = registerImages_cjy(uint16(max(dapi_raw_3D,[],3)),RotateSurEm);
%plot the register result
%method 1
% figure();
% imshow(RotateSurEm);
% hold on;
% contour(MOVINGREG.RegisteredImage,'r');
% hold off;
%method 2
figure();
[D RD] = imfuse(RotateSurEm,MOVINGREG.RegisteredImage,'falsecolor');
imshow(D);
% saveas(gcf, sprintf('ImageRegEM2_20230106HbNC13_2_GFP2focalplane.tif'));
% imwrite(D,'ImageRegEM4_20221004hbNC13_1_G2.tiff');
%apply the transform matrix to the nuclei coordinate
%https://www.mathworks.com/help/images/ref/affine2d.transformpointsforward.html
%https://www.mathworks.com/help/symbolic/rotation-matrix-and-transformation-matrix.html
tform = MOVINGREG.Transformation;
ErodeNucWeiCenScale = ErodeNucWeiCen(:,1:2);
ErodeNucWeiCenScale(:,1) = ErodeNucWeiCen(:,1).*ScaleFactor;
ErodeNucWeiCenScale(:,2) = ErodeNucWeiCen(:,2).*ScaleFactor;
U = ErodeNucWeiCenScale;
TransErodeNucWeiCen = transformPointsForward(tform,U);
%get AP axis in the midsagittal image in 20X lens
%obtain midsagittal image embryo mask
RotateMidEmFil = imgaussfilt(RotateMidEm,200);

[T EM] = graythresh(RotateMidEmFil);
EMMaskRaw = imbinarize(RotateMidEm,T);
EMMaskRaw= imfill(EMMaskRaw,'holes');
% figure();
% imshow(EMMaskRaw)
seEm = strel("disk",50);
EMMaskRawopen = imopen(EMMaskRaw,seEm);
% EMMaskRawopen_true = EMMaskRawopen;
% 
% 
%  EMMaskRawclose = imclose(EMMaskRaw,seEm);
% % % EMMaskRawclose = imfill(EMMaskRawclose,'holes');
% % % seEm2 = strel("disk",50);
% % % EMMaskRawcloseopen = imopen(EMMaskRawclose,seEm2);

% 

LLEMMaskRawopen = bwlabeln(EMMaskRawopen);
stsEMMaskRawopen = regionprops3(LLEMMaskRawopen,'Volume');
VolEMMaskRawopen = [stsEMMaskRawopen.Volume];
[M I] = max(VolEMMaskRawopen);
EMMaskRawopen_true = ismember(LLEMMaskRawopen,I);

  bw2 = ~bwareaopen(~EMMaskRawopen_true, 15);%10 15
   D = bwdist(~bw2);
   D = -D;
   %local minimum mark
   LocalMinMask = imextendedmin(D,2);%2
   %varified local minima location 
   % imshowpair(D,LocalMinMask,'blend');
   % Modify the distance transform so it only has minima at the desired locations,
   D2 = imimposemin(D,LocalMinMask);
   L = watershed(D2);
   L2 =  watershed(D);
   
   L2(~bw2) = 0;
   % imshow(label2rgb(L));
   %ovelap image
   % figure();
   % imshow(labeloverlay(AdjEmSingleOri,L,'Transparency',0.6));
   %the logical ture of the watershed segmentation located in the region of
   %the raw nuclei mask
   ZerosCross = EMMaskRawopen_true & (L==0);
   %assign the position find above to 0 in the mask image 
   EMMaskRawopen_true(ZerosCross == 1) = 0;
%    RefineNucMask = RawNucMask;

figure();
imshow(RotateMidEm);
hold on
[B,L] = bwboundaries(EMMaskRawopen_true);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
   hold on
end
hold off
%%
load('EM9_20250424_hb_nc14_GCN4_TR');
nx = TransErodeNucWeiCen(:,1);
ny = TransErodeNucWeiCen(:,2);
emmask = EMMaskRawopen_true;
flipAP =1;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn] = convertXY2AP(NorNucMeanRNAClusterIndex, nx, ny, emmask, fig,flipAP);
%[ventralData dorsalData APlength] = convertXY2AP(NucleiIntensi


nx = TransErodeNucWeiCen(:,1);
ny = TransErodeNucWeiCen(:,2);
emmask = EMMaskRawopen_true;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn] = convertXY2AP(RNANumForNucTrue, nx, ny, emmask, fig,flipAP);

figure();
scatter(NormAPn,NorNucMeanRNAClusterIndex,50,'filled');
xlabel('EL');
ylabel('Normalized mRNA-to-mRNA distance (a.u.)');
title('total  mRNA');

figure();
scatter(NormAPn,NorNucMean_On_RNAClusterIndex,50,'filled');
title('On  mRNA');
xlabel('EL');
ylabel('Normalized mRNA-to-mRNA distance (a.u.)');

figure();
scatter(NormAPn,NorNucMean_Off_RNAClusterIndex,50,'filled');
title('Off  mRNA');
xlabel('EL');
ylabel('Normalized mRNA-to-mRNA distance (a.u.)');

figure();
scatter(RNANumForNucTrue,NorNucMeanRNAClusterIndex,50,'filled');
xlabel('mRNA counts');
ylabel('Normalized mRNA-to-mRNA distance (a.u.)');
title('total  mRNA');



%%
load('EM9_20250424_hb_nc14_GCN4_TR');
nx = TransErodeNucWeiCen(:,1);
ny = TransErodeNucWeiCen(:,2);
emmask = EMMaskRawopen_true;
flipAP =1;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn] = convertXY2AP(NucMeanRNADistance, nx, ny, emmask, fig,flipAP);
%[ventralData dorsalData APlength] = convertXY2AP(NucleiIntensi


nx = TransErodeNucWeiCen(:,1);
ny = TransErodeNucWeiCen(:,2);
emmask = EMMaskRawopen_true;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn2] = convertXY2AP(RNANumForNucTrue, nx, ny, emmask, fig,flipAP);

figure();
scatter(NormAPn,NucMeanRNADistance,50,'filled');
xlabel('EL');
ylabel('mRNA-to-mRNA distance (nm)');
title('total  mRNA');

figure();
scatter(NormAPn,NucMeanRNA_On_Distance,50,'filled');
title('On  mRNA');
xlabel('EL');
ylabel('mRNA-to-mRNA distance (nm)');

figure();
scatter(NormAPn,NucMeanRNA_Off_Distance,50,'filled');
title('Off  mRNA');
xlabel('EL');
ylabel('mRNA-to-mRNA distance (nm)');

figure();
scatter(RNANumForNucTrue,NucMeanRNADistance,50,'filled');
xlabel('mRNA counts');
ylabel('Normalized mRNA-to-mRNA distance (nm)');
title('total  mRNA');
%%
figure();
scatter(NormAPn,1./NorNucMeanRNAClusterIndex,50,'filled');
xlabel('EL');
ylabel('Clustering index (a.u.)');
title('total  mRNA');

figure();
scatter(NormAPn,1./NorNucMean_On_RNAClusterIndex,50,'filled');
title('On  mRNA');
xlabel('EL');
ylabel('Clustering index (a.u.)');

figure();
scatter(NormAPn,1./NorNucMean_Off_RNAClusterIndex,50,'filled');
title('Off  mRNA');
xlabel('EL');
ylabel('Clustering index (a.u.)');

figure();
scatter(RNANumForNucTrue,1./NorNucMeanRNAClusterIndex,50,'filled');
xlabel('mRNA counts');
ylabel('Clustering index (a.u.)');
title('total  mRNA');