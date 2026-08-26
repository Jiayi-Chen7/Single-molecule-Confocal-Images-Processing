clear all
close all
clc

%%
%initial parameters
%folder
% in_folder = 'E:/CJY/20220314/EM1/';
% mid_sur_in_folder = 'E:/CJY/20220314/EM1/MID SUR/';
% out_folder =  'E:/CJY/20220314/EM1/output/';
%the image names of zstack
dapi_channel = '20230221hb_Em1.lif_Em1_20230221_z*_ch00.tif';
% gfp_channel = '20221114_hb_G2_S7.lif_Em5Large_z*_ch01.tif';
mRNA_channel = '20230221hb_Em1.lif_Em1_20230221_z*_ch02.tif';
% image list
imlist_dapi = dir(dapi_channel);
% imlist_gfp = dir(gfp_channel);
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
  
%register parametes
HighMag = 63;
LowMag = 20;
ZoomFactor = 1.3;
RegScaleFactor = LowMag./(HighMag.*ZoomFactor);
TrueNucVol = 20000;% 3D volume threshold to remove the yolk nuclei


%% read in raw data

dapi_raw_3D = zeros(size_in_merge_image(1),size_in_merge_image(2),numel(imlist_dapi));

mRNA_raw_3D =zeros(size_in_merge_image(1),size_in_merge_image(2),5);

% ZstackRange = linspace(1,size(mRNA_raw_3D,3),size(mRNA_raw_3D,3));
ZstackRange = [28 29 30 31 32];
for image_I = 1:5
    image_I2 = ZstackRange(image_I);
%     dapi_raw_3D(:,:,image_I) = imread([imlist_dapi(image_I2).name]);
%     gfp_raw_3D(:,:,image_I) = imread([imlist_gfp(image_I2).name]);
    mRNA_raw_3D(:,:,image_I) = imread([imlist_mRNA(image_I2).name]);
end

parfor image_I = 1:length(imlist_dapi)
    dapi_raw_3D(:,:,image_I) = imread([imlist_dapi(image_I).name]);
end
RefineZ = (numel(imlist_dapi));
DapiRawMask = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
NCEst = 14;
grayth = 1;
graythAdj = 0;

tic 
fig = 0;
parfor ImageMaskIndex = 1:RefineZ
EmSingleOri = dapi_raw_3D(:,:,ImageMaskIndex);
% [Mask T] = Nuclei_seg_cjy(EmSingleOri, ExVal, IxVal,low_th,cir_th,ImageMaskIndex,NCEst);
% [Mask] = Nuclei_seg_cjy(EmSingleOri,5, 25,[],0.6,ImageMaskIndex,NCEst,500,5,fig,5);
[Mask] = Nuclei_seg_cjy(EmSingleOri,2,[],[],0.7,ImageMaskIndex,NCEst,[],12,fig,12,grayth,graythAdj,6000,1);
DapiRawMask(:,:,ImageMaskIndex) = Mask;
end
toc
BWDapiRawMask = bwlabeln(DapiRawMask);
StsDapiRawMask = regionprops3(BWDapiRawMask,dapi_raw_3D(:,:,1:RefineZ),'Volume','MeanIntensity','Centroid');
NucVol = [StsDapiRawMask.Volume];
NucMeanInt = [StsDapiRawMask.MeanIntensity];
NucCenRaw = [StsDapiRawMask.Centroid];
% figure();
% histogram(NucVol);
%the threshold is estimate from the prior knowledge of the size of a nuclei
%the diameter of ~ 5um: 5 ./ 0.069 um/pixel ~ 72
%should be observed in at least 5 layers of the z stacks
TrueNucVolLB = 20000;%20000
TrueNucVolHB = 80000;%20000
NucIntenThLB = 9000;
NucIntenThHB = 16000;
TrueNucVol2 = 0;%20000
NucIntenTh2 = 0;
TrueNucInd = find(NucVol > TrueNucVol2& NucMeanInt>NucIntenTh2);
Nuc_true = ismember(BWDapiRawMask,TrueNucInd);
LLNuc_true = bwlabeln(Nuc_true);
stsNuc_true = regionprops3(LLNuc_true,dapi_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid');
NuccVol = [stsNuc_true.Volume];
NucWeiCen = [stsNuc_true.WeightedCentroid];
% length(NuccVol); 401
%shrink the muclei mask for more accurate nuclei assignment of the mRNA and
%GFP spot
ShrinkNucPara = 5;
se = strel('cube',ShrinkNucPara);
ErodeNuc_true = imerode(Nuc_true,se);
LLErodeNuc_true = bwlabeln(ErodeNuc_true);
stsLLErodeNuc_true = regionprops3(LLErodeNuc_true,dapi_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid');
ErodeNucVol = [stsLLErodeNuc_true.Volume];
ErodeNucWeiCen = [stsLLErodeNuc_true.WeightedCentroid];

%%
mRNA_raw_MIP = max(mRNA_raw_3D,[],3);
LowHighRNA = stretchlim(uint16(mRNA_raw_MIP));
AdjmRNA_raw_MIP = imadjust(uint16(mRNA_raw_MIP),LowHighRNA,[]);
figure();
imshow(AdjmRNA_raw_MIP);
% saveas(gcf, sprintf('MIPRNATest_I2.fig'));
ExRNACore = 1.2;
IxRNACore = 2.2;
Ex = fspecial('gaussian',21,ExRNACore);%10
Ix = fspecial('gaussian',21,IxRNACore);%20
outE = imfilter(single(AdjmRNA_raw_MIP),Ex,'replicate'); 
outI = imfilter(single(AdjmRNA_raw_MIP),Ix,'replicate'); 
outims = outE - outI;  

% 
figure();
imshow(outims,[]);


Gauoutims = imgaussfilt(single(outims),1.5);
% figure();
% imshow(Gauoutims(1536:2048,1536:2048),[]);
LowHighRNA = stretchlim(uint16(Gauoutims));

AdjOutims = imadjust(uint16(Gauoutims),LowHighRNA,[]);
% GauFil = imgaussfilt(single(SinglemRNA),1.5);


NumOfCanRNAMask = zeros(1,200);
ThVal = linspace(0,1,200);
for i = 1:200
    CanRawRNAMask = imbinarize(AdjOutims,ThVal(i));
    BWCanRawRNAMask = bwlabel(CanRawRNAMask);
    StsCanRawRNAMask = regionprops(BWCanRawRNAMask,'Area');
    CanRNAArea = [StsCanRawRNAMask.Area];
    NumOfCanRNAMask(i) = length(CanRNAArea);
end
% figure();
% axis square
% plot(NumOfCanRNAMask,'LineWidth',2) 
ElbowPointIndexRelative = GetElbowIndex(NumOfCanRNAMask);
% ThVal(22)
RawRNAMask = imbinarize(AdjOutims,ThVal(ElbowPointIndexRelative));
LLRawRNAMask = bwlabeln(RawRNAMask);
stsRawRNAMask = regionprops3(LLRawRNAMask,'volume');
Vol = [stsRawRNAMask.Volume];
TrueRNAInd = find(Vol > 10);
RNA_trueMask = ismember(LLRawRNAMask,TrueRNAInd);


OutimsAfterMask = double(AdjOutims).*RNA_trueMask;
RawRNAMask2 = imregionalmax(OutimsAfterMask);
LLRawRNAMask2 = bwlabeln(RawRNAMask2);
stsRawRNAMask2 = regionprops3(LLRawRNAMask2,'Centroid');

WCentroid = [stsRawRNAMask2.Centroid];
figure();
imshow(AdjmRNA_raw_MIP);
hold on
scatter(WCentroid(:,1),WCentroid(:,2),14,'r','filled');
hold off   

XAxisEdge = linspace(1,2048,10);
YAxisEdge = linspace(1,2048,10);
SampleInd = (length(XAxisEdge)-1).*(length(YAxisEdge)-1);
% % SampleRNAClusterIndex = zeros(SampleInd,1);
% % SampleCenPosition = zeros(SampleInd,2);
% % Samplesize = zeros(SampleInd,1);
% % RefMat = reshape([1:1:SampleInd],[(length(XAxisEdge)-1),(length(YAxisEdge)-1)])
% % A = reshape(RefMat,[1,SampleInd])
SampleRNAClusterIndex = zeros((length(XAxisEdge)-1),(length(YAxisEdge)-1));
SampleCenPositionX = zeros((length(XAxisEdge)-1),(length(YAxisEdge)-1));
SampleCenPositionY = zeros((length(XAxisEdge)-1),(length(YAxisEdge)-1));
Samplesize =  zeros((length(XAxisEdge)-1),(length(YAxisEdge)-1));
SimIndex = zeros((length(XAxisEdge)-1),(length(YAxisEdge)-1));

for i = 1:(length(XAxisEdge)-1)
   for j = 1:(length(YAxisEdge)-1)
%     i = 1; j = 2;
    SampleRNAInd = find(WCentroid(:,1)>=XAxisEdge(i)&WCentroid(:,1)<XAxisEdge(i+1)&WCentroid(:,2)>=YAxisEdge(j)&WCentroid(:,2)<YAxisEdge(j+1));
    ind = i*j
    Samplesize(i,j) = length(SampleRNAInd);
    SampleCenPositionX(i,j) = XAxisEdge(i)+(XAxisEdge(i+1)-XAxisEdge(i))./2;
     SampleCenPositionY(i,j) = YAxisEdge(j)+(YAxisEdge(j+1)-YAxisEdge(j))./2;
     SampleArea = (XAxisEdge(i+1)-XAxisEdge(i)).*(YAxisEdge(j+1)-YAxisEdge(j));
     PoiClusterDis = 0.5.*((SampleArea./length(SampleRNAInd)).^(0.5));
%      [d,pass] = CalRNAMutualDis([WCentroid(SampleRNAInd,:),WCentroid(SampleRNAInd,:),2);
SimdMat1 = zeros(100,1);
SimdMat2 = zeros(100,1);
for k = 1:100
     num_points = length(SampleRNAInd);
     x_min = 1;
     x_max = XAxisEdge(i+1)-XAxisEdge(i);
     y_min = 1;
     y_max = YAxisEdge(j+1)-YAxisEdge(j);
     x = x_min + (x_max - x_min) * rand(num_points, 1);
     y = y_min + (y_max - y_min) * rand(num_points, 1);
     [dsim] = NearestNeighborDistance([x,y]);
     SimdMat1(k) = mean(dsim);
end

 for k = 1:100
     num_points = length(SampleRNAInd);
     x_min = 1;
     x_max = XAxisEdge(i+1)-XAxisEdge(i);
     y_min = 1;
     y_max = YAxisEdge(j+1)-YAxisEdge(j);
     x = x_min + (x_max - x_min) * rand(num_points, 1);
     y = y_min + (y_max - y_min) * rand(num_points, 1);
     [dsim] = NearestNeighborDistance([x,y]);
     SimdMat2(k) = mean(dsim);
     end
SimIndex(i,j) = mean(SimdMat2)./mean(SimdMat1);
    [min_d] =  NearestNeighborDistance(WCentroid(SampleRNAInd,:));
     MeanRNAAdjDis = mean(min_d);
%    SampleRNAClusterIndex(i,j) = MeanRNAAdjDis./mean(SimdMat1);
SampleRNAClusterIndex(i,j) = mean(SimdMat1)./mean(SimdMat1);
   end
end
Samplesize = reshape(Samplesize,[SampleInd,1]);
SampleCenPositionX = reshape(SampleCenPositionX,[SampleInd,1]);
SampleCenPositionY = reshape(SampleCenPositionY,[SampleInd,1]);
SampleRNAClusterIndex = reshape(SampleRNAClusterIndex,[SampleInd,1]);
SimIndex = reshape(SimIndex,[SampleInd,1]);
%%
MidEm = imread(['20230221hb_Em1.lif_Em1-Mid_ch00.tif']);
SurEm = imread(['20230221hb_Em1.lif_Em1-Sur_ch00.tif']);
RotateMidEm = rot90(MidEm,-1);
RotateSurEm = rot90(SurEm,-1);
ScaleFactor = 20./(63*1.3);
%shrink the nuclei image for accurate feature extraction
MIPNuc = max(dapi_raw_3D(:,:,1:RefineZ),[],3);
% % % MIPErodeNuc_true = max(ErodeNuc_true,[],3);
% % % AA = MIPNuc.*MIPErodeNuc_true;
AA = MIPNuc;
ScaleMIPNuc = imresize(AA,ScaleFactor);
[MOVINGREG,ErrorMteric] = registerImages_cjy(uint16(ScaleMIPNuc),RotateSurEm);
tform = MOVINGREG.Transformation;
ErodeNucWeiCenScale = zeros(length(SampleCenPositionX),2);
ErodeNucWeiCenScale(:,1) = SampleCenPositionX.*ScaleFactor;
ErodeNucWeiCenScale(:,2) = SampleCenPositionY.*ScaleFactor;
U = ErodeNucWeiCenScale;
TransErodeNucWeiCen = transformPointsForward(tform,U);
RotateMidEmFil = imgaussfilt(RotateMidEm,5);
[T EM] = graythresh(RotateMidEmFil);
EMMaskRaw = imbinarize(RotateMidEm,T);
EMMaskRaw= imfill(EMMaskRaw,'holes');
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
seEm = strel("disk",60);
EMMaskRawopen = imopen(EMMaskRaw,seEm);
LLEMMaskRawopen = bwlabeln(EMMaskRawopen);
stsEMMaskRawopen = regionprops3(LLEMMaskRawopen,'Volume');
VolEMMaskRawopen = [stsEMMaskRawopen.Volume];
[M I] = max(VolEMMaskRawopen);
EMMaskRawopen_true = ismember(LLEMMaskRawopen,I);
fig = 1;
% intens = 0.5.*RNANumForNuc;
nx = TransErodeNucWeiCen(:,1);
ny = TransErodeNucWeiCen(:,2);
emmask = EMMaskRawopen_true;

flipAP =0;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn] = convertXY2AP(SampleRNAClusterIndex, nx, ny, emmask, fig,flipAP);
%[ventralData dorsalData APlength] = convertXY2AP(NucleiIntensi
NonZeroInd = find(Samplesize>0);
NonZeroSamplesize = Samplesize(NonZeroInd);

plotInd = find(Samplesize>(mean(NonZeroSamplesize) - std(NonZeroSamplesize)));
figure();
yyaxis left
scatter(NormAPn(plotInd),SampleRNAClusterIndex(plotInd),50,'filled');
ylabel("Cluster Index");
yyaxis right
scatter(NormAPn(plotInd),Samplesize(plotInd),50,'filled');
ylabel("Sample size");
ant = find(Samplesize>(mean(NonZeroSamplesize) -std(NonZeroSamplesize))&NormAPn<0.38);
pos = find(Samplesize>(mean(NonZeroSamplesize) -std(NonZeroSamplesize))&NormAPn>=0.38&NormAPn<=0.45);
xlim([0.15,0.65]);
xticks([0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]);
xlabel("EL");
box on;
set(gcf,'position',[250 250 600 600]);
axis square
set (gca,'linewidth',2,'fontsize',24);
[h p] = ttest2(SampleRNAClusterIndex(ant),SampleRNAClusterIndex(pos))

figure();
scatter(Samplesize,SampleRNAClusterIndex);
figure();
scatter(Samplesize(plotInd),SampleRNAClusterIndex(plotInd));
hold on
scatter(Samplesize,SimIndex)
xlabel('Sample Size');
ylabel('Cluster Index');
legend('Exp/Sim','Sim1/Sim2');
title('The number of XY steps: 80');
[p r] = corr(Samplesize(plotInd),SampleRNAClusterIndex(plotInd));

plotInd = find(Samplesize>0);
figure();
yyaxis left
scatter(NormAPn(plotInd),SampleRNAClusterIndex(plotInd),50,'filled');
ylabel("Cluster Index");
yyaxis right
scatter(NormAPn(plotInd),Samplesize(plotInd),50,'filled');
ylabel("Sample size");
xlim([0.15,0.65]);
xticks([0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]);
xlabel("EL");
box on;
set(gcf,'position',[250 250 600 600]);
axis square
set (gca,'linewidth',2,'fontsize',24);
ant = find(NormAPn<0.38);
pos = find(NormAPn>=0.38&NormAPn<0.5);
ter = find(NormAPn>=0.5);

[h p] = ttest2(SampleRNAClusterIndex(ter),SampleRNAClusterIndex(pos))
save Em1_0221_SpatialCluster.mat SampleRNAClusterIndex NormAPn Samplesize

rng('default')  % For reproducibility
x1 = SampleRNAClusterIndex(ant);
x2 = SampleRNAClusterIndex(pos);
x3 = SampleRNAClusterIndex(ter);
x = [x1; x2; x3];

g1 = repmat({'Ant'},length(x1),1);
g2 = repmat({'Pos'},length(x2),1);
g3 = repmat({'Ter'},length(x3),1);
g = [g1; g2; g3];
boxplot(x,g)