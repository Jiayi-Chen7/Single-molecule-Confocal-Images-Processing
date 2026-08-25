clear all
close all
clc

%%
%input filename
dapi_channel = 'BcdE1S2Em2.lif_BcdE1S2Em2_z*_ch00.tif';
gfp_channel = 'BcdE1S2Em2.lif_BcdE1S2Em2_z*_ch02.tif';
mRNA_channel = 'BcdE1S2Em2.lif_BcdE1S2Em2_z*_ch01.tif';
% image list
imlist_dapi = dir(dapi_channel);
imlist_gfp = dir(gfp_channel);
imlist_mRNA = dir(mRNA_channel);
% nuclei segmentation parameters
cir_th = 0.75;
%the size of the image in each z step
size_in_merge_image = [2048,2048];
size_in_plane_image = [1024,1024];
% count scale considering the size of nucleus volume unit between nc 13 and
% nc 14
CountScale = 1;
%register parametes
HighMag = 63;
LowMag = 20;
ZoomFactor = 1.3;
RegScaleFactor = LowMag./(HighMag.*ZoomFactor);



%% import raw data

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
%a optional z step range 
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
for ImageMaskIndex = 1:RefineZ
EmSingleOri = dapi_raw_3D(:,:,ImageMaskIndex);
[Mask] = NucleiSeg(EmSingleOri,2,[],200,0.3,ImageMaskIndex,NCEst,[],12,fig,12,grayth,graythAdj,3500);
DapiRawMask(:,:,ImageMaskIndex) = Mask;
end
toc

% % % 
%for visulization
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


BWDapiRawMask = bwlabeln(DapiRawMask);
StsDapiRawMask = regionprops3(BWDapiRawMask,dapi_raw_3D(:,:,1:RefineZ),'Volume','MeanIntensity','Centroid');
NucVol = [StsDapiRawMask.Volume];
NucMeanInt = [StsDapiRawMask.MeanIntensity];
NucCenRaw = [StsDapiRawMask.Centroid];
%for visulization
% figure();
% histogram(NucMeanInt);

TrueNucVol2 = 0;
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
% % % % % % 
% for j = 1:RefineZ
% figure();
% [B,L] = bwboundaries(Nuc_true(:,:,j));
% imshow(dapi_raw_3D(:,:,j),[]);
% hold on
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
%    hold on
% end
%      hold off    
% end
% %      saveas(gcf, sprintf('RefineNucleiWithMask%d.tif',j));
% % end

%%
%%RNA mask processing 
%test to obtain the absolute single mRNA intensity threshold and the size
%of the transcription foci
%segmentation of transcription foci in nuclei
ExRNACore = 1.2;
IxRNACore = 2.2;

ExFociCore = 4.5;
InFociCore = 15;%15
seval = 4;
MIPRNA = max(mRNA_raw_3D,[],3);
MaxNuc_true = max(Nuc_true,[],3);
[MaxTranscriptionFoci,MaxTranscriptionFociMask] = MaxFoci(MIPRNA,ExFociCore,InFociCore,seval);
% se = strel('disk',3);
% MaxTranscriptionFociMask2 = imopen(MaxTranscriptionFociMask,se);

% figure();
% imshow(MIPRNA,[]);
% hold on
% [B,L] = bwboundaries(MaxTranscriptionFociMask);
% % [BN,LN] = bwboundaries(MaxNuc_true(1:500,1:500));
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5);
%    hold on
% end
% % for k = 1:length(BN)
% %    boundary = BN{k};
% %    plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 1.5);
% %    hold on
% % end
%      hold off   
% exportgraphics(gcf,"FociMIPRaw.tiff")

for img_ind = 1:RefineZ
SinglemRNA = mRNA_raw_3D(:,:,img_ind);
SinglePlaneNuc = Nuc_true(:,:,img_ind);
mRNAFociRawMask(:,:,img_ind) = SinglePlaneFoci_Ori(SinglemRNA,SinglePlaneNuc,ExFociCore,InFociCore,img_ind,MaxTranscriptionFociMask);
end
Nr = 1;
mRNAFociMask = MultiLayerSpotIdentify(mRNAFociRawMask,Nr);


% % for visulization
% for img_ind = 1:RefineZ
% figure();
% imshow(mRNA_raw_3D(:,:,img_ind),[]);
% hold on
% [B,L] = bwboundaries(mRNAFociMask(:,:,img_ind));
% [BF,LF] = bwboundaries(MaxTranscriptionFociMask);
% [BN,LN] = bwboundaries(Nuc_true(:,:,img_ind));
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5);
%    hold on
% end
% for k = 1:length(BF)
%    boundary = BF{k};
%    plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 1.5);
%    hold on
% end
% for k = 1:length(BN)
%    boundary = BN{k};
%    plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth',2);
%    hold on
% end
% hold off      
% % % % % % saveas(gcf, sprintf('FociMask%d.tif',img_ind));
% end

%preliminary segemntation of single-molecule mRNA in cytoplasm
ExRNACore = 1.2;
IxRNACore = 2.2;
tic  
for RNAStackIndex = 1:RefineZ
SinglemRNA = mRNA_raw_3D(:,:,RNAStackIndex);
SinglePlaneNuc = Nuc_true(:,:,RNAStackIndex);
mRNAFociRawMaskPlane = mRNAFociMask(:,:,RNAStackIndex);
% figure();
% imshow(SinglemRNA(1536:2048,1536:2048),[]);
[CytoSingleRNAMask,GauFil] = SinglePlaneRNAMaskLocalMax2(SinglemRNA,SinglePlaneNuc,ExRNACore,IxRNACore,RNAStackIndex,mRNAFociRawMaskPlane);
mRNARawMask(:,:,RNAStackIndex) = CytoSingleRNAMask;
mRNAGauFil(:,:,RNAStackIndex) = GauFil;
end
toc

%%
%%the true cytoplasmic mRNA should be observed in three consecutive  z step
%%layers
tic
Nr = 1;
[MultiLayerMask,mask_out] = MultiLayerSpotIdentify(mRNARawMask,Nr);
toc
%%
%further refinment of cytoplasmic RNA masks
CytoRNAImg = double(mRNA_raw_3D(:,:,1:RefineZ)).*MultiLayerMask;
IntLower = 2500;
CytoRNAImg_Hmax = imhmax(CytoRNAImg,IntLower);
% CytoSingleRNAMask1 = imregionalmax(CytoRNAImg,26);

CytoSingleRNAMask2 = imregionalmax(CytoRNAImg_Hmax,26);
SE = strel('cube',2);
CytoSingleRNAMask3 = imdilate(CytoSingleRNAMask2,SE);
%%

LLMultiLayerMask = bwlabeln(CytoSingleRNAMask3);
stsMultiLayerMask = regionprops3(LLMultiLayerMask,mRNA_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid','MeanIntensity','MaxIntensity');
VolRNA = [stsMultiLayerMask.Volume];
WCentroid = [stsMultiLayerMask.WeightedCentroid];
MeanIntensity =[stsMultiLayerMask.MeanIntensity];
GauMaxIntensity = [stsMultiLayerMask.MaxIntensity];

%%%for visulization and manual inspection
zrange = 57:60;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(WCentroid(:,3)>=55 & WCentroid(:,3)<=62);
WCentroidZRange = WCentroid(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
saveas(gcf, sprintf('RNARawMaskBcdE1S2Em2_57_60.fig'));

zrange = 20:24;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(WCentroid(:,3)>=18 & WCentroid(:,3)<=26);
WCentroidZRange = WCentroid(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
zoom(6)
saveas(gcf, sprintf('RNARawMaskBcdE1S2Em2_20_24.fig'));

zrange = 5:10;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(WCentroid(:,3)>=3 & WCentroid(:,3)<=12);
WCentroidZRange = WCentroid(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
zoom(6)
saveas(gcf, sprintf('RNARawMaskBcdE1S2Em2_5_10.fig'));

zrange = 45:49;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(WCentroid(:,3)>=43 & WCentroid(:,3)<=51);
WCentroidZRange = WCentroid(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
zoom(6)
saveas(gcf, sprintf('RNARawMaskBcdE1S2Em2_45_49.fig'));

%%
%%gaussian fitting of the candidate RNA spots in cytoplasmic
lb = [0,0,0,0,0,0];
ub = [inf,inf,inf,inf,inf,inf];
XYrange = 3;
Zrange = 2;
sigmax0 = 1.5;   %%% Initial value of sigma_x
sigmay0 = 1.5;   %%% Initial value of sigma_y
gau2D = @(x,xdata) x(1)*exp(-((xdata(:,1) - x(2)).^2./(2.*x(3).^2)+(xdata(:,2) - x(4)).^2./(2.*x(5).^2)))+x(6);
para = zeros(length(VolRNA),6);
TotSpotVal = zeros(length(VolRNA),1);
resnorm = zeros(length(VolRNA),1);
exitflag = zeros(length(VolRNA),1);
dim0 = size(mRNA_raw_3D);
Zposition = zeros(length(VolRNA),1);
options = optimset('Display','off');
Xboundary = size_in_merge_image(1);
Yboundary = size_in_merge_image(2);
tic 
for i = 1:length(VolRNA)
LinIndex = find(LLMultiLayerMask == i);
[X0,Y0,Z0] = ind2sub(dim0,LinIndex);
XPeak = round(mean(X0));
YPeak = round(mean(Y0));
Zpeak = round(mean(Z0));
% Xmin = min(X0);
% Xmax = max(X0);
% Ymin = min(Y0);
% Ymax = max(Y0);
% Zmin = min(Z0);
% Zmax = max(Z0);
Xmin = max((XPeak-XYrange),1);
Xmax = min((XPeak+XYrange),Xboundary);
Ymin = max((YPeak-XYrange),1);
Ymax = min((YPeak+XYrange),Yboundary);
Zmin = max((Zpeak-Zrange),1);
Zmax = min((Zpeak+Zrange),RefineZ);
[X,Y] = ndgrid(Xmin:Xmax,Ymin:Ymax);
xdata = [X(:),Y(:)];
szz = length(Zmin:Zmax);
TotIntensityInZPlane = zeros(1,szz);
SpotVal = mRNA_raw_3D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
% TotSpotVal(i,:) = sum(sum(sum(SpotVal)));
for indz = 1:szz
   IntVal =  max(max(SpotVal(:,:,indz)));
   TotIntensityInZPlane(indz) = IntVal;
end
[V I] = max(TotIntensityInZPlane);
% % %
GauFilMaxSpotIntensity = SpotVal(:,:,I);
% %  GauFilMaxSpotIntensity = max(SpotVal,[],3);
% figure();
% imshow(GauFilMaxSpotIntensity,[]);
% exportgraphics(gcf,"SingleRNAPoint_30000.tiff")

 ydata = GauFilMaxSpotIntensity(:);
 Amp = (max(ydata) - min(ydata));
para0 = [Amp, XPeak, sigmax0, YPeak, sigmay0,min(ydata)];
[paraoutput,resnormoutput,~,exitflagoutput,~] = lsqcurvefit(gau2D,para0,xdata,ydata,lb,ub,options);
para(i,:) = paraoutput;
resnorm(i,:) = resnormoutput;
exitflag(i,:) = exitflagoutput;
Zposition(i,:) = Zmin+I;
end
toc
save 2DGaussianFittingBcdE1S2Em2_GCN4.mat para resnorm exitflag
%load('2DGaussianFittingRaw0915Em3hb_NC14_1_GFP_2_Alfinal.mat'); 


% % % %plot for single spot raw data and fitting result
% % % xa = xdata(:,1);
% % % ya = xdata(:,2);
% % % value = ydata;
% % % figure();
% % % x = linspace(Xmin,Xmax,40);
% % % y = linspace(Ymin,Ymax,40);
% % % [x,y] = meshgrid(Xmin:0.1:Xmax,Ymin:0.1:Ymax);
% % % a = paraoutput(1);
% % % b = paraoutput(2);
% % % c = paraoutput(3);
% % % d = paraoutput(4);
% % % e = paraoutput(5);
% % % f = paraoutput(6);
% % % val_cal = a*exp(-((x - b).^2./(2.*c.^2)+(y - d).^2./(2.*e.^2)))+f;
% % % mesh(x,y,val_cal);
% % % hold on 
% % % scatter3(xdata(:,1),xdata(:,2),ydata,'filled','k');
% % % colorbar
% % % xlabel('X axis');
% % % set(gca,'fontname','Arial');
% % % ylabel('Y axis');
% % % set(gca,'fontname','Arial');
% % % zlabel('Intensity (a.u.)');
% % % set(gca,'fontname','Arial');
% % % set (gca,'linewidth',2,'fontsize',12);
% % % % set(gcf,'position',[100 100 300 300]);
% % % hold off
% % % %end of the plot

%6 columes correspond to the xyz axis in X63 lens, the integral intensity,
% the RNA count derived from the typical intensity of single mRNA and the
% corresponding nuclei indices

%criteria reference: Heng Xu, et al. Nature Method, 2015: code: spfilter.m
PSFStdRatio = para(:,3)./para(:,5);

% GoodRNASpotInd = find(0.4 <para(:,3)& para(:,3)<3 & 0.4<para(:,5) & para(:,5)<3 & (exitflag == 1|exitflag == 3)&(sqrt(1-(para(:,3)./para(:,5)).^2) <= 0.9)&para(:,1)<100000 & para(:,1)>4000);
GoodRNASpotInd = find(0.4 <para(:,3)& para(:,3)<3 & 0.4<para(:,5) & para(:,5)<3 & (exitflag == 1|exitflag == 3)&(PSFStdRatio>0.1 & PSFStdRatio<10)&para(:,1)<100000 & para(:,1)>IntLower);

GoodPara = para(GoodRNASpotInd,:);
GoodRNAMask = ismember(LLMultiLayerMask,GoodRNASpotInd);
LLGoodRNAMask = bwlabeln(GoodRNAMask);
stsGoodRNA = regionprops3(LLGoodRNAMask,mRNA_raw_3D(:,:,1:RefineZ),'WeightedCentroid');
GoodRNAWeiCen = [stsGoodRNA.WeightedCentroid];
%there may be nan value in the calculated weighted centroid
NANInd = isnan(GoodRNAWeiCen(:,1));
PassInd = find(NANInd < 1);
GoodRNAWeiCenRemoveNAN = GoodRNAWeiCen(PassInd,:);

MatrixOfRNASpot = zeros(length(PassInd),6);
MatrixOfRNASpot(:,1:3) = GoodRNAWeiCenRemoveNAN;
paraGood = GoodPara(PassInd,:);
TotRNASpotIntensity = 2.*pi.*paraGood(:,1).*paraGood(:,3).*paraGood(:,5);
MatrixOfRNASpot(:,4) = TotRNASpotIntensity;
% figure();
% histogram(TotRNASpotIntensity);
% [Val2, Edges2] = histcounts(TotRNASpotIntensity,'Normalization','probability');
[Val2, Edges2] = histcounts(TotRNASpotIntensity);
edges2 = Edges2(2:end) - (Edges2(2)-Edges2(1))/2;
%fit 3 terms gaussian function to obtain the typical intensity of the
%single cytoplamic mRNA
lb2 = [0,0,0,0,0];
ub2 = [inf,inf,inf,inf,inf];
Amp2 = max(Val2).*0.9;
Peak2 = Amp2.*0.6;
Peak3 = Peak2.*0.5;
% Peak4 = Amp2.*0.001;
para02 = [Amp2, 75000, 40000, Peak2, Peak3];
ThreeTermGau = @(x,xdata) x(1)*exp(-((xdata - x(2)).^2./(2.*x(3).^2)))+x(4)*exp(-((xdata - 2.*x(2)).^2./(4.*x(3).^2)))+x(5)*exp(-((xdata - 3.*x(2)).^2./(6.*x(3).^2)));
[paraoutputThreeTerm,resnormoutputThreeTerm,~,exitflagoutputThreeTerm,~] = lsqcurvefit(ThreeTermGau,para02,edges2,Val2,lb2,ub2,options);

%plot
[x, y] = prepareCurveData(edges2,Val2);

GauTerm1 =  paraoutputThreeTerm(1)*exp(-((x - paraoutputThreeTerm(2)).^2./(2.*paraoutputThreeTerm(3).^2)));
GauTerm2 =  paraoutputThreeTerm(4)*exp(-((x - 2.*paraoutputThreeTerm(2)).^2./(4.*paraoutputThreeTerm(3).^2)));
GauTerm3 = paraoutputThreeTerm(5)*exp(-((x - 3.*paraoutputThreeTerm(2)).^2./(6.*paraoutputThreeTerm(3).^2)));
% GauTerm4 = paraoutputThreeTerm(6)*exp(-((x - 4.*paraoutputThreeTerm(2)).^2./(8.*paraoutputThreeTerm(3).^2)));
Gautot = paraoutputThreeTerm(1)*exp(-((x - paraoutputThreeTerm(2)).^2./(2.*paraoutputThreeTerm(3).^2)))+paraoutputThreeTerm(4)*exp(-((x - 2.*paraoutputThreeTerm(2)).^2./(4.*paraoutputThreeTerm(3).^2)))+paraoutputThreeTerm(5)*exp(-((x - 3.*paraoutputThreeTerm(2)).^2./(6.*paraoutputThreeTerm(3).^2)));
figure();
plot(x,GauTerm1,'LineWidth',3);
hold on
plot(x,GauTerm2,'LineWidth',3);
hold on
plot(x,GauTerm3,'LineWidth',3);
hold on
plot(x,Gautot, 'LineWidth',3);
hold on
scatter(x,y,30,"black","filled");
hold off
xlabel("Intensity (a.u.)");
set (gca,'linewidth',2,'fontsize',12,'FontName','Arial');
ylabel("The number of spots");
set (gca,'linewidth',2,'fontsize',12,'FontName','Arial');
box on;
axis square
set(gcf,'position',[100 100 500 500]);
xlim([0 10.^6]);
lgd = legend('Term 1','Term 2','Term 3','Sum','Raw data points');
lgd.FontSize = 12;
lgd.FontName = 'Arial';

saveas(gcf, sprintf('ThreeTermGaufitBcdE1S2Em2_GCN4.fig'));

TypicalIntensity = paraoutputThreeTerm(2);
MatrixOfRNASpot(:,5) = ceil(TotRNASpotIntensity./ TypicalIntensity);
%one term gaussian fitting 
%try
% MatrixOfRNASpot(:,5) = ceil(TotRNASpotIntensity./ 5.2e+04);
%  save GoodRNAMask.mat GoodRNAMask
% save RNAIntensity_Em1_1114_hb_nc13_2.mat MultiLayerMask mRNA_raw_3D
%%
%assign the mRNA spot in cytoplasm to the nearest nuclei
%http://www.qhull.org/

zrange = 57:60;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(MatrixOfRNASpot(:,3)>=55 & MatrixOfRNASpot(:,3)<=62);
WCentroidZRange = MatrixOfRNASpot(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
saveas(gcf, sprintf('RNAMasBcdE1S2Em2_57_60.fig'));

zrange = 20:24;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(MatrixOfRNASpot(:,3)>=18 & MatrixOfRNASpot(:,3)<=26);
WCentroidZRange = MatrixOfRNASpot(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
zoom(6)
saveas(gcf, sprintf('RNAMaskBcdE1S2Em2_20_24.fig'));

zrange = 5:10;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(MatrixOfRNASpot(:,3)>=3 & MatrixOfRNASpot(:,3)<=12);
WCentroidZRange = MatrixOfRNASpot(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
zoom(6)
saveas(gcf, sprintf('RNAMaskBcdE1S2Em2_5_10.fig'));

zrange = 45:49;
input_MIP = max(mRNA_raw_3D(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjRNARaw1 = imadjust(uint16(input_MIP),[0 0.2],[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(MatrixOfRNASpot(:,3)>=43 & MatrixOfRNASpot(:,3)<=51);
WCentroidZRange = MatrixOfRNASpot(IndZrange,:);
figure();
imshow(AdjRNARaw1);
hold on
scatter(WCentroidZRange(:,1),WCentroidZRange(:,2),20,'r','filled');
hold off
zoom(6)
saveas(gcf, sprintf('RNAMaskBcdE1S1Em3_45_49.fig'));
%%


xyPixelSize_nm = 69;
zStep_nm = 250;   
coordinateScale_nm = [xyPixelSize_nm, xyPixelSize_nm, zStep_nm];
NucCentroid_nm = ErodeNucWeiCen .* coordinateScale_nm;
RNACentroid_nm = GoodRNAWeiCenRemoveNAN .* coordinateScale_nm;
TNuc = delaunayn(NucCentroid_nm);
%k the indices of the closest points in P to the query points in PQ measured in Euclidean distance
[kRNA, distRNA2Nuc_nm] = dsearchn( ...
    NucCentroid_nm,TNuc,RNACentroid_nm);
MatrixOfRNASpot(:,6) = kRNA;
RNASpotSelectInd = find(MatrixOfRNASpot(:,3)>=Zstart & MatrixOfRNASpot(:,3)<= Zend);
MatrixOfRNASpotSelZ = MatrixOfRNASpot(RNASpotSelectInd,:);
% save EM2_20230106HbNC13_2_GFP2MatrixOfRNASpot.mat MatrixOfRNASpot
RNANumForNucSelZ = zeros(length(ErodeNucWeiCen),1);
RNANumForNuc = zeros(length(ErodeNucWeiCen),1);
parfor i = 1:length(ErodeNucWeiCen)
    RNA2NucIndSelZ = find(MatrixOfRNASpotSelZ(:,6) == i);
    RNANumForNucSelZ(i,:) = sum(MatrixOfRNASpotSelZ(RNA2NucIndSelZ,5));
    RNA2NucInd = find(MatrixOfRNASpot(:,6) == i);
    RNANumForNuc(i,:) = sum(MatrixOfRNASpot(RNA2NucInd,5));
end

%%
%register the data
%Nuclei channel registration to obtain the trasnfrom matrix
%register the 3D image to the
%%surface image
MidEm = imread(['BcdE1S2.lif_Em2 nc14 Mid_ch00.tif']);
SurEm = imread(['BcdE1S2.lif_Em2 nc14 Sur_ch00.tif']);

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
[MOVINGREG] = registerImagesNoRotation(ScaleMIPNuc,RotateSurEm);
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
imshow(RotateMidEm,[]);
hold on
[B,L] = bwboundaries(EMMaskRawopen_true);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
   hold on
end
hold off

%%
fig = 1;

intens = CountScale.*RNANumForNuc;
nx = TransErodeNucWeiCen(:,1);
ny = TransErodeNucWeiCen(:,2);
emmask = EMMaskRawopen_true;
% save ConvertAPDataEM2_20230106HbNC13_2_GFP2focalplane.mat intens nx ny emmask
% 
% intens = ConvertAPData.intens;
% nx = ConvertAPData.nx;
% ny = ConvertAPData.ny;
% emmask = ConvertAPData.emmask;
flipAP =1;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn] = convertXY2AP(intens, nx, ny, emmask, fig,flipAP);
%[ventralData dorsalData APlength] = convertXY2AP(NucleiIntensity, NucleiPosition(:,1), NucleiPosition(:,2), embryomask, fig);

% mRNACount = dorsalData(:,3);

figure();
scatter(NormAPn,intens,50,'filled');
xlim([0.15,0.65]);
ylim([0 300]);
xticks([0.2,0.3,0.4,0.5,0.6]);
yticks([0,100,200,300,400,500]);
xlabel("EL");
ylabel("mRNA counts");
box on;
set(gcf,'position',[250 250 600 600]);
axis square
set (gca,'linewidth',2,'fontsize',24);
saveas(gcf, sprintf('BcdE1S2Em2_RNACount.fig'));

%%
%GFP foci mask 
%normalize the GFP intensity 
gfp_raw_3DNor = gfp_raw_3D.*(~Nuc_true);
%Normalize the GFP in each single plane
% for  GFPStackIndex = 1:RefineZ
%     Plane = gfp_raw_3D(:,:,GFPStackIndex);
%     gfp_raw_3DNor(:,:,GFPStackIndex) = ((Plane - min(min(Plane )))./(max(max(Plane)) - min(min(Plane )))).*65535;
% end
% %global GFP normalization in the whole z stack
% gfp_raw_3DNor = ((gfp_raw_3D - min(min(min(gfp_raw_3D ))))./(max(max(max(gfp_raw_3D))) - min(min(min(gfp_raw_3D ))))).*65535;

%%
%obtain the GFP mask
tic
ExGFPCore = 1.2;
IxGFPCore = 2.2;
%
for GFPStackIndex = 1:RefineZ
SingleGFP = gfp_raw_3DNor(:,:,GFPStackIndex);
SinglePlaneNuc = Nuc_true(:,:,GFPStackIndex);
mRNAFociRawMaskPlane = mRNAFociMask(:,:,GFPStackIndex);
[CytoSingleGFPMask,GauGFPFil] = SinglePlaneRNAMaskLocalMax2(SingleGFP,SinglePlaneNuc,ExRNACore,IxRNACore,GFPStackIndex,mRNAFociRawMaskPlane);

% [CytoSingleGFPMask,GauGFPFil] = SinglePlaneGFPMask2(SingleGFP,SinglePlaneNuc,ExGFPCore,IxGFPCore,GFPStackIndex,[]);
GFPRawMask(:,:,GFPStackIndex) = CytoSingleGFPMask;
GFPGauFil(:,:,GFPStackIndex) = GauGFPFil;
end
toc
% % % %Multilayer criterion
tic
Nr = 1;
[MultiLayerMaskGFPRaw,mask_outGFP] = MultiLayerSpotIdentify(GFPRawMask,Nr);
toc

CytoGFPImg = double(gfp_raw_3D(:,:,1:RefineZ)).*MultiLayerMaskGFPRaw;

LLCytoSingleGFPMask = bwlabeln(MultiLayerMaskGFPRaw);
stsCytoSingleGFPMask = regionprops3(LLCytoSingleGFPMask,CytoGFPImg,'Volume','MaxIntensity','MeanIntensity');
GFPfociVolTest = [stsCytoSingleGFPMask.Volume];
GFPFociMaxIntTest = [stsCytoSingleGFPMask.MaxIntensity];
GFPFociMeanIntTest = [stsCytoSingleGFPMask.MeanIntensity];
%manual inspection
GFPFociMaxIntTest = GFPFociMaxIntTest(GFPFociMaxIntTest>2000);
%optional, or the threshold value can be set manually after inspection
[mu_log, sigma_log, mu_gauss, sigma_gauss, pi_log, pi_gauss, gamma_log, threshold] = mixed_em_gmm(GFPFociMaxIntTest, 'InitMethod', 'kmeans');




CytoGFPImg_Hmax = imhmax(CytoGFPImg,round(threshold));
CytoSingleGFPMask1 = imregionalmax(CytoGFPImg,26);

CytoSingleGFPMask2 = imregionalmax(CytoGFPImg_Hmax,26);
SE = strel('cube',2);
CytoSingleGFPMask3 = imdilate(CytoSingleGFPMask2,SE);
% se = strel('cube',2);
% MultiLayerMaskGFPRawOpen = imopen(MultiLayerMaskGFPRaw,se);
% LLMultiLayerMaskGFPRaw = bwlabeln(MultiLayerMaskGFPRawOpen);
LLCytoSingleGFPMask3 = bwlabeln(CytoSingleGFPMask3);
stsCytoSingleGFPMask3 = regionprops3(LLCytoSingleGFPMask3,gfp_raw_3DNor(:,:,1:RefineZ),'Volume','MaxIntensity');
GFPfociVol = [stsCytoSingleGFPMask3.Volume];
GFPFociMaxInt = [stsCytoSingleGFPMask3.MaxIntensity];
GFPFociInd = find(GFPfociVol>0);
GFPFociMaxInt2 = GFPFociMaxInt(GFPFociInd);
% % figure();
% % histogram(GFPFociMaxInt);
% % hold on
% % histogram(GFPFociMaxInt2);
CytoSingleGFPMask5 = ismember(LLCytoSingleGFPMask3,GFPFociInd);
MultiLayerMaskGFP = CytoSingleGFPMask5;
% CytoGFPImg = double(gfp_raw_3DNor(:,:,1:RefineZ)).*CytoSingleGFPMask5;
% CytoGFPImg_Hmax = imhmax(CytoGFPImg,1000);
% CytoSingleGFPMask1 = imregionalmax(CytoGFPImg,26);
% 
% CytoSingleGFPMask2 = imregionalmax(CytoGFPImg_Hmax,26);
% SE = strel('cube',3);
% CytoSingleGFPMask3 = imdilate(CytoSingleGFPMask2,SE);
% CytoSingleGFPMask4 = imdilate(MultiLayerMaskGFPRaw,SE);


%%

save MultiLayerMaskGFPBcdE1S2Em2_GCN4.mat MultiLayerMaskGFP -v7.3
% % % %%14:18, 12-20

% IndZrange = find(WCentroidGFP(:,3)>=12 & WCentroidGFP(:,3)<=20);
% WCentroidZRangegfp = WCentroidGFP(IndZrange,:);
% figure();
% imshow(AdjRNARaw1);
% hold on
% scatter(CytoRNACen(:,1),CytoRNACen(:,2),5,'r','filled');
% hold off
% exportgraphics(gcf,"GFP3D14_18.tiff")
% 
% input_MIP2 = max(gfp_raw_3DNor(:,:,zrange),[],3);
% LowHighRNARaw1 = stretchlim(uint16(input_MIP2));
% AdjRNARaw2 = imadjust(uint16(input_MIP2),[0 0.15],[]);
% figure();
% imshow(AdjRNARaw2);
% hold on
% BWMax = max(CytoSingleGFPMask5(:,:,zrange),[],3);
% [B,L] = bwboundaries(BWMax);
% 
% for k = 1:length(B)
%    boundary = B{k};
%    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
%    hold on
% end
% % % 
% % % hold off  
%%

%obtain the intensity and centroid of the GFP mask
LLMultiLayerMaskGFP = bwlabeln(MultiLayerMaskGFP);
stsMultiLayerMaskGFP = regionprops3(LLMultiLayerMaskGFP,gfp_raw_3DNor(:,:,1:RefineZ),'Volume','WeightedCentroid','VoxelValues');
TotGFPMaskWeiCenRaw = [stsMultiLayerMaskGFP.WeightedCentroid];
 TotVoxelValuesRaw = [stsMultiLayerMaskGFP.VoxelValues];
 TotGFPSpotIntensityRaw = zeros(length(TotGFPMaskWeiCenRaw),1);
for i = 1:length(TotGFPMaskWeiCenRaw)
    TotGFPSpotIntensityRaw(i) = sum([TotVoxelValuesRaw{i}]);
end

NANInd = isnan(TotGFPMaskWeiCenRaw);
SumNANInd = sum(NANInd,2);
PassInd = find(SumNANInd<1);
TotGFPMaskWeiCenRaw = TotGFPMaskWeiCenRaw(PassInd,:);
TotGFPSpotIntensityRaw = TotGFPSpotIntensityRaw(PassInd);
%verify the distribution of the intensity of the GFP spots
%  figure();
%  histogram(TotGFPSpotIntensity);
%  [Val2GFP, Edges2GFP] = histcounts(TotGFPSpotIntensity,'Normalization','probability');
% edges2GFP = Edges2GFP(2:end) - (Edges2GFP(2)-Edges2GFP(1))/2;
% [x, y] = prepareCurveData(edges2GFP,Val2GFP);
% figure();
% scatter(x,y,15,"black","filled");
%%
%plot
zrange = 57:60;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCenRaw(:,3)>=55 & TotGFPMaskWeiCenRaw(:,3)<=62);
WCentroidZRangegfp = TotGFPMaskWeiCenRaw(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMaskRaw_BcdE1S1Em3_GCN4_57_60.fig'));

zrange = 20:24;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCenRaw(:,3)>=18 & TotGFPMaskWeiCenRaw(:,3)<=26);
WCentroidZRangegfp = TotGFPMaskWeiCenRaw(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMaskRaw_BcdE1S2Em2_GCN4_20_24.fig'));

zrange = 5:10;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCenRaw(:,3)>=3 & TotGFPMaskWeiCenRaw(:,3)<=12);
WCentroidZRangegfp = TotGFPMaskWeiCenRaw(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMaskRaw_BcdE1S2Em2_GCN4_5_10.fig'));

zrange = 45:49;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCenRaw(:,3)>=43 & TotGFPMaskWeiCenRaw(:,3)<=51);
WCentroidZRangegfp = TotGFPMaskWeiCenRaw(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMaskRaw_BcdE1S2Em2_GCN4_45_49.fig'));
%%
LLMultiLayerMaskGFP = bwlabeln(MultiLayerMaskGFP);
stsMultiLayerMaskgfp = regionprops3(LLMultiLayerMaskGFP,gfp_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid','MeanIntensity','MaxIntensity');
VolGFP = [stsMultiLayerMaskgfp.Volume];
% WCentroidGFP = [stsMultiLayerMaskgfp.WeightedCentroid];
% MeanIntensityGFP =[stsMultiLayerMaskgfp.MeanIntensity];
% GauMaxIntensityGFP = [stsMultiLayerMaskgfp.MaxIntensity];

lb = [0,0,0,0,0,0];
ub = [inf,inf,inf,inf,inf,inf];
XYrange = 3;
Zrange = 2;
sigmax0 = 1.5;   %%% Initial value of sigma_x
sigmay0 = 1.5;   %%% Initial value of sigma_y
gau2D = @(x,xdata) x(1)*exp(-((xdata(:,1) - x(2)).^2./(2.*x(3).^2)+(xdata(:,2) - x(4)).^2./(2.*x(5).^2)))+x(6);
paraGFP = zeros(length(VolGFP),6);
% TotSpotVal = zeros(length(VolGFP),1);
resnormGFP = zeros(length(VolGFP),1);
exitflagGFP = zeros(length(VolGFP),1);
dim0 = size(gfp_raw_3D);
ZpositionGFP = zeros(length(VolGFP),1);
options = optimset('Display','off');
Xboundary = size_in_merge_image(1);
Yboundary = size_in_merge_image(2);
tic 
for i = 1:length(VolGFP)
LinIndex = find(LLMultiLayerMaskGFP == i);
[X0,Y0,Z0] = ind2sub(dim0,LinIndex);
XPeak = round(mean(X0));
YPeak = round(mean(Y0));
Zpeak = round(mean(Z0));
% Xmin = min(X0);
% Xmax = max(X0);
% Ymin = min(Y0);
% Ymax = max(Y0);
% Zmin = min(Z0);
% Zmax = max(Z0);
Xmin = max((XPeak-XYrange),1);
Xmax = min((XPeak+XYrange),Xboundary);
Ymin = max((YPeak-XYrange),1);
Ymax = min((YPeak+XYrange),Yboundary);
Zmin = max((Zpeak-Zrange),1);
Zmax = min((Zpeak+Zrange),RefineZ);
[X,Y] = ndgrid(Xmin:Xmax,Ymin:Ymax);
xdata = [X(:),Y(:)];
szz = length(Zmin:Zmax);
TotIntensityInZPlane = zeros(1,szz);
SpotVal = gfp_raw_3D(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
% TotSpotVal(i,:) = sum(sum(sum(SpotVal)));
for indz = 1:szz
   IntVal =  max(max(SpotVal(:,:,indz)));
   TotIntensityInZPlane(indz) = IntVal;
end
[V I] = max(TotIntensityInZPlane);
% % %
GauFilMaxSpotIntensity = SpotVal(:,:,I);
% %  GauFilMaxSpotIntensity = max(SpotVal,[],3);
% figure();
% imshow(GauFilMaxSpotIntensity,[]);
% exportgraphics(gcf,"SingleRNAPoint_30000.tiff")

 ydata = GauFilMaxSpotIntensity(:);
 Amp = (max(ydata) - min(ydata));
para0 = [Amp, XPeak, sigmax0, YPeak, sigmay0,min(ydata)];
[paraoutput,resnormoutput,~,exitflagoutput,~] = lsqcurvefit(gau2D,para0,xdata,ydata,lb,ub,options);
paraGFP(i,:) = paraoutput;
resnormGFP(i,:) = resnormoutput;
exitflagGFP(i,:) = exitflagoutput;
ZpositionGFP(i,:) = Zmin+I;
end
toc
save 2D_GFP_GaussianFittingBcdE1S2Em2_GCN4.mat paraGFP resnormGFP exitflagGFP


PSFStdRatioGFP = paraGFP(:,3)./paraGFP(:,5);

% GoodRNASpotInd = find(0.4 <para(:,3)& para(:,3)<3 & 0.4<para(:,5) & para(:,5)<3 & (exitflag == 1|exitflag == 3)&(sqrt(1-(para(:,3)./para(:,5)).^2) <= 0.9)&para(:,1)<100000 & para(:,1)>4000);
GoodGFPSpotInd = find(0.4 <paraGFP(:,3) & 0.4<paraGFP(:,5) & (exitflagGFP == 1|exitflagGFP == 3)&(PSFStdRatioGFP>0.1 & PSFStdRatioGFP<10)&paraGFP(:,1)<100000 & paraGFP(:,1)>1000);

GoodParaGFP = paraGFP(GoodGFPSpotInd,:);
GoodGFPMask = ismember(LLMultiLayerMaskGFP,GoodGFPSpotInd);
LLGoodGFPMask = bwlabeln(GoodGFPMask);
stsGoodGFP = regionprops3(LLGoodGFPMask,gfp_raw_3D(:,:,1:RefineZ),'WeightedCentroid');
GoodGFPWeiCen = [stsGoodGFP.WeightedCentroid];
%there may be nan value in the calculated weighted centroid
NANInd = isnan(GoodGFPWeiCen(:,1));
PassInd = find(NANInd < 1);
GoodGFPWeiCenRemoveNAN = GoodGFPWeiCen(PassInd,:);

% MatrixOfGFPSpot = zeros(length(PassInd),5);
% MatrixOfGFPSpot(:,1:3) = GoodGFPWeiCenRemoveNAN;
paraGoodGFP = GoodParaGFP(PassInd,:);
TotGFPSpotIntensity = 2.*pi.*paraGoodGFP(:,1).*paraGoodGFP(:,3).*paraGoodGFP(:,5);
% MatrixOfGFPSpot(:,4) = TotGFPSpotIntensity;
TotGFPMaskWeiCen = GoodGFPWeiCenRemoveNAN;

%%
%plot
zrange = 57:60;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCen(:,3)>=55 & TotGFPMaskWeiCen(:,3)<=62);
WCentroidZRangegfp = TotGFPMaskWeiCen(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMask_BcdE1S2Em2_GCN4_57_60.fig'));

zrange = 20:24;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCen(:,3)>=18 & TotGFPMaskWeiCen(:,3)<=26);
WCentroidZRangegfp = TotGFPMaskWeiCen(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMask_BcdE1S1Em3_GCN4_20_24.fig'));

zrange = 5:10;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCen(:,3)>=3 & TotGFPMaskWeiCen(:,3)<=12);
WCentroidZRangegfp = TotGFPMaskWeiCen(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMaskBcdE1S2Em2_GCN4_5_10.fig'));

zrange = 45:49;
input_MIP = max(gfp_raw_3DNor(:,:,zrange),[],3);
LowHighRNARaw1 = stretchlim(uint16(input_MIP));
AdjGFPRaw1 = imadjust(uint16(input_MIP),LowHighRNARaw1,[]);
% CytoRNALL = bwlabeln(CytoSingleRNAMask1(120:270,220:370,zrange));
% CytoRNALLsts = regionprops3(CytoRNALL,mRNA_raw_3D(120:270,220:370,zrange),'WeightedCentroid');
% CytoRNACen = [CytoRNALLsts.WeightedCentroid];
IndZrange = find(TotGFPMaskWeiCen(:,3)>=43 & TotGFPMaskWeiCen(:,3)<=51);
WCentroidZRangegfp = TotGFPMaskWeiCen(IndZrange,:);
figure();
imshow(AdjGFPRaw1);
hold on
scatter(WCentroidZRangegfp(:,1),WCentroidZRangegfp(:,2),5,'r','filled');
hold off
saveas(gcf, sprintf('GFPMask_BcdE1S2Em2_GCN4_45_49.fig'));
%%
%Find the nearest RNA spot to the GFP
% P = GoodRNAWeiCen;
xyPixelSize_nm = 69;
zStep_nm = 250;   
coordinateScale_nm = [xyPixelSize_nm, xyPixelSize_nm, zStep_nm];

P = TotGFPMaskWeiCen.*coordinateScale_nm;
T = delaunayn(P);
PQ = GoodRNAWeiCenRemoveNAN.*coordinateScale_nm;
%k the indices of the closest points in P to the query points in PQ measured in Euclidean distance
[kRNA2GFPGFP distRNA2GFP] = dsearchn(P,T,PQ);
%determine colocalization
%300 nm as the co-localization distance threshold
ColThre = 300;
Temp = kRNA2GFPGFP(distRNA2GFP < ColThre);
TempIso = kRNA2GFPGFP(distRNA2GFP >= ColThre);
OnRNAIndUnique = find(distRNA2GFP < ColThre);
OffRNAIndUnique = find(distRNA2GFP >= ColThre);
OnMatrixOfRNASpot = MatrixOfRNASpot(OnRNAIndUnique,:);
OFFMatrixOfRNASpot = MatrixOfRNASpot(OffRNAIndUnique,:);


%%
[SelectOnInd,SelectOnInda,SelectOnIndb] = intersect(RNASpotSelectInd,OnRNAIndUnique,'stable');
OnMatrixOfRNASpotSelZ = MatrixOfRNASpot(SelectOnInd,:);

OnRNANumForNucSelZ = zeros(length(ErodeNucWeiCen),1);
OnRNANumForNuc = zeros(length(ErodeNucWeiCen),1);
for i = 1:length(ErodeNucWeiCen)
    OnRNA2NucIndSelZ = find(OnMatrixOfRNASpotSelZ(:,6) == i);
    OnRNANumForNucSelZ(i,:) = sum(OnMatrixOfRNASpotSelZ(OnRNA2NucIndSelZ,5));
    OnRNA2NucInd = find(OnMatrixOfRNASpot(:,6) == i);
    OnRNANumForNuc(i,:) = sum(OnMatrixOfRNASpot(OnRNA2NucInd,5));
end

SingleOnGFPInt = TotGFPSpotIntensity(Temp);
SingleOnGFPCent = TotGFPMaskWeiCen(Temp,:);
SingleOnGFPCorresponseRNACount = MatrixOfRNASpot(OnRNAIndUnique,5);
SingleOnGFPTE = SingleOnGFPInt./SingleOnGFPCorresponseRNACount;
SingleGFPTETot = zeros(size(SingleOnGFPCent,1),size(SingleOnGFPCent,2)+3);
SingleGFPTETot(:,1:3) = SingleOnGFPCent;
SingleGFPTETot(:,4) = SingleOnGFPInt;
SingleGFPTETot(:,5) = SingleOnGFPCorresponseRNACount;
SingleGFPTETot(:,6) = SingleOnGFPTE;
SingleGFPTESelZInd = find(SingleGFPTETot(:,3)>=Zstart & SingleGFPTETot(:,3)<=Zend);
SingleGFPTESelZ = SingleGFPTETot(SingleGFPTESelZInd,:);

[SelOFFRNAInd,SelOFFRNAInda,SelOFFRNAIndb] = intersect(RNASpotSelectInd,OffRNAIndUnique,'stable');
OFFMatrixOfRNASpotSelZ = MatrixOfRNASpot(SelOFFRNAInd,:);
OnRNARatioPerNuc = OnRNANumForNuc./RNANumForNuc;
OnRNARatioPerNucSelZ = OnRNANumForNucSelZ./RNANumForNucSelZ;
fig = 1;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn_ON] = convertXY2AP(OnRNARatioPerNuc, nx, ny, emmask, fig ,flipAP);
 figure();
scatter(NormAPn_ON,OnRNARatioPerNuc,50,'filled');
xlim([0.15,0.65]);
ylim([0 1]);
xticks([0.2,0.3,0.4,0.5,0.6]);
yticks([0,0.2,0.4,0.6,0.8,1]);
xlabel("EL");
ylabel("The under translation mRNA ratio");
box on;
set(gcf,'position',[250 250 600 600]);
axis square
set (gca,'linewidth',2,'fontsize',24);
saveas(gcf, sprintf('OnRNAPerNucBcdE1S2Em2_GCN4.fig'));


RNANumForNucTrue = intens;
OnRNANumForNucTrue = CountScale.*OnRNANumForNuc;
RNANumForNucTrueSelZ = CountScale.*RNANumForNucSelZ;
OnRNANumForNucTrueSelZ = CountScale.*OnRNANumForNucSelZ;
%GFP indices
%Co-localized GFP 
P = GoodRNAWeiCenRemoveNAN.*coordinateScale_nm;
%  NANInd = isnan(P);
%  SumNANInd = sum(NANInd,2);
% PassInd = find(SumNANInd<1);
% P = P(PassInd,:);
T = delaunayn(P);
PQ = TotGFPMaskWeiCen.*coordinateScale_nm;
%k the indices of the closest points in P to the query points in PQ measured in Euclidean distance
[kGFP2RNA distGFP2RNA] = dsearchn(P,T,PQ);
%determine colocalization
%300 nm as the co-localization distance threshold

OnGFPIndUnique = find(distGFP2RNA < ColThre);
OffGFPIndUnique = find(distGFP2RNA >= ColThre);

Temp2 = kGFP2RNA(distGFP2RNA < ColThre);
TempIso2 = kGFP2RNA(distGFP2RNA >= ColThre);

ColGFPIntensity = TotGFPSpotIntensity(OnGFPIndUnique);
ColGFPWeiCen = TotGFPMaskWeiCen(OnGFPIndUnique,:);
TotGFPMatrix = zeros(size(TotGFPMaskWeiCen,1),size(TotGFPMaskWeiCen,2)+1);
TotGFPMatrix(:,1:3) = TotGFPMaskWeiCen;
TotGFPMatrix(:,4) = TotGFPSpotIntensity;
GFPCentroid_nm = TotGFPMatrix(:,1:3) .* coordinateScale_nm;
[kGFPtoNuc, distGFPtoNuc_nm] = dsearchn( ...
    NucCentroid_nm,TNuc,GFPCentroid_nm);
TotGFPMatrix(:,5) = kGFPtoNuc;

% ColGFPMatrixPre = TotGFPMatrixPre(OnGFPIndUnique,:);
% TotGFPMatrix = TotGFPMatrixPre(GFPSelectInd,:);
ColGFPMatrix = TotGFPMatrix(OnGFPIndUnique,:);

SingleOnRNAMat = MatrixOfRNASpot(Temp2,:);

SingleRNATE = TotGFPMatrix(OnGFPIndUnique,4)./MatrixOfRNASpot(Temp2,5);
SingleRNATETot = zeros(size(SingleOnRNAMat,1),size(SingleOnRNAMat,2)+3);
SingleRNATETot(:,1:3) = SingleOnRNAMat(:,1:3);
SingleRNATETot(:,4) = TotGFPMatrix(OnGFPIndUnique,4);
SingleRNATETot(:,5) = MatrixOfRNASpot(Temp2,5);
SingleRNATETot(:,6) = SingleRNATE;
SingleRNATESelZInd = find(SingleRNATETot(:,3)>=Zstart & SingleRNATETot(:,3)<=Zend);
SingleRNATESelZ = SingleRNATETot(SingleRNATESelZInd,:);
%Colocalized GFP per nuclei
ColGFPCentroid_nm = ColGFPMatrix(:,1:3) .* coordinateScale_nm;
[kColGFPtoNuc, distColGFPtoNuc_nm] = dsearchn( ...
    NucCentroid_nm,TNuc,ColGFPCentroid_nm);
ColGFPMatrix(:,5) = kColGFPtoNuc;
OnGFPSpotSelectInd = find(ColGFPMatrix(:,3)>=Zstart & ColGFPMatrix(:,3)<= Zend);
ColGFPMatrixSelZ = ColGFPMatrix(OnGFPSpotSelectInd,:);

%k the indices of the closest points in P to the query points in PQ measured in Euclidean distance

ColGFPForNucSelZ = zeros(length(ErodeNucWeiCen),1);
ColGFPForNuc = zeros(length(ErodeNucWeiCen),1);
for i = 1:length(ErodeNucWeiCen)
    CoGFP2NucIndSelZ = find(ColGFPMatrixSelZ(:,5) == i);
    ColGFPForNucSelZ(i,:) = sum(ColGFPMatrixSelZ(CoGFP2NucIndSelZ,4));
    CoGFP2NucInd = find(ColGFPMatrix(:,5) == i);
    ColGFPForNuc(i,:) = sum(ColGFPMatrix(CoGFP2NucInd,4));
end
fig = 0;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn] = convertXY2AP(ColGFPForNuc, nx, ny, emmask, fig,flipAP);
ColGFPForNucTrue = CountScale.*ColGFPForNuc;
ColGFPForNucTrueSelZ = CountScale.*ColGFPForNucSelZ;
Alpha_ON = ColGFPForNucTrue./OnRNANumForNucTrue;
Alpha_ONSelZ = ColGFPForNucTrueSelZ./OnRNANumForNucTrueSelZ;
Alpha = ColGFPForNucTrue./RNANumForNucTrue;
AlphaSelZ = ColGFPForNucTrueSelZ./RNANumForNucTrueSelZ;

figure();
scatter(NormAPn,Alpha_ON,'filled');
colorbar
colormap jet
% ColGFPForNucTrue = 0.5.*ColGFPForNuc;
figure();
scatter(NormAPn,ColGFPForNucTrue,20,'filled');
xlabel("EL");
ylabel("The protein production");
box on;
set(gcf,'position',[250 250 600 600]);
axis square
set (gca,'linewidth',2,'fontsize',24);
%%
%map the single molecule TE to the AP axis
%RNA-based TE
RNATEWetCen = SingleRNATETot(:,1:2);
RNATEWetCenScale(:,1) = RNATEWetCen(:,1).*ScaleFactor;
RNATEWetCenScale(:,2) = RNATEWetCen(:,2).*ScaleFactor;
U = RNATEWetCenScale;
TransRNATEWetCenScale = transformPointsForward(tform,U);
nx_TERNA = TransRNATEWetCenScale(:,1);
ny_TERNA = TransRNATEWetCenScale(:,2);
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn_TERNA] = convertXY2AP(SingleRNATETot(:,6), nx_TERNA, ny_TERNA, emmask, fig,flipAP);

%GFP-based TE
GFPTEWetCen = SingleGFPTETot(:,1:2);
GFPTEWetCenScale(:,1) = GFPTEWetCen(:,1).*ScaleFactor;
GFPTEWetCenScale(:,2) = GFPTEWetCen(:,2).*ScaleFactor;
U = GFPTEWetCenScale;
TransGFPTEWetCenScale = transformPointsForward(tform,U);
nx_TEGFP = TransGFPTEWetCenScale(:,1);
ny_TEGFP = TransGFPTEWetCenScale(:,2);
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn_TEGFP] = convertXY2AP(SingleGFPTETot(:,6), nx_TEGFP, ny_TEGFP, emmask, fig,flipAP);
save BcdE1S2Em2_GCN4_TEAP.mat NormAPn_TERNA NormAPn_TEGFP
%%


save BcdE1S2Em2_GCN4_TR.mat RNANumForNucTrue OnRNANumForNucTrue OnRNARatioPerNuc ColGFPForNucTrue Alpha_ON Alpha Alpha_ONSelZ AlphaSelZ NormAPn RNANumForNucTrueSelZ OnRNANumForNucTrueSelZ OnRNARatioPerNucSelZ ColGFPForNucTrueSelZ
save RNAMat_BcdE1S2Em2_GCN4.mat MatrixOfRNASpot
save OnRNAMat_BcdE1S2Em2_nc14_GCN4.mat OnMatrixOfRNASpot OnMatrixOfRNASpotSelZ
save OFFRNAMat_BcdE1S2Em2_nc14_GCN4.mat OFFMatrixOfRNASpot OFFMatrixOfRNASpotSelZ
% save TotGFPMat_EM2_20230221hbNC13_1_G2_10_selZ.mat TotGFPMatrix
save ColGFPMat_BcdE1S2Em2_GCN4.mat ColGFPMatrix ColGFPMatrixSelZ
save SingleGFPTESelZ_BcdE1S2Em2_GCN4.mat SingleGFPTETot SingleGFPTESelZ
save SingleRNATESelZ_BcdE1S2Em2_GCN4.mat SingleRNATETot SingleRNATESelZ

%%
IsoGFPMatrix = TotGFPMatrix(OffGFPIndUnique,:);
IsoGFPCountPerNuc = zeros(length(ErodeNucWeiCen),1);
TotGFPCountPerNuc = zeros(length(ErodeNucWeiCen),1);
for i = 1:length(ErodeNucWeiCen)
    IsoGFPCountInd = find(IsoGFPMatrix(:,5) == i);
    IsoGFPCountPerNuc(i,:) = length(IsoGFPCountInd);
    TotFPCountInd = find(TotGFPMatrix(:,5) == i);
    TotGFPCountPerNuc(i,:) = length(TotFPCountInd);
end
IsoGFPRatioPerNuc = ...
    IsoGFPCountPerNuc ./ TotGFPCountPerNuc;
[ventralData, dorsalData, APlength, AP_x, AP_y,NormAPn_GFPIso] = convertXY2AP(IsoGFPRatioPerNuc, nx, ny, emmask, fig ,flipAP);
save GFPIsoRatioAlongAP_BcdE1S2Em2.mat NormAPn_GFPIso IsoGFPRatioPerNuc


%%
%select the region for Pon alpha_on quantitative comparing along
%apical-basal axis
% zstepRange = 1:1:RefineZ;
% SRefineZ = RefineZ.*250;
% zstepRange = 1:1:RefineZ;
% zstepRangeReal =  250:250:SRefineZ;

zstepRange = 1:1:RefineZ;
RNACountAtEachZ = zeros(length(zstepRange)-1,1);
ONRNACountAtEachZ = zeros(length(zstepRange)-1,1);
NucMeanInt = zeros(length(zstepRange)-1,1);
NucMeanIntOri = zeros(length(zstepRange)-1,1);

OnRNANumPerZPerNuc = zeros(length(zstepRange)-1,length(ErodeNucWeiCen));
RNANumPerZPerNuc = zeros(length(zstepRange)-1,length(ErodeNucWeiCen));

for i = 1:(length(zstepRange)-1)   
ind = find(MatrixOfRNASpot(:,3) < zstepRangeReal(i+1) & MatrixOfRNASpot(:,3) >=zstepRangeReal(i));
RNACountAtEachZ(i) = sum(MatrixOfRNASpot(ind,5)); 
ind2 = find(OnMatrixOfRNASpot(:,3) < zstepRangeReal(i+1) & OnMatrixOfRNASpot(:,3) >=zstepRangeReal(i));
ONRNACountAtEachZ(i) = sum(OnMatrixOfRNASpot(ind2,5));
LLDAPIMaskZ = bwlabeln(ErodeNuc_true(:,:,zstepRange(i):zstepRange(i+1)));
stsDAPIMaskZ = regionprops3(LLDAPIMaskZ,dapi_raw_3D(:,:,zstepRange(i):zstepRange(i+1)),'WeightedCentroid','MeanIntensity');
MeanIntOfEachNuc = [stsDAPIMaskZ.MeanIntensity];
NucMeanInt(i) = sum(MeanIntOfEachNuc);
%
LLDAPIMaskOriZ = bwlabeln(Nuc_true(:,:,zstepRange(i):zstepRange(i+1)));
stsDAPIMaskOriZ = regionprops3(LLDAPIMaskOriZ,dapi_raw_3D(:,:,zstepRange(i):zstepRange(i+1)),'WeightedCentroid','MeanIntensity');
MeanIntOfEachNucOri = [stsDAPIMaskOriZ.MeanIntensity];
NucMeanIntOri(i) = sum(MeanIntOfEachNucOri);

for j = 1:length(ErodeNucWeiCen)
    RNA2NucIndPerZ = find(MatrixOfRNASpot(:,6) == j & MatrixOfRNASpot(:,3) < zstepRangeReal(i+1) & MatrixOfRNASpot(:,3) >=zstepRangeReal(i));
    RNANumPerZPerNuc(i,j) = sum(MatrixOfRNASpot(RNA2NucIndPerZ,5));
    OnRNA2NucIndPerZ = find(OnMatrixOfRNASpot(:,6) == j & OnMatrixOfRNASpot(:,3) < zstepRangeReal(i+1) & OnMatrixOfRNASpot(:,3) >=zstepRangeReal(i));
    OnRNANumPerZPerNuc(i,j) = sum(OnMatrixOfRNASpot(OnRNA2NucIndPerZ,5));
end
end


%GFP colocalized
clear i
for i = 1:(length(zstepRange)-1)   
ind = find(ColGFPMatrix(:,3) < zstepRangeReal(i+1) & ColGFPMatrix(:,3) >=zstepRangeReal(i));
ColGFPIntAtEachZ(i) = sum(ColGFPMatrix(ind,4)); 
end
%RNATEAtEachZ
clear i
for i = 1:(length(zstepRange)-1)   
ind = find(SingleRNATETot(:,3) < zstepRangeReal(i+1) & SingleRNATETot(:,3) >=zstepRangeReal(i));
RNATEAtEachZ(i) = sum(SingleRNATETot(ind,6)); 
end
%GFPTEAtEachZ
clear i
for i = 1:(length(zstepRange)-1)   
ind = find(SingleGFPTETot(:,3) < zstepRangeReal(i+1) & SingleGFPTETot(:,3) >=zstepRangeReal(i));
GFPTEIntAtEachZ(i) = sum(SingleGFPTETot(ind,6)); 
end 

save SignalAlongABAXis_BcdE1S2Em2_GCN4.mat NucMeanInt NucMeanIntOri ColGFPIntAtEachZ RNATEAtEachZ GFPTEIntAtEachZ RNACountAtEachZ ONRNACountAtEachZ RNANumPerZPerNuc OnRNANumPerZPerNuc
% figure();