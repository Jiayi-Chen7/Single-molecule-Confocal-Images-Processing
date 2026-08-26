clear all
close all
clc

%%
%initial parameters
%folder
in_folder = 'E:/CJY/20220314/EM1/';
mid_sur_in_folder = 'E:/CJY/20220314/EM1/MID SUR/';
out_folder =  'E:/CJY/20220314/EM1/output/';
%the image names of zstack
dapi_channel = '20220314thomasaqualEM1_EM1_z*_ch00.tif';
gfp_channel = '20220314thomasaqualEM1_EM1_z*_ch01.tif';
mRNA_channel = '20220314thomasaqualEM1_EM1_z*_ch02.tif';
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
RefineZ = (numel(imlist_dapi)-5);
DapiRawMask = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
GFPRawMask = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
mRNARawMask =zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
mRNAFociRawMask = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
mRNAGauFil = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
GFPGauFil = zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
thresholdvalue = zeros(1,RefineZ);
%estimate the nuclei cycle of the embryo
NCEst = 14;

tic
fig = 0;
parfor ImageMaskIndex = 1:RefineZ
EmSingleOri = dapi_raw_3D(:,:,ImageMaskIndex);
% [Mask T] = Nuclei_seg_cjy(EmSingleOri, ExVal, IxVal,low_th,cir_th,Nuc_th,ImageMaskIndex,NCEst);
[Mask] = Nuclei_seg_cjy(EmSingleOri, [], [],[],[],ImageMaskIndex,NCEst,[],[],fig);
DapiRawMask(:,:,ImageMaskIndex) = Mask;
end
toc



% % % for j = 1:RefineZ
% % % figure();
% % % [B,L] = bwboundaries(DapiRawMask(:,:,j));
% % % imshow(dapi_raw_3D(:,:,j),[]);
% % % hold on
% % % for k = 1:length(B)
% % %    boundary = B{k};
% % %    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
% % %    hold on
% % % end
% % %      hold off      
% % % %      saveas(gcf, sprintf('RawNucleiWithMask%d.tif',j));
% % % end

%remove yolk nuclei and refine
BWDapiRawMask = bwlabeln(DapiRawMask);
StsDapiRawMask = regionprops3(BWDapiRawMask,dapi_raw_3D(:,:,1:RefineZ),'Volume','MeanIntensity');
NucVol = [StsDapiRawMask.Volume];
NucMeanInt = [StsDapiRawMask.MeanIntensity];
% % figure();
% % histogram(NucVol);
%the threshold is estimate from the prior knowledge of the size of a nuclei
%the diameter of ~ 5um: 5 ./ 0.069 um/pixel ~ 72
%should be observed in at least 5 layers of the z stacks
TrueNucVol = 20000;
TrueNucInd = find(NucVol > TrueNucVol& NucMeanInt>30000);
Nuc_true = ismember(BWDapiRawMask,TrueNucInd);
LLNuc_true = bwlabeln(Nuc_true);
stsNuc_true = regionprops3(LLNuc_true,dapi_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid');
NuccVol = [stsNuc_true.Volume];
NucWeiCen = [stsNuc_true.WeightedCentroid];
length(NucVol);
%shrink the muclei mask for more accurate nuclei assignment of the mRNA and
%GFP spot
se = strel('cube',3);
ErodeNuc_true = imerode(Nuc_true,se);
LLErodeNuc_true = bwlabeln(ErodeNuc_true);
stsLLErodeNuc_true = regionprops3(LLErodeNuc_true,dapi_raw_3D(:,:,1:RefineZ),'Volume','WeightedCentroid');
ErodeNucVol = [stsLLErodeNuc_true.Volume];
ErodeNucWeiCen = [stsLLErodeNuc_true.WeightedCentroid];
% % % parfor j =  1:RefineZ
% % % figure();
% % % [B,L] = bwboundaries(Nuc_true(:,:,j));
% % % imshow(dapi_raw_3D(:,:,j),[]);
% % % hold on
% % % for k = 1:length(B)
% % %    boundary = B{k};
% % %    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
% % %    hold on
% % % end
% % %      hold off      
% % %      saveas(gcf, sprintf('RefineNucleiWithMask%d.tif',j));
% % % end
%%ISOSURFACE
% figure();
% imshow(bw_true);
% figure();
% isosurface(DapiRawMask(1:1024,1:1024,1:28));
% colormap("hsv");
% 

%%
%%RNA mask processing 
%test to obtain the absolute single mRNA intensity threshold and the size
%of the transcription foci
ExRNACore = 1.2;
IxRNACore = 2.2;

ExFociCore = 1.2;
InFociCore = 15;%15
seval = 5;
MIPRNA = max(mRNA_raw_3D,[],3);
MaxNuc_true = max(Nuc_true,[],3);
[MaxTranscriptionFoci,MaxTranscriptionFociMask] = MaxFoci(MIPRNA,ExFociCore,InFociCore,seval);
% se = strel('disk',3);
% MaxTranscriptionFociMask2 = imopen(MaxTranscriptionFociMask,se);

% % % figure();
% % % imshow(MIPRNA,[]);
% % % hold on
% % % [B,L] = bwboundaries(MaxTranscriptionFociMask);
% % % [BN,LN] = bwboundaries(MaxNuc_true);
% % % for k = 1:length(B)
% % %    boundary = B{k};
% % %    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
% % %    hold on
% % % end
% % % for k = 1:length(BN)
% % %    boundary = BN{k};
% % %    plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 0.5);
% % %    hold on
% % % end
% % %      hold off   

parfor img_ind = 1:RefineZ
SinglemRNA = mRNA_raw_3D(:,:,img_ind);
SinglePlaneNuc = Nuc_true(:,:,img_ind);
mRNAFociRawMask(:,:,img_ind) = SinglePlaneFoci(SinglemRNA,SinglePlaneNuc,ExFociCore,InFociCore,img_ind,MaxTranscriptionFociMask);
end
Nr = 1;
mRNAFociMask = MultiLayerSpotIdentify(mRNAFociRawMask,Nr);



% % % parfor img_ind = 1:RefineZ
% % % figure();
% % % imshow(mRNA_raw_3D(:,:,img_ind),[]);
% % % hold on
% % % [B,L] = bwboundaries(mRNAFociMask(:,:,img_ind));
% % % [BF,LF] = bwboundaries(MaxTranscriptionFociMask);
% % % [BN,LN] = bwboundaries(Nuc_true(:,:,img_ind));
% % % for k = 1:length(B)
% % %    boundary = B{k};
% % %    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
% % %    hold on
% % % end
% % % for k = 1:length(BF)
% % %    boundary = BF{k};
% % %    plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 0.7);
% % %    hold on
% % % end
% % % for k = 1:length(BN)
% % %    boundary = BN{k};
% % %    plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 0.5);
% % %    hold on
% % % end
% % % hold off      
% % % saveas(gcf, sprintf('FociMask%d.tif',img_ind));
% % % end

tic
parfor RNAStackIndex = 1:RefineZ
SinglemRNA = mRNA_raw_3D(:,:,RNAStackIndex);
SinglePlaneNuc = Nuc_true(:,:,RNAStackIndex);
mRNAFociRawMaskPlane = mRNAFociRawMask(:,:,RNAStackIndex);
[CytoSingleRNAMask,GauFil] = SinglePlaneRNAMask_cjy(SinglemRNA,SinglePlaneNuc,ExRNACore,IxRNACore,RNAStackIndex,mRNAFociRawMaskPlane);
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
RawRNAAfterMask = mRNA_raw_3D(:,:,1:RefineZ).*MultiLayerMask;

LocalMaxRNAMask =zeros(size_in_merge_image(1),size_in_merge_image(2),RefineZ);
tic 
parfor RNAStackIndex = 1:RefineZ
PreLocalMaxRNAMask = imregionalmax(RawRNAAfterMask(:,:,RNAStackIndex));
% LocalMaxRNAMask(:,:,RNAStackIndex) = bwmorph(PreLocalMaxRNAMask,'shrink',Inf);
LocalMaxRNAMask(:,:,RNAStackIndex) = PreLocalMaxRNAMask;
end
toc
% save MultiLayerMask.mat MultiLayerMask -v7.3
LLMultiLayerMask = bwlabeln(LocalMaxRNAMask);
stsMultiLayerMask = regionprops3(LLMultiLayerMask,mRNAGauFil,'Volume','WeightedCentroid','MeanIntensity','MaxIntensity');
VolRNA = [stsMultiLayerMask.Volume];
WCentroid = [stsMultiLayerMask.WeightedCentroid];
MeanIntensity =[stsMultiLayerMask.MeanIntensity];
GauMaxIntensity = [stsMultiLayerMask.MaxIntensity];
%%
%%gaussian fitting of the candidate RNA spots in cytoplasmic
lb = [0,0,0,0,0,0];
ub = [inf,inf,inf,inf,inf,inf];
XYrange = 5;
Zrange = 2;
sigmax0 = 1.5;   %%% Initial value of sigma_x
sigmay0 = 1.5;   %%% Initial value of sigma_y
gau2D = @(x,xdata) x(1)*exp(-((xdata(:,1) - x(2)).^2./(2.*x(3).^2)+(xdata(:,2) - x(4)).^2./(2.*x(5).^2)))+x(6);
para = zeros(length(VolRNA),6);
TotSpotVal = zeros(length(VolRNA),1);
resnorm = zeros(length(VolRNA),1);
exitflag = zeros(length(VolRNA),1);
dim0 = size(mRNA_raw_3D(:,:,1:28));
Zposition = zeros(length(VolRNA),1);
options = optimset('Display','off');
Xboundary = size_in_merge_image(1);
Yboundary = size_in_merge_image(2);
tic 
parfor i = 1:length(VolRNA)
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
%  GauFilMaxSpotIntensity = max(SpotVal,[],3);
% figure();
% imshow(GauFilMaxSpotIntensity,[]);
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
save 2DGaussianFittingRaw.mat para resnorm exitflag
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
% % % hold off
% % % %end of the plot

%6 columes correspond to the xyz axis in X63 lens, the integral intensity,
% the RNA count derived from the typical intensity of single mRNA and the
% corresponding nuclei indices

%criteria reference: Heng Xu, et al. Nature Method, 2015: code: spfilter.m
GoodRNASpotInd = find(0.4 <para(:,3)& para(:,3)<3 & 0.4<para(:,5) & para(:,5)<3 & (exitflag == 1|exitflag == 3)&(sqrt(1-(para(:,3)./para(:,5)).^2) <= 0.9)&para(:,1)<100000);
MatrixOfRNASpot = zeros(length(GoodRNASpotInd),6);

paraGood = para(GoodRNASpotInd,:);
GoodRNAMask = ismember(LLMultiLayerMask,GoodRNASpotInd);
LLGoodRNAMask = bwlabeln(GoodRNAMask);
stsGoodRNA = regionprops3(LLGoodRNAMask,mRNA_raw_3D(:,:,1:RefineZ),'WeightedCentroid');
GoodRNAWeiCen = [stsGoodRNA.WeightedCentroid];
MatrixOfRNASpot(:,1:3) = GoodRNAWeiCen;
TotRNASpotIntensity = 2.*pi.*paraGood(:,1).*paraGood(:,3).*paraGood(:,5);
MatrixOfRNASpot(:,4) = TotRNASpotIntensity;
% figure();
% histogram(TotRNASpotIntensity);
[Val2, Edges2] = histcounts(TotRNASpotIntensity);
edges2 = Edges2(2:end) - (Edges2(2)-Edges2(1))/2;
%fit 3 terms gaussian function to obtain the typical intensity of the
%single cytoplamic mRNA
[RNAfitresult, RNAgof] = RNAGaussianFit(edges2, Val2);
%plot
[x, y] = prepareCurveData(edges2,Val2);
a1 = RNAfitresult.a1;
b1 = RNAfitresult.b1;
c1 = RNAfitresult.c1;
a2 = RNAfitresult.a2;
b2 = RNAfitresult.b2;
c2 = RNAfitresult.c2;
a3 = RNAfitresult.a3;
b3 = RNAfitresult.b3;
c3 = RNAfitresult.c3;
x = edges2;
GauTerm1 =  a1.*exp(-((x-b1)./c1).^2);
GauTerm2 =  a2.*exp(-((x-b2)./c2).^2);
GauTerm3 = a3.*exp(-((x-b3)./c3).^2);
Gautot = a1.*exp(-((x-b1)./c1).^2)+a2.*exp(-((x-b2)./c2).^2)+ a3.*exp(-((x-b3)./c3).^2);
figure();
plot(x,GauTerm3,'LineWidth',2);
hold on
plot(x,GauTerm1,'LineWidth',2);
hold on
plot(x,GauTerm2,'LineWidth',2);
hold on
plot(x,Gautot, 'LineWidth',2);
hold on
scatter(x,y,15,"black","filled");
hold off
saveas(gcf, sprintf('ThreeTermGaufitEM1.fig'));

 [RNAfitresult2, RNAgof2] = RNAGaussianFitTwo(edges2, Val2);
 [x2, y2] = prepareCurveData(edges2,Val2);
a11 = RNAfitresult2.a1;
b11 = RNAfitresult2.b1;
c11 = RNAfitresult2.c1;
a21 = RNAfitresult2.a2;
b21 = RNAfitresult2.b2;
c21 = RNAfitresult2.c2;
x = edges2;
GauTerm12 =  a11.*exp(-((x-b11)./c11).^2);
GauTerm22 =  a21.*exp(-((x-b21)./c21).^2);

Gautot2 = a11.*exp(-((x-b11)./c11).^2)+a21.*exp(-((x-b21)./c21).^2);
figure();
plot(x2,GauTerm12,'LineWidth',2);
hold on
plot(x2,GauTerm22,'LineWidth',2);
hold on
plot(x2,Gautot2, 'LineWidth',2);
hold on
scatter(x2,y2,15,"black","filled");
hold off
saveas(gcf, sprintf('TwoTermGaufitEM1.fig'));
TypicalIntensity = b3;
MatrixOfRNASpot(:,5) = round(TotRNASpotIntensity./TypicalIntensity);
%  save GoodRNAMask.mat GoodRNAMask
%%
%assign the mRNA spot in cytoplasm to the nearest nuclei
%http://www.qhull.org/

P = ErodeNucWeiCen;
T = delaunayn(P);
PQ = GoodRNAWeiCen;
%k the indices of the closest points in P to the query points in PQ measured in Euclidean distance
[kRNA distRNA] = dsearchn(P,T,PQ);

MatrixOfRNASpot(:,6) = kRNA;
RNANumForNuc = zeros(length(ErodeNucVol),1);
parfor i = 1:length(ErodeNucVol)
    RNA2NucInd = find(MatrixOfRNASpot(:,6) == i);
    RNANumForNuc(i,:) = sum(MatrixOfRNASpot(RNA2NucInd,5));
end
% mean(RNANumForNuc);
% figure();
% histogram(RNANumForNuc);
%%
%register the data
%Nuclei channel registration to obtain the trasnfrom matrix
%register the 3D image to the
%%surface image
MidEm = imread(['20220314ThomasAqual-1_em1-2_ch00.tif']);
SurEm = imread(['20220314ThomasAqual-1_em1_ch00.tif']);
RotateMidEm = rot90(MidEm,-1);
RotateSurEm = rot90(SurEm,-1);
ScaleFactor = 20./(63*1.3);
%shrink the nuclei image for accurate feature extraction
MIPNuc = max(dapi_raw_3D(:,:,1:RefineZ),[],3);
MIPErodeNuc_true = max(ErodeNuc_true,[],3);
AA = MIPNuc.*MIPErodeNuc_true;
ScaleMIPNuc = imresize(AA,ScaleFactor);
%shoule transform the ScaleMIPNuc to uint16(ScaleMIPNuc)
[MOVINGREG,ErrorMteric] = registerImages_cjy(uint16(ScaleMIPNuc),RotateSurEm);
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
saveas(gcf, sprintf('ImageReg.tif'));
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
RotateMidEmFil = imgaussfilt(RotateMidEm,50);
[T EM] = graythresh(RotateMidEmFil);
EMMaskRaw = imbinarize(RotateMidEm,T);
seEm = strel("disk",5);
EMMaskRawclose = imclose(EMMaskRaw,seEm);
EMMaskRawclose = imfill(EMMaskRawclose,'holes');
seEm2 = strel("disk",50);
EMMaskRawcloseopen = imopen(EMMaskRawclose,seEm2);
% % % % figure();
% % % % imshow(RotateMidEm);
% % % % hold on
% % % % [B,L] = bwboundaries(EMMaskRawcloseopen);
% % % % for k = 1:length(B)
% % % %    boundary = B{k};
% % % %    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 0.5);
% % % %    hold on
% % % % end
% % % % hold off

%%
fig = 1;
intens = RNANumForNuc;
nx = TransErodeNucWeiCen(:,1);
ny = TransErodeNucWeiCen(:,2);
emmask = EMMaskRawcloseopen;
save ConvertAPData.mat intens nx ny emmask
% 
% intens = ConvertAPData.intens;
% nx = ConvertAPData.nx;
% ny = ConvertAPData.ny;
% emmask = ConvertAPData.emmask;
[ventralData, dorsalData, APlength, AP_x, AP_y] = convertXY2AP(intens, nx, ny, emmask, fig);
%[ventralData dorsalData APlength] = convertXY2AP(NucleiIntensity, NucleiPosition(:,1), NucleiPosition(:,2), embryomask, fig);
APPos = ventralData(:,1);
mRNACount = mean(ventralData(:,3)+dorsalData(:,3));
figure();
scatter(APPos,mRNACount);
saveas(gcf, sprintf('APMmRNACount.fig'));