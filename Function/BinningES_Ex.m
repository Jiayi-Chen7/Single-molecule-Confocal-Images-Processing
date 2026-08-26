function BinningResults=BinningES_Ex(InputData,binpoint,numbin,outlierflag,nboots)
%UNTITLED2 Summary of this function goes here
% Bin the sorted data and each bin is equally spaced 
% also store the number of points in each bin 
% estimate the errors of the std by boottrapping 
% add outlierflag to get rid of outliner data
%   Detailed explanation goes here

if ~exist('numbin') || isempty(numbin), numbin = 100; end
if ~exist('outlierflag') || isempty(outlierflag), outlierflag = 0; end
if ~exist('nboots') || isempty(nboots), nboots = 50; end


[~, n]=size(InputData);
xmin=min(binpoint);
xmax=max(binpoint);
xstep=(xmax-xmin)/numbin;


BinningResults=zeros(numbin,3*n+1);

for i=1:numbin
    
x1=xmin+(i-1)*xstep;
x2=x1+xstep;

Data1=InputData(InputData(:,1)>=x1 & InputData(:,1)<x2,:);


 BinningResults(i,1:n)=mean(Data1,1);
 
 if  outlierflag
 for j=2:n
  Data1(abs(Data1(:,j)- BinningResults(i,j))>2*std(Data1(:,j),1),:)=[]   ;
 end


 BinningResults(i,1:n)=mean(Data1,1);
 end
 

for j=1:n
 BinningResults(i,n+j)=std(Data1(:,j),0,1);
end


    
   




numpnts=size(Data1,1);


for j=1:n
 if numpnts>2
     [bootsstd bootsam]= bootstrp(nboots,@std,Data1(:,j));
     bootsave=zeros(nboots,1);
     for k=1:nboots
     Datasam=Data1(bootsam(:,k),j);
     bootsave(k,1)=mean(Datasam);
     end    
     bootsRstd=bootsstd./bootsave;
     stdbootrp=std(bootsRstd);
 else
     stdbootrp=NaN;
 end
    BinningResults(i,2*n+j)=stdbootrp;
end

BinningResults(i,3*n+1)=numpnts;


end






