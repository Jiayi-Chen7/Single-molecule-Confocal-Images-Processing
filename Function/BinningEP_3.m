function BinningResults=BinningEP_3(InputData,dataperbin,outlierflag,nboots)
%UNTITLED2 Summary of this function goes here
% Bin the sorted data and each bin is equally populated
%   Detailed explanation goes here
% Sorting data

if ~exist('outlierflag') || isempty(outlierflag), outlierflag = 0; end
if ~exist('nboots') || isempty(nboots), nboots = 50; end


[~,Idx]=sort(InputData(:,1));
InputData=InputData(Idx,:);


[numpnt n]=size(InputData);
%numbin=floor(numpnt/dataperbin);
numbin=round(numpnt/dataperbin);

BinningResults=zeros(numbin,3*n+1);

for i=1:numbin
    
 if i~=numbin  
 Data1=InputData(dataperbin*(i-1)+1:dataperbin*i,:); 
 else
 Data1=InputData(dataperbin*(i-1)+1:numpnt,:); 
 end
 NANInd = isnan(Data1);
 infInd = isinf(Data1)
PassInd = find(sum(NANInd,2) < 1 & sum(infInd,2) < 1);
Data1 = Data1(PassInd,:); 

BinningResults(i,1:n)=nanmean(Data1,1);


if  outlierflag
    
 for j=2:n
  Data1(abs(Data1(:,j)- BinningResults(i,j))>3*nanstd(Data1(:,j),1),:)=[]   ;
 end
 
 BinningResults(i,1:n) = nanmean(Data1,1);
 end

for j=1:n
 BinningResults(i,n+j)=nanstd(Data1(:,j),0,1);
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



