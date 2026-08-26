function [ElbowPointIndexRelative] = GetElbowIndex(InPutVector)
%reference:
%modified from Jona's answer on Stackoverflow
%https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve?utm_source=chatgpt.com
[PeakPointVal, IndPeakPoint] = max(InPutVector);
NumOfPoint = length(InPutVector(IndPeakPoint:end));
AllCoord = [1:NumOfPoint;InPutVector(IndPeakPoint:end)]';
%the first point
FirstPoint = AllCoord(1,:);
%get the vector betweeen first and last point - this is the line
LineVec = AllCoord(end,:) - FirstPoint;
%Normalize the line vectorto obtain the unit vector of the line
NorLineVec = LineVec / sqrt(sum(LineVec.^2));

%find the distance from each point to the line
%vector between all points and first point
VectorFromFirst = bsxfun(@minus,AllCoord, FirstPoint);
%obtain the vector parallel to the line for each point vector on curve
%obtain the scalar by point product of the unit vectot in the
%direction of the line vector and the vectors to be projected to the line
%vector
Scalar = dot(VectorFromFirst,repmat(NorLineVec,NumOfPoint,1),2);
VectorFromFirstParallel = Scalar .* NorLineVec;
VecToLine = VectorFromFirst - VectorFromFirstParallel;

DistanceToLine = sqrt(sum(VecToLine.^2,2));
% figure();
% plot(DistanceToLine);

[maxval IndexOfMax] = max(DistanceToLine);

TureIndexOfMax = IndexOfMax+(IndPeakPoint - 1);
ElbowPointIndexRelative = TureIndexOfMax;
end
