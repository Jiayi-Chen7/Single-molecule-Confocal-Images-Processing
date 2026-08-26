function [maxval] = PointToLineDistance(InPutVector,AP_x,AP_y)

NumOfPoint = length(InPutVector);
AllCoord = [InPutVector]';
%the first point
FirstPoint = [AP_x(1),AP_y(1)]
endPoint = [AP_x(2),AP_y(2)]
%get the vrctor berweeen first and last point - this is the line
LineVec = endPoint - FirstPoint
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
VecToLine = VectorFromFirst - VectorFromFirstParallel

% DistanceToLine = sqrt(sum(VecToLine.^2,2));
DistanceToLine = norm(VecToLine);
% figure();
% plot(DistanceToLine);

[maxval IndexOfMax] = max(DistanceToLine);

end
