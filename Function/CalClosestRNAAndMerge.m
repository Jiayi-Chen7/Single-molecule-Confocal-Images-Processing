function [newx] = CalClosestRNAAndMerge(x,xi,DistanceRange,RemoveThre)
%%
%the definition of the x, tri and xi are the same as the dsearchn in matlab
%DistanceRange determine the number of mRNA in the neighborhood to
%calculate the RNA-RNA distance distribution
%xi: query point
newx = zeros(size(x));
    d = zeros(size(xi,1),1);
    pass = zeros(size(xi,1),1);
    for i = 1:size(xi,1) 
        yi = repmat(xi(i,:),size(x,1),1);
        d_temp = sqrt(sum((x-yi).^2,2));
        [Sort_d_temp xiInd] = sort(d_temp);
         
        if length(d_temp) < DistanceRange
            TrueDistanceRange = length(d_temp);
            pass(i) = 1;
         
        else
            TrueDistanceRange = DistanceRange;
       
         end
      k = find(Sort_d_temp(2:TrueDistanceRange)<RemoveThre);
      if isempty(k) == 1
          newx(i,:) = xi(i);
      else
      outputxiInd = xiInd(k);
       outputxi = xi(outputxiInd,:);
       CombinSet = [xi(i,:);outputxi];
       newx(i,:) = mean(CombinSet);
      end
    end
end