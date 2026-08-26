function [min_d] =  NearestNeighborDistance(X)
%%
%from New Being

  d = pdist(X);
d = squareform(d);
d(d==0) = NaN;
% min_d = min(d);
%cjy
[M I] = min(d,[],'all','omitnan');
min_d = M;
   
       
end