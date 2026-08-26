function [d] = CalRNANearestRNA(x,xi)
%%
%the definition of the x, tri and xi are the same as the dsearchn in matlab
%DistanceRange determine the number of mRNA in the neighborhood to
%calculate the RNA-RNA distance distribution
%xi: query point

    d = zeros(size(xi,1),1);
    pass = zeros(size(xi,1),1);
    for i = 1:size(xi,1) 
        yi = repmat(xi(i,:),size(x,1),1);
        d_temp = sqrt(sum((x-yi).^2,2));
        Sort_d_temp = sort(d_temp);
        % 
        % if length(d_temp) < DistanceRange
        %     TrueDistanceRange = length(d_temp);
        %     pass(i) = 1;
        % 
        % else
        %     TrueDistanceRange = DistanceRange;
        
        % d(i) = sum(Sort_d_temp(2:TrueDistanceRange));
     
        d(i) = Sort_d_temp(2);
    end
   
       
end