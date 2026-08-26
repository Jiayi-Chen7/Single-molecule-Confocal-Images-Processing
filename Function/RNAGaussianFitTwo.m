function [fitresult, gof] = RNAGaussianFitTwo(edges2, Val2)
%CREATEFIT(EDGES2,VAL2)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : edges2
%      Y Output: Val2
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 07-Apr-2022 15:16:15 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( edges2, Val2 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [1461 105000 57037.7008546008 807.841902699953 195000 69221.4560801877];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end


