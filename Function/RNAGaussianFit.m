function [RNAfitresult, RNAgof] = RNAGaussianFit(edges2, Val2)
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

%  由 MATLAB 于 07-Apr-2022 13:56:22 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( edges2, Val2 );

% Set up fittype and options.
ft = fittype( 'gauss3' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [1461 105000 38025.1339030672 1065.84190269995 165000 41859.5992728578 671.663712736221 55000 76093.1397077353];

% Fit model to data.
[RNAfitresult, RNAgof] = fit( xData, yData, ft, opts );

end


