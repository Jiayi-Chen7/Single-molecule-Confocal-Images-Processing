function [Err] = TranslationRate(para,xdata,ydata)
% optimization :~
% pop=lsqcurvefit(@(p,x) kveri(p),p0,[0],[0])
%
%%
%   youtput = (1./(1+(para(1)./xdata).^para(2))).*(1./(1+(xdata./para(3)).^para(4)));
% youtput = (1./(1+(30./xdata).^para(1))).*(1./(1+(xdata./300).^para(2)));
%   youtputpre = para(1).*exp(-para(2).*(log(78.7523./xdata))).*xdata2;
%   youtputpre = 0.0872.*(1-(1.0864.*(log(78.7523./(xdata))).^2))+para(1);
%  youtputpre = para(1).*exp(-0.7868.*(log(78.7523./xdata))).*(0.0872.*(1-(1.0864.*(log(78.7523./(xdata))).^2))+0.5322).*xdata;
%  youtputpre = para(1).*exp(0.7868.*(log(78.7523./xdata))).*(0.0872.*(1-(1.0864.*(log(78.7523./(xdata))).^2))+0.5322).*xdata;
%   youtputpre = para(1).*exp(-para(2).*(log(para(3)./xdata))).*xdata;
youtputpre = para(1).*((83.7995./xdata).^(-para(2))).*(0.0507.*(1-(2.3810.*(log(83.7995./(xdata))).^2))+0.3529).*xdata;
% youtputpre =  para(1).*(1-( 0.2474.*(log( 78.7523./(xdata))).^2));
%   PosInd = find(youtputpre>0);
%   youtput = youtputpre(PosInd);
%   youtputinter = interp1(xdata(PosInd),youtput,xdata,'linear','extrap');
%   Err = youtputpre - ydata;
  youtputpre2 = 0.5.*youtputpre+abs(0.5.*youtputpre);
  Err = youtputpre2 - ydata;
%   youtput = para(1).*xdata.*(log(xdata./para(2))).^2;
end