function [fitresult, gof, E0] = gaussFit(x, y)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Jun-2019 12:13:57


%% Fit: 'gaussian'.
[xData, yData] = prepareCurveData( x, y );
indE0 = y == max(y);
E0 = x(indE0);
limX = round(0.98*indE0);

% Set up fittype and options.
ft = fittype( 'gauss1' );
excludedPoints = excludedata( xData, yData, 'Indices', 1:limX );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [1437948 E0 2.71819914850558];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'Gaussian fit' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'y vs. x', 'Excluded y vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel y
% grid on


