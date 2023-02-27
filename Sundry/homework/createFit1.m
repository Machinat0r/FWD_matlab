%CREATEFIT1(B,U)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : B
%      Y Output: U
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 09-Dec-2022 19:08:23 自动生成
close all

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( B, MFS_07e(:,end) );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
figure( 'Name', 'untitled fit 1' );
s1 = scatter(xData,yData,100,'p','filled','MarkerFaceColor','#7E2F8E');hold on;
x = 0:1:180;
plot(x,fitresult.p1*x+fitresult.p2,'color','#8470FF');hold on;
% h = plot( fitresult, xData, yData,'color','#8470FF');
legend('Measuring Result', 'Transfer Function', 'Location', 'NorthWest');
text(141,28.1,{['SSE:',num2str(gof.sse)];['R:',num2str(gof.rsquare)];['RMSE:',num2str(gof.rmse)]}); hold on;
% Label axes
xlabel( 'B [nT]');
ylabel( 'V_{ind} [mV]');
grid on
box on

