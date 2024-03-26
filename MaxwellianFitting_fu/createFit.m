clear;clc
close all
B = [0,8.60000000000000,34.8500000000000,54.7000000000000,73.1500000000000,91.2000000000000,135.600000000000,180.050000000000];
I = [0,5,20,30,40,50,75,100];
%% Plot1-Fit
[xData, yData] = prepareCurveData( I, B );

% Set up fittype and options.
ft = fittype( {'x'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a'} );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
% Plot fit with data.
figure(1);
x = linspace(0,100,1000);
% h = plot(x,1.807*x,'color','#4DBEEE','LineWidth',1.5);hold on;
% h = scatter(xData,yData,'k');
h = plot( fitresult,'k',xData, yData);
% plot(xData, yData,'k') 
legend( h, 'B vs. I', 'B = 1.807I', 'Location', 'NorthEast', 'Interpreter', 'none' );
str = {'SSE:4','R-square:0.9998','RMSE:0.8165'};
text(5,175,str);
% Label axes
xlabel( 'I [mA]', 'Interpreter', 'none' );
ylabel( 'B [nT]', 'Interpreter', 'none' );
grid on

%% Plot2-Residuals
figure(2)
% h = plot( fitresult, xData, yData, 'residuals' );
err = 100*(1.807*xData-yData)./xData;
errorbar(xData,0*xData,err)
% legend( h, 'untitled fit 1 - residuals', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'I [mA]', 'Interpreter', 'none' );
ylabel( 'dB/B [%]', 'Interpreter', 'none' );
grid on


