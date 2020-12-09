function testLifeCycleHypothesis()
% Tests the Life Cycle-Permanent Income Hypothesis
%
% Tests the theory laid out in Hall JPE 1978. Download the quarterly real
% consumption from FRED and fits a AR(1) model to see if it is reasonable.
%
% Part of HW 2 in 202 A with Chris Ackerman, Ekaterina Gurkova, and Ben Prie

% Ali Haider Ismail, 2020

%% Setup
close all;
status = license('test', 'datafeed_toolbox');
if ~status
  error('Datafeed toolbox does not exist. Cannot download data dynamically from Fred');
end
startDate = '01/01/1950';
endDate = '12/31/2019';
series = 'PCECC96';
url = 'https://fred.stlouisfed.org/';
c = fred(url);
modelOrder = 1;

%% Download the data
% Dwonload Quarterly Real PCE from fred within appropriate dates
rawData = fetch(c, series, startDate, endDate);
% Clean up the data that was downloaded
dates = datetime(rawData.Data(:, 1), 'ConvertFrom', 'datenum');
data = timetable(rawData.Data(:, 2), log(rawData.Data(:, 2)), 'RowTimes', dates);
data.Properties.VariableNames = {'RPCE', 'lnRPCE'};
Y = data{:, 2};

%% Fit the model
% Set model modelOrder AR lags, 0 multiplicative components, and 0 MA lags
Mdl = arima(modelOrder, 0, 0); 
EstMdl = estimate(Mdl, Y);

%% Test that the model fits well
% Estimate the R squared
[~, ~, E] = calculateR2(EstMdl, Y);

% Obtain a plot of the autocorrelation in the data
autocorr(Y);
exportgraphics(gcf, 'data-autocorrelation-plot.pdf');
close;

% Check that there is autocorrelation in the residuals
estimate(Mdl, E);
autocorr(E);
exportgraphics(gcf, 'residual-autocorrelation-plot.pdf');
close;

%% Appendix 1 - what would a true AR(1) look like
testMdl = arima('Constant', EstMdl.Constant, 'AR', EstMdl.AR);
testMdl.Variance = EstMdl.Variance;
[testY, testE] = simulate(testMdl, 1000);
autocorr(testY)
exportgraphics(gcf, 'data-simulated-autocorrelation-plot.pdf');
close;

autocorr(testE)
exportgraphics(gcf, 'residual-simulated-autocorrelation-plot.pdf');
close;

%% Appendix 2 - does using AR(2) fit better?
% Based on Hall's paper, if the true model is AR(1), AR(2) should not do better
modelOrder = 2;
Mdl = arima(modelOrder, 0, 0); 
EstMdl = estimate(Mdl, Y);
calculateR2(EstMdl, Y);

%% Run the Box-Ljyung test
% see Matlab documentation and Lee's notes for more info
stdE = E/sqrt(EstMdl.Variance); % Standardized residuals
lags = 10;
dof = lags - modelOrder; % One autoregressive parameter
[~, pValue] = lbqtest(stdE, 'Lags', lags, 'DOF', dof); 
fprintf(['The pValue of whether to reject the null hypothesis that there\n' ...
    'is no autocorrelation for 10 lags in the residuals is %f\n'], pValue);

end

function [Rsq, Rsqadj, E] = calculateR2(EstMdl, Y)
% Calculates the R^2 and adjusted R^2
%
% Inputs:
%   EstMdl - Arima model
%     Estimated Model output from estimate function
%   Y - numeric column data or numeric matrix
%     Response data
%
% Source:
%   https://stackoverflow.com/a/56497638/5101261

% Ali Haider Ismail, 2020

%% Set up
n = length(Y);

%% Get residuals
E = infer(EstMdl, Y);

%% Compute statistics
SSquares = dot(E,E);
Stotal = dot(Y - mean(Y), Y - mean(Y));
Rsq = 1 - SSquares/Stotal;
Rsqadj = 1 - (1-Rsq)*(n-1)/(n-2);

fprintf('R-Squared is %f, Adjusted R-Squared is %f\n', Rsq, Rsqadj);

end

