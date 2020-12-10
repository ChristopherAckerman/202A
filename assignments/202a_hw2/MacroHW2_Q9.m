readtable('simulated_log_linear_economy.csv');

%%drop first 100 obs
ans(1:101,:)=[]; 
c = ans{:, 2};

%% Fit the model
% Set model 1 AR lags, 0 multiplicative components, and 0 MA lags
Mdl2 = arima(1, 0, 0); 
EstMdl2 = estimate(Mdl2, c);

%infer the residuals
res=infer(EstMdl2, c);

%plot autocorrelation in residuals;
autocorr(res);
