%% MODEL 1: Latitude, Elevation, Coast, CoastSE, CoastNO
  
Ck = [Longitude, Latitude]; % Coord's for known locations
Cv = [X_valid(:,1), X_valid(:,2)]; % Validation coord's for known locations
Cu = [X_grid(I_land, 1), X_grid(I_land, 2)]; % Coord's of unknown locations
Dkk = distance_matrix(Ck, Ck); 
Duk = distance_matrix(Cu, Ck); 
Duu = distance_matrix(Cu, Cu);
Dku = distance_matrix(Ck, Cu);
Dvk = distance_matrix(Cv, Ck);
%%
Kmax = 64;
Dmax = 0.5*max(Dkk(:));
                                        
[rhat,s2hat,m,n,d] = covest_nonparametric(Dkk,res,Kmax,Dmax); 
Nr = length(res); % number of residuals

for ii = 1:100 % Bootstrapping 100 permutations
    i = randperm(Nr);
    rp = res(i);
    [rr(ii,:), s2hat,m,n,d] = covest_nonparametric(Ck, rp, Kmax, Dmax);
end
rr_std = std(rr);
I_top = mean(rr) + std(rr)*1.96;
I_bot = mean(rr) - std(rr)*1.96;

figure,
plot(rhat, 'k-o');
hold on
plot(I_top, '--b');
plot(I_bot, '--b');

par = zeros(4,1);
[par_1, beta_Krig] = covest_ml(Dkk, Y, 'matern', par, X1, 'ml');
%%
Sigma_kk = matern_covariance(Dkk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Dkk)) ;
Sigma_uk = matern_covariance(Duk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Duk)) ;
Sigma_uu = matern_covariance(Duu, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Duu));
Sigma_ku = matern_covariance(Dku, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Dku));
Sigma_vk = matern_covariance(Dvk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Dvk)); 
 
Y_hat_2 = X1_grid * beta_Krig + Sigma_uk*((Sigma_kk) \ (Y-X1*beta_Krig));
V_y_eta = Sigma_uu - Sigma_uk * (Sigma_kk \ Sigma_ku);
V_y_mu = (X1_grid' - X1'*(Sigma_kk \ Sigma_ku))' * ((X1'*(Sigma_kk \ X1)) \ (X1_grid' - X1'*(Sigma_kk \ Sigma_ku)));
V_y_Kriging = diag(V_y_eta + V_y_mu);
%%
Kriging = nan(sz);
Kriging(I_land) = Y_hat_2;
figure,
imagesc(Kriging); 
colorbar;

Y_krig = X_true*beta_Krig;
res_krig_sum = sum((Y_valid-Y_krig).^2);
res_ols_sum =  sum((Y_valid-Y_hat).^2);

V_Kriging = nan(sz);
V_Kriging(I_land) = sqrt(V_y_Kriging);
figure, 
imagesc(V_Kriging);
colorbar;

%
Y_ML_val =  X_true*beta_Krig + Sigma_vk*((Sigma_kk) \ (Y-X1*beta_Krig));
%

%% REML

[par_2, beta_Krig_2] = covest_ml(Dkk, Y, 'matern', par, X1, 'reml');

Sigma_kk = matern_covariance(Dkk, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Dkk));
Sigma_uk = matern_covariance(Duk, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Duk));
Sigma_uu = matern_covariance(Duu, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Duu));
Sigma_ku = matern_covariance(Dku, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Dku));

 
Y_hat_2 = X1_grid * beta_Krig_2 + Sigma_uk*((Sigma_kk) \ (Y-X1*beta_Krig_2));

V_y_eta = Sigma_uu - Sigma_uk * (Sigma_kk \ Sigma_ku);

V_y_mu = (X1_grid' - X1'*(Sigma_kk \ Sigma_ku))' * ((X1'*(Sigma_kk \ X1)) \ (X1_grid' - X1'*(Sigma_kk \ Sigma_ku)));

V_y_Kriging_2 = diag(V_y_eta + V_y_mu);
%%
Kriging_2 = nan(sz);
Kriging_2(I_land) = Y_hat_2;
figure,
imagesc(Kriging_2); 
colorbar;
%
V_Kriging_2 = nan(sz);
V_Kriging_2(I_land) = sqrt(V_y_Kriging_2);
figure, 
imagesc(V_Kriging_2);
colorbar;
%%
 
Y_REML_val = X_true*beta_Krig_2 + Sigma_vk*((Sigma_kk) \ (Y-X1*beta_Krig_2));
res_ML = Y_valid - Y_ML_val;
res_REML = Y_valid - Y_REML_val;
Y_OLS_val = X_true*Beta;
res_OLS = Y_valid - Y_OLS_val;
RS_ML = sum((res_ML).^2);
RS_REML = sum((res_REML).^2);
RS_OLS = sum((res_OLS).^2);
%%
RMSE_ML = sqrt(mean((res_ML).^2));
RMSE_REML = sqrt(mean((res_REML).^2));
RMSE_OLS = RMSE(1);
