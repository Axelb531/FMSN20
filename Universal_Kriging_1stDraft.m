%% MODEL 1: Latitude, Elevation, Coast, CoastSE, CoastNO
  
Ck = [Longitude, Latitude]; % Coord's for known locations
Cv = [X_valid(:,1), X_valid(:,2)]; % Validation coord's for known locations
Cu = [X_grid(I_land, 1), X_grid(I_land, 2)]; % Coord's of unknown locations

% Distance matrices
Dkk = distance_matrix(Ck, Ck); 
Duk = distance_matrix(Cu, Ck); 
Duu = distance_matrix(Cu, Cu);
Dku = distance_matrix(Ck, Cu);
Dvk = distance_matrix(Cv, Ck);
%%
Kmax = 64; % Bins
Dmax = 0.5*max(Dkk(:));
                                        
[rhat,s2hat,m,n,d] = covest_nonparametric(Dkk,res,Kmax,Dmax); 
Nr = length(res); % Number of residuals

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
%% ML
[par_1, beta_Krig] = covest_ml(Dkk, Y, 'matern', par, X1, 'ml');
Sigma_kk = matern_covariance(Dkk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Dkk));
Sigma_uk = matern_covariance(Duk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Duk));
Sigma_uu = matern_covariance(Duu, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Duu));
Sigma_ku = matern_covariance(Dku, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Dku));
Sigma_vk = matern_covariance(Dvk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Dvk)); 

% Validation 
Y_hat_2 = X1_grid * beta_Krig + Sigma_uk*((Sigma_kk) \ (Y-X1*beta_Krig));
V_y_eta = Sigma_uu - Sigma_uk * (Sigma_kk \ Sigma_ku);
V_y_mu = (X1_grid' - X1'*(Sigma_kk \ Sigma_ku))' * ((X1'*(Sigma_kk \ X1)) \ (X1_grid' - X1'*(Sigma_kk \ Sigma_ku)));
V_y_Kriging = diag(V_y_eta + V_y_mu);
Y_ML_val =  X_true*beta_Krig + Sigma_vk*((Sigma_kk) \ (Y-X1*beta_Krig)); 
res_ML = Y_valid - Y_ML_val;
RS_ML = sum((res_ML).^2);
RMSE_ML = sqrt( mean( (Y_valid - Y_ML_val).^2 ) );
%% REML
[par_2, beta_Krig_2] = covest_ml(Dkk, Y, 'matern', par, X1, 'reml');
Sigma_kk_REML = matern_covariance(Dkk, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Dkk));
Sigma_uk_REML  = matern_covariance(Duk, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Duk));
Sigma_uu_REML  = matern_covariance(Duu, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Duu));
Sigma_ku_REML  = matern_covariance(Dku, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Dku));
Sigma_vk_REML  = matern_covariance(Dvk, par_2(1), par_2(2), par_2(3)) + par_2(4).*eye(size(Dvk));

% Validation
Y_hat_3 = X1_grid * beta_Krig_2 + Sigma_uk_REML *((Sigma_kk_REML ) \ (Y-X1*beta_Krig_2));
V_y_eta_2 = Sigma_uu_REML  - Sigma_uk_REML  * (Sigma_kk_REML  \ Sigma_ku_REML);
V_y_mu_2 = (X1_grid' - X1'*(Sigma_kk_REML  \ Sigma_ku_REML ))' * ((X1'*(Sigma_kk_REML  \ X1)) \ (X1_grid' - X1'*(Sigma_kk_REML  \ Sigma_ku_REML )));
V_y_Kriging_2 = diag(V_y_eta_2 + V_y_mu_2);
Y_REML_val = X_true*beta_Krig_2 + Sigma_vk_REML*((Sigma_kk_REML) \ (Y-X1*beta_Krig_2));
res_REML = Y_valid - Y_REML_val;
RMSE_REML = sqrt(mean((res_REML).^2)); 
RS_REML = sum((res_REML).^2);
%% OLS RMSE (used for comparison reasons)
Y_OLS_val = X_true*Beta;
res_OLS = Y_valid - Y_OLS_val;
RS_OLS = sum((res_OLS).^2);
RMSE_OLS = RMSE(7);
%% Plotting Y_hat using ML
Kriging = nan(sz);
Kriging(I_land) = Y_hat_2;
figure,
imagesc([4.9 30.9], [71.1 55.4], Kriging, 'alphadata', I_img)
axis xy; 
hold on;
plot(Border(:,1),Border(:,2))
colorbar;
% Plotting variance
V_Kriging = nan(sz);
V_Kriging(I_land) = sqrt(V_y_Kriging);
figure, 
imagesc([4.9 30.9], [71.1 55.4], V_Kriging, 'alphadata', I_img)
axis xy; 
hold on;
plot(Border(:,1),Border(:,2))
colorbar; 
%% Plotting Y_hat using REML
Kriging_2 = nan(sz);
Kriging_2(I_land) = Y_hat_3;
figure,
imagesc([4.9 30.9], [71.1 55.4], Kriging_2, 'alphadata', I_img)
axis xy; 
hold on;
plot(Border(:,1),Border(:,2))
colorbar;
% Plotting variance
V_Kriging_2 = nan(sz);
V_Kriging_2(I_land) = sqrt(V_y_Kriging_2);
figure, 
imagesc([4.9 30.9], [71.1 55.4], V_Kriging_2, 'alphadata', I_img)
axis xy; 
hold on;
plot(Border(:,1),Border(:,2))
colorbar;
%% Comparing different covariance functions
covfunc = 'spherical';
if strcmpi(covfunc,'gaussian')
    par_temp = zeros(3,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, 'gaussian', par_temp, X1, 'ml');
    Sigma_vk_covf = gaussian_covariance(Dvk, par_comp(1), par_comp(2)).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp(1) = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );
elseif strcmpi(covfunc,'exponential')
    par_temp = zeros(3,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, covfunc, par_temp, X1, 'ml');
    Sigma_vk_covf = exponential_covariance(Dvk, par_comp(1), par_comp(2)).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp(2) = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );
elseif strcmpi(covfunc,'matern')
    par_temp = zeros(4,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, covfunc, par_temp, X1, 'ml');
    Sigma_vk_covf = matern_covariance(Dvk, par_comp(1), par_comp(2), par_comp(3)) + par_comp(4).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp(3) = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );
elseif strcmpi(covfunc,'spherical')
    par_temp = zeros(3,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, covfunc, par_temp, X1, 'ml');
    Sigma_vk_covf = spherical_covariance(Dvk, par_comp(1), par_comp(2)).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp(4) = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );    
end
%% Plotting RMSE using different covariance functions
figure(2)
b = bar(RMSE_comp,'FaceColor', colors(5,:));
xlabel('Covariance function');
ylabel('Root Mean Square Error of Residues');
grid on;