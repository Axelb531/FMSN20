function [RMSE_comp] = Comp_covf(covfunc)
%COMP_COVF Summary of this function goes here
%   Detailed explanation goes here
load HA1_temp_2015_06
Longitude = X(:,1);
Latitude = X(:,2);
Elevation = X(:,3);
Coast = X(:,4);
Coast_SE = X(:,5);
Coast_NO = X(:,6);
n_x = size(X,1); 
X1 = [ones(n_x,1), Latitude , Elevation, Coast, Coast_SE, Coast_NO];
Ck = [Longitude, Latitude]; % Coord's for known locations
Cv = [X_valid(:,1), X_valid(:,2)]; % Validation coord's for known locations
Dkk = distance_matrix(Ck, Ck); 
Dvk = distance_matrix(Cv, Ck);

if strcmpi(covfunc,'gaussian')
   % covfunc = 'gaussian';
    par_temp = zeros(3,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, 'gaussian', par_temp, X1, 'ml');
    Sigma_vk_covf = gaussian_covariance(Dvk, par_comp(1), par_comp(2)).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );
elseif strcmpi(covfunc,'exponential')
    par_temp = zeros(3,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, covfunc, par_temp, X1, 'ml');
    Sigma_vk_covf = exponential_covariance(Dvk, par_comp(1), par_comp(2)).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );
elseif strcmpi(covfunc,'matern')
    par_temp = zeros(4,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, covfunc, par_temp, X1, 'ml');
    Sigma_vk_covf = matern_covariance(Dvk, par_comp(1), par_comp(2), par_comp(3)) + par_comp(4).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );
elseif strcmpi(covfunc,'spherical')
    par_temp = zeros(4,1); 
    [par_comp, beta_Krig_comp] = covest_ml(Dkk, Y, covfunc, par_temp, X1, 'ml');
    Sigma_vk_covf = spherical_covariance(Dvk, par_comp(1), par_comp(2), par_comp(3)) + par_comp(4).*eye(size(Dvk)); 
    Y_ML_temp =  X_true*beta_Krig_comp + Sigma_vk_covf*((Sigma_kk) \ (Y-X1*beta_Krig_comp)); 
    RMSE_comp = sqrt( mean( (Y_valid - Y_ML_temp).^2 ) );    
end

   


