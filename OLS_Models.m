%% MODEL 1: latitude, elevation, coast, coastSE, coastNO
x_lat = X(:,2);
x_elevation = X(:,3);
x_coast = X(:,4);
x_coastSE = X(:,5);
x_coastNO = X(:,6);
X_1 = [ones(425,1), x_lat, x_elevation, x_coast, x_coastSE, x_coastNO];
beta_hat_1 = (X_1'*X_1) \ (X_1'*Y);
X_valid_1 = X_valid;
X_valid_1(:,1) = [];
sigma2_hat_1 = sum(( (Y - X_1*beta_hat_1).^2) /  (length(Y) - size(X,2) - 1));
Y_hat_1 = X_valid_1 .* beta_hat_1; % Nugget?
%accuracy = sum((Y_valid-Y_hat_1).^2);
%sigma2_hat_1 = res_1'*res_1 / (length(Y_hat_1) - size(X,2) - 1);

% Computing the covariance matrix
Sigma = sigma2_hat_1 * inv(X_1'*X_1);
% The uncertainty in beta_hat_1 is
X_1_transpose = X_1';
inv_X_1 = inv(X_1_transpose*X_1);
var_beta_hat_1 = sigma2_hat_1 * inv_X_1; 
Y_hat_1_nugg =  X_valid_1 * beta_hat_1 + randn(106,1)*sigma2_hat_1; % 106 = size(X_valid_1 * beta_hat_1)
% RMSE 
RMSE_1 = sqrt(mean( (Y_valid - Y_hat_1).^2) );

%% MODEL 2 latitude, elevation
X_2 = [x_lat, x_elevation];
beta_hat_2 = (X_2.'*X_2) \ (X_2.'*Y);
X_valid_2 = X_valid;
X_valid_2(:,1) = [];
X_valid_2(:,3) = [];
X_valid_2(:,3) = [];
X_valid_2(:,3) = [];
 
Y_hat_2 = X_valid_2 * beta_hat_2;
res_2 = Y_valid - Y_hat_2;
sigma2_hat_2 = res_2.'*res_2 / (length(Y_hat_2) - size(X,2) - 4);

% The uncertainty in beta_hat_1 is
X_2_transpose = X_2';
inv_X_2 = inv(X_2_transpose*X_2);
var_beta_hat_2 = sigma2_hat_2 * inv_X_2; 

% RMSE 
RMSE_2 = sqrt(mean( (res_2).^2) );

%% MODEL 3 latitude, elevation, coast
X_3 = [x_lat, x_elevation, x_coast];
beta_hat_3 = (X_3.'*X_3) \ (X_3.'*Y);
X_valid_3 = X_valid;
X_valid_3(:,1) = [];
X_valid_3(:,4) = [];
X_valid_3(:,4) = [];
 
Y_hat_3 = X_valid_3 * beta_hat_3;
res_3 = Y_valid - Y_hat_3;
sigma2_hat_3 = res_3.'*res_3 / (length(Y_hat_3) - size(X,2) - 3);

% The uncertainty in beta_hat_1 is
X_3_transpose = X_3';
inv_X_3 = inv(X_3_transpose*X_3);
var_beta_hat_3 = sigma2_hat_3 * inv_X_3; 

% RMSE 
RMSE_3 = sqrt(mean( (res_3).^2) );

%% MODEL 4 latitude, elevation, coastSE
X_4 = [x_lat, x_elevation, x_coastSE];
beta_hat_4 = (X_4.'*X_4) \ (X_4.'*Y);
X_valid_4 = X_valid;
X_valid_4(:,1) = [];
X_valid_4(:,3) = [];
X_valid_4(:,4) = [];

 
Y_hat_4 = X_valid_4 * beta_hat_4;
res_4 = Y_valid - Y_hat_4;
sigma2_hat_4 = res_4.'*res_4 / (length(Y_hat_4) - size(X,2) - 3);

% The uncertainty in beta_hat_1 is
X_4_transpose = X_4';
inv_X_4 = inv(X_4_transpose*X_4);
var_beta_hat_4 = sigma2_hat_4 * inv_X_4; 

% RMSE 
RMSE_4 = sqrt(mean( (res_4).^2) );

%% MODEL 5 latitude, elevation, coastNO
X_5 = [x_lat, x_elevation, x_coastNO];
beta_hat_5 = (X_5.'*X_5) \ (X_5.'*Y);
X_valid_5 = X_valid;
X_valid_5(:,1) = [];
X_valid_5(:,3) = [];
X_valid_5(:,3) = [];
 
Y_hat_5 = X_valid_5 * beta_hat_5;
res_5 = Y_valid - Y_hat_5;
sigma2_hat_5 = res_5.'*res_5 / (length(Y_hat_5) - size(X,2) - 3);

% The uncertainty in beta_hat_1 is
X_5_transpose = X_5';
inv_X_5 = inv(X_5_transpose*X_5);
var_beta_hat_5 = sigma2_hat_5 * inv_X_5; 

% RMSE 
RMSE_5 = sqrt(mean( (res_5).^2) );

%% MODEL 6 latitude, elevation, coastSE, coastNO
X_6 = [x_lat, x_elevation, x_coastSE, x_coastNO];
beta_hat_6 = (X_6.'*X_6) \ (X_6.'*Y);
X_valid_6 = X_valid;
X_valid_6(:,1) = []; 
X_valid_6(:,3) = [];
 
Y_hat_6 = X_valid_6 * beta_hat_6;
res_6 = Y_valid - Y_hat_6;
sigma2_hat_6 = res_6.'*res_6 / (length(Y_hat_6) - size(X,2) - 2);

% The uncertainty in beta_hat_1 is
X_6_transpose = X_6';
inv_X_6 = inv(X_6_transpose*X_6);
var_beta_hat_6 = sigma2_hat_6 * inv_X_6; 

% RMSE 
RMSE_6 = sqrt(mean( (res_6).^2) );

%% MODEL 7 elevation, coastSE, coastNO
X_7 = [x_elevation, x_coastSE, x_coastNO];
beta_hat_7 = (X_7.'*X_7) \ (X_7.'*Y);
X_valid_7 = X_valid;
X_valid_7(:,1) = []; 
X_valid_7(:,1) = []; 
X_valid_7(:,2) = [];
 
Y_hat_7 = X_valid_7 * beta_hat_7;
res_7 = Y_valid - Y_hat_7;
sigma2_hat_7 = res_7.'*res_7 / (length(Y_hat_7) - size(X,2) - 3);

% The uncertainty in beta_hat_1 is
X_7_transpose = X_7';
inv_X_7 = inv(X_7_transpose*X_7);
var_beta_hat_7 = sigma2_hat_7 * inv_X_7; 

% RMSE 
RMSE_7 = sqrt(mean( (res_7).^2) );