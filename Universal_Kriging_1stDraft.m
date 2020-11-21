%% MODEL 1: Latitude, Elevation, Coast, CoastSE, CoastNO
% Creating dataframe for model
X1 = [ ones(n_x,1), Latitude , Elevation, Coast, Coast_SE, Coast_NO];
Y_k = Y; % Observed locations
dist_ = distance_matrix(X1);
sigma2 = Sigma2_1;
kappa = 1;
nu = 1;
Sigma = matern_covariance(dist, sigma2, kappa, nu);
