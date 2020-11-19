function [Beta, res, sigma2, Cov] = OLS_func(X,y)
% Amount of Data-points
N = length(y);
% Number of columns 
p = size(X,2);
% Determine Beta-coeff
Beta = X\y;
% Determine R-matrix
R = X'*X;
% Calculating residues
res = y-X*Beta;
% Calculating Variance of residues
sigma2 = norm(y-X*Beta)/(N-p);
% And finally covariance matrix
Cov = sigma2 * inv(R);