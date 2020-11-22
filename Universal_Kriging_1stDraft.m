%% MODEL 1: Latitude, Elevation, Coast, CoastSE, CoastNO
 

Ck = [Longitude, Latitude]; % Coord's for known locations
Cv = [X_valid(:,1), X_valid(:,2)]; % Validation coord's for known locations
Cu = [X_grid(:,1), X_grid(:,2)]; % Coord's of unknown locations
Dkk = distance_matrix(Ck); 
Duk = distance_matrix(Cu, Ck);
Duu = distance_matrix(Cu);
Dku = distance_matrix(Ck, Cu);

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
  
Y_hat_2 = X1_grid * beta_Krig + Sigma_uk*((Sigma_kk) \ (Y-X1*beta_Krig));

V_y_eta = Sigma_uu - Sigma_uk * (Sigma_kk \ Sigma_ku);
V_y_mu = X_true' - (X1'*(Sigma_kk \ Sigma_ku))' * (X1'*(Sigma_kk \ X1)) \ (X_true' - X1'*(Sigma_kk \ Sigma_ku));

Y_image_2 = reshape(Y_hat_2, sz);
figure,
imagesc(Y_image_2); 
colorbar;

Y_krig = X_true*beta_Krig;
res_krig_sum = sum((Y_valid-Y_krig).^2);
res_ols_sum =  sum((Y_valid-Y_hat).^2);