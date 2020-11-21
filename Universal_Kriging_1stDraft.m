%% MODEL 1: Latitude, Elevation, Coast, CoastSE, CoastNO
 
Ck = [Longitude, Latitude]; % Coord's for known locations
Cv = [X_valid(:,1), X_valid(:,2)]; % Validation coord's for known locations
Cu = [X_grid(:,1), X_grid(:,2)]; % Coord's of unknown locations
Dk = distance_matrix(Ck);
Dvk = distance_matrix(Cv, Ck);
Duk = distance_matrix(Cu, Ck);

Kmax = 64;
Dmax = 0.5*max(Dk(:));

[rhat,s2hat,m,n,d] = covest_nonparametric(Dk,res,Kmax,Dmax); 
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
[par_1, beta_test] = covest_ml(Dk, Y, 'matern', par, X1, 'ml');
%%
Sigma_kk = matern_covariance(Dk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Dk)) ;
Sigma_uk = matern_covariance(Duk, par_1(1), par_1(2), par_1(3)) + par_1(4).*eye(size(Duk)) ;

Y_hat_2 = X1_grid * beta_test + Sigma_uk*((Sigma_kk) \ (Y-X1*beta_test));

Y_image_2 = reshape(Y_hat_2, sz);
figure,
imagesc(Y_image_2); 
colorbar;

Y_krig = X_true*beta_test;
res_krig_sum = sum((Y_valid-Y_krig).^2);
res_ols_sum =  sum((Y_valid-Y_hat).^2);