%% MODEL 1: Latitude, Elevation, Coast, CoastSE, CoastNO
% Creating dataframe for model
X1 = [ones(n_x,1), Latitude , Elevation, Coast, Coast_SE, Coast_NO];

% Computing regression Parameters
[Beta_1, resid_1, Sigma2_1, Cov_1] = OLS_func(X1,Y);

% RMSE 
RMSE(1) = sqrt(mean((resid_1).^2) );

%% MODEL 2: Latitude, Elevation
X2 = [ones(n_x,1), Latitude, Elevation];
[Beta_2, resid_2, Sigma2_2, Cov_2] = OLS_func(X2,Y);
RMSE(2) = sqrt(mean((resid_2).^2) );

%% MODEL 3: Latitude, Elevation, Coast
X3 = [ones(n_x,1), Latitude, Elevation, Coast];
[Beta_3, resid_3, Sigma2_3, Cov_3] = OLS_func(X3,Y);
RMSE(3) = sqrt(mean((resid_3).^2) );

%% MODEL 4: Latitude, elevation, coastSE
X4 = [ones(n_x,1), Latitude, Elevation, Coast_SE];
[Beta_4, resid_4, Sigma2_4, Cov_4] = OLS_func(X4,Y);
RMSE(4) = sqrt(mean((resid_4).^2) );

%% MODEL 5: Latitude, elevation, coastNO
X5 = [ones(n_x,1), Latitude, Elevation, Coast_NO];
[Beta_5, resid_5, Sigma2_5, Cov_5] = OLS_func(X5,Y);
RMSE(5) = sqrt(mean((resid_5).^2) );

%% MODEL 6: Latitude, elevation, coastSE, coastNO
X6 = [ones(n_x,1), Latitude, Elevation,Coast_SE, Coast_NO];
[Beta_6, resid_6, Sigma2_6, Cov_6] = OLS_func(X6,Y);
RMSE(6) = sqrt(mean((resid_6).^2) );

%% MODEL 7: Elevation, coastSE, coastNO
X7 = [ones(n_x,1), Elevation, Coast_SE, Coast_NO];
[Beta_7, resid_7, Sigma2_7, Cov_7] = OLS_func(X7,Y);
RMSE(7) = sqrt(mean((resid_7).^2) );

%% Plotting RMSE
figure(2)
b = bar(RMSE,'FaceColor', colors(5,:));
xlabel('Model');
ylabel('Root Mean Square Error of Residues');
grid on;

% Model 1 is best fit

Beta = Beta_1;
Sigma2 = Sigma2_1;
Cov = Cov_1;
X1_grid = [ones(40820,1),X_grid(:,2:6)]; 
res = resid_1;
%% Validation of model
n_y = size(Y_valid, 1);
X_true = [ones(n_y, 1),X_valid(:,2:6)];
Y_hat = X_true*Beta;
sigma_2 = norm(Y_valid-Y_hat)^2/(n_y-size(X_true,2));
RMSE(7) = sqrt(mean((Y_valid-Y_hat).^2) );
%Confidence interval
ci_h = Y_valid+1.96*sqrt(var(Y_valid)/n_y);
ci_l = Y_valid-1.96*sqrt(var(Y_valid)/n_y);
figure(3);
plot(Y_valid, 'x');
hold on;
plot(Y_hat, '*');
%ylim([0,15]);
plot(ci_h ,'--b');
plot(ci_l,'--b');
legend;
%plotting Y-hat on map:
Y_hat_land = X1_grid*Beta;
Y_image = reshape(Y_hat_land, sz);
figure(4);
imagesc(Y_image); 
colorbar;
%%
Vbeta = Sigma2 * inv(X1' * X1);
Vmu = sum((X1_grid*Vbeta).*X1_grid,2);
V_y = Sigma2 + Vmu; 
V_y_image = reshape(V_y, sz);
V_std = sqrt(V_y_image);
figure(5);
imagesc(V_std); 
colorbar;
