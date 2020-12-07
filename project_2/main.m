%% Data Init. Always Run
load( 'AL_data_2000_BP.mat' );
spde.C = speye(prod(sz));
spde.G = igmrfprec(sz,1);
spde.G2 = spde.G*spde.G;
I_obs = reshape(sum(A,1),sz);
I_grid = reshape(sum(A_grid,1),sz);
q_beta = 1e-6;
% Car or Sar?
alpha = 1;
% Selection of model (Car:1, Sar:2)
%% Specify which betas we want  
B = [B(:,1), B(:,3), B(:,4), B(:,6)];
B_grid = [B_grid(:,1), B_grid(:,3), B_grid(:,4), B_grid(:,6)];

%% Evaluating Xmode
global x_mode;
x_mode = [];
par = fminsearch( @(theta) GMRF_negloglike(theta, Y, A, B, spde, q_beta, alpha), [0 0]);
E_xy = x_mode;
%%
n_B = size(B,2);
tau = exp(par(1));
kappa2 = exp(par(2));
if alpha==1
  Q_x = tau*(kappa2*spde.C + spde.G);
elseif alpha==2
  Q_x = tau*(kappa2*kappa2.*spde.C + 2*kappa2*spde.G + spde.G^2);
end
Q_beta = q_beta*eye(n_B);
Qtilde = blkdiag(Q_x,Q_beta);
Atilde = [A B];
[~, ~, Q_xy] = GMRF_taylor(E_xy, Y, Atilde, Qtilde);
E_beta = x_mode(end - n_B + 1:end);
%% Simulation 1000 samples of the approximate posterior
Rxy = chol(Q_xy);
x_samp = repmat(x_mode(:),[1,1000]) + Rxy\randn(size(Rxy,1),1000);
% Applying diffrent results
Atilde_full= [A_grid B_grid];
Atilde_spatial= [A_grid zeros(size(B_grid))];
Atilde_mean= [zeros(size(A_grid)) B_grid];

E_zy = Atilde_full*x_mode;

E_zy_spat = Atilde_spatial*x_mode;

E_zy_mean = Atilde_mean*x_mode;

z = Atilde_full*x_samp;
z_spat = Atilde_spatial*x_samp;
z_mean = Atilde_mean*x_samp;
% Variance of reconstruction in correct dimension
V_zy  = var(z,0,2);
% Mean of the variance
mean_Vzy = mean(V_zy);
% Variance of mean and spatial component

V_zy_spat = var(z_spat,0,2);
V_zy_mean = var(z_mean,0,2);

CI95_zy = 1.96*sqrt(V_zy);

CI95_zy_spat = 1.96*sqrt(V_zy_spat);
CI95_zy_mean = 1.96*sqrt(V_zy_mean);

%% Plotting
figure(3)
imagesc(longitude, latitude([2 1]),reshape(A_grid'*logsig(E_zy), sz),'alphadata', I_grid);
title('Reconstruction (Probability)')
axis xy tight
colorbar
%ax = gca;
%ax.TitleFontSizeMultiplier = 3;
set(gca,'FontSize',30)
set(gca,'TitleFontSizeMultiplier',1.5)
%%
% Mean & Spatial
% Logsig for transformation to probability
figure(4);
subplot(1,2,1)
imagesc(longitude, latitude([2 1]),reshape(A_grid'*logsig(E_zy_spat),sz),'alphadata', I_grid);
axis xy tight
title('Spatial Component')
colorbar
set(gca,'FontSize',30)
subplot(1,2,2)
imagesc(longitude, latitude([2 1]), reshape(A_grid'*logsig(E_zy_mean),sz),'alphadata', I_grid);
axis xy tight
title('Mean Component')
colorbar
set(gca,'FontSize',30)
%%
% Plotting CI for A_tilde_full, spatial and mean
figure(5);
sgtitle('CI95 - Reconstruction (Probability)');
set(gca,'FontSize',30)
subplot(1,2,1)
imagesc(longitude, latitude([2 1]),reshape(A_grid'*logsig(E_zy-CI95_zy),sz),'alphadata', I_grid);
axis xy tight
title('95% CI - Lower limit')
colorbar
set(gca,'FontSize',30)
subplot(1,2,2)
imagesc(longitude, latitude([2 1]), reshape(A_grid'*logsig(E_zy+CI95_zy),sz),'alphadata', I_grid);
axis xy tight
title('95% CI - Upper limit')
colorbar
set(gca,'FontSize',30)

%%


figure(6);
sgtitle('CI95 - Spatial Component (Probability)');
subplot(1,2,1)
imagesc(longitude, latitude([2 1]),reshape(A_grid'*logsig(E_zy_spat-CI95_zy_spat),sz),'alphadata', I_grid);
axis xy tight
title('95% CI - Lower limit')
colorbar
set(gca,'FontSize',30)
subplot(1,2,2)
imagesc(longitude, latitude([2 1]), reshape(A_grid'*logsig(E_zy_spat+CI95_zy_spat),sz),'alphadata', I_grid);
axis xy tight
title('95% CI - Upper limit')
set(gca,'FontSize',30)
colorbar

%%
figure(7);
sgtitle('CI95 - Mean Component (Probability)');
set(gca,'FontSize',30)
subplot(1,2,1)
imagesc(longitude, latitude([2 1]),reshape(A_grid'*logsig(E_zy_mean-CI95_zy_mean),sz),'alphadata', I_grid);
axis xy tight
title('95% CI - Lower limit (Probability) ')
colorbar
set(gca,'FontSize',30)
subplot(1,2,2)
imagesc(longitude, latitude([2 1]),reshape(A_grid'*logsig(E_zy_mean+CI95_zy_mean),sz),'alphadata', I_grid);
axis xy tight
title('95% CI - Upper limit (Probability)')
colorbar
set(gca,'FontSize',30)