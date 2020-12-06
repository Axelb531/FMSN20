%% Data Init. Always Run
load( 'AL_data_2000_BP.mat' );
spde.C = speye(prod(sz));
spde.G = igmrfprec(sz,1);
spde.G2 = spde.G*spde.G;
I_obs = reshape(sum(A,1),sz);
I_grid = reshape(sum(A_grid,1),sz);
q_beta = 1e-6;
% Selection of model (Car:1, Sar:2)
alpha = 2;
%% Specify which betas we want  
B = [B(:,1) B(:,2) B(:,3) B(:,4) B(:,5) B(:,6)];
B_grid = [B_grid(:,1) B_grid(:,2) B_grid(:,3) B_grid(:,4) B_grid(:,5) B_grid(:,6)];
%% Evaluating Xmode
global x_mode;
x_mode = [];
par = fminsearch( @(theta) GMRF_negloglike(theta, Y, A, B, spde), [0 0]);
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
%


V_z  = var(z,0,2);
V_z_mean = mean(V_z);
% Vzy = Vx + Vbeta;
%Vout = exp(V_zy);

%% Plotting
figure(3)
imagesc(reshape(A_grid'*logsig(E_zy), sz),'alphadata', I_grid);
title('Estimated counts')
colorbar
%%
% Mean & Spatial
figure(4);
subplot(1,2,1)
imagesc(reshape(A_grid'*logsig(E_zy_spat),sz),'alphadata', I_grid);
title('Spatial Component')
colorbar
subplot(1,2,2)
imagesc( reshape(A_grid'*logsig(E_zy_mean),sz),'alphadata', I_grid)
title('Mean Component')
colorbar

