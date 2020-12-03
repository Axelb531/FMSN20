%% Init. Always Run
load( 'AL_data_2000_BP.mat' );
spde.C = speye(prod(sz));
spde.G = igmrfprec(sz,1);
spde.G2 = spde.G*spde.G;
I_obs = reshape(sum(A,1),sz);
I_grid = reshape(sum(A_grid,1),sz);
q_beta = 1e-6;
%% Specify B 
B;
%% Creaate X-mode
global x_mode;
x_mode = [];
par = fminsearch( @(theta) GMRF_negloglike(theta, Y, A, B, spde), [0 0]);
E_xy = x_mode;

%%
n_B = size(B,2);
tau = exp(par(1));
kappa2 = exp(par(2));

% CAR(1)
Q_x = tau*(kappa2*spde.C + spde.G);
Qbeta = q_beta*speye(n_B);
% 
Qtilde = blkdiag(Q_x,Qbeta);
Atilde = [A B];
[~, ~, Q_xy] = GMRF_taylor(E_xy, Y, Atilde, Qtilde);
%% Simulation
E_xy = A_grid * x_mode(1:end - n_B);
E_beta = x_mode(end - n_B + 1:end);
E_Bbeta = B_grid * E_beta;
E_zy = E_xy + E_Bbeta;
E_out = exp(E_zy);

%reuse taylor expansion to compute posterior precision

e = [zeros(size(Q_xy,1)-n_B, n_B); eye(n_B)];
V_beta0 = e'*(Q_xy\e);
beta_int_size = sqrt(diag(V_beta0)) * 1.96;

for i = 1:n_B
  fprintf(1, 'Beta %d: %11.4f +- %.4f\n', i, E_beta(i), beta_int_size(i));
end

% Plotting V_beta in boxplot:
boxplot(E_beta);

% Simulation over 1000 samples of the approximate posterior

Rxy = chol(Q_xy);
x_samp = repmat(x_mode(:),[1,1000]) + Rxy\randn(size(Rxy,1),1000);

%% Test
%Rxy \ randn(size(Rxy,1), 1000));
meanx_samp = mean(x_samp(1:end-n_B,:));
z_samp = Atilde*x_samp;
Vzy = var(z_samp,0,2);
Vx  = 1 / (size(x_samp,2) - 1) * sum(x_samp(1:end-n_B,:).^2,2);
Vbeta = sum((B_grid * V_beta0) .* B_grid, 2);
% Vzy = Vx + Vbeta;
Vout = exp(Vzy);

%% Plotting
% Complete
figure(3)
sigmoid_E_zy = exp(E_zy)./(1+exp(E_zy));
imagesc(reshape(A_grid'*sigmoid_E_zy, sz));
title('Estimated counts')
colorbar
%%
% Mean
figure(4);
subplot(1,2,1)
imagesc(reshape(A_grid'*E_Bbeta,sz));
title('Mean component')
colorbar
subplot(1,2,2)
imagesc(reshape(sqrt(A_grid'*Vbeta),sz));
title('Standard deviation')
colorbar

