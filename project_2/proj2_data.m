%% load data
load( 'AL_data_2000_BP.mat' )
alpha = 2;
% B vector containing standard data on human land use
B = [B(:,1)  B(:,2)  B(:,3)  B(:,4) B(:,5) B(:,6), B(:,7)];
% B vector with logitransformed human land use
%B = [B(:,1)  B(:,3)  B(:,4)  B(:,6)];
%% create precision matrix
spde.C = speye(prod(sz));
spde.G = igmrfprec(sz,1);
spde.G2 = spde.G*spde.G;

%% Initial plots
%indicators of observed grid cells
I_obs = reshape(sum(A,1),sz);
I_grid = reshape(sum(A_grid,1),sz);

%plot covariates
figure(1)
for i=2:size(B_grid,2)
  subplot(2,3,i-1)
  imagesc(longitude, latitude([2 1]), reshape(A_grid'*B_grid(:,i),sz), ...
    'alphadata', I_grid)
  axis xy tight
  title(names{i})
  colorbar
  set(gca,'FontSize',30);
end

%plot observations and their covariates
figure(2)
subplot(3,3,1)
imagesc(longitude, latitude([2 1]), reshape(A'*Y,sz), 'alphadata', I_obs)
axis xy tight
title('Pollen data')
set(gca,'FontSize',25);
colorbar
for i=2:size(B,2)-1
  subplot(3,3,i)
  imagesc(longitude, latitude([2 1]), reshape(A'*B(:,i),sz), ...
    'alphadata', I_obs)
  axis xy tight
  title(names{i})
  colorbar
  set(gca,'FontSize',25);
end

%% Fit model

%we need a global variable for x_mode to reuse it
%between optimisation calls
global x_mode;
x_mode = [];
%attempt to estimate parameters (optim in optim...)
%subset to only observed points here
% [0 0] - two parameters to estimate, lets start with 0 in log-scale.

par = fminsearch( @(theta) GMRF_negloglike(theta, Y, A, B, spde), [0 0]);
%par2 = fminsearch( @(theta) GMRF_negloglike(theta, Y, A, B2, spde, alpha), [0 0]);
%conditional mean is now given be the mode
E_xy = x_mode;

%use the taylor expansion to compute posterior precision
%you need to reuse some of the code from GMRF_negloglike_NG 
%to create inputs for this function call
N = size(B,2);
tau = exp(par(1));
kappa2 = exp(par(2));
 % Q_x = exp(par(1))*exp((par(2))*spde.C + spde.G); <- Förra (fel på
 % parenteserna)

Q_x = tau*(kappa2*spde.C + spde.G);
qbeta=1e-6;
Qbeta = qbeta*speye(N);
Qall = blkdiag(Q_x,Qbeta);
Atilde = [A B];

[~, ~, Q_xy] = GMRF_taylor(E_xy, Y, Atilde, Qall);
%% Calculation of relevant Betas
e = [zeros(size(Q_xy,1)-N, N); eye(N)];
V_beta0 = e'*(Q_xy\e);
beta_CI95 =  sqrt(diag(V_beta0)) * 1.96;
E_beta = x_mode(end - N + 1:end);

for i = 1:N
 print = ['Beta ', num2str(i), ' : ', num2str(E_beta(i)), ' +- ', num2str(beta_CI95(i))];
 disp(print);
end

%CI of beta6 is further from zero -> B(1:6) is more appropriate

%CI_U = E_beta(2:end) + beta_i(2:end);
%CI_L = E_beta(2:end) - beta_i(2:end);

%Plotted using matlab r2014b, broken on newer versions.
%figure(6);
%candle(CI_U,CI_U,CI_L,CI_L 'b');
%xlim([0.5,5.5]);
%line([0 7],[0 0],'LineWidth',2);
%ylim([-2 2]);
%grid on;

