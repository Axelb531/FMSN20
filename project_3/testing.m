load fmri.mat

%size of data¨
sz = size(img);
nbr_class = 4;
%Option 1: regress onto indicator functions
beta = X\colstack(img)';
%reshape back to an image where the beta-coefs are in each "color"-layer
beta = reshape(beta', sz(1), sz(2), []);
%and treat the beta:s as the image
%Perform a PCA on the regression coefficients to find important components
[y_beta, ~, P_beta] = pca(colstack(beta(:,:,3:end))); % PCA med alla regressionskoefficienter?
y_beta = reshape(y_beta, sz(1), sz(2), []);
%% KMEANS on regression
xc_beta = colstack(y_beta(:,:,1:end));
[cl_beta,theta_beta] = kmeans(xc_beta, nbr_class, inf, 0);
figure, 
imagesc(reshape(cl_beta, 87, 102))

%% KMEANS on timedata
[y,~,P_time] = pca(colstack(img));
y = reshape(y, sz(1), sz(2), []);
xc_time = colstack(y); 
[cl_time,theta_time] = kmeans(xc_time, nbr_class, inf, 0);
figure, 
imagesc(reshape(cl_time, 87, 102))
figure, 
imagesc(beta(:,:,1));
%% GMM on Regression model
[theta_beta,prior]=normmix_gibbs(xc_beta,nbr_class);
[cl,cl_ind,~]=normmix_classify(xc_beta,theta_beta,prior);
figure,
imagesc(reshape(cl, sz(1:2)));
%% creating alpha posteriors from GMM
alpha_posterior = log(prior(2:nbr_class)/prior(1));
cl_img = reshape(cl_ind,[sz(1:2),nbr_class ]);

%% Simulating MRF
N1 = [0 1 0; 1 0 1; 0 1 0]; 
N2 = [1 1 1; 1 0 1; 1 1 1]; 
N3 = [1 0 1; 0 0 0; 1 0 1];

y_m = y_beta(:,:,1:end);
y_c = colstack(y_m);
iter = 1000;
Burnin = 350;
beta = 0; 
% Suggested prior
beta_prior = 1/10;

%Plog = zeros(sz(1),sz(2),nbr_class);
Plog = zeros(length(iter),1);
z_z = cl_img;
z_sum = zeros(size(z_z));
for i=1:iter
    alpha_beta_post=mrf_gaussian_post([0 alpha_posterior],theta_beta,y_m);
    z_z= mrf_sim(z_z, N2, alpha_beta_post, beta, 1);
    [~,Mz_b]=mrf_sim(z_z,N2,alpha_beta_post,beta,0);
    Plog(i) = sum(log(Mz_b(logical(z_z))));
    
    [alpha_posterior, beta, acc_beta] = gibbs_alpha_beta(alpha_posterior, beta, z_z, N2, beta_prior);
for k = 1:nbr_class
        ind_beta = find(z_z(:,:,k));
        [theta_beta{k}.mu, theta_beta{k}.Sigma] = gibbs_mu_sigma(y_c(ind_beta,:));
       
end
 
 figure(6)
 image(rgbimage(z_z))
 title('MRK-sim, Beta ')
 drawnow    
     if i > Burnin   
         z_sum = z_sum + z_z;
     end
end

z_mean = z_sum./(iter-Burnin);

%%
figure(7) 
plot(Plog);
title('Pseduo-likelihood regression')
figure(8)
[~,z_plot] = max(z_mean,[],3);
imagesc(z_plot);
