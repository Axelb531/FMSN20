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
[y_beta, ~, P_beta] = pca(colstack(beta));
y_beta = reshape(y_beta, sz(1), sz(2), []);
%% PCA on regression
xc = colstack(y_beta(:,:,3:11));

[cl2,theta] = kmeans(xc, nbr_class, inf, 0);
figure, 
imagesc(reshape(cl2, 87, 102))

%% PCA on timedata
[y,V,P] = pca(colstack(img));
y = reshape(y, sz(1), sz(2), []);
xc_2 = colstack(y); 
[cl_2,theta] = kmeans(xc_2, nbr_class, inf, 0);
figure, 
imagesc(reshape(cl_2, 87, 102))

%% GMM on Regression model

[theta,prior]=normmix_gibbs(xc,nbr_class,[200, 100],0);
[cl,cl_ind,p]=normmix_classify(xc,theta,prior);
figure,
y=reshape(p,[87, 102 ,nbr_class]);
image(rgbimage(y))




%% creating alpha posteriors from GMM

alpha_posterior = log(prior(2:nbr_class)/prior(1));
cl_img= reshape(cl_ind,[sz(1:2),nbr_class ]);


%% Simulationg MRF
N1 = [0 1 0; 1 0 1; 0 1 0];
N2 = [1 1 1; 1 0 1; 1 1 1];
N3 = [1 0 1; 0 0 0; 1 0 1];

%%

y_m = y_beta(:,:,3:11);
y_c = colstack(y_m);
iter = 1000;
Burnin = 250;
beta = 0; 
beta_prior = 1/10;

Plog = zeros(sz(1),sz(2),nbr_class);


z = cl_img;


z_sum = zeros(size(z));
for i=1:iter
    alpha_beta_post=mrf_gaussian_post([0 alpha_posterior],theta,y_m);
    z= mrf_sim(z, N1, alpha_beta_post, beta, 1);
    [~,Mz_b]=mrf_sim(z,N1,alpha_beta_post,beta,0);
    Plog(i) = sum(log(Mz_b(logical(z))));
    
    [alpha_posterior, beta, acc_beta] = gibbs_alpha_beta(alpha_posterior, beta, z, N1, beta_prior);
for k = 1:nbr_class
        ind_beta = find(z(:,:,k));
        [theta_beta{k}.mu, theta_beta{k}.Sigma] = gibbs_mu_sigma(y_c(ind_beta,:));
       
end
 
 figure(6)
 image(rgbimage(z_sum))
 title('MRK-sim, Beta ')
 drawnow
     
     
     if i > Burnin   
         z_sum = z_sum + z;
     end
       
end
figure(7) 
 
plot(Plog); 
title('Pseduo likelyhood regression')
