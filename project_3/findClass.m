function c = findClass(y_beta, y_time, sz, nbr_class, N, drawflag)
%% KMEANS on regression
xc_beta = colstack(y_beta(:,:,1:end));
[cl_beta,~] = kmeans(xc_beta, nbr_class, inf, 0);
figure(1);
sgtitle('KMEANS classification on regression parameters');
subplot(2,2,(nbr_class-2)); 
imagesc(reshape(cl_beta, 87, 102));
title(string((nbr_class)) + ' classes ');
%% KMEANS on timedata
xc_time = colstack(y_time); 
[cl_time,~] = kmeans(xc_time, nbr_class, inf, 0);
figure(2);
sgtitle('KMEANS classification on timeseries');
subplot(2,2,(nbr_class-2)); 
imagesc(reshape(cl_time, 87, 102))
title(string((nbr_class)) + ' classes ');
%% GMM on Regression model
[theta_beta,prior]=normmix_gibbs(xc_beta,nbr_class);
[cl,cl_ind,~]=normmix_classify(xc_beta,theta_beta,prior);
figure(3);
sgtitle('Gaussian Mixture Model classification');
subplot(2,2,(nbr_class-2)   ); 
imagesc(reshape(cl, sz(1:2)));
title(string((nbr_class)) + ' classes ');
%% creating alpha posteriors from GMM
alpha_posterior = log(prior(2:nbr_class)/prior(1));
cl_img = reshape(cl_ind,[sz(1:2),nbr_class ]);

%% Simulating MRF

y_m = y_beta(:,:,1:end);
y_c = colstack(y_m);
iter = 1000;
Burnin = 350;
alpha_posterior = ones(iter+1,nbr_class-1)/nbr_class;
%Same beta for each class
%beta = zeros(iter+1,1) ; 
% Diffrent beta for each class
beta = zeros(iter+1,nbr_class);
% Suggested prior
beta_prior = 1/10*eye(nbr_class);

%Plog = zeros(sz(1),sz(2),nbr_class);
Plog = zeros(length(iter),1);
z_z = cl_img;
z_sum = zeros(size(z_z));
for i=1:iter
    alpha_beta_post=mrf_gaussian_post([0 alpha_posterior(i,:)],theta_beta,y_m);
    z_z= mrf_sim(z_z, N, alpha_beta_post, beta(i,:), 1);
    [~,Mz_b]=mrf_sim(z_z,N,alpha_beta_post, beta(i,:) ,0);
    Plog(i) = sum(log(Mz_b(logical(z_z))));
    
    [alpha_posterior(i+1,:), beta(i +1,:), ~] = gibbs_alpha_beta(alpha_posterior(i,:), beta(i,:), z_z, N, beta_prior);
for k = 1:nbr_class
        ind_beta = find(z_z(:,:,k));
        [theta_beta{k}.mu, theta_beta{k}.Sigma] = gibbs_mu_sigma(y_c(ind_beta,:));
       
end
 if drawflag == 1
 figure,
 image(rgbimage(z_z))
 title('MRK-sim, Beta ')
 drawnow    
 end
     if i > Burnin   
         z_sum = z_sum + z_z;
     end
end

z_mean = z_sum./(iter-Burnin);
% Choosing plot number based on N
if (sum(find(N,2)) == 6)
    plt = [4, 5];
    tit = string('4 Neighbours');
elseif(sum(find(N,2)) == 3)
    plt = [6,7];
    tit = string('8 Neighbours');
end

figure(plt(1)),
subplot(2,2,(nbr_class-2)); 
plot(Plog);
sgtitle('Pseduo-likelihood of diffrent number of classes - '+ tit );
title(string((nbr_class)) + ' classes ');
figure(plt(2)),
[~,z_plot] = max(z_mean,[],3);
subplot(2,2,(nbr_class-2)); 
imagesc(z_plot);
sgtitle('GMRF classification - '+tit);
title(string((nbr_class)) + ' classes ');

end
