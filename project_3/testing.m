load fmri.mat

%size of data
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

figure,  
[cl2,theta] = kmeans(xc, nbr_class, inf, 2);
figure, 
imagesc(reshape(cl2, 87, 102))

%% PCA on img
[y,V,P] = pca(colstack(img));
y = reshape(y, sz(1), sz(2), []);
xc_2 = colstack(y);
figure,  
[cl_2,theta] = kmeans(xc_2, nbr_class, inf, 2);
figure, 
imagesc(reshape(cl_2, 87, 102))

%% 

figure,
[theta,prior]=normmix_gibbs(xc,nbr_class,[200, 100],2)
%%
figure,
[cl,cl_ind,p]=normmix_classify(xc,theta,prior)
figure,
y=reshape(p,[87, 102 ,nbr_class]);
image(rgbimage(y))

%%
N1 = [0 1 0; 1 0 1; 0 1 0];
N2 = [1 1 1; 1 0 1; 1 1 1];
N3 = [0 0 0 1 1; 0 0 1 1 1; 0 1 0 1 0; 1 1 1 0 0; 1 1 0 0 0];
%alpha_post = mrf_gaussian_post(alpha,theta,y)
  
x = zeros(87,102,3); % M-by-N-by-P array of zeros
iter = 100;
Plog = zeros(iter,1);
for i = 1:iter
x=mrf_sim(x,N3,0,1,1); 
[~,Mz,Mf]=mrf_sim(x,N3,0,0.1,0);
Plog(i) = sum(log(Mz(logical(x))));
imagesc(x)
drawnow
end
