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

%Perform a PCA on the time-series data to find important components
[y,~,P_time] = pca(colstack(img));
y = reshape(y, sz(1), sz(2), []);


N4 = [0 1 0; 1 0 1; 0 1 0]; 
N8 = [1 1 1; 1 0 1; 1 1 1]; 




% Classifications with diffrent amount of classes
for i = 3:6
    findClass(y_beta, y, sz, i,N4, 0);
    findClass(y_beta, y, sz, i,N8, 0);
end
