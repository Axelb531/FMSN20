%% load data
load HA1_temp_2015_06
%points inside of Sweden (i.e. not nan)
I_land = ~any(isnan(X_grid),2);
%reshape grid to images
X_img = reshape(X_grid, sz(1), sz(2), []);
I_img = reshape(I_land, sz);

%% Data Visualization

longitude = reshape(X_grid(:,1),sz);
latitude = reshape(X_grid(:,2),sz);
elevation = reshape(X_grid(:,3),sz);
coast = reshape(X_grid(:,4),sz);
coast_SE = reshape(X_grid(:,5),sz);
coast_NO = reshape(X_grid(:,6),sz);

figure(1);
%Test
title('Data');
subplot(3,2,1);
imagesc(longitude);
colorbar;
subplot(3,2,2);
imagesc(latitude);
colorbar;
subplot(3,2,3);
imagesc(elevation);
colorbar;
subplot(3,2,4);
imagesc(coast);
colorbar;
subplot(3,2,5);
imagesc(coast_SE);
colorbar;
subplot(3,2,6);
imagesc(coast_NO);
colorbar;


%plotting diffrent correlations
figure(2), subplot(3,2,1);
title('Kovarians');
scatter(X(:,1),Y,20,'filled');
xlabel('longitude');
subplot(3,2,2);
scatter(X(:,2),Y,20,'filled');
xlabel('latitude');
subplot(3,2,3);
scatter(X(:,3),Y,20,'filled');
xlabel('elevation');
subplot(3,2,4);
scatter(X(:,4),Y,20,'filled');
xlabel('coast');
subplot(3,2,5);
scatter(X(:,5),Y,20,'filled');
xlabel('coast_SE');
subplot(3,2,6);
scatter(X(:,6),Y,20,'filled');
xlabel('coast_NO');

%% BHAT
Y_nb = Y-(X(:,4).*Y);
cd = X(:,4);
figure(3), scatter(X_coast,Y_nb, 20, 'filled');
bhat = (cd.'*cd) \ (cd.'*Y_nb);
