%%  Data Init. Always run. 
load HA1_temp_2015_06
%points inside of Sweden (i.e. not nan)
I_land = ~any(isnan(X_grid),2);
%reshape grid to images
X_img = reshape(X_grid, sz(1), sz(2), []);
I_img = reshape(I_land, sz);
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250];[0.4660, 0.6740, 0.1880];[0.6350, 0.0780, 0.1840];[0.4940, 0.1840, 0.5560]];

Longitude = reshape(X_grid(:,1),sz);
Latitude = reshape(X_grid(:,2),sz);
Elevation = reshape(X_grid(:,3),sz);
Coast = reshape(X_grid(:,4),sz);
Coast_SE = reshape(X_grid(:,5),sz);
Coast_NO = reshape(X_grid(:,6),sz);

%% Data Visualization

% Plots covariates
figure(1);
subplot(3,2,1),imagesc(Longitude),colorbar,title('Longitude');
subplot(3,2,2),imagesc(Latitude),colorbar,title('Latitude');
subplot(3,2,3),imagesc(Elevation),colorbar,title('Elevation');
subplot(3,2,4),imagesc(Coast),colorbar,title('Distance to coast');
subplot(3,2,5),imagesc(Coast_SE),colorbar,title('Distance to Swedish coast');
subplot(3,2,6);imagesc(Coast_NO),colorbar,title('Distance to Norwegian coast');


%plotting diffrent correlations
figure(2), subplot(3,2,1);

title('Kovarians');
scatter(X(:,1),Y,20,colors(1,:), 'filled');
grid on;
xlabel('Longitude');
subplot(3,2,2);
scatter(X(:,2),Y,20,colors(2,:),'filled');
grid on;
xlabel('Latitude');
subplot(3,2,3);
scatter(X(:,3),Y,20,colors(3,:),'filled');
grid on;
xlabel('Elevation');
subplot(3,2,4);
scatter(X(:,4),Y,20,colors(4,:),'filled');
grid on;
xlabel('Distance to coast');
subplot(3,2,5);
scatter(X(:,5),Y,20,colors(5,:),'filled');
grid on;
xlabel('Distance to Swedish coast');
subplot(3,2,6);
scatter(X(:,6),Y,20,colors(6,:),'filled');
grid on;
xlabel('Distance to Norwegian coast');

%% OLS 
OLS_Models;
