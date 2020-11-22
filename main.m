%%  Data Init. Always run. 
load HA1_temp_2015_06
%points inside of Sweden (i.e. not nan)
I_land = ~any(isnan(X_grid),2);
%reshape grid to images
X_img = reshape(X_grid, sz(1), sz(2), []);
I_img = reshape(I_land, sz);
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250];[0.4660, 0.6740, 0.1880];[0.6350, 0.0780, 0.1840];[0.4940, 0.1840, 0.5560]];
n_x = size(X,1); 
% Data
Longitude = X(:,1);
Latitude = X(:,2);
Elevation = X(:,3);
Coast = X(:,4);
Coast_SE = X(:,5);
Coast_NO = X(:,6);
% IMG-data
Longitude_img = reshape(X_grid(:,1),sz);
Latitude_img = reshape(X_grid(:,2),sz);
Elevation_img = reshape(X_grid(:,3),sz);
Coast_img = reshape(X_grid(:,4),sz);
Coast_SE_img = reshape(X_grid(:,5),sz);
Coast_NO_img = reshape(X_grid(:,6),sz);

%% Data Visualization

% Plots covariates
figure(1);
for i= 1:6
subplot(2,3,i)
  imagesc([4.9 30.9], [71.1 55.4], X_img(:,:,i), 'alphadata', I_img)
  axis xy; 
  hold on;
  plot(Border(:,1),Border(:,2))
  scatter(X(:,1), X(:,2), 25, X(:,i),...
    'filled','markeredgecolor','k')
  colorbar
  hold off
  title(names{i})
end

% Plotting diffrent correlations
figure(2), 
for i =1:6 
subplot(3,2,i);
scatter(X(:,i),Y,20,colors(i,:), 'filled');
grid on;
xlabel(names(i));
end


%% OLS 
OLS_Models;

%% Kriging
Universal_Kriging_1stDraft;
