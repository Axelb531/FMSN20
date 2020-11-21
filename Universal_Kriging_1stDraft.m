%% MODEL 1: Latitude, Elevation, Coast, CoastSE, CoastNO
% Creating dataframe for model
X1 = [ ones(n_x,1), Latitude , Elevation, Coast, Coast_SE, Coast_NO];
Y_k = Y; % Observed locations

