%% Load data from file %%
calib_data = readmatrix('calibrationmodel.csv'); %turns my inputted data into a 2 x (number of positions) matrix
siize = size(calib_data);
N = max(siize); %the number of positions

%% Calculate the mean and standard deviation of the data %%
az = calib_data(:,1);
el = calib_data(:,2);

mu_az = mean(az);
mu_el = mean(el);
std_az = std(az);
std_el = std(el);


%% Calculate the uncertainty of the values using the normal distribution %%
coeffs_az = lsq_fit(mu_az, az);
coeffs_el = lsq_fit(mu_el, el);

uncertainty_az = std_az/sqrt(N); 
uncertainty_el = std_el/sqrt(N);

%% Plot residuals %%

residual_az = az - mu_az;
plot(residual_az);

residual_el = el - mu_el;
plot(residual_el);

%% Calculate the coefficients of the normal distribution using least squares %%
% The coefficients of the normal distribution are calculated using least squares,
% where x is the deviation of each data point from the mean, 
% A is the matrix for the linear least squares problem,
% and coeffs are the least squares coefficients.

function coeffs = lsq_fit(mu, data)

x = data - mu; %residual
A = [x, ones(size(x))];

end