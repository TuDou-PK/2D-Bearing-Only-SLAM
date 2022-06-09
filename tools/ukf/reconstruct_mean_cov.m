#this function reconstruct mean and covariance given sigma points and weights
# inputs: 
#   sigmaP: is a matrix 3x(2n+1) containing the sigma points (one per column)
#   weightsM and weightsC: the weights of the sigma points
#
# outputs: 
#  mu: is the mean of the robot pose (its size determine the number of sigma points)
#  sigma: is the mean of the previously estimated robot pose (3x3 matrix)

function [mu, sigma] = reconstruct_mean_cov(sigmaP, wM, wC)

	state_dim = size(sigmaP,1);
	num_of_sigma_points = size(sigmaP,2);

	%initialize mean
	mu = zeros(state_dim,1);
	%populate mean
	for i=1:num_of_sigma_points
		mu += wM(i)*sigmaP(:,i);
	endfor
	
	%initialize covariance
	sigma = zeros(state_dim,state_dim);
	%populate covariance
	for i=1:num_of_sigma_points
		delta = (sigmaP(:,i) - mu);
		sigma += (delta)*(delta)' * wC(i);
	endfor
	

end
