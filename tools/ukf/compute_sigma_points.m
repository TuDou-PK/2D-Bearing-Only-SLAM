#this function extract the sigma points
# inputs: 
#  mu: is the mean of the robot pose (its size determine the number of sigma points)
#  sigma: is the mean of the previously estimated robot pose (3x3 matrix)
#
# outputs 
#   sigmaP: is a matrix 3x(2n+1) containing the sigma points (one per column)
#   weightsM and weightsC: the weights of the sigma points

function [sigmaP, weightsM, weightsC] = compute_sigma_points(mu, sigma)

	state_dim = size(mu,1); 
	num_of_sigma_points = 2*state_dim+1;

	% From theory
	k = 0;
	alpha = 1e-3;
	lambda = alpha^2 * ( state_dim + k );% - d;
	bet = 2;

	%init the sigma points matrix and its weight vectors
	sigmaP = zeros(size(mu,1),num_of_sigma_points);
	weightsM = zeros(num_of_sigma_points,1);
	weightsC = zeros(num_of_sigma_points,1);

	%% STEP 1: compute weights
	% the first weights can be computed as
	weightsM(1) = lambda/(state_dim+lambda);
	weightsC(1) = weightsM(1) + (1 + alpha^2 + bet);

	
	%% STEP 2: compute sigma points
	% the first sigma point is the mean
	sigmaP(:,1) = mu;


	%compute the other weights
	weight_value = .5/(state_dim + lambda);
	weightsM(2:end) = repmat(weight_value, num_of_sigma_points-1,1);
	weightsC(2:end) = weightsM(2:end);

	% chol extract the Cholesky Decomposition
	% with the arg "lower" it returns the lower
	% triangular matrix L, such that L'*L = A
	[L] = chol(sigma*(state_dim+lambda) , "lower");

	point_idx = 2;
	for i=1:state_dim
		L_i = L(:,i);
		%half of the remaining 2*state_dim sigma points
		sigmaP(:,point_idx) = mu + L_i;
		point_idx++;
		sigmaP(:,point_idx) = mu - L_i;
		point_idx++;
		
	end

endfunction
