function [logp, dlogp, d2logp]= GMRF_taylor_Be(x_0, y, A, Q)
% GMRF_TAYLOR_BE  Taylor expansion of the conditional for non-Gaussian observations
%
% [logp, dlogp, d2logp] = GMRF_taylor_Be(x_0, y, A, Q)
%
% x_0 = value at which to compute taylor expansion
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% Q = the precision matrix, sparse N-by-N
%
% Function should return taylor expansion, gradient and Hessian including 
% both the data likelihood and the prior on x (i.e. what is used to finding 
% the mode).
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_negloglike_skeleton.m 4454 2011-10-02 18:29:12Z johanl $
 
%Observation model
%Y_i ~ Be(p_i)  with p_i = exp(z_i)/(exp(z_i)+1)

%compute log observations, and derivatives
z = A*x_0;

binom = 1 / (factorial(y*(factorial(1-y))));
f = log(binom) + z*y - log(e^z + 1); 

%compute the function, log p(x|y,theta)
logp = (x_0'*Q*x_0)/2 - sum(f);

if nargout>1
  %compute derivatives (if needed, i.e. nargout>1)
  d_f = y - e^z / (e^z + 1); % first derivative wrt z_i
  dlogp = Q*x_0 - A'*d_f;
end

if nargout>2
  %compute hessian (if needed, i.e. nargout>2)
  d2_f = -e^z / (e^z + 1)^2; % second derivative wrt z_i
  n = size(A,1);
  d2logp = Q - A'*spdiags(d2_f,0,n,n)*A;
end