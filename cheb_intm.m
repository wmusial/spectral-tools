function m = cheb_intm(N, xmin, xmax)
% CHEB_INTM Computes integration matrices for use with CHEB_INT
%
%   M = CHEB_INT(N, XMIN, XMAX)
%    
%      N, XMIN, XMAX should match values used in CHEB_GRID.
%      Returns an integration matrix for use with CHEB_INT.
%      
%  
%      Adapted from:
%           Trefethen "Spectral Methods in Matlab"  Chapter 6
%
% @author  Wojciech Musial wmusial@mit.edu
% @created 10/19/2013
% 

  if N < 2
    error('can only integrate non-singleton dimensions')
  end

  % compute differentiation matrices
  M = N-1;
  x = cos(pi*(0:M)./M)';
  c = [2; ones(M-1,1); 2].*(-1).^(0:M)';
  X = repmat(x, 1, M+1);
  dX = X-X';
  D = (c*(1./c)')./(dX+(eye(M+1)));
  D = D - diag(sum(D'));

  % compute the integration matrices
  m = zeros(N);

  i = 2:N;
  m(i,i) = inv(D(i,i));

  % fix normalization
  normalization = (xmin-xmax)./2;
  m = m*normalization;
end
