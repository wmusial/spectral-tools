function x = cheb_grid(N, xmin, xmax)
% CHEB_DER Creates a Chebyshev grid.
%
%   X = CHEB_GRID(N, XMIN, XMAX)
%
%      Creates a Chebyshev grid [XMIN, XMAX) for later use
%      with CHEB_DER and CHEB_INT.
%
% @author  Wojciech Musial wmusial@mit.edu
% @created 10/19/2013
% 
  x = (xmin - xmax)./2.*(cos(pi*(2*(0:N-1))./(2*(N-1)))-1) + xmin;
end
