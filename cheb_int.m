function w = cheb_int(g, cheb_m)
% CHEB_INT Numerically integrates data using the fft-based Chebyshev 
% integration scheme.
%
%   IG = CHEB_INT(G, CHEB_M)
%
%      G specifies a vector of data to be integrated.
%      The vector dimension of  G must conform to CHEB_GRID.
%      CHEB_M is a forward/backward integration matrix obtained 
%      from cheb_intm.
% 
%      Adapted from:
%           Trefethen "Spectral Methods in Matlab"  Chapter 6
%
% @author  Wojciech Musial wmusial@mit.edu
% @created 10/19/2013
%
  g = g(:);
  w = cheb_m*g;

end
