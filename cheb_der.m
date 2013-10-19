function [dg ddg] = cheb_der(g, xmin, xmax, varargin)
% CHEB_DER Numerically approximates first and second derivatives using 
% the fft-based Chebyshev differentiation scheme.
%
%   [DG, DDG] = CHEB_DER(G, XMIN, XMAX, DIM=1)
%          DG = CHEB_DER(G, XMIN, XMAX, DIM=1)
%
%      G specifies data to differentiate along dimension DIM.
%      The differentiation dimension DIM of G must conform to 
%      CHEB_GRID with XMIN and XMAX.
% 
%      Second-derivative computations are carried out only if 
%      CHEB_DER is called with both output arguments.
%  
%      Adapted from:
%           Trefethen "Spectral Methods in Matlab"  Chapter 8
%
% @author  Wojciech Musial wmusial@mit.edu
% @created 10/19/2013
% 
  if nargout > 2
    error('can only compute derivatives of order 1 and 2');
  end

  ord2 = (nargout == 2);

  if length(varargin) == 0
    dim = 1;
    if numel(g) == max(size(g))
      g = g(:);
    end
  else
    dim = varargin{1};
  end

  if size(g, dim) == 1
    error('attempting to differentiate along a singleton dimension!');
  end

  % shift the differentiation dimension to dim 1
  g = shiftdim(g,dim-1); 

  % store the original size of g
  gs = size(g);

  % reshape g for convenience
  g = reshape(g, gs(1), prod(gs)./gs(1));

  % first dimension
  M = size(g, 1);
  N = M-1;

  % replicate the first dimension of g into V
  V = zeros(2.*N, size(g, 2));
  V(1:M,:) = g;
  V(M+1:end,:) = g(M-1:-1:2,:);

  % compute the fourier transform 
  V = fft(V, [], 1);

  % multiply the fourier transform by k
  k = [0:N-1 N 1-N:-1]';
  k = repmat(k, [1 size(g,2)]);
  dV  = (1i .* k) .* V;

  % take the inverse fft
  dV = ifft(dV, [], 1);

  if ord2
    % multiply the fourier transform by k
    ddV = (-k.^2) .* V;
    % take the inverse fft
    ddV = ifft(ddV, [], 1);
    % truncate the replicated first dimension
    ddV = ddV(1:M, :);
  end

  % truncate the replicated first dimension
  V  =  V(1:M, :);
  dV = dV(1:M, :);
  k  =  k(N:end, :);

  
  % compute the algebraic derivative
  x = cos(pi .* [0:N]' ./ N);
  x = repmat(x, [1 size(g,2)]);

  dg = -dV ./ sqrt(1-x.^2);

  % fill the first and last point.
  k = [0:N]';
  k = repmat(k, [1 size(g,2)]);
  k2 = k.^2;

  V(1,:)   = V(1,:) ./ 2;
  V(end,:) = V(end,:) ./ 2;

  dg(1,:)   = sum(V .* k2, 1) ./ N;
  dg(end,:) = sum(-V .* (-1).^k .* k2, 1) ./ N;

  dg = reshape(dg, gs);
  dg = shiftdim(dg, length(gs) + 1 - dim) .* 2 ./ (xmin - xmax);

  if isreal(g)
    dg = real(dg);
  end

  if ord2
    ddg = -dV .* x ./ (1-x.^2).^(3/2) + ddV ./ (1 - x.^2);

    % fill the first and last point.
    k = [0:N]';
    k = repmat(k, [1 size(g,2)]);

    V(1,:)   = V(1,:) ./ 2;
    V(end,:) = V(end,:) ./ 2;

    ddg(1,:)   = sum(V .* k2 .* (k2-1) ./ 3, 1) ./ N;
    ddg(end,:) = sum(V.* (-1).^k .* k2 .* (k2-1) ./ 3, 1) ./ N;

    ddg = reshape(ddg, gs);
    ddg = shiftdim(ddg, length(gs) + 1 - dim) .* (2 ./ (xmin - xmax)).^2;

    if isreal(g)
      ddg = real(ddg);
    end
  end
end
