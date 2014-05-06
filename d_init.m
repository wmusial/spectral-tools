function d_init(N, D)
% D_INIT precomputes spectral derivative matrices.
%
%     N number of grid points
%     D physical distance between adjacent grid cells

  tic;

  fprintf('computing spectral derivative matrices...  ');

  % TODO can upgrade on this
  % implement fft-based differentiation -- trefethen pg 6


  % the derivative code assumes Nx = Ny for now. 
  % Looks like that might in fact stay this way forever.
  Nx = N;
  Ny = N;
  Dx = D;
  Dy = D;

  global d1x d1y;
  %global d2x d2y;

  if size(d1x) ~= [Nx, Nx] | size(d1y) ~= [Ny, Ny]
    d1x = zeros(Nx, Nx);
    d1y = zeros(Ny, Ny);
    d2x = zeros(Nx, Nx);
    d2y = zeros(Ny, Ny);

    % x derivatives
    for x = 1:Nx
      for y=1:Nx
        s1 = 0;
        s2 = 0;
        for k = -Nx/2 : Nx/2
          s1 = s1 + 2. * pi * i / Nx / Nx * k * exp(2 * pi * i * k * (x - y) / Nx);
          %s2 = s2 + (2. * pi * i / Nx * k)^2 / Nx * exp(2 * pi * i * k * (x - y) / Nx);
        end
        d1x(x,y) = real(s1)/Dx;
        %d2x(x,y) = real(s2)/Dx/Dx;
      end
    end

    % y derivatives
    for x = 1:Ny
      for y=1:Ny
        s1 = 0;
        s2 = 0;
        for k = -Ny/2 : Ny/2
          s1 = s1 + 2. * pi * i / Ny / Ny * k * exp(2 * pi * i * k * (x - y) / Ny);
          %s2 = s2 + (2. * pi * i / Ny * k)^2 / Ny * exp(2 * pi * i * k * (x - y) / Ny);
        end
        d1y(x,y) = real(s1)/Dy;
        %d2y(x,y) = real(s2)/Dy/Dy;
      end
    end

  end

  tder = toc;
  fprintf('done in %.2f sec\n', toc);

