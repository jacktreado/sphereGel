function [kbins, sk, Sxy, Sxz, Syz, lambda, V] = fourierSpaceCorrelation(binaryLattice,L)
%% FUNCTION to compute fourier space correlations based on binary lattice

% lattice size (NOTE: rows = y, cols = x for images)
[ny, nx, nz] = size(binaryLattice);

% box lengths
Lx = L(1);
Ly = L(2);
Lz = L(3);


fprintf('\t\t -- Getting correlations: computing FFTs...\n');
F_microCT = fftn(binaryLattice);
F_microCT = fftshift(F_microCT);
S = (F_microCT.*conj(F_microCT))./(Lx*Ly*Lz);

% get wavevector values (diff for even and odd)

% kx
Fs = (2.0*pi)/Lx;
if mod(nx,2) == 0
    kx = Fs*(-nx/2:(nx/2 - 1));
else
    kx = Fs*(-floor(nx/2):(round(nx/2) - 1));
end
kx0 = find(kx == 0);

% ky
Fs = (2.0*pi)/Ly;
if mod(ny,2) == 0
    ky = Fs*(-ny/2:(ny/2 - 1));
else
    ky = Fs*(-floor(ny/2):(round(ny/2) - 1));
end
ky0 = find(ky == 0);

% kz
Fs = (2.0*pi)/Lz;
if mod(nz,2) == 0
    kz = Fs*(-nz/2:(nz/2 - 1));
else
    kz = Fs*(-floor(nz/2):(round(nz/2) - 1));
end
kz0 = find(kz == 0);

% subtract out contribution from k = 0
S(ky0,kx0,kz0) = 0;

% compute structure factor by radial average
fprintf('\t\t -- Computing isotropic structure factor...\n');

% loop over distances in each image, compute radial correlations
kmin        = 0;
kmax        = mean([max(kx) max(ky) max(kz)]);
dk          = 1.5*(2.0*pi)/mean([nx ny nz]);
binedges    = (kmin:dk:kmax)';
nkbins      = length(binedges) - 1;
leftBins    = binedges(1:end-1);
rightBins   = binedges(2:end);
kbins       = 0.5*(leftBins + rightBins);
sk          = zeros(nkbins,1);

% compute knorm lattice
fprintf('\t\t -- Computing kcube matrix...\n');
kyx                         = (ky').^2 + kx.^2;
kzpage                      = reshape(kz,1,1,nz);
knormsq                     = kyx + kzpage.^2;
knorm                       = sqrt(knormsq);

% loop over bins
fprintf('\t\t -- Looping over bins and calculating S(k)...\n');
for bb = 1:nkbins
    % get location of k values in 3D matrix
    ktmpmin = leftBins(bb);
    ktmpmax = rightBins(bb);
    kinds = knorm(:) < ktmpmax & knorm(:) > ktmpmin;

    % take standard mean
    sk(bb) = mean(S(kinds),'all');
end

% Gyration tensor (see ctnm limit def on wikipedia)
G = zeros(3);
Sint = 0.0;

% loop over all points in S
NP = nx*ny*nz;
nxy = nx*ny;
fprintf('\t\t -- Getting gyration tensor G...\n');
for pp = 1:NP
    % get indices based on pp
    pm1 = pp-1;
    rowi = mod(pm1,ny)+1;
    pgei = floor(pm1/nxy)+1;
    coli = mod(floor(pm1/ny),nx)+1;
    
    % get wavevector
    ktmp = [kx(coli), ky(rowi), kz(pgei)];

    % add to gyration tensor
    for ii = 1:3
        for jj = ii:3
            G(ii,jj) = G(ii,jj) + S(pp)*ktmp(ii)*ktmp(jj);
        end
    end
    
    % add to integral over all S
    Sint = Sint + S(pp);
end

% symmetrize
G(2,1) = G(1,2);
G(3,2) = G(2,3);
G(3,1) = G(1,3);

% normalize
G = G./Sint;

fprintf('\t\t -- Getting eigenvalues of G...\n');

% get principal axes from gyration tensor
[V,D] = eig(G);
lambda = diag(D)';

% get 2D correlators
fprintf('\t\t -- Smoothing S, getting 2D correlators...\n');
Ssmooth = smooth3(S,'box',3);

% yz: average over columns (X)
Syz = mean(Ssmooth,2);
Syz = reshape(Syz,ny,nz);

% yz: average over rows (Y)
Sxz = mean(Ssmooth,1);
Sxz = reshape(Sxz,nx,nz);

% yz: average over pages (Z)
Sxy = mean(Ssmooth,3);

end