function cahnHilliardStrain(NT, NPRINT, NSKIPSTRAIN, Lx0, Ly0, Lz0, dx, dy, dz, phi0, seed, ftype)
%% FUNCTION TO SIMULATION CAHN-HILLIARD EQUATION 
% in 3D with variable seeds, parameters and initial conditions

% set random seed
rng(seed);

% system boundary info
NDIM    = 3;

% simulation duration info
dt      = 0.01;

% get max possible matrix size
Lxmax = 2*Lx0;
Lymax = 2*Ly0;
Lzmax = 2*Lz0;

% initialize Lx, Ly, Lz
Lx = Lx0;
Ly = Ly0;
Lz = Lz0;

% use fixed vectors
rx = 1:Lx;
ry = 1:Ly;
rz = 1:Lz;

% concentration field (pad to max)
phi             = zeros(Lymax,Lxmax,Lzmax);
phi0mat         = repmat(phi0,Lymax,Lxmax,Lzmax);
phi(ry,rx,rz)   = phi0mat(ry,rx,rz) + 0.1*randn(Ly,Lx,Lz);
psi             = phi.^3 - phi;

% FFT phi ONCE (is updated in Euler, so no need to constantly FFT back and
% forth)
fphi    = fftn(phi(ry,rx,rz));

% track time
t = 0.0;
frame = 1;

%% Wavevector matrices

% initial wavevectors
kx = karray(Lx);
ky = karray(Ly);
kz = karray(Lz);

% matrices
[KX, KY, KZ] = meshgrid(kx,ky,kz);
K2 = KX.^2 + KY.^2 + KZ.^2;
K4 = K2.^2;


%% Loop over time, update based on pseudospectral semi-implicit method (see notes)

% time loop
for tt = 1:NT
    % strain box
    if mod(tt,NSKIPSTRAIN) == 0
        % increment box length
        if Lx + dx <= Lxmax
            Lx = Lx + dx;
        end
        if Ly + dy <= Lymax
            Ly = Ly + dy;
        end
        if Lz + dz <= Lzmax
            Lz = Lz + dz;
        end
        
        % reset range
        rxnew = 1:Lx;
        rynew = 1:Ly;
        rznew = 1:Lz;
        
        % wavevectors
        kx = karray(Lx);
        ky = karray(Ly);
        kz = karray(Lz);

        % matrices
        [KX, KY, KZ] = meshgrid(kx,ky,kz);
        K2 = KX.^2 + KY.^2 + KZ.^2;
        K4 = K2.^2;
        
        % linearly interpolate
    end
    % Print to console and figure
    if mod(tt,NPRINT) == 0
        % print to console
        fprintf('   ** t = %d: Simulating CH eq in d = %d, with Lx=%d, Ly=%d, Lz=%d, phi0=%0.5g\n',tt,NDIM,Lx,Ly,Lz,phi0);
        
        % print configuration
        fstr = [ftype '_' num2str(frame) '.pos'];
        frame = frame + 1;
        plot2File(fstr,phi,t,dt,Lx,Ly,Lz);
    end
    
    % FFT phi and psi
    fpsi = fftn(psi(ry,rx,rz));
    
    % update phi based on semi-implicit scheme in Fourier space
    fphi = (fphi - dt.*K2.*fpsi)./(1 + K4.*dt);
    
    % IFFT back to update psi
    phi(ry,rx,rz) = ifftn(fphi);
    psi(ry,rx,rz) = phi(ry,rx,rz).^3 - phi(ry,rx,rz);
    
    % update time
    t = t + dt;
end

end



%% function to plot configuration to file
function plot2File(fstr,phi,t,dt,Lx,Ly,Lz)
% open file
fid = fopen(fstr,'w');

% print
fprintf(fid,'%d \n',Lx);
fprintf(fid,'%d \n',Ly);
fprintf(fid,'%d \n',Lz);
fprintf(fid,'%f \n',t);
fprintf(fid,'%f \n',dt);
for zz = 1:Lz
    for yy = 1:Ly
        for xx = 1:Lx
            fprintf(fid,'%f ',phi(yy,xx,zz));
        end
        fprintf(fid,'\n');
    end
end

% close
fclose(fid);
end


%% Function to compute wavevector arrays given a certain number of grid points
function k = karray(L)
k = zeros(1,L);
Fs = 2.0*pi/L;
k(1:(L/2)) = Fs*(0:(L/2-1));
k((L/2)+1:end) = Fs*((L/2:(L-1))-L);
end