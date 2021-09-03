function cahnHilliard(NT, NPRINT, Lx, Ly, Lz, phi0, seed, ftype)
%% FUNCTION TO SIMULATION CAHN-HILLIARD EQUATION 
% in 3D with variable seeds, parameters and initial conditions

% set random seed
rng(seed);

% system boundary info
NDIM    = 3;

% simulation duration info
dt      = 0.01;

% concentration field
phi     = phi0 + 0.1*randn(Ly,Lx,Lz);
psi     = phi.^3 - phi;

% FFT phi ONCE (is updated in Euler, so no need to constantly FFT back and
% forth)
fphi    = fftn(phi);

% track time
t = 0.0;
frame = 1;

%% Wavevector matrices

% kx
kx = zeros(1,Lx);
Fs = 2.0*pi/Lx;
kx(1:(Lx/2)) = Fs*(0:(Lx/2-1));
kx((Lx/2)+1:end) = Fs*((Lx/2:(Lx-1))-Lx);

% ky
ky = zeros(1,Ly);
Fs = 2.0*pi/Ly;
ky(1:(Ly/2)) = Fs*(0:(Ly/2-1));
ky((Ly/2)+1:end) = Fs*((Ly/2:(Ly-1))-Ly);

% kz
kz = zeros(1,Lz);
Fs = 2.0*pi/Lz;
kz(1:(Lz/2)) = Fs*(0:(Lz/2-1));
kz((Lz/2)+1:end) = Fs*((Lz/2:(Lz-1))-Lz);

% matrices
[KX, KY, KZ] = meshgrid(kx,ky,kz);
K2 = KX.^2 + KY.^2 + KZ.^2;
K4 = K2.^2;


%% Loop over time, update based on pseudospectral semi-implicit method (see notes)

% time loop
for tt = 1:NT
    % Print to console and figure
    if mod(tt,NPRINT) == 0
        % print to console
        fprintf('   ** t = %d: Simulating CH eq in d = %d, with Lx=%d, Ly=%d, Lz=%d, phi0 = %f\n',tt,NDIM,Lx,Ly,Lz,phi0);
        
        % print configuration
        fstr = [ftype '_' num2str(frame) '.pos'];
        frame = frame + 1;
        plot2File(fstr,phi,t,dt,Lx,Ly,Lz);
    end
    
    % FFT phi and psi
    fpsi = fftn(psi);
    
    % update phi based on semi-implicit scheme in Fourier space
    fphi = (fphi - dt.*K2.*fpsi)./(1 + K4.*dt);
    
    % IFFT back to update psi
    phi = ifftn(fphi);
    psi = phi.^3 - phi;
    
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