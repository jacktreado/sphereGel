function cahnHilliardUniExt(NT, NPRINT, NSKIPSTRAIN, Lx, Ly, Lz0, phi0, seed, ftype)
%% FUNCTION TO SIMULATION CAHN-HILLIARD EQUATION 
% in 3D with variable seeds, parameters and initial conditions
% -- Need to figure out odd integers

% set random seed
rng(seed);

% system boundary info
NDIM    = 3;

% simulation duration info
dt      = 0.01;

% get max possible matrix size
Lzmax = 2*Lz0;

% initialize Lx, Ly, Lz
Lz = Lz0;
rz = 1:Lz;

% concentration field (pad to max)
phi             = zeros(Ly,Lx,Lzmax);
phi0mat         = repmat(phi0,Ly,Lx,Lzmax);
phi(:,:,rz)     = phi0mat(:,:,rz) + 0.1*randn(Ly,Lx,Lz);
psi             = phi.^3 - phi;

% FFT phi ONCE (is updated in Euler, so no need to constantly FFT back and
% forth)
fphi    = fftn(phi(:,:,rz));

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
    if mod(tt,NSKIPSTRAIN) == 0 && Lz + 1 <= Lzmax
        % linearly interpolate
        sc = (Lz - (1:Lz-1))./Lz;
        fwd = 2:Lz;
        bwd = 1:Lz-1;
        sc = reshape(sc,1,1,Lz-1);
        fwd = reshape(fwd,1,1,Lz-1);
        bwd = reshape(bwd,1,1,Lz-1);
        for xx = 1:Lx
            for yy = 1:Ly
                phi(yy,xx,Lz+1) = phi(yy,xx,Lz);
                phi(yy,xx,fwd) = sc.*(phi(yy,xx,fwd) - phi(yy,xx,bwd)) + phi(yy,xx,bwd);
            end
        end
        Lz = Lz + 1;
        
        % reset range
        rz = 1:Lz;
        
        % wavevectors
        kz = karray(Lz);

        % matrices
        [KX, KY, KZ] = meshgrid(kx,ky,kz);
        K2 = KX.^2 + KY.^2 + KZ.^2;
        K4 = K2.^2;
        
        % update fourier version of phi
        fphi = fftn(phi(:,:,rz));
    end
    % Print to console and figure
    if mod(tt,NPRINT) == 0
        % print to console
        fprintf('   ** t = %d: Simulating CH eq in d = %d, with Lx=%d, Ly=%d, Lz=%d, phi0=%0.5g\n',tt,NDIM,Lx,Ly,Lz,phi0);
        
        % print configuration
        fstr = [ftype '_' num2str(frame) '.pos'];
        frame = frame + 1;
        plot2File(fstr,phi(:,:,rz),t,dt,Lx,Ly,Lz);
    end
    
    % FFT phi and psi
    fpsi = fftn(psi(:,:,rz));
    
    % update phi based on semi-implicit scheme in Fourier space
    fphi = (fphi - dt.*K2.*fpsi)./(1 + K4.*dt);
    
    % IFFT back to update psi
    phi(:,:,rz) = ifftn(fphi);
    psi(:,:,rz) = phi(:,:,rz).^3 - phi(:,:,rz);
    
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