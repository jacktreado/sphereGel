function cahnHilliardUniExtQD(NT, NPRINT, NSKIPSTRAIN, A, Lx, Ly, Lz0, phi0, seed, ftype)
%% FUNCTION to model phase separation with quenched disorder (with strength A)
% in 3D with variable seeds, parameters and initial conditions

% set random seed
rng(seed);

% simulation duration info
dt      = 0.01;

% get max possible matrix size
Lzmax = 2*Lz0;

% initialize Lx, Ly, Lz
Lz = Lz0;
rz = 1:Lz;

% initialize quenched-disorder
Du = A*rand(Ly,Lx,Lzmax);
Dr = A*rand(Ly,Lx,Lzmax);

% concentration field (pad to max)
phi             = zeros(Ly,Lx,Lzmax);
phi0mat         = repmat(phi0,Ly,Lx,Lzmax);
phi(:,:,rz)     = phi0mat(:,:,rz) + 0.1*randn(Ly,Lx,Lz);
psi             = (1 + Du).*phi.^3 - (1 + Dr).*phi;

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
                % extend phi
                phi(yy,xx,Lz+1) = phi(yy,xx,Lz);
                phi(yy,xx,fwd) = sc.*(phi(yy,xx,fwd) - phi(yy,xx,bwd)) + phi(yy,xx,bwd);
                
                % also extend quenched disorder
                Du(yy,xx,Lz+1) = Du(yy,xx,Lz);
                Du(yy,xx,fwd) = sc.*(Du(yy,xx,fwd) - Du(yy,xx,bwd)) + Du(yy,xx,bwd);
                
                Dr(yy,xx,Lz+1) = Dr(yy,xx,Lz);
                Dr(yy,xx,fwd) = sc.*(Dr(yy,xx,fwd) - Dr(yy,xx,bwd)) + Dr(yy,xx,bwd);
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
        fprintf('   ** t = %d: Simulating CH eq with Lz=%d, phi(1,1,Lz)=%0.5g, Dr(1,1,Lz)=%0.5g\n',tt,Lz,phi(1,1,Lz),Dr(1,1,Lz));
        
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
    psi(:,:,rz) = (1 + Du(:,:,rz)).*phi(:,:,rz).^3 - (1 + Dr(:,:,rz)).*phi(:,:,rz);
    
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
if mod(L,2) == 0
    H = L/2;
else
    H = ceil(L/2);
end
k = zeros(1,L);
Fs = 2.0*pi/L;
k(1:H) = Fs*(0:(H-1));
k(H+1:end) = Fs*((H:(L-1))-L);
end