%% Plot contacts and structure during sim
% 
% TO-DO:
%   -- Compute various structural quantities for Jammed spheres
%   -- Figure out how to best compare to plant samples ... always scale by
%   xi?
%   -- Turn into function to put onto cluster, part of workflow
%   -- Add larger systems to cluster (~2000 particles)

clear;
close all;
clc;

% simulation info
Nstr        = '400';
drstr       = '0.15';
dphistr     = '0.0005';
seedstr     = '1';

% sim pattern
fprefix     = ['sgel_N' Nstr '_dr' drstr '_dphi' dphistr '_dlz'];
fpattern    = ['xyz/' fprefix '*_seed' seedstr '.xyz'];

% file list
flist       = dir(fpattern);
NF          = length(flist);
if NF == 0
    fprintf('\t** No files found with pattern %s, ending.\n',fpattern);
    error('No Files found');
else
    fprintf('-- Found %d files with pattern, processing.\n\n',NF);
end

% Grid size (unit of length is mean particle radius)
NGPERRADIUS = 8;
g = 1.0/NGPERRADIUS;

%% Get parameters

% max number of different parameters
MAXNUMDLZSIMS   = 1e3;
MAXNUML2SIMS    = 1e3;

% save last stable phi
dlz             = zeros(MAXNUMDLZSIMS,1);
l2              = zeros(MAXNUML2SIMS,1);
params          = zeros(NF,2);
inds            = zeros(NF,2);
last            = zeros(2,1);

% get parameters
for ff = 1:NF
    % file info
    fname       = flist(ff).name;
    floc        = flist(ff).folder;
    fstr        = [floc '/' fname];
    
    % print simulation name to file
    fprintf('\t** processing file %s\n',fname);
    
    % get parameters, determine simulation ind
    ptmp        = sscanf(fname,[fprefix '%f_l2%f']);
    dlztmp      = ptmp(1);
    l2tmp       = ptmp(2);
    
    dlzdiff     = abs(dlztmp - dlz) < 1e-8;
    dlzcheck    = sum(dlzdiff);
    if dlzcheck == 0
        last(1) = last(1) + 1;
        dlz(last(1)) = dlztmp;
        dlzind = last(1);
    else
        dlzind = find(dlzdiff);
    end
    
    l2diff      = abs(l2tmp - l2) < 1e-8;
    l2check     = sum(l2diff);
    if l2check == 0
        last(2) = last(2) + 1;
        l2(last(2)) = l2tmp;
        l2ind = last(2);
    else
        l2ind = find(l2diff);
    end
    
    % save param values and indices
    params(ff,:)    = [dlztmp l2tmp];
    inds(ff,:)      = [dlzind l2ind];
end


% delete extra entries
dlz(last(1)+1:end)          = [];
l2(last(2)+1:end)           = [];

% number of parameters
NDLZ                        = last(1);
NL2                         = last(2);


%% Loop over files, get z and structure over time

% memory allocation
lastphi         = zeros(NDLZ,NL2);
phiImgList      = cell(NF,1);
phi0List        = cell(NF,1);
skList          = cell(NF,1);
corr2DList      = cell(NF,1);
lambdaList      = cell(NF,1);

% open figure windows
figure(1), clf, hold on, box on;

% marker character cells
lstyles         = {'-','--','-.',':'};
NLS             = length(lstyles);
clr             = jet(NL2);

% loop over simulations, read in data
for ff = 1:NF
    % file info
    fname       = flist(ff).name;
    floc        = flist(ff).folder;
    fstr        = [floc '/' fname];
    
    % print simulation name to file
    fprintf('\t** processing file %s\n',fname);
    
    % get parameters, determine simulation ind
    ptmp        = params(ff,:);
    dlztmp      = ptmp(1);
    l2tmp       = ptmp(2);
    
    dlzind      = inds(ff,1);
    l2ind       = inds(ff,2);
    
    % load xyz data
    posdata     = readSGelXYZ(fstr);
    N           = posdata.N;
    NFRAMES     = posdata.NFRAMES;
    radii       = posdata.radii;
    L           = posdata.L;

    % get packing fraction over time
    pvols       = (4.0/3.0)*pi*radii.^3;
    boxvols     = L(:,1).*L(:,2).*L(:,3);
    pvsum       = sum(pvols,2);
    phi         = pvsum./boxvols;
    phi0List{ff} = phi;
    
    % load contact data
    cmfstr      = [floc '/' fname(1:end-4) '.cm'];
    fid         = fopen(cmfstr);

    % contact matrix
    fprintf('\t** reading in contact matrix...');
    
    NPAIRS      = 0.5*N*(N-1);
    frmt        = repmat('%f ',1,NPAIRS);
    data        = textscan(fid,frmt,NFRAMES);
    cm          = cell2mat(data);
    
    % close file, print to console
    fclose(fid);
    fprintf('done!\n\n');
    
    % plot z vs phi
    Nc          = sum(cm,2);
    z           = 2.0*Nc./N;
    
    % create marker string for plotting
    lind = mod(dlzind-1,NLS) + 1;
    lc = lstyles{lind};
    
    figure(1),
    plot(phi(2:end),z,[lc 'o'],'linewidth',1.5,'color',clr(l2ind,:),'markersize',6);
    
    % save last phi before loosing rigidity
    lastRigidZInd = find(z(1:end-1) > 6 & z(2:end) < 6);
    lastRigidZInd = lastRigidZInd(1);
    lastphi(dlzind,l2ind) = phi(lastRigidZInd);
    
    
    % populate lattice to turn into image
    fprintf('\t** populating polydisperse spheres onto lattice for image...\n');
    
    % save positions
    xpos            = posdata.xpos;
    ypos            = posdata.ypos;
    zpos            = posdata.zpos;
    
    % create tmp cells to store correlation/image info
    phiImg      = zeros(lastRigidZInd,1);
    sk          = cell(lastRigidZInd,2);
    corr2D      = cell(lastRigidZInd,3);
    lambda      = zeros(lastRigidZInd,6);
    
    % loop over frames with z > 6
    for ii = 1:lastRigidZInd
        % get positions in frame
        x = xpos(ii,:)';
        y = ypos(ii,:)';
        z = zpos(ii,:)';
        
        % get box size in frame
        Lx = L(ii,1);
        Ly = L(ii,2);
        Lz = L(ii,3);
        
        % get grid size
        gx = 0:g:(Lx-g);
        gy = 0:g:(Ly-g);
        gz = 0:g:(Lz-g);
        
        nx = length(gx);
        ny = length(gy);
        nz = length(gz);
        
        % get grid of all points
        [GRIDX, GRIDY, GRIDZ] = meshgrid(gx,gy,gz);
        
        % fill image matrix
        binaryLattice = zeros(nx,ny,nz);
        for nn = 1:N
            % print to console
            fprintf('(%d,%d)  ** ** In frame %d / %d, Adding sphere %d / %d to image \n',ff,NF,ii,lastRigidZInd,nn,N);

            % particle positions
            xi = xpos(nn);
            yi = ypos(nn);
            zi = zpos(nn);

            % get distances
            dx = xi - GRIDX(:);
            dx = dx - Lx*round(dx./Lx);

            dy = yi - GRIDY(:);
            dy = dy - Ly*round(dy./Ly);

            dz = zi - GRIDZ(:);
            dz = dz - Lz*round(dz./Lz);

            dr2 = dx.*dx + dy.*dy + dz.*dz;
            sphereInds = dr2 < radii(ii,nn)*radii(ii,nn);
            NSPHERINDS = sum(sphereInds);
            binaryLattice(sphereInds) = ones(NSPHERINDS,1);
        end
        phiImg(ii) = mean(binaryLattice,'all');
        
        % Use external function to plot correlation 
        fprintf('\t ** frame %d, computing correlation functions...',ii);
        [sk{ii,1}, sk{ii,2}, corr2D{ii,1}, corr2D{ii,2}, corr2D{ii,3}, lambda(ii,:)] = fourierSpaceCorrelation(binaryLattice,L(ff,:));
        fprintf('\t done computing correlation functions.\n\n');
    end
    fprintf('\t ** ...done populating spheres onto lattice.\n\n');
    
    % save correlation info
    phiImgList{ff} = phiImg;
    skList{ff} = sk;
    corr2DList{ff} = corr2D;
    lambdaList{ff} = lambda;
end

%% Save all data

% save('sgel_processed_data.mat','lambdaList','corr2DList','skList','phiImgList','phi0List','lastphi','params','inds','flist','g','NF');
    

%% Plotting

% load data
load('sgel_processed_data.mat');
NF = length(flist);

% number of color bins
nclrbins = 50;

% contact info
figure(1),
xlabel('$\phi_0$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;


% plot 1D S(k) over time while sim is rigid
figure(10), clf, hold on, box on;
skInfo = skList{15};
nftmp = length(skInfo);
clr = jet(nftmp);
for ff = 1:nftmp
    ktmp = skInfo{ff,1};
    sktmp = skInfo{ff,2};
    plot(ktmp,sktmp,'-','linewidth',2,'color',clr(ff,:));
end
xlabel('$k$','Interpreter','latex');
ylabel('$S(k)$','Interpreter','latex');
ax = gca;
ax.FontSize = 18;
ax.XScale = 'log';
ax.YScale = 'log';


% Make 3D plots of structural features right at contact
phiEnd = zeros(NF,1);
xi = zeros(NF,1);
Kg = zeros(NF,1);
b = zeros(NF,1);
c = zeros(NF,1);
kap2 = zeros(NF,1);
for ff = 1:NF
    lambdatmp = lambdaList{ff};
    l1 = lambdatmp(end,1);
    l2 = lambdatmp(end,2);
    l3 = lambdatmp(end,3);
    
    Kg(ff) = sqrt(l1 + l2 + l3);
    b(ff) = l3 - 0.5*(l1 + l2);
    c(ff) = l2 - l1;
    kap2(ff) = (b(ff)^2 + 0.75*c(ff)^2)/(Kg(ff)^4);
    
    phiImgTmp = phi0List{ff};
    phiEnd(ff) = phiImgTmp(end);
    
    skInfoTmp = skList{ff};
    kbintmp = skInfoTmp{end,1};
    sktmp = skInfoTmp{end,2};
    k1 = sum(kbintmp.*sktmp)./sum(sktmp);
    xi(ff) = (2.0*pi)/k1;
end

figure(11), clf, hold on, box on;

kap2min = log10(0.95*min(kap2));
kap2max = log10(1.05*max(kap2));
kap2clrbins = logspace(kap2min,kap2max,nclrbins+1);
clr = jet(nclrbins);
for ff = 1:NF
    kap2tmp = kap2(ff);
    kap2ind = kap2tmp > kap2clrbins(1:end-1) & kap2tmp < kap2clrbins(2:end);
    plot(params(ff,1),params(ff,2),'ko','markersize',10,'markerfacecolor',clr(kap2ind,:));
end
xlabel('$\delta L_z$','Interpreter','latex');
ylabel('$l_2$','Interpreter','latex');
ax = gca;
ax.FontSize = 18;
ax.YScale = 'log';



figure(12), clf, hold on, box on;

Kgmin = log10(0.95*min(Kg));
Kgmax = log10(1.05*max(Kg));
Kgclrbins = logspace(Kgmin,Kgmax,nclrbins+1);
clr = jet(nclrbins);
for ff = 1:NF
    Kgtmp = Kg(ff);
    Kgind = Kgtmp > Kgclrbins(1:end-1) & Kgtmp < Kgclrbins(2:end);
    plot(params(ff,1),params(ff,2),'ko','markersize',10,'markerfacecolor',clr(Kgind,:));
end
xlabel('$\delta L_z$','Interpreter','latex');
ylabel('$l_2$','Interpreter','latex');
ax = gca;
ax.FontSize = 18;
ax.YScale = 'log';




% make big 4D plots similar to plant scans
figure(21), clf, hold on, box on;
ax21 = gca;
figure(22), clf, hold on, box on;
ax22 = gca;
figure(23), clf, hold on, box on;
ax23 = gca;
figure(24), clf, hold on, box on;
ax24 = gca;
figure(25), clf, hold on, box on;
ax25 = gca;
clr = jet(nclrbins);
for ff = 1:NF
    % stuff to compute
    lambdatmp = lambdaList{ff};
    skInfoTmp = skList{ff};

    % get phi colors
    phitmp = phi0List{ff};
    phimin = 0.95*min(phitmp);
    phimax = 1.05*max(phitmp);
    phibins = linspace(phimin,phimax,nclrbins+1);
    
    % stuff to save
    NFRAMES = size(lambdatmp,1);
    
    Kgtmp = zeros(NFRAMES,1);
    btmp = zeros(NFRAMES,1);
    ctmp = zeros(NFRAMES,1);
    kap2tmp = zeros(NFRAMES,1);
    xitmp = zeros(NFRAMES,1);
    
    % loop over frames, plot structural quantities
    for ii = 1:NFRAMES
        l1 = lambdatmp(ii,1);
        l2 = lambdatmp(ii,2);
        l3 = lambdatmp(ii,3);

        Kgtmp(ii) = sqrt(l1 + l2 + l3);
        btmp(ii) = l3 - 0.5*(l1 + l2);
        ctmp(ii) = l2 - l1;
        kap2tmp(ii) = (btmp(ii)^2 + 0.75*ctmp(ii)^2)/(Kgtmp(ii)^4);
        
        kbintmp = skInfoTmp{ii,1};
        sktmp = skInfoTmp{ii,2};
        kbintmp = kbintmp(~isnan(sktmp));
        sktmp = sktmp(~isnan(sktmp));
        k1 = sum(kbintmp.*sktmp)./sum(sktmp);
        xitmp(ii) = (2.0*pi)/k1;
        
        % get bin for phi color
        phi = phitmp(ii);
        phiind = phi > phibins(1:end-1) & phi < phibins(2:end);
        
        % plot Kg, b, k2, phi (color)
        plot3(ax21,Kgtmp(ii)*xitmp(ii),btmp(ii)*(xitmp(ii)^2),kap2tmp(ii),'ko','markersize',10,'markerfacecolor',clr(phiind,:));
        
        % plot c, b, k2, phi (color)
        plot3(ax22,ctmp(ii)*(xitmp(ii)^2),btmp(ii)*(xitmp(ii)^2),kap2tmp(ii),'ko','markersize',10,'markerfacecolor',clr(phiind,:));
        
        % plot xi, phi, kap2
        plot3(ax25,xitmp(ii),phi,kap2tmp(ii),'ko','markersize',10,'markerfacecolor','b');
        
        fprintf('** In sim %d / %d, frame %d / %d\n',ff,NF,ii,NFRAMES);
    end
    
    % plot 2D things
    plot(ax23,xitmp,Kgtmp,'ks','markersize',10);
    plot(ax24,xitmp,kap2tmp,'ks','markersize',10);
end

% label axes
xlabel(ax21,'$K_g \xi$','Interpreter','latex');
ylabel(ax21,'$b \xi^2$','Interpreter','latex');
zlabel(ax21,'$\kappa^2$','Interpreter','latex');
ax21.FontSize = 24;
ax21.XScale = 'log';
ax21.YScale = 'log';
ax21.ZScale = 'log';

% label axes
xlabel(ax25,'$\xi$','Interpreter','latex');
ylabel(ax25,'$\phi$','Interpreter','latex');
zlabel(ax25,'$\kappa^2$','Interpreter','latex');
ax25.FontSize = 24;
ax25.ZScale = 'log';


figure(22),
xlabel('$c \xi^2$','Interpreter','latex');
ylabel('$b \xi^2$','Interpreter','latex');
zlabel('$\kappa^2$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;
ax.XScale = 'log';
ax.YScale = 'log';
ax.ZScale = 'log';


figure(23),
xlabel('$\xi$','Interpreter','latex');
ylabel('$K_g$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;

figure(24),
xlabel('$\xi$','Interpreter','latex');
ylabel('$\kappa^2$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;
