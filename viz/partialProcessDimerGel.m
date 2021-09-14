function partialProcessDimerGel(floc,fpattern,savestr)
%% FUNCTION to process data from sphere gel .xyz files
% NOTE: need to redo analysis of data on cluster (08/27), does not include
% phiImg

% Grid size (unit of length is mean particle radius)
NGPERRADIUS = 8;
g = 1.0/NGPERRADIUS;

% get list of files with pattern
simlist = dir([floc '/' fpattern '*.xyz']);
NS = length(simlist);
if NS == 0
    error('partialProcessSphereGel:noInputFilesFound','No input files found with string %s in loc %s, ending here.\n',fpattern,floc);
else
    fprintf('NF = %d files found using pattern %s, processing...\n',NS,fpattern);
end

% number of snapshots
NSNAPS = 10;

%% Loop over files, extract structural features

% files to skip
fskip = false(NS,1);

% data to save
fnameList = cell(NS,1);
paramList = zeros(NS,4);   % dg, del, l2, seed
N_LIST = zeros(NS,1);

% save things only at each snapshot
SNAP_LIST = zeros(NS,NSNAPS);
phiImgList = zeros(NS,NSNAPS);
KgList = zeros(NS,NSNAPS);
xiList = zeros(NS,NSNAPS);
k2List = zeros(NS,NSNAPS);
LList = cell(NS,NSNAPS);
radiiList = cell(NS,NSNAPS);
xposList = cell(NS,NSNAPS);
yposList = cell(NS,NSNAPS);
zposList = cell(NS,NSNAPS);
skList = cell(NS,NSNAPS);
kbinList = cell(NS,NSNAPS);
SiList = cell(NS,NSNAPS);

% save everything
zList = cell(NS,1);
SnormList = cell(NS,1);
SshearList = cell(NS,1);


% loop over sim files
for ss = 1:NS
    % get file info
    fsize       = simlist(ss).bytes;
    fname       = simlist(ss).name;
    ffldr       = simlist(ss).folder;

    % check input is not empty
    if fsize == 0
        fprintf('Input file %s empty, skipping.\n',fname);
        fskip(ss) = true;
        continue;
    else
        fprintf('File %s is %0.4g MB, processing...\n',fname,fsize/1e6);
    end
    
    % save file name
    fnameList{ss} = fname;
    
    % access pos data 
    posdata     = readDGelXYZ([ffldr '/' fname]);
    N           = posdata.N;
    L           = posdata.L;
    
    % save positions & radii
    xpos            = posdata.xpos;
    ypos            = posdata.ypos;
    zpos            = posdata.zpos;
    radii           = posdata.radii;
    
    % save z and S
    Snorm = posdata.Snorm;
    Sshear = posdata.Sshear;
    z = mean(posdata.z,2);
    
    % save all info
    SnormList{ss} = Snorm;
    SshearList{ss} = Sshear;
    zList{ss} = z;

    % save last index before loosing rigidity (add 1, missing frame in cm
    % list)
    if min(z) > 6
        lastRigidZInd = posdata.NFRAMES;
        NRIGID = lastRigidZInd;
    else
        lastRigidZInd = find(z(1:end-1) > 6 & z(2:end) < 6);
        if isempty(lastRigidZInd)
            fprintf('Could not find lastRigidZInd, skipping.\n');
            fskip(ss) = true;
            continue;
        else
            NRIGID = lastRigidZInd(1);
            fprintf('Found %s rigid frames (including start), processing...\n',NRIGID);
        end
    end
    
    % save number of particles and frames
    N_LIST(ss) = N;
    
    % get frame ids of snapshots
    snapids = round(linspace(1,NRIGID,NSNAPS));
    SNAP_LIST(ss,:) = snapids';
    
    % save into cells
    for pp = 1:NSNAPS
        % particle info
        xposList{ss,pp} = xpos(snapids(pp),:);
        yposList{ss,pp} = xpos(snapids(pp),:);
        zposList{ss,pp} = xpos(snapids(pp),:);
        radiiList{ss,pp} = radii(snapids(pp),:);
        
        % box info
        LList{ss,pp} = L(snapids(pp),:);
    end
    
    % extract parameters from file name
    paramidx    = strfind(fname,'_dg');
    paramstr    = fname(paramidx:end);
    params      = sscanf(paramstr,'_dg%f_del%f_l2%f_seed%f.xyz');
    dg          = params(1);
    del         = params(2);
    l2          = params(3);
    seed        = params(4);
    
    % save
    paramList(ss,:) = [dg, del, l2, seed];
    
    
    % -- Construct sphere image, get structural features

    % populate lattice to turn into image
    fprintf('\t** Populating spheres onto lattice for image...\n');

    % loop over frames with z > 6
    for kk = 1:NSNAPS
        % ii is frame
        ii = snapids(kk);
        
        % get positions in frame
        x = xpos(ii,:)';
        y = ypos(ii,:)';
        z = zpos(ii,:)';
        r = radii(ii,:)';
        r2 = r.*r;

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
        
         % print to console
        fprintf('** ** (%d/%d, %d/%d) On frame %d, adding sphere %d / %d to image \n',ss,NS,kk,NSNAPS,ii);

        % fill image matrix
        binaryLattice = zeros(ny,nx,nz);
        for nn = 1:N
            % particle positions
            xi = x(nn);
            yi = y(nn);
            zi = z(nn);

            % get distances
            dx = xi - GRIDX(:);
            dx = dx - Lx*round(dx./Lx);

            dy = yi - GRIDY(:);
            dy = dy - Ly*round(dy./Ly);

            dz = zi - GRIDZ(:);
            dz = dz - Lz*round(dz./Lz);

            dr2 = dx.*dx + dy.*dy + dz.*dz;
            sphereInds = dr2 < r2(nn);
            NSPHERINDS = sum(sphereInds);
            binaryLattice(sphereInds) = ones(NSPHERINDS,1);
        end
        phiImgList(ss,kk) = mean(binaryLattice,'all');

        % Use external function to plot correlation 
        fprintf('\t ** frame %d, computing correlation functions...',ii);
        [kbinList{ss,kk}, skList{ss,kk}, ~, ~, ~, lambda, ~] = fourierSpaceCorrelation(binaryLattice,L(ii,:));
        fprintf('\t done computing correlation functions.\n\n');
        
        % sort eigenvalues
        lambda = sort(lambda);
        
        % first moment of s(k)
        xiList(ss,kk) = sum(skList{ss,kk}.*kbinList{ss,kk})./sum(skList{ss,kk});
        
        % correlation anistropy
        Kgtmp = sqrt(lambda(1) + lambda(2) + lambda(3));
        btmp = lambda(3) - 0.5*(lambda(1) + lambda(2));
        ctmp = lambda(2) - lambda(1);
        k2tmp = (btmp^2 + 0.75*ctmp^2)/Kgtmp^4;
        
        % save info
        KgList(ss,kk) = Kgtmp;
        k2List(ss,kk) = k2tmp;
        
        % compute nematic order
        SiList{ss,kk} = nematicOrder(x,y,z,r,Lx,Ly,Lz,4.0);
    end
    fprintf('\t ** ...done populating spheres onto lattice.\n\n');
    
    % save matfile as you go
    saveStruct.fnameList = fnameList(~fskip(1:ss));
    saveStruct.paramList = paramList(~fskip(1:ss),:);
    saveStruct.N_LIST = N_LIST(~fskip(1:ss));
    saveStruct.SNAP_LIST = SNAP_LIST(~fskip(1:ss),:);
    saveStruct.phiImgList = phiImgList(~fskip(1:ss),:);
    saveStruct.LList = LList(~fskip(1:ss),:);
    saveStruct.radiiList = radiiList(~fskip(1:ss),:);
    saveStruct.xposList = xposList(~fskip(1:ss),:);
    saveStruct.yposList = yposList(~fskip(1:ss),:);
    saveStruct.zposList = zposList(~fskip(1:ss),:);
    saveStruct.zList = zList(~fskip(1:ss),:);
    saveStruct.skList = skList(~fskip(1:ss),:);
    saveStruct.kbinList = kbinList(~fskip(1:ss),:);
    saveStruct.KgList = KgList(~fskip(1:ss),:);
    saveStruct.k2List = k2List(~fskip(1:ss),:);
    saveStruct.xiList = xiList(~fskip(1:ss),:);
    saveStruct.SiList = SiList(~fskip(1:ss),:);
    saveStruct.SnormList = SnormList(~fskip(1:ss),:);
    saveStruct.NS = sum(~fskip(1:ss));
    save(savestr,'-struct','saveStruct');
end

% remove extra files
fnameList(fskip)    = [];
paramList(fskip,:)  = [];
N_LIST(fskip)       = [];
SNAP_LIST(fskip,:)  = [];
phiImgList(fskip,:) = [];
LList(fskip,:)      = [];
radiiList(fskip,:)  = [];
xposList(fskip,:)   = [];
yposList(fskip,:)   = [];
zposList(fskip,:)   = [];
zList(fskip,:)      = [];
skList(fskip,:)     = [];
kbinList(fskip,:)   = [];
KgList(fskip,:)     = [];
k2List(fskip,:)     = [];
xiList(fskip,:)     = [];
SiList(fskip,:)     = [];
SnormList(fskip,:)  = [];
SshearList(fskip,:) = [];

NS = sum(~fskip);


%% Save data, end function

fprintf('Saving data to savestr %s, ending.\n',savestr);
save(savestr,'fnameList','paramList','N_LIST','SNAP_LIST','kbinList','LList','radiiList','phiImgList','savestr','NS','NSNAPS',...
    'xposList','yposList','zposList','zList','skList','KgList','k2List','xiList','g','floc','fpattern',...
    'SiList','SnormList','SshearList');


end



%% Function to compute nematic order

function S = nematicOrder(x,y,z,r,Lx,Ly,Lz,cutoff)
% number of dimers
N = length(x)/2;

% local nematic order parameter
S = zeros(N,1);

% compute all directors
ex = zeros(N,1);
ey = zeros(N,1);
ez = zeros(N,1);
cx = zeros(N,1);
cy = zeros(N,1);
cz = zeros(N,1);
for nn = 1:N
    % distances
    dx = x(2*nn) - x(2*nn-1);
    dx = dx - Lx*round(dx/Lx);
    
    dy = y(2*nn) - y(2*nn-1);
    dy = dy - Ly*round(dy/Ly);
    
    dz = z(2*nn) - z(2*nn-1);
    dz = dz - Lz*round(dz/Lz);
    
    % normalize
    dnrm = sqrt(dx*dx + dy*dy + dz*dz);
    
    % unit vectors
    ex(nn) = dx/dnrm;
    ey(nn) = dy/dnrm;
    ez(nn) = dz/dnrm;
    
    % also save centers of mass
    cx(nn) = x(2*nn-1) + 0.5*dx;
    cy(nn) = y(2*nn-1) + 0.5*dy;
    cz(nn) = z(2*nn-1) + 0.5*dz;
end

% compute center-to-center distances
dcx = cx - cx';
dcx = dcx - Lx*round(dcx./Lx);

dcy = cy - cy';
dcy = dcy - Ly*round(dcy./Ly);

dcz = cz - cz';
dcz = dcz - Lz*round(dcz./Lz);

dc2 = dcx.^2 + dcy.^2 + dcz.^2;

% loop back over dimers, average directors of neighbors based on cutoff
for nn = 1:N
    neighborInds = dc2(:,nn) < (2.0*cutoff*r(nn))^2;
    neighborInds(nn) = 0;
    
    if sum(neighborInds) > 0
        Sij = 1.5*abs(ex(nn).*ex(neighborInds) + ey(nn).*ey(neighborInds)) - 0.5;
        S(nn) = sum(Sij)/sum(neighborInds);
    else
        S(nn) = 0.0;
    end
end

end