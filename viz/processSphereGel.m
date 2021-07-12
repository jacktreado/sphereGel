function processSphereGel(floc,fpattern,savestr)
%% FUNCTION to process data from sphere gel .xyz files

% Grid size (unit of length is mean particle radius)
NGPERRADIUS = 8;
g = 1.0/NGPERRADIUS;

% get list of files with pattern
simlist = dir([floc '/' fpattern '*.xyz']);
NS = length(simlist);
if NS == 0
    error('processSphereGel:noInputFilesFound','No input files found with string %s in loc %s, ending here.\n',fpattern,floc);
else
    fprintf('NF = %d files found using pattern %s, processing...\n',NS,fpattern);
end

%% Loop over files, extract structural features

% files to skip
fskip = false(NS,1);

% data to save
fnameList = cell(NS,1);
paramList = zeros(NS,3);   % dlz, l2, seed
NFRAME_LIST = zeros(NS,1);
NRIGID_LIST = zeros(NS,1);
N_LIST = zeros(NS,1);

phiList = cell(NS,2);
LList = cell(NS,1);
radiiList = cell(NS,1);
xposList = cell(NS,1);
yposList = cell(NS,1);
zposList = cell(NS,1);
cmList = cell(NS,1);
zList = cell(NS,1);

skList = cell(NS,1);
corr2DList = cell(NS,1);
lambdaList = cell(NS,1);
evList = cell(NS,1);

% phi binning: assume 0.9 -> 0.1
phimax = 0.9;
phimin = 0.1;
dphi = 0.01;
phic = phimax:-dphi:phimin;
NBINS = length(phic);
phie = (phimax+0.5*dphi):-dphi:(phimin-0.5*dphi);

bStats = nan(NBINS,NS);
cStats = nan(NBINS,NS);
KgStats = nan(NBINS,NS);
k1Stats = nan(NBINS,NS);
k2Stats = nan(NBINS,NS);

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
    posdata     = readSGelXYZ([ffldr '/' fname]);
    N           = posdata.N;
    NFRAMES     = posdata.NFRAMES;
    radii       = posdata.radii;
    L           = posdata.L;
    
    % save positions
    xpos            = posdata.xpos;
    ypos            = posdata.ypos;
    zpos            = posdata.zpos;
    
    % load contact data
    cmfstr      = [ffldr '/' fname(1:end-4) '.cm'];
    % check contact file exists
    if ~exist(cmfstr,'file')
        fprintf('Contact file %s does not exist, skipping.\n',cmfstr);
        fskip(ff) = true;
        continue;
    else
        fprintf('\t\t CM File %s exists, processing...\n',cmfstr);
        cmfinfo     = dir(cmfstr);
        cmfsize     = cmfinfo.bytes;
        cmfname     = cmfinfo.name;
        if fsize == 0
            erfprintfror('CM file %s empty, skipping.\n',cmfname);
            fskip(ff) = true;
            continue;
        else
            fprintf('\t\t CM File %s is %0.4g MB, processing...\n',cmfname,cmfsize/1e6);
        end
    end
    
    % contact matrix
    % NOTE: 1 less frame in CM than in xyz
    fprintf('\t** Reading in contact matrix...');
    fid         = fopen(cmfstr);
    NPAIRS      = 0.5*N*(N-1);
    frmt        = repmat('%f ',1,NPAIRS);
    data        = textscan(fid,frmt,NFRAMES-1);
    cm          = cell2mat(data);
    
    % close file, print to console
    fclose(fid);
    fprintf('done!\n\n');
    
    % get mean contact number z
    Nc          = sum(cm,2);
    z           = 2.0*Nc./N;

    % save last index before loosing rigidity (add 1, missing frame in cm
    % list)
    lastRigidZInd   = find(z(1:end-1) > 6 & z(2:end) < 6);
    NRIGID          = lastRigidZInd(1)+1;
    
    % save number of particles and frames
    N_LIST(ss) = N;
    NFRAME_LIST(ss) = NFRAMES;
    NRIGID_LIST(ss) = NRIGID;
    
    % save rigid part of contact matrix (get 2 - NRIGID, first frame is
    % random init condition)
    cmList{ss} = cm(1:(NRIGID-1),:);
    zList{ss} = z;
    
    % save positions in rigid network
    xposList{ss} = xpos(2:NRIGID,:);
    yposList{ss} = ypos(2:NRIGID,:);
    zposList{ss} = zpos(2:NRIGID,:);
    radiiList{ss} = radii(2:NRIGID,:);
    LList{ss} = L(1:NRIGID,:);
    
    % extract parameters from file name
    paramidx    = strfind(fname,'dlz');
    paramstr    = fname(paramidx:end);
    params      = sscanf(paramstr,'dlz%f_l2%f_seed%f.xyz');
    dlz         = params(1);
    l2          = params(2);
    seed        = params(3);
    
    % save
    paramList(ss,:) = [dlz, l2, seed];

    % get apparent packing fraction over time
    pvols       = (4.0/3.0)*pi*radii.^3;
    boxvols     = L(:,1).*L(:,2).*L(:,3);
    pvsum       = sum(pvols,2);
    phi         = pvsum./boxvols;
    
    % save
    phiList{ss,1} = phi;
    
    
    % -- Construct sphere image, get structural features

    % populate lattice to turn into image
    fprintf('\t** Populating spheres onto lattice for image...\n');

    % create tmp cells to store correlation/image info
    phiImg          = zeros(NRIGID-1,1);
    sk              = cell(NRIGID-1,2);
    corr2D          = cell(NRIGID-1,3);
    lambda          = zeros(NRIGID-1,3);
    eV              = cell(NRIGID-1,1);

    % loop over frames with z > 6
    for ii = 1:NRIGID-1
        % get positions in frame
        x = xpos(ii+1,:)';
        y = ypos(ii+1,:)';
        z = zpos(ii+1,:)';

        % get box size in frame
        Lx = L(ii+1,1);
        Ly = L(ii+1,2);
        Lz = L(ii+1,3);

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
            fprintf('** ** (%d/%d) On frame %d / %d (note: skipping frame 1), Adding sphere %d / %d to image \n',ss,NS,ii+1,NRIGID,nn,N);

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
            sphereInds = dr2 < radii(ii+1,nn)*radii(ii+1,nn);
            NSPHERINDS = sum(sphereInds);
            binaryLattice(sphereInds) = ones(NSPHERINDS,1);
        end
        phiImg(ii) = mean(binaryLattice,'all');

        % Use external function to plot correlation 
        fprintf('\t ** frame %d (note: skipping random frame 1), computing correlation functions...',ii + 1);
        [sk{ii,1}, sk{ii,2}, corr2D{ii,1}, corr2D{ii,2}, corr2D{ii,3}, lambda(ii,:), eV{ii}] = fourierSpaceCorrelation(binaryLattice,L(ii+1,:));
        fprintf('\t done computing correlation functions.\n\n');
        
        % tmp correlation length
        kbtmp = sk{ii,1};
        sktmp = sk{ii,2};
        k1tmp = nansum(kbtmp.*sktmp)/nansum(sktmp);
        
        % tmp structural stats (note eig are sorted in ascending order)
        lb3 = lambda(ii,3);
        lb2 = lambda(ii,2);
        lb1 = lambda(ii,1);
                
        btmp = lb3 - 0.5*(lb2 + lb1);
        ctmp = lb2 - lb1;
        Kgtmp = sqrt(lb1 + lb2 + lb3);
        k2tmp = (btmp^2 + 0.75*ctmp^2)/(Kgtmp^4);
        
        % add to stats
        phitmp = phi(ii+1);
        phiind = phitmp < phie(1:end-1) & phitmp > phie(2:end);
        
        bStats(phiind,ss) = btmp;
        cStats(phiind,ss) = ctmp;
        KgStats(phiind,ss) = Kgtmp;
        k1Stats(phiind,ss) = k1tmp;
        k2Stats(phiind,ss) = k2tmp;
    end
    fprintf('\t ** ...done populating spheres onto lattice.\n\n');

    % save 
    phiList{ss,2} = phiImg;
    skList{ss} = sk;
    corr2DList{ss} = corr2D;
    lambdaList{ss} = lambda;
    evList{ss} = eV;
end

% remove extra files
fnameList(fskip)    = [];
paramList(fskip,:)  = [];
NFRAME_LIST(fskip)  = [];
NRIGID_LIST(fskip)  = [];
N_LIST(fskip)       = [];
phiList(fskip,:)    = [];
LList(fskip)        = [];
radiiList(fskip)    = [];
xposList(fskip)     = [];
yposList(fskip)     = [];
zposList(fskip)     = [];
cmList(fskip)       = [];
zList(fskip)        = [];
skList(fskip)       = [];
corr2DList(fskip)   = [];
lambdaList(fskip)   = [];
evList(fskip)       = [];


%% Save data, end function

fprintf('Saving data to savestr %s, ending.\n',savestr);
save(savestr,'fnameList','paramList','N_LIST','NFRAME_LIST','NRIGID_LIST','phiList','LList','radiiList',...
    'xposList','yposList','zposList','cmList','zList','skList','corr2DList','lambdaList','evList','g',...
    'bStats','cStats','KgStats','k1Stats','k2Stats','phic','floc','fpattern','savestr','NS');


end