function partialProcessSphereGel(floc,fpattern,savestr)
%% FUNCTION to process data from sphere gel .xyz files

% Grid size (unit of length is mean particle radius)
NGPERRADIUS = 6;
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
NSNAPS = 3;

%% Loop over files, extract structural features

% files to skip
fskip = false(NS,1);

% data to save
fnameList = cell(NS,1);
paramList = zeros(NS,4);   % dg, del, l2, seed
N_LIST = zeros(NS,1);

% save things only at each snapshot
SNAP_LIST = zeros(NS,NSNAPS);
phiList = zeros(NS,NSNAPS);
phiImgList = zeros(NS,NSNAPS);
KgList = zeros(NS,NSNAPS);
xiList = zeros(NS,NSNAPS);
k2List = zeros(NS,NSNAPS);

LList = cell(NS,NSNAPS);
radiiList = cell(NS,NSNAPS);
xposList = cell(NS,NSNAPS);
yposList = cell(NS,NSNAPS);
zposList = cell(NS,NSNAPS);
zList = zeros(NS,NSNAPS);
skList = cell(NS,NSNAPS);
kbinList = cell(NS,NSNAPS);

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
    L           = posdata.L;
    
    % save positions & radii
    xpos            = posdata.xpos;
    ypos            = posdata.ypos;
    zpos            = posdata.zpos;
    radii           = posdata.radii;
    
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
    
    % get apparent packing fraction over time
    pvols       = (4.0/3.0)*pi*radii.^3;
    boxvols     = L(:,1).*L(:,2).*L(:,3);
    pvsum       = sum(pvols,2);
    phi         = pvsum./boxvols;
    
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
    
    % clear memory by removing cm
    clear('cm');

    % save last index before loosing rigidity (add 1, missing frame in cm
    % list)
    lastRigidZInd = find(z(1:end-1) > 6 & z(2:end) < 6);
    if isempty(lastRigidZInd)
        fprintf('Could not find lastRigidZInd, skipping.\n');
        fskip(ss) = true;
        continue;
    else
        NRIGID = lastRigidZInd(1)+1;
        fprintf('Found %s rigid frames (including start), processing...\n',NRIGID);
    end
    
    % save number of particles and frames
    N_LIST(ss) = N;
    
    % get frame ids of snapshots
    snapids = round(linspace(2,NRIGID,NSNAPS));
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
        
        % contact info
        zList(ss,pp) = z(snapids(pp),:);
        
        % packing fraction
        phiList(ss,pp) = phi(snapids(pp));
    end
    
    % extract parameters from file name
    paramidx    = strfind(fname,'dg');
    paramstr    = fname(paramidx:end);
    params      = sscanf(paramstr,'dg%f_del%f_l2%f_seed%f.xyz');
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
        phiImgList(ss,kk) = mean(binaryLattice,'all');

        % Use external function to plot correlation 
        fprintf('\t ** frame %d (note: skipping random frame 1), computing correlation functions...',ii + 1);
        [kbinList{ss,kk}, skList{ss,kk}, ~, ~, ~, lambda, ~] = fourierSpaceCorrelation(binaryLattice,L(ii+1,:));
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
    end
    fprintf('\t ** ...done populating spheres onto lattice.\n\n');
end

% remove extra files
fnameList(fskip)    = [];
paramList(fskip,:)  = [];
N_LIST(fskip)       = [];
SNAP_LIST(fskip,:)  = [];
phiList(fskip,:)    = [];
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

NS = sum(~fskip);


%% Save data, end function

fprintf('Saving data to savestr %s, ending.\n',savestr);
save(savestr,'fnameList','paramList','N_LIST','SNAP_LIST','kbinList','phiList','LList','radiiList',...
    'xposList','yposList','zposList','zList','skList','KgList','k2List','xiList','g','floc','fpattern','savestr','NS','NSNAPS');


end