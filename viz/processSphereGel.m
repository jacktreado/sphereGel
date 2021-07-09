function processSphereGel(fstr,savestr)
%% FUNCTION to process data from sphere gel .xyz files

% Grid size (unit of length is mean particle radius)
NGPERRADIUS = 8;
g = 1.0/NGPERRADIUS;

% check input exists
if ~exist(fstr,'file')
    error('processSphereGel:inputFileDoesNotExist','Input file %s does not exist, ending here.\n',fstr);
else
    fprintf('String %s exists, processing...\n',fstr);
end

% file info
finfo       = dir(fstr);
fsize       = finfo.bytes;
fname       = finfo.name;
floc        = finfo.folder;

% check input is not empty
if fsize == 0
    error('processSphereGel:inputFileEmpty','Input file %s empty, ending here.\n',fname);
else
    fprintf('File %s is %0.4g MB, processing...\n',fname,fsize/1e6);
end


% extract parameters from file name
paramidx    = strfind(fname,'dlz');
paramstr    = fname(paramidx:end);
params      = sscanf(paramstr,'dlz%f_l2%f_seed%f.xyz');
dlz         = params(1);
l2          = params(2);
seed        = params(3);


% access pos data 
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

% load contact data
cmfstr      = [floc '/' fname(1:end-4) '.cm'];
fid         = fopen(cmfstr);

% contact matrix
fprintf('\t** Reading in contact matrix...');

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

% save last index before loosing rigidity
lastRigidZInd = find(z(1:end-1) > 6 & z(2:end) < 6);
lastRigidZInd = lastRigidZInd(1);


%% Construct sphere image, get structural features

% populate lattice to turn into image
fprintf('\t** Populating spheres onto lattice for image...\n');

% save positions
xpos            = posdata.xpos;
ypos            = posdata.ypos;
zpos            = posdata.zpos;

% create tmp cells to store correlation/image info
phiImg          = zeros(lastRigidZInd,1);
sk              = cell(lastRigidZInd,2);
corr2D          = cell(lastRigidZInd,3);
lambda          = zeros(lastRigidZInd,6);

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
        fprintf('** ** On frame %d / %d, Adding sphere %d / %d to image \n',ii,lastRigidZInd,nn,N);

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
        sphereInds = dr2 < radii(ii,nn)*radii(ii,nn);
        NSPHERINDS = sum(sphereInds);
        binaryLattice(sphereInds) = ones(NSPHERINDS,1);
    end
    phiImg(ii) = mean(binaryLattice,'all');

    % Use external function to plot correlation 
    fprintf('\t ** frame %d, computing correlation functions...',ii);
    [sk{ii,1}, sk{ii,2}, corr2D{ii,1}, corr2D{ii,2}, corr2D{ii,3}, lambda(ii,:)] = fourierSpaceCorrelation(binaryLattice,L(ii,:));
    fprintf('\t done computing correlation functions.\n\n');
end
fprintf('\t ** ...done populating spheres onto lattice.\n\n');


%% Save data, end function

fprintf('Saving data to savestr %s, ending.\n',savestr);
save(savestr,'fname','floc','dlz','l2','seed','N','NFRAMES','radii','L',...
    'xpos','ypos','zpos','phiImg','sk','corr2D','lambda','g','phi','cm','z');


end