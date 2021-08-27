function data = readDGelXYZ(fstr)
%% read in xyz position data from sphere gel xyz file

% open file
fid = fopen(fstr);

% read number of particles
ndata = textscan(fid,'%f',1);
N = ndata{1};

% print to console
finfo = dir(fstr);
fsize = finfo.bytes;
fprintf('reading in data from %s, N = %d, size = %0.4g MB\n',fstr,N,fsize/1e6);

% stuff to save
MAXNFRAMES = 1e5;
L = zeros(MAXNFRAMES,3);
xpos = zeros(MAXNFRAMES,N);
ypos = zeros(MAXNFRAMES,N);
zpos = zeros(MAXNFRAMES,N);
radii = zeros(MAXNFRAMES,N);
z = zeros(MAXNFRAMES,N);
frame = 1;

% loop over file, read in data
while ~feof(fid)
    % print to console
    fprintf('\t ** on frame %d\n',frame);
    
    % read box info
    ldata = textscan(fid,'Lattice="%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f"',1);
    tline = fgetl(fid);
    L(frame,1) = ldata{1};
    L(frame,2) = ldata{2};
    L(frame,3) = ldata{3};
    
    % read position info
    posdata = textscan(fid,'%*c %f %f %f %f %f',N);
    xpos(frame,:) = posdata{1}';
    ypos(frame,:) = posdata{2}';
    zpos(frame,:) = posdata{3}';
    radii(frame,:) = posdata{4}';
    z(frame,:) = posdata{5}';
    
    % read next cell
    tline = fgetl(fid);
    tline = fgetl(fid);
    if tline == -1
        fprintf('\t ** End of file found, NFRAMES = %d\n',frame);
        break;
    end
    
    % increment frame label
    frame = frame + 1;
end

% save number of frames
NFRAMES         = frame;

% delete extra
L(NFRAMES+1:end,:) = [];
xpos(NFRAMES+1:end,:) = [];
ypos(NFRAMES+1:end,:) = [];
zpos(NFRAMES+1:end,:) = [];
radii(NFRAMES+1:end,:) = [];
z(NFRAMES+1:end,:) = [];

% save data in struct
data            = struct('N',N);
data.NFRAMES    = NFRAMES;
data.L          = L;
data.xpos       = xpos;
data.ypos       = ypos;
data.zpos       = zpos;
data.radii      = radii;
data.z          = z;

% close fid
fclose(fid);

end