function processCahnHilliardUniExt(simloc,simfrmt,savestr)
%% FUNCTION to analyze structural features of cahn hilliard sims, make movies

% check that directory exists
if ~exist(simloc,'dir')
    error('process_cahnHilliardUniExt:dirDoesNotExist','Directory %s not found, ending.',simloc);
else
    fprintf('Looking in directory %s for data\n',simloc);
end

% get list of files
searchstr = [simloc '/' simfrmt '*.pos'];
simlist = dir(searchstr);
NSIMS = length(simlist);
if NSIMS == 0
    error('process_cahnHilliardUniExt:noFilesFound','No files found matching format %s, ending\n',simfrmt);
else
    fprintf('NSIMS = %d found in %s, processing ...\n',NSIMS,simloc);
end

% sort
frameID = zeros(NSIMS,1);
LList = zeros(NSIMS,3);
for ss = 1:NSIMS
    fname = simlist(ss).name;
    fstr = [simloc '/' fname];
    
    % read data into phi 
    fid = fopen(fstr);
    Lx = textscan(fid,'%f',1); Lx = Lx{1};
    Ly = textscan(fid,'%f',1); Ly = Ly{1};
    Lz = textscan(fid,'%f',1); Lz = Lz{1};
    LList(ss,:) = [Lx,Ly,Lz];
    fclose(fid);
    
    frameID(ss) = sscanf(fname,[simfrmt '_%f.pos']);
end
[frameID,srtidx] = sort(frameID);
simlist = simlist(srtidx);
LList = LList(srtidx,:);
% Lplot = LList(end,:);

%% Loop over files, make movie, save structural features

% save
tList = zeros(NSIMS,1);
skList = cell(NSIMS,1);
kbinList = cell(NSIMS,1);
phiImgList = zeros(NSIMS,1);
xiList = zeros(NSIMS,1);
kvarList = zeros(NSIMS,1);
KgList = zeros(NSIMS,1);
k2List = zeros(NSIMS,1);

% loop over frames
for ss = 1:NSIMS
    % file info
    fname = simlist(ss).name;
    fstr = [simloc '/' fname];
    finfo = dir(fstr);
    fsize = finfo.bytes;
    
    % read data into phi 
    fprintf('\n\nFrame = %d, frameID = %d, fname = %s, size = %0.4g MB\n',ss,frameID(ss),fname,fsize/1e6);
    fid = fopen(fstr);
    Lx = textscan(fid,'%f',1); Lx = Lx{1};
    Ly = textscan(fid,'%f',1); Ly = Ly{1};
    Lz = textscan(fid,'%f',1); Lz = Lz{1};
    t = textscan(fid,'%f',1); t = t{1};
    dt = textscan(fid,'%f',1); dt = dt{1};
    psiin = textscan(fid,repmat('%f ',1,Lx),Ly*Lz);
    psiin = cell2mat(psiin);
    psi = zeros(Ly,Lx,Lz);
    last = 1;
    for zz = 1:Lz
        next = last + Ly - 1;
        psi(:,:,zz) = psiin(last:next,:);
        last = next + 1;
    end
    tList(ss) = t;
        
%     % get grid points
%     [PX, PY, PZ] = meshgrid(1:Lx,1:Ly,1:Lz);
% 
%     % isovalue
%     iv = 0;
% 
%     fprintf('** Generating isosurface in frame %d, t = %0.4g...\n',frameID(ss),t);
%     fv = isosurface(PX,PY,PZ,psi,iv);
%     vpos = fv.vertices;
%     finfo = fv.faces;
% 
%     % plot isosurface
%     fprintf('** Plotting isosurface....not showing figure, but writing to %s...\n',mvstr);
%     f = figure('visible','off');
%     clf, hold on, box on;
%     f.Color = [1 1 1];
%     patch('Faces',finfo,'Vertices',vpos,'FaceColor','b','EdgeColor','none');
%     patch(isocaps(PX, PY, PZ, psi, iv), 'FaceColor', 'interp', 'EdgeColor', 'none');
%     view(3);
%     axis vis3d equal;
%     colormap('jet');
%     ax = gca;
%     ax.Clipping = 'off';
%     ax.XLim = [1 Lplot(1)];
%     ax.YLim = [1 Lplot(2)];
%     ax.ZLim = [1 Lplot(3)];
%     camlight; lighting phong
% 
%     % write to video
%     frame = getframe(f);
%     img = frame2im(frame);
%     [imind,cm] = rgb2ind(img,256);
%     if ss == 1
%         imwrite(imind,cm,mvstr,'gif','Loopcount',inf); 
%     else
%         imwrite(imind,cm,mvstr,'gif','WriteMode','append'); 
%     end
    
    % get structural information
    microCT = psi > mean(psi(:));
    phiImgList(ss,1) = mean(microCT,'all');
    
    fprintf('** ff = %d, frameID = %d: Getting correlations, computing FFTs...\n',ss,frameID(ss));
    F_microCT = fftn(microCT);
    F_microCT = fftshift(F_microCT);
    S = (F_microCT.*conj(F_microCT))./(Lx*Ly*Lz);

    % get wavevector values (diff for even and odd)

    % kx
    Fs = (2.0*pi)/Lx;
    if mod(Lx,2) == 0
        kx = Fs*(-Lx/2:(Lx/2 - 1));
    else
        kx = Fs*(-floor(Lx/2):(round(Lx/2) - 1));
    end
    kx0 = find(kx == 0);

    % ky
    Fs = (2.0*pi)/Ly;
    if mod(Ly,2) == 0
        ky = Fs*(-Ly/2:(Ly/2 - 1));
    else
        ky = Fs*(-floor(Ly/2):(round(Ly/2) - 1));
    end
    ky0 = find(ky == 0);

    % kz
    Fs = (2.0*pi)/Lz;
    if mod(Lz,2) == 0
        kz = Fs*(-Lz/2:(Lz/2 - 1));
    else
        kz = Fs*(-floor(Lz/2):(round(Lz/2) - 1));
    end
    kz0 = find(kz == 0);

    % subtract out contribution from k = 0
    S(ky0,kx0,kz0) = 0;


    % compute structure factor by radial average
    fprintf('** ff = %d, frameID = %d: Computing isotropic structure factor...\n',ss,frameID(ss));

    % loop over distances in each image, compute radial correlations
    kmin        = 0;
    kmax        = mean([max(kx) max(ky) max(kz)]);
    dk          = 1.5*(2.0*pi)/mean([Lx Ly Lz]);
    binedges    = (kmin:dk:kmax)';
    nkbins      = length(binedges) - 1;
    leftBins    = binedges(1:end-1);
    rightBins   = binedges(2:end);
    kbin        = 0.5*(leftBins + rightBins);
    sk          = zeros(nkbins,1);

    % compute knorm lattice
    kyx                         = (ky').^2 + kx.^2;
    kzpage                      = reshape(kz,1,1,Lz);
    knormsq                     = kyx + kzpage.^2;
    knorm                       = sqrt(knormsq);

    % clear excess data
    clear('knormsq');
    clear('kzpage');
    clear('kyx');

    % loop over bins
    for bb = 1:nkbins
        % get location of k values in 3D matrix
        ktmpmin = leftBins(bb);
        ktmpmax = rightBins(bb);
        kinds = knorm(:) < ktmpmax & knorm(:) > ktmpmin;

        % take standard mean
        sk(bb) = mean(S(kinds),'all');
    end
    skList{ss} = sk;
    kbinList{ss} = kbin;
    k1tmp = sum(sk.*kbin)/sum(sk);
    k2tmp = sum(sk.*kbin.^2)/sum(sk);
    kvarList(ss) = k2tmp - k1tmp^2;
    xiList(ss) = 2.0*pi/k1tmp;

    fprintf('** ff = %d, frameID = %d: Computing Gyration Tensor...\n',ss,frameID(ss));

    % Gyration tensor (see ctnm limit def on wikipedia)
    G = zeros(3);
    Sint = 0.0;

    % loop over all points in S
    NP = Lx*Ly*Lz;
    nxy = Lx*Ly;
    for pp = 1:NP
        % get indices based on pp
        pm1 = pp-1;
        rowi = mod(pm1,Ly)+1;
        pgei = floor(pm1/nxy)+1;
        coli = mod(floor(pm1/Ly),Lx)+1;

        % get wavevector
        ktmp = [kx(coli), ky(rowi), kz(pgei)];

        % add to gyration tensor
        for ii = 1:3
            for jj = ii:3
                G(ii,jj) = G(ii,jj) + S(pp)*ktmp(ii)*ktmp(jj);
            end
        end

        % add to integral over all S
        Sint = Sint + S(pp);
    end

    % symmetrize
    G(2,1) = G(1,2);
    G(3,2) = G(2,3);
    G(3,1) = G(1,3);

    % normalize
    G = G./Sint;
    
    % get eigenvalues
    lambda = eig(G);
    lambda = sort(lambda);
    
    % anaylyze
    Kgtmp = sqrt(sum(lambda));
    btmp = lambda(3) - 0.5*(lambda(1) + lambda(2));
    ctmp = lambda(2) - lambda(1);
    k2List(ss) = ((btmp^2) + 0.75*(ctmp^2))/Kgtmp^4;
    KgList(ss) = Kgtmp;
end

%% Save in matfile

save(savestr,'tList','LList','simlist','skList','kbinList','phiImgList',...
    'xiList','kvarList','KgList','k2List');

end