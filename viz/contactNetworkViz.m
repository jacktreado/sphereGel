%% Plot contact network during sphere gel sim

clear;
close all;
clc;

% simulation info
Nstr        = '400';
drstr       = '0.15';
dphistr     = '0.0005';
dlzstr      = '1.1';
l2str       = '0.2';
seedstr     = '1';

Ninput      = str2double(Nstr);
dr          = str2double(drstr);
dphi        = str2double(dphistr);
dlz         = str2double(dlzstr);
l2          = str2double(l2str);
seed        = str2double(seedstr);

% simulation string
simstr      = ['sgel_N' Nstr '_dr' drstr '_dphi' dphistr '_dlz' dlzstr '_l2' l2str '_seed' seedstr];

% file strings
fstr        = ['xyz/' simstr '.xyz'];
finfo       = dir(fstr);
fsz         = finfo.bytes;
fname       = finfo.name;
cmfstr      = ['xyz/' simstr '.cm'];

% print file info to console
if fsz == 0
    fprintf('-- Choosen file empty, ending with error.\n');
    error('fstr empty');
else
    fprintf('-- Reading in file %s, size = %0.3f MB\n',fname,fsz/1e6);
end

% parse data from xyz file
posdata     = readSGelXYZ(fstr);

N           = posdata.N;
NFRAMES     = posdata.NFRAMES;
xpos        = posdata.xpos;
ypos        = posdata.ypos;
zpos        = posdata.zpos;
radii       = posdata.radii;
L           = posdata.L;

% get packing fraction over time
pvols       = (4.0/3.0)*pi*radii.^3;
boxvols     = L(:,1).*L(:,2).*L(:,3);
pvsum       = sum(pvols,2);
phi         = pvsum./boxvols;


% option to make a movie
makeAMovie = 1;


%% get contact matrix, plot number of contacts over time

% open cm
fid = fopen(cmfstr);

% contact matrix
NPAIRS = 0.5*N*(N-1);
frmt = repmat('%f ',1,NPAIRS);
data = textscan(fid,frmt,NFRAMES);
cm = cell2mat(data);

fclose(fid);

%% Plots

% number of contacts over time
Nc = sum(cm,2);

figure(1), clf, hold on, box on;
plot(phi(2:end),Nc,'ko-','linewidth',2);
xlabel('$\phi$','Interpreter','latex');
ylabel('$N_c$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;


% contacts / particle over time

z = 2.0*Nc./N;

figure(2), clf, hold on, box on;
plot(phi(2:end),z,'ko-','linewidth',2);
xlabel('$\phi$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;



%% Dry contacts from simulation

% movie information
if makeAMovie == 1
    mvstr = ['movies/' simstr '.mp4'];
    vobj = VideoWriter(mvstr,'MPEG-4');
    vobj.FrameRate = 10;
    open(vobj);
end

FSTART = 2;
FEND = NFRAMES;
FSTEP = 1;

for ff = FSTART:FSTEP:FEND
    % open figure window
    figure(10), clf, hold on, box on
    
    % print to console
    fprintf('viewing cm for ff = %d\n',ff);
    
    Lx = L(ff,1);
    Ly = L(ff,2);
    Lz = L(ff,3);
    
    % loop over contacts, draw line
    cc = 1;
    for nn = 1:N
        for mm = nn+1:N
            if cm(ff-1,cc) == 1
                dz = zpos(ff,mm) - zpos(ff,nn);
                dz = dz - Lz*round(dz/Lz);
                    
                dy = ypos(ff,mm) - ypos(ff,nn);
                dy = dy - Ly*round(dy/Ly);
                        
                dx = xpos(ff,mm) - xpos(ff,nn);
                dx = dx - Lx*round(dx/Lx);
                        
                % draw line between contacts
                xl = [xpos(ff,nn), xpos(ff,nn) + dx];
                yl = [ypos(ff,nn), ypos(ff,nn) + dy];
                zl = [zpos(ff,nn), zpos(ff,nn) + dz];
                line(xl,yl,zl,'linewidth',1.5,'color','r');
            end
            cc = cc + 1;
        end
        plot3(xpos(ff,nn),ypos(ff,nn),zpos(ff,nn),'bo','markersize',4,'MarkerFaceColor','b');
    end
    
    % draw boundaries
    plot3([0 Lx Lx 0 0], [0 0 Ly Ly 0], [0 0 0 0 0], 'k-', 'linewidth', 1.4);
    plot3([0 Lx Lx 0 0], [0 0 Ly Ly 0], [Lz Lz Lz Lz Lz], 'k-', 'linewidth', 1.4);
    plot3([0 0 Lx Lx], [0 0 0 0], [0 Lz Lz 0], 'k-', 'linewidth', 1.4);
    plot3([0 0 Lx Lx], [Ly Ly Ly Ly], [0 Lz Lz 0], 'k-', 'linewidth', 1.4);
    axis('equal');
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.ZTick = [];
    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    ax.FontSize = 22;
    
    ax.XLim = Lx*[-0.25 1.25];
    ax.YLim = Ly*[-0.25 1.25];
    ax.ZLim = Lz*[-0.25 1.25];
    
    view([45 36]);
    
    % write to movie file
    if makeAMovie == 1
        drawnow;
        currframe = getframe(gcf);
        writeVideo(vobj,currframe);
    end
end

% close video file
if makeAMovie == 1
    close(vobj);
end
