% Changkun Li, June 4, 2019
% Visualization of FTE on Magnetopause

%% Read in the names of all 3D files
clear;
! ls 3d*.dat > tmp.txt
! wc -l tmp.txt > numline.txt
f = fopen('numline.txt');
N_output = fscanf(f,'%i');
fclose(f);
ThreeD_namelist = cell(N_output,1);
f = fopen('tmp.txt');
for i = 1:N_output
   ThreeD_namelist{i} = fgetl(f); 
end
fclose(f);
!rm -f tmp.txt
!rm -f numline.txt

%% Plot out magnetopause for each 3D file and save them into .png format
% png2avi.m can be used to generate avi from png,
% /opt/local/bin/convert(part of ImageMagick package) can be used to
% generate pdf from png.

screen=get(0,'ScreenSize');
W=screen(3);H=screen(4);    % W here will be overwritten by z-component of magnetopause normal: [V, W] = differentiate(fitresult, yq, zq);
w = 0.6*W;                  % width of the idl window
h = 0.8*H;                  % hight of the idl window

! rm -rf Two_D_magnetopause
! mkdir Two_D_magnetopause

for i = 1:N_output
    [data, connect, header, N_ele, N_con, time, VAR] = read_3d(ThreeD_namelist{i});
    valueSet = 1:max(size(VAR)); 
    Dict = containers.Map(VAR,valueSet); % Similar to python dictionary

    % Courtesy to Hongyang (https://github.com/henry2004y/)
    % Start-------------------------------------------
    x = data(:,Dict('X [R]'));
    y = data(:,Dict('Y [R]'));
    z = data(:,Dict('Z [R]'));
    status = data(:,Dict('Status'));
    
    x3_1 = x(status==3 & x>0);
    y3_1 = y(status==3 & x>0);
    z3_1 = z(status==3 & x>0);
    
    % Pick the compact boundary points coords.
    k31  = boundary(x3_1,y3_1,z3_1,1);
    bc31 = unique(k31);

    x3bc = x3_1(bc31);
    y3bc = y3_1(bc31);
    z3bc = z3_1(bc31);

    % Find the outer boundary points
    rThres = 1.1;
    xThres = 0.7;
    
    mapindex_ = (x3bc.^2 + y3bc.^2 + z3bc.^2) > rThres^2 & x3bc > xThres;
    x3bc = x3bc(mapindex_);
    y3bc = y3bc(mapindex_);
    z3bc = z3bc(mapindex_);
    
    figure(1); scatter3(x3bc,y3bc,z3bc,'.'); axis equal
    
    % Fit the closed field line boundary with hypersurface fit

    [fitresult,gof] = surface_fit(x3bc,y3bc,z3bc);

    % Generate mesh points from fitted surface
    ymin = -1.1+1/15; ymax = 1.1-1/15; zmin = -0.7+1/15; zmax = 1-1/15;
    dy = 1/30; dz = dy;
    [yq,zq] = ndgrid(ymin:dy:ymax,zmin:dz:zmax);

    xq = fitresult(yq,zq);

    % Transform into local coordinate system

    % Calculate the normal direction to the fitted surface
    [V, W] = differentiate(fitresult, yq, zq);

    U = -ones(size(V));

    % [U,V,W]
    figure(1); hold on
    quiver3(xq,yq,zq,U,V,W,10,'color','r')
    hold off

    % get the three local directions
    % dipole-direction unit vector
    unitDipole = [0 0 1];
    % Initialize local vectors
    dL = Inf(3,size(xq,1),size(xq,2));
    dM = dL; dN = dL;

    % This part could potentially be optimized!
    for ix=1:size(xq,1)
        for iy=1:size(xq,2)
            dN(:,ix,iy) = -[U(ix,iy) V(ix,iy) W(ix,iy)];
            dM(:,ix,iy) = cross(dN(:,ix,iy),unitDipole);
            dL(:,ix,iy) = cross(dM(:,ix,iy),dN(:,ix,iy));
      
            % Normalization
            dL(:,ix,iy) = dL(:,ix,iy) / norm(dL(:,ix,iy));
            dM(:,ix,iy) = dM(:,ix,iy) / norm(dM(:,ix,iy));
            dN(:,ix,iy) = dN(:,ix,iy) / norm(dN(:,ix,iy));
        end
    end
    
    bx = data(:,Dict('B_x [nT]'));
    by = data(:,Dict('B_y [nT]'));
    bz = data(:,Dict('B_z [nT]'));
    
    Coef = 1.1;  % Hongyang's default, should change this value later
    bxv= griddata(x, y, z, bx, Coef*xq, Coef*yq, Coef*zq);
    byv= griddata(x, y, z, by, Coef*xq, Coef*yq, Coef*zq); 
    bzv= griddata(x, y, z, bz, Coef*xq, Coef*yq, Coef*zq);
    
    % Transform vectors into local coordinate system
    bL = Inf(size(xq)); bM = bL; bN = bL;
    
    % This could potentially be improved!
    for ix=1:size(xq,1)
        for iy=1:size(xq,2)
            bL(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dL(:,ix,iy);
            bM(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dM(:,ix,iy);
            bN(ix,iy) = [bxv(ix,iy) byv(ix,iy) bzv(ix,iy)]*dN(:,ix,iy);
        end
    end
    
    handle = figure('Position',[0,0,w,h]);
    contourf(yq,zq,bN,50,'Linestyle','none'); colorbar;
    axis tight equal
    set(gca,'Xdir','reverse')
    caxis([-40 40])
    xlabel('y [R_G]'); ylabel('z [R_G]');
    title(['B_N [nT], ' time]);
    
    filename = sprintf('%4.4d.png', i);
    saveas(handle,['Two_D_magnetopause/' filename]);
    close(handle);
    % End--------------------------------------------
    
end

! cd Two_D_magnetopause; /opt/local/bin/convert *.png B_N.pdf

writerObj = VideoWriter('Two_D_magnetopause/B_N.avi');
writerObj.FrameRate=4;
open(writerObj);

for K = 1 : N_output
  filename = sprintf('%4.4d.png', K);
  thisimage = imread(['Two_D_magnetopause/' filename]);
  writeVideo(writerObj, thisimage);
end
close(writerObj);

%% Function definition

function [data, connect, header, N_ele, N_con, time, VAR] = read_3d(filename)
    f = fopen(filename);
    header = cell(23,1);
    for j = 1:23 % number of lines of file header
        header{j} = fgetl(f);
    end
    k = find(header{22} == '"');
    time = header{22}(k(1)+1:k(2)-1);
    a = strfind(header{3},'N= ');
    b = strfind(header{3},'E= ');
    c = strfind(header{3},'F=');
    N_ele = str2double(header{3}(a+2:b-2)); % number of elements
    N_con = str2double(header{3}(b+2:c-2)); % number of connectivity list
    data = fscanf(f,'%f');
    fclose(f);
    
    total = max(size(data));
    connect = data((total - N_con*8+1):total);
    connect = reshape(connect,[8,N_con]);
    connect = connect';
    data = data(1:(total - N_con*8));
    data = reshape(data,[(total - N_con*8)/N_ele,N_ele]);
    data = data';
    
    k = find(header{2} == '"');
    VAR = cell(max(size(k))/2,1);
    for i = 1:max(size(k))/2
        VAR{i} = header{2}((k(2*i - 1)+1):(k(2*i)-1));
    end
end

function [fitresult,gof] = surface_fit(x3bc,y3bc,z3bc,varargin)
    %SURFACE_FIT Fit the closed field line boundary with hypersurface fit
    %
    % INPUT:
    % x3bc,y3bc,z3bc: 1D array of coordinates
    % TypeFit: fit model type (see MATLAB document)
    % DoPlot : logical, display figure or not
    %
    % OUTPUT:
    % fitresult: fit model
    % gof: goodess of fit
    %
    % Hongyang Zhou, hyzhou@umich.edu 06/29/2018

    if nargin<3
       error('Not enough input arguments.')
    end

    optargs = {'poly55' false}; % default parameters
    optargs(1:nargin-3) = varargin;
    [TypeFit,DoPlot] = optargs{:};


    % Set up fittype and options.
    ft = fittype( TypeFit );

    % Fit model to data.
    [fitresult, gof] = fit( [y3bc, z3bc], x3bc, ft );

    if DoPlot
       % Plot fit with data.
       figure(1); %hold on
       h = plot( fitresult );
       legend( h, 'poly5, x=x(y,z)', 'Location', 'NorthEast' );
       % Label axes
       xlabel('x [R_G]')
       ylabel('y [R_G]')
       zlabel('z [R_G]')
       grid on; axis equal

       xx = get(h, 'XData');
       yy = get(h, 'YData');
       zz = get(h, 'Zdata');
       set(h, 'XData', zz, 'YData', xx, 'ZData', yy);

       hold on;
       scatter3(x3bc,y3bc,z3bc,20,'r','filled'); hold off
       %axis tight
       xlim([-2 0]); ylim([-1.5 1.5]); zlim([-0.6 0.8]);
       set(gca,'FontSize',14,'LineWidth',1.2)
    end

end



