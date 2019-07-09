function region = sat_generator(x1,y1,z1,x2,y2,z2,nx,ny,nz)
    region = cell(nx*ny*nz,1);
    x = linspace(x1,x2,nx);
    y = linspace(y1,y2,ny);
    z = linspace(z1,z2,nz);
    [X,Y,Z] = meshgrid(x,y,z);
    for i = 1:nx*ny*nz
       region{i} = [X(i) Y(i) Z(i)]; 
    end
end