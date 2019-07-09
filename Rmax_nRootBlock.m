ncell = 4; % number of cell per block
n_u = 6; % number of cell underneath the surface

nRootBlock = 10:35;
R_max = (5/4).^(nRootBlock * ncell/n_u - 1);
plot(nRootBlock,R_max);

ylabel('R_{max}');
xlabel('nRootBlock');