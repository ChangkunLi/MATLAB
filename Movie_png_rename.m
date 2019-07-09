% This script is designed to solve the problem rising from multiple IDL
% animate _data, for instance, you want 0001-0501.png, but if you have to
% stop IDL and re-animate_data, you will probably get
% 0001-0230.png,0001-0271.png, the goal of this script is to change the
% name of 0001-0271.png to 231-501.png

displacement = 252;
for K = 1:249
    filename_old = sprintf('%4.4d.png', K);
    filename_new = sprintf('%4.4d.png', K+displacement);
    movefile(filename_old,filename_new);
end