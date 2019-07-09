%% This script requires MATLAB to be opened from a terminal, otherwise you can not use IDL inside MATLAB & all the system PATH parameters will be lost
% To start MATLAB from terminal, simply open up a terminal ans type the
% following commands:
%  > matlab &
%  > disown %1 
%  You are all set!

%% Calling IDL::animate_data from matlab
clear;

IDL_var = input("IDL variable to be plotted, exemplary input:    'By',    'Rho',    'Uz'\n ");
IDL_cut = input("Which cut would you like? exemplary input:  'x',     'y',    'z'\n ");

! touch animation.pro
f = fopen('animation.pro','w');

if strcmp(IDL_cut,'x')
    fprintf(f,'%s\n','filename=''x=*outs''');
elseif strcmp(IDL_cut,'y')
    fprintf(f,'%s\n','filename=''y=*outs''');
elseif strcmp(IDL_cut,'z')
    fprintf(f,'%s\n','filename=''z=*outs''');
end

fprintf(f,'%s\n','npict = 1');
fprintf(f,'%s\n','dpict = 1');
fprintf(f,'%s\n','npictmax = 2001');
fprintf(f,'%s\n','read_data');

if strcmp(IDL_cut,'y')
    if strcmp(IDL_var,'By')
        fprintf(f,'%s\n','func=''{by}<50>(-50) bx;bz''');
    elseif strcmp(IDL_var,'Rho')
        fprintf(f,'%s\n','func=''{rho}<150 bx;bz''');
    else
        fprintf(f,'%s\n',['func='''  IDL_var ' bx;bz''']);
    end
elseif strcmp(IDL_cut,'z')
    if strcmp(IDL_var,'By')
        fprintf(f,'%s\n','func=''{by}<50>(-50) bx;by''');
    elseif strcmp(IDL_var,'Rho')
        fprintf(f,'%s\n','func=''{rho}<150 bx;by''');
    else
        fprintf(f,'%s\n',['func='''  IDL_var ' bx;by''']);
    end
else
    if strcmp(IDL_var,'By')
        fprintf(f,'%s\n','func=''{by}<50>(-50)''');
    elseif strcmp(IDL_var,'Rho')
        fprintf(f,'%s\n','func=''{rho}<150''');
    else
        fprintf(f,'%s\n',['func='''  IDL_var '''']);
    end
end

screen=get(0,'ScreenSize');
W=screen(3);H=screen(4);
w = 0.6*W;              % width of the idl window
h = 0.8*H;              % hight of the idl window

fprintf(f,'window,xsize=%i,ysize=%i\n',w,h);
fprintf(f,'%s\n','plotmode = ''contbar streamoverbody''');

x_range = [-6,3];
y_range = [-3,3];

rangeX = ['!x.range = [' num2str(x_range(1)) ',' num2str(x_range(2)) ']'];
rangeY = ['!y.range = [' num2str(y_range(1)) ',' num2str(y_range(2)) ']'];
fprintf(f,'%s\n',rangeX);
fprintf(f,'%s\n',rangeY);
fprintf(f,'%s\n','savemovie=''png''');
fprintf(f,'%s\n','animate_data');
fprintf(f,'%s\n','exit');

fclose(f);
! idl animation.pro
! rm -f animation.pro

Dir_name = ['Movie_' IDL_var '_' IDL_cut '=0'];

! touch tmp.txt
tmp = fopen('tmp.txt','w');
fprintf(tmp,'%s',Dir_name);
fclose(tmp);

PDF_name = [IDL_var '_' IDL_cut '=0.pdf'];
! touch pdf.txt
pdf = fopen('pdf.txt','w');
fprintf(pdf,'%s',PDF_name);
fclose(pdf);

! rm -rf $(cat tmp.txt)
! mv Movie $(cat tmp.txt) 
! cd $(cat tmp.txt); /opt/local/bin/convert *.png $(cat ../pdf.txt)

% This segment is copied from png2avi.m

writerObj = VideoWriter([Dir_name '/' IDL_var '_' IDL_cut '=0.avi']);
writerObj.FrameRate=10;
open(writerObj);

! cd $(cat tmp.txt); ls *.png | wc -l > ../number.txt
num = fopen('number.txt');
N_png = fscanf(f,'%i');
fclose(num);
! rm -f number.txt;
! rm -f tmp.txt
! rm -f pdf.txt

for K = 1 : N_png
    filename = sprintf('%4.4d.png', K);
    thisimage = imread([Dir_name '/' filename]);
    writeVideo(writerObj, thisimage);
end
close(writerObj);


%% Loading .png files, then find the location of satellite

A = imread([Dir_name '/0001.png']);
size_A = size(A);

%------------------------------
prompt = "Would you like to compare satellites at different locations?(yes/no)";
flag = input(prompt,'s');
if strcmp(flag,"yes")
    N = input("How many satellites would you like?      ");    
else
    N = 1;
end

clear prompt;

% Input the locations of satellites

Locations = cell(N,1);
for i = 1:N
    if i == 1
       disp("Exemplary format of input:    [x  y z]      (press ""return"" to proceed)"); 
       pause;
    end
    
    if N == 1
        disp("Plaese enter the location of satellite");
        Locations{i} = input("");
    else
        disp(strcat("Please enter the location of satellite " , num2str(i)));
        Locations{i} = input("");
    end
end
%------------------------------

color = double(A(:,:,1)) + double(A(:,:,2)) + double(A(:,:,3));
color_max = 2^8 -1; % Unit 8
% color_max = 2^16 - 1; % Unit 16

color = color ~= 0 & color ~= 3*color_max;
index = find(color);

corner_1 = min(index);
corner_2 = max(index);

column_1 = ceil(corner_1 / size_A(1));
row_1 = corner_1 - size_A(1)*(column_1 - 1);
% A(row_1-5:row_1+5,column_1-5:column_1+5,1) = color_max;
% A(row_1-5:row_1+5,column_1-5:column_1+5,2) = 0;
% A(row_1-5:row_1+5,column_1-5:column_1+5,3) = color_max;

column_2 = ceil(corner_2 / size_A(1));
row_2 = corner_2 - size_A(1)*(column_2 - 1);

for i = column_2:-1:1
   if color(row_2,i) == 0
       column_2 = i;
       break; 
   end
end

for i = column_2:-1:1
   if color(row_2,i) ~= 0
       column_2 = i;
       break; 
   end
end

% A(row_2-5:row_2+5,column_2-5:column_2+5,1) = color_max;
% A(row_2-5:row_2+5,column_2-5:column_2+5,2) = 0;
% A(row_2-5:row_2+5,column_2-5:column_2+5,3) = color_max;

clear A;
! rm -rf My_movie_1
! mkdir My_movie_1

for K = 1 : N_png
    filename = sprintf('%4.4d.png', K);
    A = imread([Dir_name '/' filename]);

    sat_size = 10;

    for i = 1:N
        if strcmp(IDL_cut,'x')
	        sat_row = row_1 + (row_2 - row_1)*(y_range(2) - Locations{i}(3)/100)/(y_range(2) - y_range(1));
            sat_column = column_1 + (column_2 - column_1)*(Locations{i}(2)/100 - x_range(1))/(x_range(2) - x_range(1));
            sat_row = round(sat_row);
            sat_column = round(sat_column);
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,1) = color_max;
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,2) = 0;
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,3) = color_max;
        elseif strcmp(IDL_cut,'y')
	        sat_row = row_1 + (row_2 - row_1)*(y_range(2) - Locations{i}(3)/100)/(y_range(2) - y_range(1));
            sat_column = column_1 + (column_2 - column_1)*(Locations{i}(1)/100 - x_range(1))/(x_range(2) - x_range(1));
            sat_row = round(sat_row);
            sat_column = round(sat_column);
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,1) = color_max;
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,2) = 0;
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,3) = color_max;
        elseif strcmp(IDL_cut,'z')
	        sat_row = row_1 + (row_2 - row_1)*(y_range(2) - Locations{i}(2)/100)/(y_range(2) - y_range(1));
            sat_column = column_1 + (column_2 - column_1)*(Locations{i}(1)/100 - x_range(1))/(x_range(2) - x_range(1));
            sat_row = round(sat_row);
            sat_column = round(sat_column);
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,1) = color_max;
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,2) = 0;
            A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,3) = color_max;
        end
    end

    figure('Position',[0,0,w,h]);
    image(A);

    % This block of code is used to eliminate the white space in MATLAB figure
    %----------------------------------------------------
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    %----------------------------------------------------
    
    saveas(gcf,['My_movie_1/' filename]);
    close gcf;
end








