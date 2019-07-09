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