clear;

global IDL_var IDL_cut x_range y_range w h 
global IDL_var_g IDL_cut_g 

screen = get(0,'ScreenSize');
W_g = screen(3); H_g = screen(4);
w_g = 0.6*W_g;              % width of the GUI window
h_g = 0.8*H_g;              % hight of the GUI window

figure('Color',[1 1 1],'Position',[0,0,w_g,h_g],'Name','User Input GUI','NumberTitle','off','MenuBar','none');

IDL_var_g = {'Rho' 'Ux' 'Uy' 'Uz' 'b1x' 'b1y' 'b1z' 'Bx' 'By' 'Bz' 'p' 'eta' 'jx' 'jy' 'jz' 'dt' 'dtblk' 'cons' 'impl' 'dx' 'hall' 'absdivb' 'pe' 'pic'};
t1 = uicontrol(gcf,'Style','text','String','IDL variable to be plotted','Position',[0.01*w_g 0.925*h_g 0.2*w_g 0.05*h_g]);
t1.FontSize = 25;
h1 = uicontrol(gcf,'style','popupmenu','Position',[0.01*w_g 0.7*h_g 0.2*w_g 0.2*h_g],'string', IDL_var_g,'callback',@h1_callback);
h1.FontSize = 20;

IDL_cut_g = {'x' 'y' 'z'};
t2 = uicontrol(gcf,'Style','text','String','IDL cut','Position',[0.25*w_g 0.925*h_g 0.1*w_g 0.05*h_g]);
t2.FontSize = 25;
h2 = uicontrol(gcf,'style','popupmenu','Position',[0.25*w_g 0.7*h_g 0.1*w_g 0.2*h_g],'string', IDL_cut_g,'callback',@h2_callback);
h2.FontSize = 20;

t3 = uicontrol(gcf,'Style','text','String','!x.range','Position',[0.4*w_g 0.925*h_g 0.1*w_g 0.05*h_g]);
t3.FontSize = 25;
h3 = uicontrol(gcf,'style','edit','string',' ','Position',[0.4*w_g 0.7*h_g 0.1*w_g 0.2*h_g],'callback',@h3_callback);
h3.FontSize = 20;

t4 = uicontrol(gcf,'Style','text','String','!y.range','Position',[0.55*w_g 0.925*h_g 0.1*w_g 0.05*h_g]);
t4.FontSize = 25;
h4 = uicontrol(gcf,'style','edit','string',' ','Position',[0.55*w_g 0.7*h_g 0.1*w_g 0.2*h_g],'callback',@h4_callback);
h4.FontSize = 20;

t5 = uicontrol(gcf,'Style','text','String','Width of IDL plot: 0~1','Position',[0.7*w_g 0.925*h_g 0.1*w_g 0.06*h_g]);
t5.FontSize = 25;
h5 = uicontrol(gcf,'style','edit','string',' ','Position',[0.7*w_g 0.7*h_g 0.1*w_g 0.2*h_g],'callback',@h5_callback);
h5.FontSize = 20;

t6 = uicontrol(gcf,'Style','text','String','Height of IDL plot: 0~1','Position',[0.85*w_g 0.925*h_g 0.1*w_g 0.06*h_g]);
t6.FontSize = 25;
h6 = uicontrol(gcf,'style','edit','string',' ','Position',[0.85*w_g 0.7*h_g 0.1*w_g 0.2*h_g],'callback',@h6_callback);
h6.FontSize = 20;

h15 = uicontrol(gcf,'style','push','string','Plot','position',[0.01*w_g 0.1*h_g 0.95*w_g 0.5*h_g],'callback',@h15_callback);
h15.FontSize = 100;

% Default setting
IDL_var = 'Rho';
IDL_cut = 'y';
x_range = [-6,3];
y_range = [-3,3];
w = 0.6*W_g;
h = 0.8*H_g;

function h1_callback(hObject,~)
    global IDL_var IDL_var_g
    
    m = get(hObject,'value');
    IDL_var = IDL_var_g{m};
end
    
function h2_callback(hObject,~)
    global IDL_cut IDL_cut_g
    
    m = get(hObject,'value');
    IDL_cut = IDL_cut_g{m};
end
    
function h3_callback(hObject,~)
    global x_range

    x_range = [-6,3];
    m = get(hObject,'string');
    eval(['x_range = ' m ';']);
end

function h4_callback(hObject,~)
    global y_range

    y_range = [-3,3];
    m = get(hObject,'string');
    eval(['y_range = ' m ';']);
end

function h5_callback(hObject,~)
    global w
    
    screen = get(0,'ScreenSize');
    W = screen(3);
    m = get(hObject,'string');
    w = str2double(m) * W;
end

function h6_callback(hObject,~)
    global h
    
    screen = get(0,'ScreenSize');
    H = screen(4);
    m = get(hObject,'string');
    h = str2double(m) * H;
end

function h15_callback(~,~)
    global  IDL_var IDL_cut x_range y_range w h
    
    %% Calling IDL::animate_data from matlab
%     clear;
% 
%     IDL_var = input("IDL variable to be plotted, exemplary input:    'By',    'Rho',    'Uz'\n ");
%     IDL_cut = input("Which cut would you like? exemplary input:  'x',     'y',    'z'\n ");

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
    fprintf(f,'%s\n','npictmax = 1');
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

%     w = 0.6*W;              % width of the idl window
%     h = 0.8*H;              % hight of the idl window

    fprintf(f,'window,xsize=%i,ysize=%i\n',w,h);
    fprintf(f,'%s\n','plotmode = ''contbar streamoverbody''');

%     x_range = [-6,3];
%     y_range = [-3,3];

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

    Dir_name = ['SatPosition_' IDL_var '_' IDL_cut '=0'];

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
%     prompt = "Would you like to compare satellites at different locations?(yes/no)";
%     flag = input(prompt,'s');
%     if strcmp(flag,"yes")
%         N = input("How many satellites would you like?      ");    
%     else
%         N = 1;
%     end
% 
%     clear prompt;
% 
%     % Input the locations of satellites
% 
%     Locations = cell(N,1);
%     for i = 1:N
%         if i == 1
%            disp("Exemplary format of input:    [x  y z]      (press ""return"" to proceed)"); 
%            pause;
%         end
% 
%         if N == 1
%             disp("Plaese enter the location of satellite");
%             Locations{i} = input("");
%         else
%             disp(strcat("Please enter the location of satellite " , num2str(i)));
%             Locations{i} = input("");
%         end
%     end

    region_1 = sat_generator(1.1,-1.5,-1.5,2.3,1.5,1.5,7,5,5);
    region_2 = sat_generator(0,-0.8,0.9,1,0.8,2,4,5,4);
    region_3 = sat_generator(0.2,-0.8,-0.9,1.2,0.8,-2,4,5,4);
    region_4 = sat_generator(1,0,-0.8,1,0,-0.8,1,1,1);
    region_5 = sat_generator(0.9,0,0,-0.9,0,0,2,1,1);
    region_6 = sat_generator(0,0.9,0,0,-0.9,0,1,2,1);
    region_7 = sat_generator(0,0,0.9,0,0,-0.9,1,1,2);
    region_8 = sat_generator(1.2,-1,-0.8,1.6,1,1.2,9,5,5);
    region_9 = sat_generator(1.15,-1,1.2,0,1,1.2,24,5,1);
    region_10= sat_generator(1.15,-1,-0.8,0,1,-0.8,24,5,1); 
    
    Locations = [region_1;region_2;region_3;region_4;region_5;region_6;region_7;region_8;region_9;region_10];
    N = max(size(Locations));
    
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
    ! rm -rf SatPosition
    ! mkdir SatPosition

    for K = 1 : N_png
        filename = sprintf('%4.4d.png', K);
        A = imread([Dir_name '/' filename]);

        sat_size = 3;

        for i = 1:N
            if strcmp(IDL_cut,'x')
                sat_row = row_1 + (row_2 - row_1)*(y_range(2) - Locations{i}(3))/(y_range(2) - y_range(1));
                sat_column = column_1 + (column_2 - column_1)*(Locations{i}(2) - x_range(1))/(x_range(2) - x_range(1));
                sat_row = round(sat_row);
                sat_column = round(sat_column);
                A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,1) = color_max;
                A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,2) = 0;
                A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,3) = color_max;
            elseif strcmp(IDL_cut,'y')
                sat_row = row_1 + (row_2 - row_1)*(y_range(2) - Locations{i}(3))/(y_range(2) - y_range(1));
                sat_column = column_1 + (column_2 - column_1)*(Locations{i}(1) - x_range(1))/(x_range(2) - x_range(1));
                sat_row = round(sat_row);
                sat_column = round(sat_column);
                A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,1) = color_max;
                A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,2) = 0;
                A(sat_row-sat_size:sat_row+sat_size,sat_column-sat_size:sat_column+sat_size,3) = color_max;
            elseif strcmp(IDL_cut,'z')
                sat_row = row_1 + (row_2 - row_1)*(y_range(2) - Locations{i}(2))/(y_range(2) - y_range(1));
                sat_column = column_1 + (column_2 - column_1)*(Locations{i}(1) - x_range(1))/(x_range(2) - x_range(1));
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

        saveas(gcf,['SatPosition/' filename]);
        close gcf;
    end

end