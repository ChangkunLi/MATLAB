clear;

global IDL_var IDL_cut x_range y_range w h Sat_namelist Sat_VAR
global IDL_var_g IDL_cut_g Sat_namelist_g Sat_VAR_g h7 h9 h11 h13

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

% Read in the names of all satellite files
! ls *.sat > tmp.txt
! wc -l tmp.txt > numline.txt
f = fopen('numline.txt');
N_output = fscanf(f,'%i');
fclose(f);
Sat_namelist_g = cell(N_output,1);
f = fopen('tmp.txt');
for i = 1:N_output
   Sat_namelist_g{i} = fgetl(f); 
end
fclose(f);
!rm -f tmp.txt
!rm -f numline.txt
Sat_namelist = {' '};

h7 = uicontrol(gcf,'style','popupmenu','Position',[0.01*w_g 0.4*h_g 0.3*w_g 0.19*h_g],'string',Sat_namelist_g);
h7.FontSize = 20;

h8 = uicontrol(gcf,'style','push','string','Add satellite','position',[0.01*w_g 0.6*h_g 0.14*w_g 0.09*h_g],'callback',@h8_callback);
h8.FontSize = 20;

h9 = uicontrol(gcf,'style','listbox','string',Sat_namelist,'position',[0.01*w_g 0.05*h_g 0.3*w_g 0.49*h_g]);
h9.FontSize = 20;

h10 = uicontrol(gcf,'style','push','string','Remove satellite','position',[0.16*w_g 0.6*h_g 0.14*w_g 0.09*h_g],'callback',@h10_callback);
h10.FontSize = 20;

% !!!!!! String array is preferred, because plot_satellite.m uses string array rather than cell array

% Sat_VAR_g = {'it' 'year' 'mo' 'dy' 'hr' 'mn' 'sc' 'msc' 't' 'X' 'Y' 'Z' 'Rho' 'Ux' 'Uy' 'Uz' 'Bx' 'By' 'Bz' 'Hyp' 'Pe' 'P' 'jx' 'jy' 'jz'};
% Sat_VAR = {' '};

Sat_VAR_g = ["it" "year" "mo" "dy" "hr" "mn" "sc" "msc" "t" "X" "Y" "Z" "Rho" "Ux" "Uy" "Uz" "Bx" "By" "Bz" "Hyp" "Pe" "P" "jx" "jy" "jz" "Bm" "Um" "jm"];
Sat_VAR = " ";

h11 = uicontrol(gcf,'style','popupmenu','Position',[0.35*w_g 0.4*h_g 0.3*w_g 0.19*h_g],'string',Sat_VAR_g);
h11.FontSize = 20;

h12 = uicontrol(gcf,'style','push','string','Add satellite VAR','position',[0.35*w_g 0.6*h_g 0.14*w_g 0.09*h_g],'callback',@h12_callback);
h12.FontSize = 20;

h13 = uicontrol(gcf,'style','listbox','string',Sat_VAR,'position',[0.35*w_g 0.05*h_g 0.3*w_g 0.49*h_g]);
h13.FontSize = 20;

h14 = uicontrol(gcf,'style','push','string','Remove satellite VAR','position',[0.51*w_g 0.6*h_g 0.14*w_g 0.09*h_g],'callback',@h14_callback);
h14.FontSize = 20;

h15 = uicontrol(gcf,'style','push','string','Plot','position',[0.7*w_g 0.1*h_g 0.25*w_g 0.5*h_g],'callback',@h15_callback);
h15.FontSize = 100;

% Default setting
IDL_var = 'Rho';
IDL_cut = 'y';
x_range = [-6,3];
y_range = [-3,3];
w = 0.6*W_g;
h = 0.8*H_g;
Sat_namelist = {' '};

% Sat_VAR = {' '};
Sat_VAR = " ";

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

function h8_callback(~,~)
    global h7 h9 Sat_namelist Sat_namelist_g
    
    Sat_namelist = [Sat_namelist Sat_namelist_g{h7.Value}];
    h9.String = Sat_namelist;
end

function h10_callback(~,~)
    global h9 Sat_namelist
    
    Size = size(Sat_namelist);
    End = Size(2);
    if h9.Value == End
        Sat_namelist = Sat_namelist(1:(h9.Value-1));
        h9.Value = h9.Value - 1;
    else
        Sat_namelist = Sat_namelist([1:(h9.Value-1) (h9.Value + 1):end]);
    end
    h9.String = Sat_namelist;
end

% function h12_callback(~,~)
%     global h11 h13 Sat_VAR Sat_VAR_g
%     
%     Sat_VAR = [Sat_VAR Sat_VAR_g{h11.Value}];
%     h13.String = Sat_VAR;
% end
% 
% function h14_callback(~,~)
%     global h13 Sat_VAR
%     
%     Size = size(Sat_VAR);
%     End = Size(2);
%     if h13.Value == End
%         Sat_VAR = Sat_VAR(1:(h13.Value-1));
%         h13.Value = h13.Value - 1;
%     else
%         Sat_VAR = Sat_VAR([1:(h13.Value-1) (h13.Value + 1):end]);
%     end
%     h13.String = Sat_VAR;
% end

function h12_callback(~,~)
    global h11 h13 Sat_VAR Sat_VAR_g
    
    Sat_VAR = [Sat_VAR Sat_VAR_g(h11.Value)];
    h13.String = Sat_VAR;
end

function h14_callback(~,~)
    global h13 Sat_VAR
    
    Size = size(Sat_VAR);
    End = Size(2);
    if h13.Value == End
        Sat_VAR = Sat_VAR(1:(h13.Value-1));
        h13.Value = h13.Value - 1;
    else
        Sat_VAR = Sat_VAR([1:(h13.Value-1) (h13.Value + 1):end]);
    end
    h13.String = Sat_VAR;
end

function h15_callback(~,~)
    global  IDL_var IDL_cut x_range y_range w h Sat_namelist Sat_VAR 
    
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
    W=screen(3);
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

    Size = size(Sat_namelist);
    N = Size(2) - 1;
    
    Locations = cell(N,1);
    for i = 1:N
        Locations{i} = location_extract(Sat_namelist{i+1});
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
    
    %% Read & Construct data

    DATA = cell(N,1);
    for i = 1:N
        DATA{i} = read_sat(Locations{i}(1),Locations{i}(2),Locations{i}(3));
    end

    %% Plot out time series for selected satellites

%     disp("Please enter the physical variable(s)");
%     disp("Exemplary format of input:    [""Bx"" ""Rho""  ""P""]      (press ""return"" to proceed)"); 
%     pause;
%     disp("Enter your choice below");
%     VAR_str = input("");
    
    VAR_str = Sat_VAR(2:end);
    keySet = {'it' 'year' 'mo' 'dy' 'hr' 'mn' 'sc' 'msc' 't' 'X' 'Y' 'Z' 'Rho' 'Ux' 'Uy' 'Uz' 'Bx' 'By' 'Bz' 'Hyp' 'Pe' 'P' 'jx' 'jy' 'jz' 'Bm' 'Um' 'jm'};
    valueSet = 1:max(size(keySet));
    Dict = containers.Map(keySet,valueSet); % Similar to python dictionary

    for i = 1:N
        DATA{i} = [DATA{i};(DATA{i}(Dict('Bx'),:).^2 + DATA{i}(Dict('By'),:).^2 + DATA{i}(Dict('Bz'),:).^2).^0.5]; % Adding Bm
        DATA{i} = [DATA{i};(DATA{i}(Dict('Ux'),:).^2 + DATA{i}(Dict('Uy'),:).^2 + DATA{i}(Dict('Uz'),:).^2).^0.5]; % Adding Um
        DATA{i} = [DATA{i};(DATA{i}(Dict('jx'),:).^2 + DATA{i}(Dict('jy'),:).^2 + DATA{i}(Dict('jz'),:).^2).^0.5]; % Adding Jm
    end
    
    unitSet = {'step','year','month','day','h','min','sec','msec','sec','R_{planet}','R_{planet}','R_{planet}',...
        'amu/cm^3','km/s','km/s','km/s','nT','nT','nT','IDK','nPa','nPa','{\mu}A/m^2','{\mu}A/m^2','{\mu}A/m^2','nT','km/s','{\mu}A/m^2'};
    Unit = containers.Map(keySet,unitSet,'UniformValues',false);

    N_var = max(size(VAR_str));

    ! rm -rf My_movie_2
    ! mkdir My_movie_2

    for K = 1 : N_png   

        figure('Position',[w,0,W-w,h]);

        for i = 1:N    % Satellite loop
            for j = 1:N_var  % Variable loop
                subplot(N_var,1,j);
                tag = ['x=' num2str(Locations{i}(1)) ' y=' num2str(Locations{i}(2)) ' z=' num2str(Locations{i}(3))];
                plot(DATA{i}(Dict('t'),:),DATA{i}(Dict(VAR_str(j)),:),'DisplayName',tag,'LineWidth',3);
                if i == 1
                    ylabel(strcat(VAR_str(j) , '(', Unit(VAR_str(j)) , ')'),'FontSize',22);
                end
                hold on;
            end
        end

        for j = 1:N_var
            subplot(N_var,1,j);
            if j == N_var  % Show legend in N_var(th) subplot, you can change this number as you wish
                legend;
            end
            X_lim = xlim;
            vl = vline(X_lim(1) + (X_lim(2) - X_lim(1))*(K-1)/(N_png-1),'k');
            set(vl,'LineWidth',3);
            ax = gca;
            ax.FontSize = 20; % The fontsize of tick label
        end

        filename = sprintf('%4.4d.png', K);
        saveas(gcf,['My_movie_2/' filename]);
        close gcf;

    end

    ! rm -rf My_movie
    ! mkdir My_movie

    for K = 1 : N_png 
        filename = sprintf('%4.4d.png', K);
        Name_1 = ['My_movie_1/' filename];
        Name_2 = ['My_movie_2/' filename];
        Name   = ['My_movie/'   filename];

        ! touch name_1.txt
        tmp = fopen('name_1.txt','w');
        fprintf(tmp,'%s',Name_1);
        fclose(tmp);

        ! touch name_2.txt
        tmp = fopen('name_2.txt','w');
        fprintf(tmp,'%s',Name_2);
        fclose(tmp);

        ! touch name.txt
        tmp = fopen('name.txt','w');
        fprintf(tmp,'%s',Name);
        fclose(tmp);

        ! /opt/local/bin/convert $(cat name_1.txt) $(cat name_2.txt) +append $(cat name.txt)
        ! rm -rf name*.txt
    end

    ! cd My_movie/; /opt/local/bin/convert *.png movie.pdf

    writerObj = VideoWriter('My_movie/movie.avi');
    writerObj.FrameRate=4;
    open(writerObj);

    for K = 1 : N_png
      filename = sprintf('%4.4d.png', K);
      thisimage = imread(['My_movie/' filename]);
      writeVideo(writerObj, thisimage);
    end
    close(writerObj);

end
