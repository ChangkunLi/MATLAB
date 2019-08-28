clear;

global w h Sat_namelist Sat_VAR N_png
global Sat_namelist_g Sat_VAR_g h7 h9 h11 h13

screen = get(0,'ScreenSize');
W_g = screen(3); H_g = screen(4);
w_g = 0.6*W_g;              % width of the GUI window
h_g = 0.8*H_g;              % hight of the GUI window

figure('Color',[1 1 1],'Position',[0,0,w_g,0.7*h_g],'Name','User Input GUI','NumberTitle','off','MenuBar','none');

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

h15 = uicontrol(gcf,'style','push','string','Plot','position',[0.7*w_g 0.05*h_g 0.25*w_g 0.49*h_g],'callback',@h15_callback);
h15.FontSize = 100;

t16 = uicontrol(gcf,'Style','text','String','Number of Frames','Position',[0.7*w_g 0.63*h_g 0.25*w_g 0.05*h_g]);
t16.FontSize = 25;
h16 = uicontrol(gcf,'style','edit','string',' ','Position',[0.7*w_g 0.56*h_g 0.25*w_g 0.05*h_g],'callback',@h16_callback);
h16.FontSize = 20;

% Default setting
w = 0.6*W_g;
h = 0.8*H_g;
Sat_namelist = {' '};

% Sat_VAR = {' '};
Sat_VAR = " ";

N_png = 501;

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
    global w h Sat_namelist Sat_VAR N_png
    
    screen=get(0,'ScreenSize');
    W=screen(3);

    Size = size(Sat_namelist);
    N = Size(2) - 1;
    
    Locations = cell(N,1);
    for i = 1:N
        Locations{i} = location_extract(Sat_namelist{i+1});
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

    ! rm -rf My_movie_sat
    ! mkdir My_movie_sat
    
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
        saveas(gcf,['My_movie_sat/' filename]);
        close gcf;

    end

    writerObj = VideoWriter('My_movie_sat/movie.avi');
    writerObj.FrameRate=4;
    open(writerObj);

    for K = 1 : N_png
      filename = sprintf('%4.4d.png', K);
      thisimage = imread(['My_movie_sat/' filename]);
      writeVideo(writerObj, thisimage);
    end
    close(writerObj);

    ! cd My_movie_sat; /opt/local/bin/convert *.png movie.pdf
end

function h16_callback(hObject,~)
    global N_png

    m = get(hObject,'string');
    N_png = str2double(m);
end
