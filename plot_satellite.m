% Changkun Li, May 31st, 2019
% The PLANETARY units are used in the satellite log file (in the format of .sat)


%% Input the number of satellites you want
clear;

prompt = "Would you like to compare satellites at different locations?(yes/no)";
flag = input(prompt,'s');
if strcmp(flag,"yes")
    N = input("How many satellites would you like?      ");    
else
    N = 1;
end

clear prompt;

%% Input the locations of satellites

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

%% Read & Construct data

DATA = cell(N,1);
for i = 1:N
    DATA{i} = read_sat(Locations{i}(1),Locations{i}(2),Locations{i}(3));
end

%% Plot out time series for selected satellites

figure;

disp("Please enter the physical variable(s)");
disp("Exemplary format of input:    [""Bx"" ""Rho""  ""P""]      (press ""return"" to proceed)"); 
pause;
disp("Enter your choice below");
VAR_str = input("");
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

for i = 1:N    % Satellite loop
    for j = 1:N_var  % Variable loop
        subplot(N_var,1,j);
        tag = ['x=' num2str(Locations{i}(1)) ' y=' num2str(Locations{i}(2)) ' z=' num2str(Locations{i}(3))];
        plot(DATA{i}(Dict('t'),:),DATA{i}(Dict(VAR_str(j)),:),'DisplayName',tag,'LineWidth',3,'Color','r');
        
        if i == 1
            ylabel(strcat(VAR_str(j) , '(', Unit(VAR_str(j)) , ')'));
        end
        hold on;
    end
end

for j = 1:N_var
    subplot(N_var,1,j);
    legend;
end


%% Function definition

function data = read_sat(x,y,z)
    filename = ['*_x=' num2str(x) '_y=' num2str(y) '_z=' num2str(z) '*.sat'];
    ! touch tmp.txt
    tmp = fopen('tmp.txt','w');
    fprintf(tmp,'%s',filename);
    fclose(tmp);
    ! ls $(cat tmp.txt) > filename.txt
    ! rm -f tmp.txt
    
    
    f = fopen('filename.txt');
    filename = fscanf(f,'%s');
    fclose(f);
    ! wc -l $(cat filename.txt) > numline.txt
    ! rm filename.txt
    f = fopen('numline.txt');
    Nt = fscanf(f,'%i');
    Nt = Nt -2; % Delete the header
    fclose(f);
    ! rm numline.txt;
    sat = fopen(filename);
    
    header_1 = fgetl(sat);
    header_2 = fgetl(sat);
    data = fscanf(sat,'%f');
    fclose(sat);
    total = max(size(data));
    data = reshape(data,[total/Nt,Nt]);
end