%% Startup

locate_satellite;

%% Read & Construct data

DATA = cell(N,1);
for i = 1:N
    DATA{i} = read_sat(Locations{i}(1),Locations{i}(2),Locations{i}(3));
end

%% Plot out time series for selected satellites

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
                ylabel(strcat(VAR_str(j) , '(', Unit(VAR_str(j)) , ')'));
            end
            hold on;
        end
    end

    for j = 1:N_var
        subplot(N_var,1,j);
        legend;
        X_lim = xlim;
        vl = vline(X_lim(1) + (X_lim(2) - X_lim(1))*(K-1)/(N_png-1),'k');
        set(vl,'LineWidth',3);
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

