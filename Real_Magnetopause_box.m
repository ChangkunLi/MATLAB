c_dir = pwd;
wk_dir = '/Users/changkul/Desktop/Utilities/Tecplot_field_line_tracer';
cd(wk_dir);
system(['python3 extract_status_box.py ''' c_dir '/3d*.plt'' -i 1']);
cd(c_dir);
cd('data_status_box');

listing = dir('3D*.plt'); % Get the name of all 3D box files
N = max(size(listing));  % Find the number of 3D box files

for i = 1:N
    f = fopen('plt2dat.mcr','w');
    fprintf(f,'%s\n','#!MC 1410');
    fprintf(f,'%s\n','$!Pick AddAtPosition');
    fprintf(f,'%s\n','  X = 0.38');
    fprintf(f,'%s\n','  Y = 2.71974358974');
    fprintf(f,'%s\n','  ConsiderStyle = Yes');
    
    filename = listing(i).name;
    filename_old = filename;
    a = strfind(filename,'.plt');
    filename(a:end) = '.dat';
    fprintf(f,'%s\n',['$!WriteDataSet  "' c_dir '/data_status_box/' filename '"']);
    fprintf(f,'%s\n','  IncludeText = No');
    fprintf(f,'%s\n','  IncludeGeom = No');
    fprintf(f,'%s\n','  IncludeCustomLabels = No');
    fprintf(f,'%s\n','  IncludeDataShareLinkage = Yes');
    fprintf(f,'%s\n','  Binary = No');
    fprintf(f,'%s\n','  UsePointFormat = Yes');
    fprintf(f,'%s\n','  Precision = 9');
    fprintf(f,'%s\n','  TecplotVersionToWrite = TecplotCurrent');
    fclose(f);  
    
    system(['tec360 -nostdaddons -b -p plt2dat.mcr ' filename_old]);
end

cd(c_dir);
! cp Magnetopause_box.m data_status_box/
cd('data_status_box');
Magnetopause_box;