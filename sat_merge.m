clear;

N = input("How many satellite files do you want to merge in one time?\n");
index = zeros(N,1);
for i = 1:N
   index(i) = input(['Input index(' num2str(i) ')\n']); 
end

Name = cell(N,1);
for i = 1:N
   Name{i} = dir(['*_n' num2str(index(i)) '.sat']); 
end

N_sat = max(size(Name{1}));

for j = 1:N_sat
    for i = 2:N
        if Name{i}(j).bytes == 0
            system(['rm ' Name{i}(j).name]);
        else
            sat_cat(Name{1}(j).name,Name{i}(j).name); 
        end
    end
end

function sat_cat(sat_1,sat_2)
    s1 = fopen(sat_1);
    s2 = fopen(sat_2);
    
    system(['wc -l ' sat_1 ' > numline.txt']);
    f_1 = fopen('numline.txt');
    Nt_1 = fscanf(f_1,'%i');
    Nt_1 = Nt_1 -2; % Delete the header
    fclose(f_1);
    ! rm numline.txt;
    
    system(['wc -l ' sat_2 ' > numline.txt']);
    f_2 = fopen('numline.txt');
    Nt_2 = fscanf(f_2,'%i');
    Nt_2 = Nt_2 -2; % Delete the header
    fclose(f_2);
    ! rm numline.txt;
    
    header_s1_1 = fgetl(s1);
    header_s1_2 = fgetl(s1);
    data_1 = fscanf(s1,'%f');
    fclose(s1);
    total_1 = max(size(data_1));
    data_1 = reshape(data_1,[total_1/Nt_1,Nt_1]);
    
    header_s2_1 = fgetl(s2);
    header_s2_2 = fgetl(s2);
    data_2 = fscanf(s2,'%f');
    fclose(s2);
    total_2 = max(size(data_2));
    data_2 = reshape(data_2,[total_2/Nt_2,Nt_2]);
    
    for i = 1:Nt_2
       if data_1(9,end) == data_2(9,i) || data_1(1,end) == data_2(1,i)
           num_deleted_line = i + 2;  % Plus 2 lines of header
           break;
       end
    end
    
    s1 = fopen(sat_1,'a');
    s2 = fopen(sat_2);
    
    for i = 1:Nt_2+2  % Plus 2 lines of header
        line = fgetl(s2);
        
        if i > num_deleted_line
            fprintf(s1,"%s\n",line);
        end
    end
    
    fclose(s1);
    fclose(s2);
    
    system(['rm ' sat_2]);
end
