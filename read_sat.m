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