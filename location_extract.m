% Extract satellite location from filename

function location = location_extract(sat_filename)
    a = strfind(sat_filename,'x=');
    b = strfind(sat_filename,'y=');
    c = strfind(sat_filename,'z=');
    d = strfind(sat_filename,'_n');
    location = [str2double(sat_filename(a+2:b-2)) str2double(sat_filename(b+2:c-2)) str2double(sat_filename(c+2:d-1))];
end