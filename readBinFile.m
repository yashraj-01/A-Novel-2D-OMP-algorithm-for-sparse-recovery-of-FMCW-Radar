function [res, rows, cols] = readBinFile(fileName)
    fileID = fopen(fileName,'r');
    file = fread(fileID,'double');
    fclose(fileID);

    rows = file(1);
    cols = file(2);

    res_vec = zeros(rows*cols, 1);
    index = 1;
    for i=3:2:length(file)
        res_vec(index) = complex(file(i), file(i+1));
        index = index+1;
    end
    
    res = reshape(res_vec, rows, cols);
end