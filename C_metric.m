clear
clc
fclose all;

num = 2;
al{1} = 'C:\Users\Desktop\SCMOEA\';
al{2} = 'C:\Users\Desktop\IMOMAD\';
sc_1 = zeros(100, 1);
sc_2 = zeros(100, 1); 

for idx = 1:100
    ds = sprintf('Mk%d/', idx);
    name = 'res';
    txt = '.txt';
    k = 1;
    x_1 = []; y_1 = [];
    x_2 = []; y_2 = [];

    for file = 1:10 
        RealPath_1 = [al{1}, ds, name, num2str(file), txt];
        RealPath_2 = [al{2}, ds, name, num2str(file), txt];

        [x, y] = DataRead(RealPath_1);
        x_1 = [x_1, x];
        y_1 = [y_1, y];

        [x, y] = DataRead(RealPath_2);
        x_2 = [x_2, x];
        y_2 = [y_2, y];
    end

    A_1 = [x_1; y_1]'; 
    A_2 = [x_2; y_2]';
  
    dominated_count_1 = 0;
    for i = 1:size(A_1, 1)
        for j = 1:size(A_2, 1)
            if all(A_1(i, :) <= A_2(j, :)) && any(A_1(i, :) < A_2(j, :))
                dominated_count_1 = dominated_count_1 + 1;
                break;
            end
        end
    end
    sc_1(idx) = dominated_count_1 / size(A_1, 1);

    dominated_count_2 = 0;
    for i = 1:size(A_2, 1)
        for j = 1:size(A_1, 1)
            if all(A_2(i, :) <= A_1(j, :)) && any(A_2(i, :) < A_1(j, :))
                dominated_count_2 = dominated_count_2 + 1;
                break;
            end
        end
    end
    sc_2(idx) = dominated_count_2 / size(A_2, 1);
end

respath = 'C.txt';
fout = fopen(respath, 'w');
fprintf(fout, '1 Set coverage: %4f \r\n', sc_1);
fprintf(fout, '2 Set coverage: %4f \r\n', sc_2);
fclose(fout);

function [x, y] = DataRead(RealPath)
fin = fopen(RealPath, 'r'); 
sizeA = [2 Inf];
A = fscanf(fin, '%f %f', sizeA);
x = A(1, :);
y = A(2, :);
fclose(fin);
end