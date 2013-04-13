BOOL_MATRIX = cell(20,20)
for i=1:20
    for j=1:20
        if Y(i,j)==0
            BOOL_MATRIX{i,j}='false';
        else
            BOOL_MATRIX{i,j}='true';
        end
    end
end