function [nvalue nr nc] = NODF_blind(matrix)

    matrix = matrix ~= 0;
    [nrows ncols] = size(matrix);

    sort_matrix;
    
    sumrows = sum(matrix,2);
    sumcols = sum(matrix,1);

    nvalue = 0;
    nr = 0;
    %Fill for columns
    for i = 1:nrows
        for j = i+1:nrows
            if( sumrows(j) < sumrows(i) && sumrows(j) > 0 )
                nr = nr + 100*sum(matrix(i, matrix(j,:)==1))/sumrows(j);
            end
        end
    end
    
    nc = 0;
    for k = 1:ncols
        for l = k+1:ncols
            if( sumcols(l) < sumcols(k) && sumcols(l) > 0 )
                nc = nc + 100*sum(matrix(matrix(:,l)==1,k))/sumcols(l);
            end
        end
    end
    
    nvalue = nr + nc;
    denom = nrows*(nrows-1)/2 + ncols*(ncols-1)/2;
    nr = nr/(nrows*(nrows-1)/2);
    nc = nc/(ncols*(ncols-1)/2);
    nvalue = nvalue / (denom);
    
    imagesc(matrix);
    
    function sort_matrix
        
        [values rowsort] = sort(sum(matrix,2),'descend');
        [values colsort] = sort(sum(matrix,1),'descend');
        
        matrix = matrix(rowsort,colsort);
    end
end
