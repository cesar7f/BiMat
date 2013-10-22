% MatrixFunctions - Static class that contains a set of useful matrix
% manipulation functions
%
% NullModels Methods:
%    SORT_MATRIX - Sort a matrix
%    SORT_MATRIX - Randomly permute a matrix
%    GET_DEGREES - Get degrees (sum) of both rows and column nodes of a bipartite matrix
%    GET_FILL - Get the connectance
%    NON_ZERO_MATRIX - Delete empty rows and columns of a bipartite matrix
%    TYPE_MATRIX - Return a matrix with only certain type of interactions
%    TYPE_MATRIX_NON_ZERO - Return a matrix with only certain type of interactions and delete all empty rows and columns.
%    BIPARTITE_TO_UNIPARTITE - Return the unipartite version of matrix.
%    NESTED_MATRIX - Gest an square perfect nested matrix
%    GET_BIGGEST_EIGENVALUE - Get the biggest eigenvalue of a bipartite matrix
%    RANDOM_SUBMATRIX - Get a random submatrix of the bipartite matrix
%    MATRIX_UNION - Create a bipartite matrix by the union of two bipartite matrices.
%
% See also:
%    
classdef MatrixFunctions < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = MatrixFunctions

        end
    end 
   
    methods(Static)
        
        function sorted_matrix = SORT_MATRIX(matrix)
        % SORT_MATRIX - Sort a matrix
        %   sorted_matrix = SORT_MATRIX(matrix) Return the matrix with a
        %   permutation in which sums in rows and columns are in decreasing
        %   order.
        %
            [~, ir] = sort(sum(matrix,2),'descend');
            [~, ic] = sort(sum(matrix,1),'descend');
           
            sorted_matrix = matrix(ir,ic);
            
        end
        
        function [matrix index_row index_col] = RANDOM_SORT(matrix)
        % SORT_MATRIX - Randomly permute a matrix
        %   matrix_permutation = RANDOM_SORT(matrix) Return the matrix in
        %   which columns and rows were permuted randomly
        %   
        %  [matrix_permutation index_row index_col] = RANDOM_SORT(matrix)
        %  In addition to return a random permutation of the matrix, we get
        %  the indexes of the permutations for both rows (index_row) and
        %  columns (index_col)
            [n_rows n_cols] = size(matrix);
            
            index_row = randperm(n_rows);
            index_col = randperm(n_cols);
            matrix = matrix(index_row,index_col);
            
        end
        
        
        function [k d] = GET_DEGREES(matrix)
        % GET_DEGREES - Get degrees (sum) of both rows and column nodes of
        % a bipartite matrix
        %   [k d] = GET_DEGREES(matrix) Get the degrees of both rows (k)
        %   and columns (d). k is the same size than the number of rows in
        %   the matrix. Similarly d is a vector of the same size than the
        %   number of columns in the matrix
        %   
            k = sum(matrix,2);
            d = sum(matrix,1);
            
        end
        
        function p = GET_FILL(matrix)
        % GET_FILL - Get the connectance
        %   p = GET_FILL(matrix) Get the connectance (or fill) p of a matrix      
            p = sum(sum(matrix))/numel(matrix);
        end
        
        function [matrix ir ic] = NON_ZERO_MATRIX(matrix)
        % NON_ZERO_MATRIX - Delete empty rows and columns of a matrix
        %   matrix = NON_ZERO_MATRIX(matrix) Return the original matrix
        %   with no empty rows and columns
        %
        %   [matrix ir ic] = NON_ZERO_MATRIX(matrix) Return the original matrix
        %   with no empty rows and columns. In addition the function return
        %   the indices of the original matrix that will continue to exist
        %   in the new matrix. ir and ic represent rows and column indices,
        %   respectivally.
            ir = find(sum(matrix,2)>0);
            ic = find(sum(matrix,1)>0);
            
            matrix = matrix(ir,ic);
         
        end
        
        function has_empty_nodes = HAS_EMPTY_NODES(matrix)
        % HAS_EMPTY_NODES - Check if empty rows/columns exist in the matrix
        % which will translate to have nodes without interactions
        %
        %   has_empty_nodes = HAS_EMPTY_NODES(matrix) Return true(1) if emty
        %   row/columns exist in the matrix, 0 otherwise.
        
            has_empty_nodes = sum(sum(matrix,1)==0)||sum(sum(matrix,2)==0);
                
        end
        
        function [matrix ir ic] = TYPE_MATRIX(matrix, n_type)
        % TYPE_MATRIX - Return a matrix with only certain type of
        % interactions
        %   matrix = TYPE_MATRIX(matrix,n_type) Return a matrix of the same size
        %   as the original matrix that will only contain the elements for
        %   which matrix(i,j) == n_type. In order words, this a filter for
        %   matrix values.
            [n_rows n_cols] = size(matrix);
            ir = 1:n_rows; ic = 1:n_cols;
            %ir = find(sum(matrix==n_type,2)>0);
            %ic = find(sum(matrix==n_type,1)>0);
            
            %matrix = matrix(ir,ic);
            matrix = (matrix==n_type)*n_type;
            
        end
        
        function [matrix ir ic] = TYPE_MATRIX_NON_ZERO(matrix, n_type)
        % TYPE_MATRIX_NON_ZERO - Return a matrix with only certain type of
        % interactions and delete all empty rows and columns.
        %   matrix = TYPE_MATRIX_NON_ZERO(matrix) Filter the matrix such that only
        %   elements for which matrix(i,j) == n_type continue to exist. In
        %   addition delete all empty rows and comuns after the filter.
        %
        %   [matrix ir ic] = TYPE_MATRIX_NON_ZERO(matrix,n_type) Filter the matrix such that only
        %   elements for which matrix(i,j) == n_type continue to exist. In
        %   addition delete all empty rows and comuns after the filter.
        %   Additionaly return the original rows (ir) and column (ic)
        %   indices of the original matrix that will continue to exist in
        %   the new matrix.
            [n_rows n_cols] = size(matrix);
            ir = 1:n_rows; ic = 1:n_cols;
            
            [matrix ir_type  ic_type] = MatrixNull.TYPE_MATRIX(matrix, n_type);
            [matrix ir_zero  ic_zero] = MatrixNull.NON_ZERO_MATRIX(matrix);
            
            ir = ir(ir_type(ir_zero));
            ic = ic(ic_type(ic_zero));
            
        end
        
        function A = BIPARTITE_TO_UNIPARTITE(matrix)
        % BIPARTITE_TO_UNIPARTITE - Return the unipartite version of
        % matrix.
        %   A = BIPARTITE_TO_UNIPARTITE(matrix) Return the unipartite
        %   version of the bipartite matrix
            [n_rows, n_cols] = size(matrix);
            A = [zeros(n_rows,n_rows),matrix;matrix',zeros(n_cols,n_cols)];
            
        end
        
        function matrix = NESTED_MATRIX(n)
        % NESTED_MATRIX - Gest an square perfect nested matrix
        %   matrix = NESTED_MATRIX(n) Get an square perfect matrix of size
        %   n by n.
            matrix = flipud(tril(ones(n)));
            
        end
        
        function max_eig = GET_BIGGEST_EIGENVALUE(matrix)
        % GET_BIGGEST_EIGENVALUE - Get the biggest eigenvalue of a bipartite matrix
        %   max_eig = GET_BIGGEST_EIGENVALUE(n)  Get the biggest eigenvalue of a bipartite matrix
            max_eig = max(eig(MatrixNull.BipartiteToUnipartite(matrix)));
            
        end
        
        function submatrix = RANDOM_SUBMATRIX(matrix, n_rows, n_cols)
        % RANDOM_SUBMATRIX - Get a random submatrix of the bipartite matrix
        %   submatrix = RANDOM_SUBMATRIX(matrix, n_rows, n_cols)  Get a
        %   random submatrix of size n_rows by n_cols of the original
        %   matrix. n_rows and n_cols must be smaller or equal to the
        %   actual values in the original matrix.
            [m n] = size(matrix);
            assert(m>=n_rows); assert(n>=n_cols);
            
            submatrix = matrix(randsample(m, n_rows), randsample(n,n_cols));
            
        end
        
        function matrix = MATRIX_UNION(matrix_1,matrix_2)
        % MATRIX_UNION - Create a bipartite matrix by the union of two bipartite matrices.
        %   matrix = MATRIX_UNION(matrix_1,matrix_2) Create a bipartie
        %   matrix that is the composition of to bipartite matrices. The
        %   new size of the composed matrix is the sum of the two original
        %   matrices.
            matrix = blkdiag(matrix_1,matrix_2);
        end
        
        function matrix = SWAP_COLUMNS(matrix,i,j)
            matrix(:,[i j]) = matrix(:,[j i]);
        end
        
        function matrix = SWAP_ROWS(matrix,i,j)
            matrix([i j],:) = matrix([j i],:);
        end
        
        function matrix = MOVE_COLUMNS(matrix,m)
            
            [~,n_cols] = size(matrix);
            assert(m<=n_cols);
            
            temp = matrix(:,1:(n_cols-m));
            matrix(:,1:m) = matrix(:,n_cols-m+1:n_cols);
            matrix(:,m+1:n_cols) = temp;
            
        end
        
        function matrix = MOVE_ROWS(matrix,m)
            
            [n_rows,~] = size(matrix);
            assert(m<=n_rows);
            
            temp = matrix(1:(n_rows-m),:);
            matrix(1:m,:) = matrix(n_rows-m+1:n_rows,:);
            matrix(m+1:n_rows,:) = temp;
            
        end
        
        function matrix = MAX_NESTED_MOTIFS(m)
            assert(mod(m,4)==0);
            
            bs = m/4; %block size
            
            matrix = 0;%[zeros(bs,bs) ones(
            
        end
        
    end
end