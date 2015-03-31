% MatrixFunctions - Static class that contains a set of useful matrix
% manipulation functions
%
% MatrixFunctions Methods:
%    SORT_MATRIX - Sort a matrix
%    SORT_MATRIX - Randomly permute a matrix
%    GET_DEGREES - Get degrees (sum) of both rows and column nodes of a bipartite matrix
%    GET_FILL - Get the connectance
%    NON_ZERO_MATRIX - Delete empty rows and columns of a bipartite matrix
%    HAS_EMPTY_NODES - Check if empty rows/columns exist in the matrix
%    TYPE_MATRIX - Return a matrix with only certain type of interactions
%    TYPE_MATRIX_NON_ZERO - Return a matrix with only certain type of interactions and delete all empty rows and columns.
%    BIPARTITE_TO_UNIPARTITE - Return the unipartite version of matrix
%    NESTED_MATRIX - Get a square perfect nested matrix
% MIX_MATRICES - Mix two matrices randomly
%    BLOCK_MATRIX - Get an square matrix with n blocks
%    NEGATE_MATRIX - Negates a boolean matrix
%    GET_BIGGEST_EIGENVALUE - Get the biggest eigenvalue of a bipartite matrix
%    RANDOM_SUBMATRIX - Get a random submatrix of the bipartite matrix
%    MATRIX_UNION - Create a bipartite matrix by the union of two bipartite matrices.
%    ISOLATED_COMPONENTS - Find the components of a bipartite adjacency matrix
%
   
classdef MatrixFunctions < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = MatrixFunctions

        end
    end 
   
    methods(Static)
        
        function [sorted_matrix, ir,ic] = SORT_MATRIX(matrix)
        % SORT_MATRIX - Sort a matrix
        %
        %   matrix = SORT_MATRIX(MATRIX) Return the matrix with a
        %   permutation in which sums in rows and columns are in decreasing
        %   order.
        %
        %   [matrix, index_row, index_col] = SORT_MATRIX(MATRIX) In addition to returning
        %   the sorted version of MATRIX, it also returns the rows (index_row)
        %   and column (index_col) indices positioned in their new positions.
        
            [~, ir] = sort(sum(matrix,2),'descend');
            [~, ic] = sort(sum(matrix,1),'descend');
           
            sorted_matrix = matrix(ir,ic);
            
        end
        
        function [matrix,index_row,index_col] = RANDOM_SORT(matrix)
        % SORT_MATRIX - Randomly permute a matrix
        %
        %   matrix_permutation = RANDOM_SORT(MATRIX) Return the matrix in
        %   which columns and rows were permuted randomly
        %   
        %  [matrix_permutation index_row index_col] = RANDOM_SORT(MATRIX)
        %  In addition to return a random permutation of MATRIX, we get
        %  the indexes of the permutations for both rows (index_row) and
        %  columns (index_col)
            [n_rows,n_cols] = size(matrix);
            
            index_row = randperm(n_rows);
            index_col = randperm(n_cols);
            matrix = matrix(index_row,index_col);
            
        end
        
        function [k,d] = GET_DEGREES(matrix)
        % GET_DEGREES - Get degrees (sum) of both rows and column nodes of
        % a bipartite matrix
        %
        %   [k d] = GET_DEGREES(MATRIX) Get the degrees of both rows (k)
        %   and columns (d). k is the same size than the number of rows in
        %   the matrix. Similarly d is a vector of the same size than the
        %   number of columns in the matrix
        %   
            k = sum(matrix,2);
            d = sum(matrix,1);
            
        end
        
        function p = GET_FILL(matrix)
        % GET_FILL - Get the connectance
        %
        %   p = GET_FILL(MATRIX) Get the connectance (or fill) p of MATRIX     
        %
            p = sum(sum(matrix))/numel(matrix);
        end
        
        function [matrix,ir,ic] = NON_ZERO_MATRIX(matrix)
        % NON_ZERO_MATRIX - Delete empty rows and columns of a matrix
        %
        %   matrix = NON_ZERO_MATRIX(MATRIX) Return the original matrix
        %   with no empty rows and columns
        %
        %   [matrix ir ic] = NON_ZERO_MATRIX(MATRIX) Return the original matrix
        %   with no empty rows and columns. In addition the function return
        %   the indices of the original matrix that will continue to exist
        %   in the new matrix. ir and ic represent rows and column indices,
        %   respectivally.
        %
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
        %
        
            has_empty_nodes = sum(sum(matrix,1)==0)||sum(sum(matrix,2)==0);
                
        end
        
        function [matrix,ir,ic] = TYPE_MATRIX(matrix, n_type)
        % TYPE_MATRIX - Return a matrix with only certain type of
        % interactions
        %
        %   matrix = TYPE_MATRIX(MATRIX,n_type) Return a matrix of the same size
        %   as the original MATRIX that will only contain the elements for
        %   which matrix(i,j) == n_type. In order words, this a filter for
        %   matrix values.
        %
            [n_rows, n_cols] = size(matrix);
            ir = 1:n_rows; ic = 1:n_cols;
            %ir = find(sum(matrix==n_type,2)>0);
            %ic = find(sum(matrix==n_type,1)>0);
            
            %matrix = matrix(ir,ic);
            matrix = (matrix==n_type)*n_type;
            
        end
        
        function [matrix,ir,ic] = TYPE_MATRIX_NON_ZERO(matrix, n_type)
        % TYPE_MATRIX_NON_ZERO - Return a matrix with only certain type of
        % interactions and delete all empty rows and columns.
        %
        %   matrix = TYPE_MATRIX_NON_ZERO(MATRIX) Filter the MATRIX such that only
        %   elements for which MATRIX(i,j) == n_type continue to exist. In
        %   addition, it deletes all empty rows and comuns after the filter.
        %
        %   [matrix, ir, ic] = TYPE_MATRIX_NON_ZERO(MATRIX,n_type) Filter the matrix such that only
        %   elements for which MATRIX(i,j) == n_type continue to exist. In
        %   addition delete all empty rows and comuns after the filter.
        %   Additionaly returns the original rows (ir) and column (ic)
        %   indices of the original matrix that will continue to exist in
        %   the new matrix.
        %
            [n_rows,n_cols] = size(matrix);
            ir = 1:n_rows; ic = 1:n_cols;
            
            [matrix,ir_type,ic_type] = MatrixNull.TYPE_MATRIX(matrix, n_type);
            [matrix,ir_zero,ic_zero] = MatrixNull.NON_ZERO_MATRIX(matrix);
            
            ir = ir(ir_type(ir_zero));
            ic = ic(ic_type(ic_zero));
            
        end
        
        function A = BIPARTITE_TO_UNIPARTITE(matrix)
        % BIPARTITE_TO_UNIPARTITE - Return the unipartite version of
        % matrix.
        %   A = BIPARTITE_TO_UNIPARTITE(MATRIX) Return the unipartite
        %   version of the bipartite matrix
        %
            [n_rows, n_cols] = size(matrix);
            A = [zeros(n_rows,n_rows),matrix;matrix',zeros(n_cols,n_cols)];
            
        end
        
        function matrix = NESTED_MATRIX(n)
        % NESTED_MATRIX - Get an square perfect nested matrix
        %
        %   matrix = NESTED_MATRIX(n) Get an square perfect matrix of size
        %   n by n.
        %
            matrix = flipud(tril(ones(n)));
            
        end
        
        function matrix = CHECKERBOARD(n,m)
        % CHECKERBOARD - Get checkerboard matrix
        %
        %   matrix = CHECKERBOARD(n,m) Get an checkerboard matrix of size
        %   n,m
        %    
            matrix = zeros(n,m);
            for i = 1:n
                for j = 1:m
                    if(mod(i+j,2) == 0)
                        matrix(i,j) = 1.0;
                    end
                end
            end
        
        end
        
        function matrix = BLOCK_MATRIX(n,s)
        % BLOCK_MATRIX - Get an square matrix of n blocks
        %
        %   matrix = BLOCK_MATRIX(n,s ) Get a square matrix with n blocks
        %   of size s x s
        %
            block_matrix = ones(s,s);
            matrix = kron(eye(n),block_matrix);
        end
        
        function matrix = MIX_MATRICES(matrix_1,matrix_2, p)
        % MIX_MATRICES - Mix two matrices randomly
        %
        %   matrix = MIX_MATRICES(MATRIX_1,MATRIX_2) Get a matrix in
        %   which each element of MATRIX_1 stays with probability 0.5 or is
        %   replaced by the corresponding element in MATRIX_2 with the same
        %   probability
        %
        %   matrix = MIX_MATRICES(MATRIX_1,MATRIX_2,p) Get a matrix in
        %   which each element of MATRIX_1 stays with probability p or is
        %   replaced by the corresponding element in MATRIX_2 with
        %   probability 1-p
        %
        
            if(nargin == 2)
                p = 0.5;
            end
            
            assert(sum(size(matrix_1) == size(matrix_2))==2);
            assert(p <= 1); assert(p >= 0);
            
            rand_matrix = rand(size(matrix_1));
        
            matrix = (rand_matrix <= p).*matrix_1 + ...
                (rand_matrix >= p).*matrix_2;
            
            
        end
        
        function matrix = NEGATE_MATRIX(matrix, prob)
        % NEGATE_MATRIX - Negates a boolean matrix
        %
        %   matrix = NEGATE_MATRIX(MATRIX) Returns a matrix in which
        %   each element of a boolean matrix is neggated 
        %
        %   matrix = NEGATE_MATRIX(MATRIX,prob) Returns a matrix in which
        %   each element of a boolean matrix is neggated with probability prob 
        
            if(nargin == 1)
                prob = 1.0;
            end
        
            matrix = matrix - abs(rand(size(matrix) < prob));
        end
        
        function max_eig = GET_BIGGEST_EIGENVALUE(matrix)
        % GET_BIGGEST_EIGENVALUE - Get the biggest eigenvalue of a bipartite matrix
        %
        %   max_eig = GET_BIGGEST_EIGENVALUE(MATRIX)  Get the biggest eigenvalue of a bipartite matrix
        %
            max_eig = max(eig(MatrixNull.BipartiteToUnipartite(matrix)));
            
        end
        
        function submatrix = RANDOM_SUBMATRIX(matrix, n_rows, n_cols)
        % RANDOM_SUBMATRIX - Get a random submatrix of the bipartite matrix
        %
        %   submatrix = RANDOM_SUBMATRIX(MATRIX, n_rows, n_cols)  Get a
        %   random submatrix of size n_rows by n_cols of the original
        %   matrix. n_rows and n_cols must be smaller or equal to the
        %   actual values in the original matrix.
        %
            [m, n] = size(matrix);
            assert(m>=n_rows); assert(n>=n_cols);
            
            submatrix = matrix(randsample(m, n_rows), randsample(n,n_cols));
            
        end
        
        function matrix = MATRIX_UNION(matrix_1,matrix_2)
        % MATRIX_UNION - Create a bipartite matrix by the union of two bipartite matrices.
        %
        %   matrix = MATRIX_UNION(matrix_1,matrix_2) Create a bipartie
        %   matrix that is the composition of to bipartite matrices. The
        %   new size of the composed matrix is the sum of the two original
        %   matrices.
        %
            matrix = blkdiag(matrix_1,matrix_2);
        end
        
        function [comp_row,comp_col,n_comp] = ISOLATED_COMPONENTS(matrix)
        % ISOLATED_COMPONENTS - Find the components of a bipartite
        % adjacency matrix
        %
        %   [comp_row comp_col c] = ISOLATED_COMPONENTS(MATRIX) Find the
        %   isolated components of a bipartite adjacency matrix. It returns
        %   the number of components in the variable c, and the
        %   corresponding component index for rows and column nodes in
        %   comp_row and comp_col, respectively.
        
            [n_rows,n_cols] = size(matrix);
            
            uni_matrix = MatrixFunctions.BIPARTITE_TO_UNIPARTITE(matrix);
            comp_idx = zeros(n_rows+n_cols,1);
            
            n_comp = 0;
            
            node_idx = find(comp_idx==0,1);
            while(~isempty(node_idx))
                n_comp = n_comp + 1;
                visited = [];
                DFS_COMPONENT(node_idx);
                comp_idx(visited) = n_comp;
                node_idx = find(comp_idx==0,1);
            end
            
            comp_row = comp_idx(1:n_rows);
            comp_col = comp_idx(1+n_rows:n_rows+n_cols);
        
            function DFS_COMPONENT(idx)
           
                visited  = [visited idx];

                neighbors = find(uni_matrix(idx,:)>0);

                for i = 1:length(neighbors)   
                    if(isempty( find(visited == neighbors(i),1))) 
                        DFS_COMPONENT(neighbors(i));
                    end
                end

            end
            
        end
        
        function matrix = GET_FILTERED_MATRIX(matrix, min_value,max_value)
        % GET_FILTERED_MATRIX - Filter the matrix according to min and max
        % value
        %
        %   matrix = GET_FILTERED_MATRIX(MATRIX, min_value,max_value)
        %   Delete elements in the matrix that are outside the range
        %   (min_value, max_value) and convert the rest to 1's.
        %
            
            assert(min_value < max_value);
            matrix = 1.0 * (matrix<max_value & matrix>min_value);
            
        end
        
        function matrix_categ = CONVERT_FROM_QUANTITATIVE_TO_CATEGORICAL(matrix,margs_or_ncategs)
        % NOT TESTED FUNCTION
        
            min_value = min(matrix(:)); max_value = max(matrix(:));
            
            if(length(margs_or_ncategs) == 1)
                n_categs = margs_or_ncategs;
                margs = linspace(min_value,max_value,n_categs+1);
            else
                n_categs = length(margs_or_ncategs)-1;
                margs = margs_or_ncategs;
            end
            
            matrix_categ = zeros(size(matrix));
            
            for i = 1:n_categs
                
                matrix_categ = matrix_categ + i * MatrixFunctions.GET_FILTERED_MATRIX(matrix,margs(i),margs(i+1));
                
            end
            
        end
        
    end
end