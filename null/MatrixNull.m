classdef MatrixNull < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = MatrixNull

        end
    end 
   
    methods(Static)
        
        function [matrix fill_sites empty_sites] = NESTED_TO_RANDOM(n_size, n_links)
        %[matrix fill_sites empty_sites] = NESTED_TO_RANDOM(n_size, n_links)%
            n = n_size;
            matrix = flipud(tril(ones(n)));
            fill_sites = find(matrix-flipud(eye(n)));
            empty_sites = find(matrix==0);

            n_fill = length(fill_sites); n_empty = length(empty_sites);

            random_fill_switched = randsample(fill_sites,min(n_fill,n_links),0);
            random_empty_switched = randsample(empty_sites,min(n_empty,n_links),0);

            matrix(random_fill_switched) = 0;
            matrix(random_empty_switched) = 1;
            %imagesc(matrix);
        end
        
        function [matrix nested_sites modular_sites] = NESTED_TO_MODULAR(n_size, n_links)
           
            n = n_size;
            matrix = flipud(tril(ones(n)));
            nested_sites = [];
            modular_sites = [];
            for i = 1:n/2
                %matrix( (i-1)*n + n/2 + 1 : n*i-i+1) = 1;
                %matrix( n/2*n + (i-1)*n + 1 : n/2*n + (i-1)*n + 1 + n/2 - i ) = 1;
                %matrix( n/2*n + (i-1)*n + n/2 : n/2*n + (i-1)*n + n) = 1;
                nested_sites = [nested_sites (i-1)*n + n/2 + 1 : n*i-i+1];
                nested_sites = [nested_sites n/2*n + (i-1)*n + 1 : n/2*n + (i-1)*n + 1 + n/2 - i];
                modular_sites = [modular_sites n/2*n + (i-1)*n + n/2 + 1: n/2*n + (i-1)*n + n ];
            end
            
            n_nest = length(nested_sites); n_modul = length(modular_sites);
            
            random_nest_switched = randsample(nested_sites,min(n_nest,n_links),0);
            random_modul_switched = randsample(modular_sites,min(n_modul,n_links),0);
            
            matrix(random_nest_switched) = 0;
            matrix(random_modul_switched) = 1;
            
            %imagesc(matrix);
            
        end
        
        function matrix = NESTED_NULL(n_size,p)
           
            n = n_size;
            
            matrix = p*flipud(tril(ones(n))) + (1-p)* fliplr(tril(ones(n),-1));
            matrix = 1.0 *(rand(n,n) < matrix);
            
        end
        
        function [matrix matrix_p] = BRADFORD(n_size,p)
           
            n = n_size;
            matrix_p = p*blkdiag(ones(floor([n/2 n/2])), ones(ceil([n/2 n/2])));
            matrix_p = (1-p)*fliplr(blkdiag(ones(floor([n/2 n/2])), ones(ceil([n/2 n/2])))) + matrix_p;
            matrix = 1.0*(rand(n,n) < matrix_p);
            
        end
        
        function [matrix matrix_p] = GABRIEL_NULL(n_size, p)
           
            n = n_size;
            matrix = flipud(tril(ones(n)));
            nested_sites = [];
            modular_sites = [];
            for i = 1:n/2
                nested_sites = [nested_sites (i-1)*n + n/2 + 1 : n*i-i+1];
                nested_sites = [nested_sites n/2*n + (i-1)*n + 1 : n/2*n + (i-1)*n + 1 + n/2 - i];
                modular_sites = [modular_sites n/2*n + (i-1)*n + n/2 + 1: n/2*n + (i-1)*n + n ];
            end
            
            matrix_p = matrix;            
            matrix(nested_sites) = 1.0*(rand(size(matrix(nested_sites)))<p);
            matrix(modular_sites) = 1.0*(rand(size(matrix(modular_sites)))<(1-p));
            
            matrix_p(nested_sites) = p;
            matrix_p(modular_sites) = 1-p;
            
        end
        
        function [matrix nested_sites modular_sites] = GABRIEL_NULL_P(n_size, p)
           
            n = n_size;
            matrix = flipud(tril(ones(n)));
            nested_sites = [];
            modular_sites = [];
            for i = 1:n/2
                %matrix( (i-1)*n + n/2 + 1 : n*i-i+1) = 1;
                %matrix( n/2*n + (i-1)*n + 1 : n/2*n + (i-1)*n + 1 + n/2 - i ) = 1;
                %matrix( n/2*n + (i-1)*n + n/2 : n/2*n + (i-1)*n + n) = 1;
                nested_sites = [nested_sites (i-1)*n + n/2 + 1 : n*i-i+1];
                nested_sites = [nested_sites n/2*n + (i-1)*n + 1 : n/2*n + (i-1)*n + 1 + n/2 - i];
                modular_sites = [modular_sites n/2*n + (i-1)*n + n/2 + 1: n/2*n + (i-1)*n + n ];
            end
            
            n_nest = length(nested_sites); n_modul = length(modular_sites);
            
            %random_nest_switched = randsample(length(nested_sites),round(p*n_nest),0);
            %random_modul_switched = randsample(length(modular_sites),round(p*n_,0);
            
            %matrix(random_nest_switched) = 0;
            %matrix(random_modul_switched) = 1;
            
            
            %imagesc(matrix);
            
        end
        
        function matrix = SORT_MATRIX(matrix)
           
            [a ir] = sort(sum(matrix,2),'descend');
            [a ic] = sort(sum(matrix,1),'descend');
           
            matrix = matrix(ir,ic);
            
        end
        
        function [matrix index_row index_col] = RANDOM_SORT(matrix)
           
            [n_rows n_cols] = size(matrix);
            
            index_row = randperm(n_rows);
            index_col = randperm(n_cols);
            matrix = matrix(index_row,index_col);
            
        end
        
        function M_new = NESTED_TO_RANDOM_2(M,n_links)
        % M: matrix
        % n_links: elements to move
            if n_links >=0
                n = size(M);
                M_new = M;
                M_fill = M;
                M_empty = ~M;
                nestedM = logical(flipud(tril(ones(n)))-flipud(eye(n)));
                nestedM_diag = logical(flipud(tril(ones(n))));
                M_fill(~nestedM)=0;
                M_empty(nestedM_diag) =0;
                fill_sites_move = find(M_fill);
                empty_sites_move = find(M_empty);
                n_fill = sum(sum(M_fill)); n_empty = sum(sum(M_empty));

                random_fill_switched = randsample(fill_sites_move,min(n_fill,n_links),0);
                random_empty_switched = randsample(empty_sites_move,min(n_empty,n_links),0);

                M_new(random_fill_switched) = 0;
                M_new(random_empty_switched) = 1;
            else
                n_links = abs(n_links)
                M=~M;
                n = size(M);
                M_new = M;
                M_fill = M;
                M_empty = ~M;
                nestedM = logical(flipud(tril(ones(n)))-flipud(eye(n)));
                nestedM_diag = logical(flipud(tril(ones(n))));
                M_fill(~nestedM)=0;
                M_empty(nestedM_diag) =0;
                fill_sites_move = find(M_fill);
                empty_sites_move = find(M_empty);
                n_fill = sum(sum(M_fill)); n_empty = sum(sum(M_empty));

                random_fill_switched = randsample(fill_sites_move,min(n_fill,n_links),0);
                random_empty_switched = randsample(empty_sites_move,min(n_empty,n_links),0);

                M_new(random_fill_switched) = 0;
                M_new(random_empty_switched) = 1;
                M_new =~M_new;
            end
        end
        
        function [k d] = GET_DEGREES(matrix)
           
            k = sum(matrix,2);
            d = sum(matrix,1);
            
        end
        
        function matrix = BASCOMPTE_NULL(k_degs, d_degs)
           
            if(nargin == 1)
                [k_degs d_degs] = MatrixNull.GET_DEGREES(k_degs);
            end
            
            n_rows = length(k_degs); n_cols = length(d_degs);
            
            if(max(d_degs)> n_rows); error('Degrees on columns are bigger than the number of rows'); end;
            if(max(k_degs)> n_cols); error('Degrees on rows are bigger than the number of columns'); end;
            
            p_rows = k_degs/n_cols; p_cols = d_degs/n_rows;
            p_rows = reshape(p_rows,n_rows,1); p_cols = reshape(p_cols,1,n_cols);
            p = (repmat(p_rows,1,n_cols) + repmat(p_cols,n_rows,1))/2;
            
            matrix = rand(n_rows,n_cols)<p;
        end
        
        function matrix = SUM_MATRIX(matrices)
            
            n = length(matrices);
            matrix = 1.0*matrices{1};
            for i = 2:n
                matrix = matrix+matrices{i};
            end
            matrix = matrix/n;
            matrix = rand(size(matrix)) < matrix;
            
        end
        
        function p = FILL(matrix)
            p = sum(sum(matrix))/numel(matrix);
        end
        
        function [matrix ir ic] = NON_ZERO_MATRIX(matrix)
        
            ir = find(sum(matrix,2)>0);
            ic = find(sum(matrix,1)>0);
            
            matrix = matrix(ir,ic);
         
        end
        
        function [matrix ir ic] = TYPE_MATRIX(matrix, n_type)
           
            [n_rows n_cols] = size(matrix);
            ir = 1:n_rows; ic = 1:n_cols;
            %ir = find(sum(matrix==n_type,2)>0);
            %ic = find(sum(matrix==n_type,1)>0);
            
            %matrix = matrix(ir,ic);
            matrix = (matrix==n_type)*n_type;
            
        end
        
        function [matrix ir ic] = TYPE_MATRIX_NON_ZERO(matrix, n_type)
           
            [n_rows n_cols] = size(matrix);
            ir = 1:n_rows; ic = 1:n_cols;
            
            [matrix ir_type  ic_type] = MatrixNull.TYPE_MATRIX(matrix, n_type);
            [matrix ir_zero  ic_zero] = MatrixNull.NON_ZERO_MATRIX(matrix);
            
            ir = ir(ir_type(ir_zero));
            ic = ic(ic_type(ic_zero));
            
        end
        
        function A = BipartiteToUnipartite(matrix)
           
            [n_rows, n_cols] = size(matrix);
            A = [zeros(n_rows,n_rows),matrix;matrix',zeros(n_cols,n_cols)];
            
        end
        
        function matrix = NestedMatrix(n)
            
            matrix = flipud(tril(ones(n)));
            
        end
        
        function max_eig = GetBiggestEigenvalue(matrix)
           
            max_eig = max(eig(MatrixNull.BipartiteToUnipartite(matrix)));
            
        end
        
        function fill = GET_FILL(matrix)
            fill = sum(sum(matrix>0)) / numel(matrix);
        end
        
        function submatrix = RANDOM_SAMPLING(matrix, n_rows, n_cols)
           
            [m n] = size(matrix);
            assert(m>=n_rows); assert(n>=n_cols);
            
            submatrix = matrix(randsample(m, n_rows), randsample(n,n_cols));
            
        end
        
        function matrix = MatrixUnion(matrix_1,matrix_2)
            matrix = blkdiag(matrix_1,matrix_2);
        end
        
    end
end