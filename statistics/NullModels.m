% NullModels - Static class for creating bipartite random networks
%
% NullModels Methods:
%    EQUIPROBABLE - Create a random matrix with equiprobable probabilities
%    AVERAGE - Create a random matrix with average probability in rows and columns.
%    AVERAGE_ROWS - Create a random matrix with average probability in rows
%    AVERAGE_COLS - Create a random matrix with average probability in columns
%    FIXED - Create a random matrix that respect the exact degree sequence
%    NULL_MODEL - Create a cell of random matrices
%
% See also:
%    StatisticalTest, InternalStatistics, MetaStatistics
classdef NullModels < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = NullModels

        end
    end 
   
    methods(Static)
        
        function matrix = EQUIPROBABLE(MatrixOrSizeRows, n_cols, ones)
        % EQUIPROBABLE - Create a random matrix with equiprobable probabilities
        %
        %   random_matrix = EQUIPROBABLE(matrix) Create a random matrix with
        %   the same number of interactions (ones) than matrix, but randomly
        %   distributed across the entire matrix
        %
        %   random_matrix = EQUIPROBABLE(n_rows, n_cols, ones) Create a
        %   random matrix of size n_rows by n_cols with ones interactions
        %   that are randomly distributed across the matrix.
        
            if(nargin == 1)

                [n_rows, n_cols] = size(MatrixOrSizeRows);
                ones = sum(sum(MatrixOrSizeRows>0));

            else 
                 n_rows = MatrixOrSizeRows;
            end

            allow_empty_nodes = Options.ALLOW_ISOLATED_NODES;
            trials = Options.TRIALS_FOR_NON_EMPTY_NODES;
            while(1 && trials > 0)
                tmpr=rand(n_rows*n_cols,1);
                [~, tmps]=sort(tmpr);
                matrix=zeros(n_rows,n_cols);
                matrix(tmps(1:ones))=1;
                if(allow_empty_nodes == 1 || ~MatrixFunctions.HAS_EMPTY_NODES(matrix))
                    break;
                end
                trials = trials - 1;
                if(trials == 0)
                    fprintf(['Not possible to create a matrix with non isolated nodes.\n', ...
                    'The random matrix was created without this constraint instead.\n', ...
                    'Consider to modify Options.ALLOW_ISOLATED_NODES and/or Options.INCLUDE_EMPTY_NODES\n\n']);
                end
                
            end

        end
        
        function matrix = AVERAGE(n_rows, n_cols, probRows,probCols)
        % AVERAGE - Create a random matrix with average probability in rows
        % and columns.
        %
        %   random_matrix = AVERAGE(matrix) Create a random matrix of the same size than matrix,
        %   such that cell i,j have a probability of having an interaction of:
        %
        %      p_ij = 1/2 * ( sum(matrix(i,:)/size(matrix,2) + sum(matrix(:,j)/size(matrix,1))
        %
        %   random_matrix = AVERAGE(n_rows, n_cols, probRows,probCols) Create
        %   a random matrix of size n_rows by n_cols with the probability of having an interaction of:
        %
        %      p_ij = 1/2 * ( probRows(i) + probCols(j))
            
            if(nargin ==1)
                matrix = n_rows;
                [n_rows, n_cols] = size(matrix);
                probCols = sum(matrix)/n_rows;
                probRows = sum(matrix,2)/n_cols;
            end
            
            if(n_rows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            if(n_cols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            ss = size(probRows);
            if(ss(1) == 1); probRows = probRows'; end;
            
            pr = repmat(probRows,1,n_cols);
            pc = repmat(probCols,n_rows,1);
                      
            mat = (pr+pc)/2;
            
            allow_empty_nodes = Options.ALLOW_ISOLATED_NODES;
            trials = Options.TRIALS_FOR_NON_EMPTY_NODES;
            while(1)
                matrix = rand(n_rows,n_cols) <= mat;
                if(allow_empty_nodes == 1 || ~MatrixFunctions.HAS_EMPTY_NODES(matrix))
                    break;
                end
                trials = trials - 1;
                if(trials == 0)
                    fprintf(['Not possible to create a matrix with non isolated nodes.\n', ...
                    'The random matrix was created without this constraint instead.\n', ...
                    'Consider to modify Options.ALLOW_ISOLATED_NODES and/or Options.INCLUDE_EMPTY_NODES\n\n']);
                end
            end

        end
        
        function matrix = AVERAGE_ROWS(n_rows, n_cols, probRows)
        % AVERAGE_ROWS - Create a random matrix with average probability in
        % rows
        %
        %   random_matrix = AVERAGE_ROWS(matrix) Create a random matrix of 
        %   the same size than matrix, such that each cell i,j have a probability
        %   of having an interaction of:
        %
        %      p_ij = sum(matrix(i,:)/size(matrix,2)
        %
        %   random_matrix = AVERAGE(n_rows, n_cols, probRows) Create
        %   a random matrix of size n_rows x n_cols with the probability of having an interaction of:
        %
        %      p_ij = probRows(i) 
            
            if(nargin ==1)
                matrix = n_rows;
                [n_rows, n_cols] = size(matrix);
                probRows = sum(matrix,2)/n_cols;
            end
            
            if(n_rows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            
            allow_empty_nodes = Options.ALLOW_ISOLATED_NODES;
            trials = Options.TRIALS_FOR_NON_EMPTY_NODES;
            
            while(1)
                matrix = rand(n_rows, n_cols) < repmat(probRows,1,n_cols);
                if(allow_empty_nodes == 1 || ~MatrixFunctions.HAS_EMPTY_NODES(matrix))
                    break;
                end
                trials = trials - 1;
                if(trials == 0)
                    fprintf(['Not possible to create a matrix with non isolated nodes.\n', ...
                    'The random matrix was created without this constraint instead.\n', ...
                    'Consider to modify Options.ALLOW_ISOLATED_NODES and/or Options.INCLUDE_EMPTY_NODES\n\n']);
                end
            end            
        end
        
        function matrix = AVERAGE_COLS(n_rows, n_cols, probCols)
        % AVERAGE_COLS - Create a random matrix with average probability in
        % columns
        %
        %   random_matrix = AVERAGE_COLS(matrix) Create a random matrix of 
        %   the same size than matrix, such that each cell i,j have a probability
        %   of having an interaction of:
        %
        %      p_ij = sum(matrix(:,j)/size(matrix,1)
        %
        %   random_matrix = AVERAGE_COLS(n_rows, n_cols, probCols) Create
        %   a random matrix of size n_rows x n_cols with the probability of having an interaction of:
        %
        %      p_ij = probRows(j)
            if(nargin ==1)
                matrix = n_rows;
                [n_rows, n_cols] = size(matrix);
                probCols = sum(matrix)/n_rows;
            end
            
            if(n_cols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            allow_empty_nodes = Options.ALLOW_ISOLATED_NODES;
            trials = Options.TRIALS_FOR_NON_EMPTY_NODES;
            
            while(1)
                matrix = rand(n_rows, n_cols) < repmat(probCols,n_rows,1);
                if(allow_empty_nodes == 1 || ~MatrixFunctions.HAS_EMPTY_NODES(matrix))
                    break;
                end
                trials = trials - 1;
                if(trials == 0)
                    fprintf(['Not possible to create a matrix with non isolated nodes.\n', ...
                    'The random matrix was created without this constraint instead.\n', ...
                    'Consider to modify Options.ALLOW_ISOLATED_NODES and/or Options.INCLUDE_EMPTY_NODES\n\n']);
                end
            end   
        end
        
        function matrix = FIXED(matrix)
        % FIXED - Create a random matrix that respect the exact degree
        % sequence
        %
        %   random_matrix = FIXED(matrix) Create a random matrix that respect the exact degree
        %   sequence. In other words, each node will continue to have the
        %   exact same number of edges.
        
            assert(sum(sum(matrix==0 + matrix==1))==numel(matrix));

            num_edges = sum(matrix(:));

            
            n_swap_edges_factor = 100;
            n_swap_edges = n_swap_edges_factor * num_edges;

            [irows, icols] = ind2sub(size(matrix),find(matrix>0));
            [~,idx] = sort(irows);
            edges = [irows(idx), icols(idx)];
            for j = 1:n_swap_edges
                idx_to_swap = ceil(num_edges*rand(2,1));
                e_1 = idx_to_swap(1); e_2 = idx_to_swap(2);

                if(matrix(edges(e_1,1), edges(e_2,2)) == 1 || ...
                   matrix(edges(e_2,1), edges(e_1,2)) == 1)
                    continue;
                end

                matrix(edges(e_1,1),edges(e_1,2)) = 0;
                matrix(edges(e_2,1),edges(e_2,2)) = 0;

                edge_temp = edges(e_1,:);
                edges(e_1,2) = edges(e_2,2);
                edges(e_2,2) = edge_temp(2);

                matrix(edges(e_1,1),edges(e_1,2)) = 1;
                matrix(edges(e_2,1),edges(e_2,2)) = 1;
            end
        
        end
        
        function rmatrices = NULL_MODEL(adjmatrix,model,replicates)
        % NULL_MODEL - Create a cell of random matrices
        %
        %   rmatrices = NULL_MODEL(adjmatrix,model,replicates) Create a
        %   cell of random matrices using adjmatrix as empirical matrix
        %   from which the required properties will be extracted. model is
        %   the name of the null model function that will be used and
        %   replicates the number of random matrices that are required.
        %
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
            
            adjmatrix = adjmatrix ~= 0;
            
            n_edges = sum(adjmatrix(:));
            [n_rows, n_cols] = size(adjmatrix);
            p_cols = sum(adjmatrix,1)/n_rows;
            p_rows = sum(adjmatrix,2)/n_cols;
            %p = n_edges / (n_rows * n_cols);
            
            nullmodel = 1;
            if(strcmp(func2str(model),'NullModels.EQUIPROBABLE')); nullmodel = 1; end;
            if(strcmp(func2str(model),'NullModels.AVERAGE')); nullmodel = 2; end;
            if(strcmp(func2str(model),'NullModels.AVERAGE_ROWS')); nullmodel = 3; end;
            if(strcmp(func2str(model),'NullModels.AVERAGE_COLS')); nullmodel = 4; end;
            if(strcmp(func2str(model),'NullModels.FIXED')); nullmodel = 5; end;
            
            rmatrices = cell(1,replicates);
            
            for i = 1:replicates
                %fprintf('iteration %i\n',i);
                switch nullmodel
                    case 1
                        rmatrices{i} = NullModels.EQUIPROBABLE(n_rows,n_cols,n_edges);
                    case 2
                        rmatrices{i} = NullModels.AVERAGE(n_rows,n_cols,p_rows,p_cols);
                    case 3
                        rmatrices{i} = NullModels.AVERAGE_ROWS(n_rows,n_cols,p_rows);
                    case 4
                        rmatrices{i} = NullModels.AVERAGE_COLS(n_rows,n_cols,p_cols);
                    case 5
                        rmatrices{i} = NullModels.AVERAGE_COLS(n_rows,n_cols,p_cols);
                end
            end
            
        end
        
    end
end