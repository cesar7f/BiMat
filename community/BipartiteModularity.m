% BipartiteModularity - Modularity parent class. 
% This class is the main part of all three bipartite modularity codes.
% It contains all the shared methods of the modularity algorithms.
% The way this class works is by isolating the graph (matrix) in
% independent connected components first, and find the best modularity
% configuration for each component. However, the modularity still try
% to maximize the global modularity.
%
% BipartiteModularity Properties:
%    matrix - Boolean Bipartite adjacency matrix
%    rr - Row community matrix. Size = n_rows*N
%    tt - Column community matrix. Size = n_cols*N
%    rr_sorted - rr sorted by community size (=rr(index_rows))
%    tt_sorted - tt sorted by community size (=tt(index_cols))
%    n_rows - Number of rows
%    n_cols - Number of columns
%    n_edges - Number of edges
%    bb - Original - Null (matrix - kidj/num_edges)
%    index_rows - Register of the swaps in Rows.
%    index_cols - Register of the swaps in Cols.
%    trials - Number of trials to find the best modularity configuration
%    Qb - Standard Modularity function (1/num_edges) * Trace(rr' * bb * tt)
%    Qr - Percentage of interactions inside modules
%    N  - Number of modules
%    row_modules - Module index for rows
%    col_modules - Module index for columns
%    done - The algorithm has been performed
%    webmatrix - Same than matrix but not neccesearly boolean.
%    optimize_by_component - Otimize modularity by component, nor by the entire network.
%    print_results - Flag to indicate if result output will be generated
%
% BipartiteModularity Methods:
%    Detect - Main modularity detection method
%    ExtractCommunityMatrices - Extract all community adjacency matrices
%    ExtractCommunityIndexes - Extract row and column community indexes
%    ExtractCommunityModules - Extract all communities as Bipartite objects
%    Print - Print modularity information
%    CALCULATE_MODULARITY_MATRIX - Calculate the modularity matrix of a bipartite network
%    CALCULATE_Qb_VALUE - Calculate the standard modularity value
%    CALCULATE_Qr_VALUE - Calculate the ratio of internal vs external interactions.
%    BRIM - Perform the standar BRIM algorithm in a set of the corresponding matrices
%    ADAPTIVE_BRIM - Calculate the modularity using the Adaptive Brim algorithm
%    LP_BRIM - Calculate the modularity using the Adaptive Brim algorithm
%    LEADING_EIGENVECTOR - Calculate the modularity using the Leading eigenvector algorithm
%
% See also:
%    AdaptiveBrim, LeadingEigenvector, and LPBrim
classdef BipartiteModularity < handle

    properties
        Qb                   = 0;   %Standard Modularity function (1/num_edges) * Trace(rr' * bb * tt)
        Qr                   = 0;   %Percentage of interactions inside modules
        N                    = 0;   %Number of modules
        matrix               = [];  %Boolean Bipartite adjacency matrix
        row_modules          = [];  %Module index for rows
        col_modules          = [];  %Module index for columns
        bb                   = [];  %Original - Null (matrix - kidj/num_edges)
        n_rows               = 0;   %Number of rows
        n_cols               = 0;   %Number of columns
        n_edges              = 0;   %Number of edges
        index_rows           = [];  %Register of the swaps in Rows.
        index_cols           = [];  %Register of the swaps in Cols.
        trials               = Options.TRIALS_MODULARITY; %Number of trials to find the best modularity configuration
        done                 = 0;   %The algorithm has been performed
        optimize_by_component= 0;   %Otimize modularity by component, nor by the entire network.
        print_results        = Options.PRINT_RESULTS % Flag to indicate if result output will be generated
    end
    
    properties(Access = protected) %Used for calculating modularity in a single graph component
        webmatrix            = [];  %Same than matrix but not neccesearly boolean.
        rr                   = [];  %Row community matrix. Size = n_rows*N
        tt                   = [];  %Column community matrix. Size = n_cols*N
        rr_sorted            = [];  %$$rr$$ sorted by community size (=rr(index_rows))
        tt_sorted            = [];  %tt sorted by community size (=tt(index_cols))
        matrix_component      = [];
        bb_component          = [];
        rr_component          = [];
        tt_component          = [];
        row_modules_component = [];
        col_modules_component = [];
        n_rows_component      = [];
        n_cols_component      = [];
        n_edges_component     = [];
        N_component           = [];
    end
        
    methods(Abstract, Access = 'protected')
        
        % DetectComponent - Abstract method to be implemented in all
        %    BipartiteModularity son classes
        % See AdaptiveBrim, LeadingEigenvector, and LPBrim
        obj = DetectComponent(obj);

    end
    
    methods(Access= 'protected')
        
        
        function obj = BipartiteModularity(bipmatrix)
        % BipartiteModularity(bipmatrix) - Main constructor
        %   It contains all the shared operations of the bipartite modularity algorithm
        %   constructors. Can be called only from a son class.
        %
        % See AdaptiveBrim, LeadingEigenvector, and LPBrim
        
            obj.webmatrix = bipmatrix;
            obj.matrix = bipmatrix > 0;
            [obj.n_rows, obj.n_cols] = size(obj.matrix);
            
            obj.n_edges = sum(sum(obj.matrix));
            
        end
        
        function obj = AssignSingleModule(obj)
        % AssignSingleModule(obj)
        % Assign all the nodes to a single module.
        
            obj.N_component = 1;
            obj.rr_component = ones(obj.n_rows_component,1);
            obj.tt_component = ones(obj.n_cols_component,1);
            obj.row_modules_component = ones(1,obj.n_rows_component);
            obj.col_modules_component = ones(1,obj.n_cols_component);
            
        end
        
        function obj = SortModules(obj)
        % obj = SortModules(obj)
        % Method using for finding the adecuate sorting of rows and columns
        % according to the module that they belong to. This method is
        % mainly used for plotting modular graph and matrix plots.
            
            %Sort by module row or column size.
            if(obj.n_rows > obj.n_cols)
                [~,idx] = sort(sum(obj.rr),'descend');
            else
                [~,idx] = sort(sum(obj.tt),'descend');
            end
            
            %Initial indexing is according in ascending order 1,2,3,...
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            
            %Sort module index matrices
            obj.rr = obj.rr(:,idx);
            obj.tt = obj.tt(:,idx);
            
            %For each module, sort according to a nested configuration.
            for i = 1:obj.N
                idx_rows = obj.rr(:,i);
                idx_cols = obj.tt(:,i);
                
                if(sum(sum(obj.matrix(idx_rows==1,:),2)==0) || sum(sum(obj.matrix(:,idx_cols==1),1))==0)
                     rr_temp = obj.rr(:,i);
                    tt_temp = obj.tt(:,i);
                    for j = i:obj.N-1
                        obj.rr(:,j) = obj.rr(:,j+1);
                        obj.tt(:,j) = obj.tt(:,j+1);
                    end
                    obj.rr(:,obj.N) = rr_temp;
                    obj.tt(:,obj.N) = tt_temp;
                    break;
                end
            end
            
            sorted_matrix = obj.matrix;
            
            row_global = zeros(obj.n_rows,1); col_global = zeros(obj.n_cols,1);
            %For each module, sort according to a nested configuration.
            i_r = 0; i_c = 0;
            for i = 1:obj.N
               
                row_i = find(obj.rr(:,i));
                col_i = find(obj.tt(:,i));
                
                %[~,row_loc] = sort(sum(sorted_matrix(row_i,col_i),2),'descend');
                %[~,col_loc] = sort(sum(sorted_matrix(row_i,col_i),1),'ascend');
                
                [~,row_loc] = sort(sum(10000*sorted_matrix(row_i,col_i),2)+sum(sorted_matrix(row_i,:),2),'descend');
                [~,col_loc] = sort(sum(10000*sorted_matrix(row_i,col_i),1)+sum(sorted_matrix(:,col_i),1),'ascend');
                
                row_global(i_r+1 : i_r+length(row_i)) = row_i(row_loc);
                col_global(i_c+1 : i_c+length(col_i)) = col_i(col_loc);
                i_r = i_r+length(row_i); i_c = i_c+length(col_i);
                
            end
           
            try
                col_global = flipud(col_global);
                obj.rr_sorted = obj.rr(row_global,:);
                obj.tt_sorted = obj.tt(col_global,:);
            catch
                error('Error in BipartiteModularity. Please report to the main author');
            end
            
            obj.index_rows = row_global;
            obj.index_cols = col_global;
              
            [row,col] = find(obj.rr);[~,ix] = sort(row);obj.row_modules = col(ix);
            [row,col] = find(obj.tt);[~,ix] = sort(row);obj.col_modules = col(ix);
            
        end
    
    end
    
    methods
        
        
        function obj = Detect(obj,ntrials)
        % Detect - Main modularity detection method
        %
        %   obj = Detect(obj) - This method calculate the modularity by first dividing the
        %   network (matrix) in isolated components. The modularity is
        %   optimized in each component separatly. The modularity is
        %   optimized globally by default, for optimizing components
        %   independently of the rest change property optimize_by_component
        %   to true before calling this method.
        %
        %   obj = Detect(obj,ntrials) - Same than the previous but instead
        %   of using the default number of restarts, the method will use
        %   ntrials as this number.
            
            if(isempty(obj.matrix))
                obj.Qb = NaN;
                obj.Qr = NaN;
                return;
            end
        
            if(nargin == 2)
                obj.trials = ntrials;
            end
            
            [obj.n_rows, obj.n_cols] = size(obj.matrix);
            
            %Divide the network in isolated graph components.
            %mm = [zeros(obj.n_rows, obj.n_rows) obj.matrix; obj.matrix' zeros(obj.n_cols, obj.n_cols)];
            %S = number of components, C = component indexes.
            %[S, C] = graphconncomp(sparse(mm),'Weak', true);
            
            %Assign component indexes to rows and column nodes.
            [C_rows, C_cols, n_comp] = MatrixFunctions.ISOLATED_COMPONENTS(obj.matrix);
            %C_rows = C(1:obj.n_rows);
            %C_cols = C(1+obj.n_rows:obj.n_rows+obj.n_cols);
            
            %Single nodes will go to the same empty module
            rows_empty = setdiff(C_rows,C_cols);
            cols_empty = setdiff(C_cols,C_rows);
            rows_empty = arrayfun(@(x) find(C_rows==x), rows_empty);
            cols_empty = arrayfun(@(x) find(C_cols==x), cols_empty);
            
            %Community indexes for rows (rr) and columns(tt)
            obj.rr = zeros(obj.n_rows - length(rows_empty));
            obj.tt = zeros(obj.n_cols - length(cols_empty));
            row_modules_global = zeros(obj.n_rows,1);
            col_modules_global = zeros(obj.n_cols,1);
            obj.N = 0;

            
            %For each component in the graph find the best modularity
            for i = 1:n_comp
               
                %Locate the nodes that belongs to component i
                idx_rows = find(C_rows==i);
                idx_cols = find(C_cols==i);
                
                %These nodes are already included in the empty module
                if(isempty(idx_rows) || isempty(idx_cols))
                    continue;
                end
                
                %Extract the component matrix from the entire matrix.
                obj.matrix_component = obj.matrix(idx_rows,idx_cols);
                [obj.n_rows_component, obj.n_cols_component] = size(obj.matrix_component);
                %Calculate the modularity matrix for the extracted matrix
                if(obj.optimize_by_component)
                    obj.bb_component = BipartiteModularity.CALCULATE_MODULARITY_MATRIX(obj.matrix_component,sum(obj.matrix_component(:)));
                else
                    obj.bb_component = BipartiteModularity.CALCULATE_MODULARITY_MATRIX(obj.matrix_component,obj.n_edges);
                end
                %If the matrix is not fully connected, proced to find the
                %best modularity distribution
                if(sum(sum(obj.matrix_component))~=numel(obj.matrix_component))
                    obj.DetectComponent();
                    %Clean empty columns in the community matrices.
                    [obj.rr_component, obj.tt_component] = BipartiteModularity.CLEAN_MODULE_MATRICES(obj.rr_component,obj.tt_component);
                    obj.N_component = size(obj.rr_component,2);
                    %Assign each row and column in the matrix to its corresponding
                    %modules.
                    obj.row_modules_component = BipartiteModularity.ASSIGN_MODULES(obj.rr_component);
                    obj.col_modules_component = BipartiteModularity.ASSIGN_MODULES(obj.tt_component);
                %otherwise assign the nodes to a single module
                else
                    obj.AssignSingleModule();
                end
                
                %Expand the global community indexes by the findings in the
                %current component.
                try
                obj.rr(idx_rows, obj.N+1:obj.N+obj.N_component) = obj.rr_component;
                obj.tt(idx_cols, obj.N+1:obj.N+obj.N_component) = obj.tt_component;
                row_modules_global(idx_rows) = obj.row_modules_component+obj.N;
                col_modules_global(idx_cols) = obj.col_modules_component+obj.N;
                obj.N = obj.N + obj.N_component;
                catch err
                    error('Error in BipartiteModularity. Please report to the main author');
                end
            end

            obj.row_modules = row_modules_global;
            obj.col_modules = col_modules_global;
            
            if(~isempty(rows_empty) || ~isempty(cols_empty))
                obj.rr(rows_empty,obj.N+1) = 1;
                obj.tt(cols_empty,obj.N+1) = 1;
                obj.N = obj.N+1;
                obj.row_modules(idx_rows) = obj.N;
                obj.col_modules(idx_cols) = obj.N;
            end
            obj.rr(:,obj.N+1:end) = [];
            obj.tt(:,obj.N+1:end) = [];
            
            %[obj.n_rows obj.n_cols] = size(obj.matrix);
            %obj.n_edges = sum(obj.matrix(:));
            obj.bb = BipartiteModularity.CALCULATE_MODULARITY_MATRIX(obj.matrix,obj.n_edges);
            obj.Qb = BipartiteModularity.CALCULATE_Qb_VALUE(obj.rr,obj.bb,obj.tt,obj.n_edges);
            obj.Qr = BipartiteModularity.CALCULATE_Qr_VALUE(obj.matrix,obj.rr,obj.tt);
            obj.SortModules();
            
            obj.done = 1;
            
            if(obj.print_results)
                obj.Print();
            end
            
        end
        
       
        function matrices = ExtractCommunityMatrices(obj)
        % ExtractCommunityMatrices - Extract all community adjacency
        % matrices
        %
        %   MATRICES = ExtractCommunityMatrices(obj)
        %   Extract all the bipartite adjacency matrices corresponding to
        %   all detected communities and return them in MATRICES, which is
        %   a cell array compossed of N matrices, where N is the number of
        %   modules of the current network. The rows and columns of each
        %   matrix will correspond to each module.
        
            matrices = cell(obj.N,1);
            for i = 1:obj.N
                idx_rows = obj.rr(:,i)==1;
                idx_cols = obj.tt(:,i)==1;
                matrices{i} = obj.matrix(idx_rows,idx_cols);
            end 
        end
        
        function [module_rows, module_cols] = ExtractCommunityIndexes(obj)
        % ExtractCommunityIndexes - Extract row and column community
        % indexes
        %
        %   [module_rows module_cols] = ExtractCommunityIndexes(obj)
        %   This method will return two cell arrays with N matrices, where N is
        %   the number of modules. Each cell array will have a the index
        %   number (row/column id) of the row/columns that belongs to module
        %   i.
        
            module_rows = cell(obj.N,1);
            module_cols = cell(obj.N,1);
            for i = 1:obj.N
                
                %Find the row and column indexes that belongs to module i
                idx_rows = find(obj.rr(:,i)==1);
                idx_cols = find(obj.tt(:,i)==1);
                
                module_rows{i} = idx_rows;
                module_cols{i} = idx_cols;
                
            end
            
        end 
        
        function q = TestRowsContribution(obj,row_ids)
        %Untested function
            q = trace(obj.rr(row_ids,:)' * obj.bb(row_ids,:) * obj.tt) / obj.n_edges;
        end
        
        function networks = ExtractCommunityModules(obj)
        % ExtractCommunityModules - Extract all communities as Bipartite
        % objects
        %
        %    networks = ExtractCommunityModules(obj)
        %    This method will return a cell array with N Bipartite objects, where N is
        %    the number of modules. Each bipartite object i will be composed of only the
        %    row and column nodes that belong to the module i.
        
            networks = cell(obj.N, 1);
            for i = 1:obj.N
                idx_rows = obj.rr(:,i)==1;
                idx_cols = obj.tt(:,i)==1;
                ma = obj.webmatrix(idx_rows,idx_cols);
                networks{i} = Bipartite(ma);
            end 
        end
        
        function str = Print(obj,filename)
        % Print - Print modularity information
        %
        %   STR = Print(obj) Print the modularity information to screen and
        %   return this information to the string STR
        %
        %   STR = Print(obj, FILE) Print the modularity information to screen and
        %   text file FILE and return this information to the string STR
            
            str = 'Modularity:\n';
            str = [str, '\tUsed algorithm:             \t', sprintf('%20s',class(obj)), '\n'];
            str = [str, '\tN (Number of modules):      \t', sprintf('%20i',obj.N), '\n'];
            str = [str, '\tQb (Standard metric):       \t', sprintf('%20.4f',obj.Qb), '\n'];
            str = [str, '\tQr (Ratio of int/ext inter):\t', sprintf('%20.4f',obj.Qr), '\n'];
            
            fprintf(str);  
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
    end
    
    methods(Static)
        
        function bb = CALCULATE_MODULARITY_MATRIX(matrix,n_edges)
        % CALCULATE_MODULARITY_MATRIX - Calculate the modularity matrix of
        % a bipartite network
        %
        %   bb = CALCULATE_MODULARITY_MATRIX(matrix,n_edges)
        %   Return the modularity matrix:
        %     bb_ij = matrix_ij - k_i d_j / n_edges
        %   where k_i and d_j are the degrees of rows and columns (sum in
        %   rows and columns of matrix), and n_edges is the number of edges
        %   or total sum of matrix
        
            if all(all(matrix == 0))
                bb = matrix;
            else
                coldeg = sum(matrix, 1);
                rowdeg = sum(matrix, 2);
                %obj.n_edges = sum(rowdeg);
                bb = matrix - (1/n_edges) * rowdeg * coldeg;
            end
        end
        
        function Qb = CALCULATE_Qb_VALUE(rr,bb,tt,n_edges)
        % CALCULATE_Qb_VALUE - Calculate the standard modularity value
        %
        %   Qb = CALCULATE_Qb_VALUE(rr,bb,tt,n_edges)
        %   Calculaate the standard bipartite modularity
        %     Qb = (1/n_edges) Tr (rr' bb tt)
        %   Value in the range 0-1 with a maximal possible value of
        %   1-fill(or connectance).    
            Qb = trace(rr' * bb * tt) / n_edges;
            
        end
        
        
        function Qr= CALCULATE_Qr_VALUE(matrix,rr,tt)
        % CALCULATE_Qr_VALUE - Calculate the ratio of internal vs external
        % interactions.
        %
        %   Qr = CALCULATE_Qr_VALUE(matrix,rr,tt) Calculate the ratio of
        %   internal vs external interactions as a way of a posteriory
        %   modularity metric.
        
            Qr = 0;
            n_modules = size(rr,2);
            m_edges = sum(matrix(:));
            for i = 1:n_modules
                row_index = find(rr(:,i));
                col_index = find(tt(:,i));
                nr = length(row_index);
                nc = length(col_index);
                
                for j = 1:nr
                    for k = 1:nc
                        if(matrix(row_index(j),col_index(k)) > 0)
                            Qr = Qr + 1;
                        end
                    end
                end
            end
            
            Qr = (2*Qr/m_edges) - 1;
        end
        
        function module_idx = ASSIGN_MODULES(module_matrix)
        % ASSIGN_MODULES - Assign modules to each node given a module
        % matrix (R or T).
        %
        %   module_idx = ASSIGN_MODULES(module_matrix) Create a vector with
        %   the indices modules given a module matrix module_matrix.
        
            [a, b] = ind2sub(size(module_matrix), find(module_matrix));
            [~, sortv] = sort(a);
            module_idx = b(sortv);
            
        end
        
        function [rr, tt] = CLEAN_MODULE_MATRICES(rr,tt)
        % CLEAN_MODULE_MATRICES - Clean module matrices
        %   [rr tt] = CLEAN_MODULE_MATRICES(rr,tt) Clean empty columns and
        %   rows in the module matrices. Such that no extra space is used.
        
            
            com1 = find(any(rr));
            com2 = find(any(tt));
            
            coms = union(com1,com2);
            
            rr = rr(:,coms);
            tt = tt(:,coms);            
        end
        
        function [rr, tt, preQ] = BRIM(rr,bb,tt,n_edges)
        % BRIM - Perform the standar BRIM algorithm in a set of the
        % corresponding matrices
        %
        %   [rr tt preQ] = BRIM(rr,bb,tt,n_edges) Perform the BRIM modularity
        %   algorithm in the module matrices rr, tt using the modularity
        %   matrix bb (bb = adjacency - null model) and the corresponding
        %   number of edges. It returns the updated rr, tt and the biggest Q
        %   value found in preQ.
            
            inducingBlueFlag = 1;
            
            preQ = BipartiteModularity.CALCULATE_Qb_VALUE(rr,bb,tt,n_edges);
            while(1)
               
                %Assign blug nodes
                if(inducingBlueFlag)
                    [~, maxind] = max(bb'*rr, [], 2);
                    modules = eye(size(rr,2));
                    tt = modules(maxind, :);
                else %Assign red nodes
                    [~, maxind] = max(bb*tt, [], 2);
                    modules = eye(size(tt,2));
                    rr = modules(maxind, :);
                end
                
                newQ = BipartiteModularity.CALCULATE_Qb_VALUE(rr,bb,tt,n_edges);
                
                if(newQ <= preQ)
                    break;
                else
                    preQ = newQ;
                    inducingBlueFlag = ~inducingBlueFlag;
                end
                
            end
        end
        
        function modul = ADAPTIVE_BRIM(matrix)
        % ADAPTIVE_BRIM - Calculate the modularity using the Adaptive Brim
        % algorithm
        %
        %   modul = ADAPTIVE_BRIM(matrix) Calculate the modularity using
        %   the Adaptive Brim algorithm, plot the basic information to
        %   screen and return an AdaptiveBrim object that contains such
        %   information in modul.
            
            modul = AdaptiveBrim(matrix);
            modul.print_results = false;
            modul.Detect();
            modul.Print();
            
        end
        
        function modul = LP_BRIM(matrix)
        % LP_BRIM - Calculate the modularity using the Adaptive Brim
        % algorithm
        %
        %   modul = LP_BRIM(matrix) Calculate the modularity using
        %   the LP&BRIM algorithm, plot the basic information to
        %   screen and return an AdaptiveBrim object that contains such
        %   information in modul.    
        
            modul = LP_BRIM(matrix);
            modul.print_results = false;
            modul.Detect();
            modul.Print();
            
        end

        function modul = LEADING_EIGENVECTOR(matrix)
        % LEADING_EIGENVECTOR - Calculate the modularity using the Leading
        % eigenvector algorithm
        %
        %   modul = LEADING_EIGENVECTOR(matrix) Calculate the modularity using
        %   the LP&BRIM algorithm, plot the basic information to
        %   screen and return an AdaptiveBrim object that contains such
        %   information in modul. 
        
            modul = LeadingEigenvector(matrix);
            modul.print_results = false;
            modul.Detect();
            modul.Print();
            
        end
        
    end
end  
