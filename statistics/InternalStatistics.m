% InternalStatistics - Multi-scale statistical analysis for a bipartite complex
% network. This class performs modularity and nestedness analysis of the
% internal modules (after performing modularity algorithms) of the
% bipartite network. In addition, it can be used to calculate diversity
% index of the labels id's of rows and/or columns inside module to see if
% there exist a correlation between modularity labeling and type of nodes.
% For a nice example of how this class can be used see the next paper:
%
%   Cesar O. Flores, Sergi Valverde, and Joshua S. Weitz. 
%   Multi-scale structure and geographic drivers of cross-infection within marine bacteria and phages.
%   ISME Journal 2012
%
% InternalStatistics Properties:
%     bipweb - Bipartite network object in which the multi-scale analysis will be performed
%     row_diversity - Data structure for row labeling diversity
%     col_diversity - Data structure for col labeling diversity
%     row_class - Labeling of rows
%     col_class - Labeling of columns
%     module_networks - Cell composed of bipartite network objects which corresponds to the internal modules of biweb
%     meta_statistics - MetaStatistics object that is populated with all modules (communities) inside bipweb 
%     idx_to_focus_on - Module indexes of where to focus the analysis. If empty all modules are evaluated.
%
% InternalStatistics Methods:
%     InternalStatistics - Main Constructor
%     TestDiversityRows - Function used for calculating the diverstiy of row classes 
%     TestDiversityColumns - Function used for calculating the diverstiy of column classes
%
% See also:
%     BipartiteModularity, Nestedness
classdef InternalStatistics < handle
   
    properties(GetAccess = 'public', SetAccess = 'public')
        bipweb          = []; % Bipartite network object in which the multi-scale analysis will be performed
        row_diversity   = []; % Data structure for row labeling diversity
        col_diversity   = []; % Data structure for col labeling diversity
        row_class       = []; % Labeling of rows
        col_class       = []; % Labeling of columns
        module_networks = {}; % Cell composed of bipartite network objects which corresponds to the internal modules of biweb
        meta_statistics = {}; % MetaStatistics object that is populated with all modules (communities) inside bipweb
        idx_to_focus_on = []; % Module indexes of where to focus the analysis. If empty all modules are evaluated.
    end
    
    
    methods
        function obj = InternalStatistics(webbip)
        % InternalStatistics - Main Constructor    
            obj.bipweb = webbip;

        end

        function TestDiversityRows(obj,n_trials,row_class,index_function)
        % TestDiversityRows - Function used for calculating the diverstiy
        % of the row labeling inside internal modules. The diversity is
        % calculated also in random permutations. The comparison of the
        % diversity of the original labeling and the random permutations
        % can be used to test if a correlation may exist between module
        % id's and row labeling.  After calling this function, you can
        % refer to the results using the property InternalStatitics.row_diversity
        %
        %   obj = TestDiversityRows(obj) Calculate the statistics of
        %   diversity for row labeling using Options.REPLICATES random permutations, the
        %   original row labeling of the bipartite network object and the
        %   Simpson's diversity index.
        %
        %   obj = TestDiversityRows(obj,n_trials) Calculate the statistics of
        %   diversity for row labeling using n_trials random permutations, the
        %   original row labeling of the bipartite network object and the
        %   Simpson's diversity index.
        %
        %   obj = TestDiversityRows(obj,n_trials,rows_ids) Calculate the statistics of
        %   diversity for row labeling using n_trials random permutations, rows_ids
        %   as row labeling of the bipartite network object and the
        %   Simpson's diversity index.
        %
        %   obj = TestDiversityRows(obj,n_trials,rows_ids,index_function) Calculate the statistics of
        %   diversity for row labeling using n_trials random permutations, rows_ids
        %   as row labeling of the bipartite network object and the
        %   diversity function index_function
        %
        %  See also:
        %   Diversity.SHANNON_INDEX, Diversity.SIMPSON_INDEX
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
            
            
            if(nargin == 1)
                n_trials = Options.REPLICATES;
                row_class = obj.bipweb.row_class;
                index_function = @Diversity.SIMPSON_INDEX;
            elseif(nargin==2)
                row_class = obj.bipweb.row_class;
                index_function = @Diversity.SIMPSON_INDEX;
            end
            
            assert(length(row_class)==obj.bipweb.community.n_rows);
            
            modul = obj.bipweb.community;

            obj.row_class = row_class;
            
            [row_modules, ~] = modul.ExtractCommunityIndexes();
            
            n = length(row_class);
            for i = 1:modul.N
                
                real_value = index_function(row_class(row_modules{i}));
                
                nr = length(row_modules{i});
                random_values = zeros(n_trials,1);
                for j = 1:n_trials
                    row_idx_random = row_class(randperm(n));
                
                    random_values(j) = index_function(row_idx_random(1:nr));
                end
                sta_vals = StatisticalTest.GET_STATISTICAL_VALUES(real_value,random_values);
                
                diversity.value(i) = sta_vals.value;
                diversity.zscore(i) = sta_vals.zscore;
                diversity.percentile(i) = sta_vals.percentile;
            end
            
            obj.row_diversity = diversity;
            obj.row_diversity.diversity_index = index_function;
            obj.row_diversity.n_permutations = n_trials;
            
        end
        
        function TestDiversityColumns(obj,n_trials,col_class,index_function)
        % TestDiversityColumns - Function used for calculating the diverstiy
        % of the column labeling inside internal modules. The diversity is
        % calculated also in random permutations. The comparison of the
        % diversity of the original labeling and the random permutations
        % can be used to test if a correlation may exist between module
        % id's and column labeling. After calling this function, you can
        % refer to the results using the property InternalStatitics.col_diversity
        %
        %   obj = TestDiversityColumns(obj) Calculate the statistics of
        %   diversity for column labeling using Options.REPLICATES random permutations, the
        %   original column labeling of the bipartite network object and the
        %   Simpson's diversity index.
        %
        %   obj = TestDiversityColumns(obj,n_trials) Calculate the statistics of
        %   diversity for column labeling using n_trials random permutations, the
        %   original column labeling of the bipartite network object and the
        %   Simpson's diversity index.
        %
        %   obj = TestDiversityColumns(obj,n_trials,cols_ids) Calculate the statistics of
        %   diversity for column labeling using n_trials random permutations, cols_ids
        %   as column labeling of the bipartite network object and the
        %   Simpson's diversity index.
        %
        %   obj = TestDiversityColumns(obj,n_trials,cols_ids,index_function) Calculate the statistics of
        %   diversity for column labeling using n_trials random permutations, cols_ids
        %   as column labeling of the bipartite network object and the
        %   diversity function index_function
        %
        %  See also:
        %   Diversity.SHANNON_INDEX, Diversity.SIMPSON_INDEX
                    
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
            
            
            if(nargin == 1)
                n_trials = Options.REPLICATES;
                col_class = obj.bipweb.col_class;
                index_function = @Diversity.SIMPSON_INDEX;
            elseif(nargin==2)
                col_class = obj.bipweb.col_class;
                index_function = @Diversity.SIMPSON_INDEX;
            end
            
            assert(length(col_class)==obj.bipweb.community.n_cols);
            
            modul = obj.bipweb.community;
            
            obj.col_class = col_class;
            
            [~,col_modules] = modul.ExtractCommunityIndexes();
            
            n = length(col_class);
            for i = 1:modul.N
                
                real_value = index_function(col_class(col_modules{i}));
                %obj.realindex.cols(i) = index_function(cols_ids(col_modules{i}));
                nr = length(col_modules{i});
                random_values = zeros(n_trials,1);
                for j = 1:n_trials
                    cols_idx_random = col_class(randperm(n));
                    %randomdata(j,i) = index_function(cols_idx_random(1:nr));
                    random_values(j) = index_function(cols_idx_random(1:nr));
                end
                sta_vals = StatisticalTest.GET_STATISTICAL_VALUES(real_value,random_values);
                
                diversity.value(i) = sta_vals.value;
                diversity.zscore(i) = sta_vals.zscore;
                diversity.percentile(i) = sta_vals.percentile;
            end
            
            obj.col_diversity = diversity;
            obj.col_diversity.diversity_index = index_function;
            obj.col_diversity.n_permutations = n_trials;
            
        end
        
        function TestInternalModules(obj,n_trials,null_model)
        % TestInternalModules - Function that calculate the modularity and
        % nestedness statistics of the internal modules of the bipartite
        % object obj.bipweb
        %
        %   obj = TestInternalModules(obj) It test for modularity, and
        %   nestedness the internal modules of the object obj.bipweb using
        %   Options.REPLICATES random matrices for each module and the EQUIPROBABLE null
        %   model.
        %
        %   obj = TestInternalModules(obj,n_trials) It test for modularity, and
        %   nestedness the internal modules of the object obj.bipweb using
        %   n_trials random matrices for each module and the EQUIPROBABLE null
        %   model.
        %
        %   obj = TestInternalModules(obj,n_trials,null_model) It test for modularity, and
        %   nestedness the internal modules of the object obj.bipweb using
        %   n_trials random matrices for each module and null_model as null
        %   model for creating the random matrices
            
            
            if(nargin == 1)
                n_trials = Options.REPLICATES;
                null_model = Options.DEFAULT_NULL_MODEL;
            elseif(nargin == 2)
                null_model = Options.DEFAULT_NULL_MODEL;
            end
            
            module = obj.bipweb.community;
            obj.module_networks = module.ExtractCommunityModules();
            
            if(~isempty(obj.idx_to_focus_on))
                obj.module_networks = obj.module_networks(obj.idx_to_focus_on);
            end
            
            gp = MetaStatistics(obj.module_networks);
            gp.replicates = n_trials;
            gp.null_model = null_model;
            gp.modularity_algorithm = str2func(class(obj.bipweb.community));
            gp.nestedness_algorithm = str2func(class(obj.bipweb.nestedness));
            
            gp.do_community = 1; % Perform Modularity analysis (default)
            gp.do_nestedness = 1; % Perform Nestedness analysis
            gp.DoMetaAnalyisis(); % Perform the analysis.
            
            obj.meta_statistics = gp;
            
        end
        
        
        
        
    end
    
    methods
       
%         function value = get.module_networks(obj)
%         % Just a get function to run modularity in case some calls this
%         % object.
%             
%             
%             if(obj.bipweb.community.done == 0)
%                 obj.bipweb.community.Detect();
%             end
%             value = obj.bipweb.community.ExtractCommunityModules();
%             
%         end
        
    end
    
end
% 
% THE NEXT FUNCTIONS WERE COMMENTED BECAUSE THEY WERE NOT TESTED
% EXAUSTIVELLY. THE USER IS FREE TO INCLUDE THEM IF HE WANTS.

% function networks = ExtractCommunityModules(obj)
%             modul = obj.bipweb.community;
%             networks = cell(modul.N,1);
%             for i = 1:modul.N
%                 idx_rows = find(modul.rr(:,i)==1);
%                 idx_cols = find(modul.tt(:,i)==1);
%                 ma = modul.webmatrix(idx_rows,idx_cols);
%                 networks{i} = Bipartite(ma);
%                 networks{i}.row_labels = obj.bipweb.row_labels(idx_rows);
%                 networks{i}.col_labels = obj.bipweb.col_labels(idx_cols);
%                 if(~isempty(obj.bipweb.row_class))
%                     networks{i}.row_class = obj.bipweb.row_class(idx_rows);
%                 end
%                 if(~isempty(obj.bipweb.col_class))
%                     networks{i}.col_class = obj.bipweb.col_class(idx_cols);
%                 end
%             end 
%         end
%         
%         function CreateModuleTableProperties(obj)
%             
%             modules = obj.module_networks;
%             nn = length(modules);
%             
%             fprintf('No. \t H \t P \t S \t I \t M \t C \t Lh \t Lp\n');
%             
%             for i = 1:nn
%                 H = modules{i}.n_rows;
%                 P = modules{i}.n_cols;
%                 S = H+P;
%                 I = sum(sum(modules{i}.matrix));
%                 M = H*P;
%                 C = I/M;
%                 Lh = I/H;
%                 Lp = I/P;
%                 fprintf('%i \t %i \t %i \t %i \t %i \t %i \t %f \t %f \t %f\n', ...
%                     i, H, P, S, I, M, C, Lh, Lp);
%             end
%             
%         end
%         
%         function nodf = TestNODFModuleContributions(obj)
%             matrix = obj.bipweb.matrix;
%             rr = obj.bipweb.community.rr;
%             tt = obj.bipweb.community.tt;
%             
%             [n_rows n_cols] = size(matrix);
%             
%             rows_sum = sum(matrix,2);
%             col_sum = sum(matrix,1);
%             
%             [row,col] = find(rr);[~,ix] = sort(row);row_modules = col(ix);
%             [row,col] = find(tt);[~,ix] = sort(row);col_modules = col(ix);
%             
%             nij_row_in = 0;
%             nij_row_out = 0;
%             for i = 1:n_rows
%                 for j = i+1:n_rows
%                     cont = sum(matrix(i,:).*matrix(j,:))*(rows_sum(i)~=rows_sum(j))/min(sum(matrix(i,:)),sum(matrix(j,:)));
%                     if(~isnan(cont))
%                         if(row_modules(i) == row_modules(j))
%                             nij_row_in = nij_row_in + cont;
%                         else
%                             nij_row_out = nij_row_out + cont;
%                         end
%                     end
%                 end
%             end
%             
%             nij_col_in = 0;
%             nij_col_out = 0;
%             for i = 1:n_cols
%                 for j = i+1:n_cols
%                     cont = sum(matrix(:,i).*matrix(:,j))*(col_sum(i)~=col_sum(j))/min(sum(matrix(:,i)),sum(matrix(:,j)));
%                     if(col_modules(i) == col_modules(j))
%                         nij_col_in = nij_col_in + cont;
%                     else
%                         nij_col_out = nij_col_out + cont;
%                     end
%                 end
%             end
%             
%             denom = n_rows*(n_rows-1)/2 + n_cols*(n_cols-1)/2;
%             nodf_in = 100.0*(nij_row_in+nij_col_in)/denom;
%             nodf_out = 100.0*(nij_row_out+nij_col_out)/denom;
%             nodf = nodf_in + nodf_out;
%             
%             nodf.total = nodf; nodf.in = nodf_in; nodf.out = nodf_out;
%                 
%         end
%         function row_nets = TestRowsByLabeling(obj,rows_idx,n_trials,null_model)
%             
%             if(nargin == 1)
%                 n_trials = Options.REPLICATES;
%                 null_model = @NullModels.NULL_1;
%                 rows_idx = obj.bipweb.rows_idx;
%             elseif(nargin == 2)
%                 n_trials = Options.REPLICATES;
%                 null_model = @NullModels.NULL_1;
%             elseif(nargin == 3)
%                 null_model = @NullModels.NULL_1;
%             end
%             
%             row_nets = obj.GetNetsByRowLabeling(rows_idx);
%             
%             InternalStatistics.TestNetworks(row_nets,n_trials,null_model)
%             
%         end
%         
%         function inter_nets = TestByInteractions(obj,n_trials,null_model)
%             
%             if(nargin == 1)
%                 n_trials = Options.REPLICATES;
%                 null_model = @NullModels.NULL_1;
%             elseif(nargin == 2)
%                 n_trials = Options.REPLICATES;
%             end
%             
%             inter_nets = obj.GetNetsByInteractions();
%             
%             InternalStatistics.TestNetworks(inter_nets,n_trials,null_model)
%             
%         end
%         
%         function [nets row_idx] = GetNetsByRowLabeling(obj,row_labels)
%             
%             if(nargin==1)
%                 row_labels = obj.bipweb.rows_idx;
%             end
%             idx = unique(row_labels);
%             n_nets = length(idx);
%             nets = cell(n_nets,1);
%             row_idx = cell(n_nets,1);
%             matrix = obj.bipweb.webmatrix;
%             
%             for i = 1:n_nets
%                 idx_rows = find(row_labels==idx(i));
%                 nets{i} = Bipartite(matrix(idx_rows,:));
%                 nets{i}.rows_idx = obj.bipweb.rows_idx(idx_rows);
%                 nets{i}.row_labels = obj.bipweb.row_labels(idx_rows);
%                 nets{i}.col_labels = obj.bipweb.col_labels;
%                 row_idx{i} = idx_rows;
%             end
%             
%         end
%         
%         function [nets row_idx col_idx]= GetNetsByInteractions(obj)
%             
%             matrix = obj.bipweb.webmatrix;
%             idx = unique(matrix);
%             if(idx(1) == 0); idx(1) = []; end;
%             n_nets = length(idx);
%             nets = cell(n_nets,1);
%             row_idx = cell(n_nets,1);
%             
%             for i = 1:n_nets
%                 [matrix_ex idx_rows idx_cols] = MatrixNull.TYPE_MATRIX_NON_ZERO(matrix, idx(i));
%                 nets{i} = Bipartite(matrix_ex);
%                 nets{i}.rows_idx = obj.bipweb.rows_idx(idx_rows);                
%                 nets{i}.row_labels = obj.bipweb.row_labels(idx_rows);
%                 %nets{i}.cols_idx = obj.bipweb.cols_idx(idx_cols);
%                 nets{i}.col_labels = obj.bipweb.col_labels(idx_cols);
%                 row_idx{i} = idx_rows;
%                 col_idx{i} = idx_cols;
%             end
%             
%         end
