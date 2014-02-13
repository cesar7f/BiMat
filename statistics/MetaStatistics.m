% MetaStatistics - Statistical analysis for a set of bipartite networks.
% See the next paper for a nice example of how this class can be used chek
% the next paper:
%
%   Cesar O. Flores, Justin R. Meyer, Sergi Valverde, Lauren Farr, and Joshua S. Weitz. 
%   Statistical structure of host-phage interactions.
%   PNAS 2011
%
% MetaStatistics Properties:
%     matrices             = {}; % The matrices that will be tested
%     names                = {}; % The names of the matrix (no required)
%     qb_values            = []; % Statistics for the modularity values
%     qr_values            = []  % Statistics for the modularity interactions ratio values
%     nodf_values          = []; % Statistics for the nodf values
%     ntc_values           = []; % Statistics for the ntc values
%     n_networks           = 0;  % Number of matrices (networks)
%     do_modul             = 1;  % Flag to Perform the tests for modularity
%     do_nodf              = 1;  % Flag to Perform the tests for nodf
%     do_ntc               = 1;  % Flag to Perform the tests for temperature
%     modularity_algorithm = Options.MODULARITY_ALGORITHM; %Algorithm for modularity
%     replicates           = Options.REPLICATES; %Number of replicates for the tests
%     null_model           = Options.DEFAULT_NULL_MODEL; %Null model that will be used during the tests
%     networks             = {}; % A cell that will contain the Bipartite objects
%     clean_nulls          = 1; % Clean the random matrices (nulls) after performing each test. Useful for not saturating memory
%
% MetaStatistics Methods:
%     MetaStatistics - Main Constructor
%     DoMetaAnalyisis - Main method to perform statistical analysis in a set of matrices
%     Print - Print all analysed results
%
% See also:
%     MetaStatisticsPlotter, StatisticalTest, NullModels
classdef MetaStatistics < handle
   
    properties
        matrices             = {}; % The matrices that will be tested
        names                = {}; % The names of the matrix (no required)
        qb_values            = []; % Statistics for the modularity values
        qr_values            = []  % Statistics for the modularity interactions ratio values
        nodf_values          = []; % Statistics for the nodf values
        ntc_values           = []; % Statistics for the ntc values
        n_networks           = 0;  % Number of matrices (networks)
        do_modul             = 1;  % Flag to Perform the tests for modularity
        do_nodf              = 1;  % Flag to Perform the tests for nodf
        do_ntc               = 1;  % Flag to Perform the tests for temperature
        modularity_algorithm = Options.MODULARITY_ALGORITHM; %Algorithm for modularity
        replicates           = Options.REPLICATES; %Number of replicates for the tests
        null_model           = Options.DEFAULT_NULL_MODEL; %Null model that will be used during the tests
        networks             = {}; % A cell that will contain the Bipartite objects
        clean_nulls          = 1; % Clean the random matrices (nulls) after performing each test. Useful for not saturating memory
    end
    
    methods
        
        function obj = MetaStatistics(matrices_or_networks_or_files,str_names)
        % MetaStatistics - Main Constructor
        %
        %   obj = MetaStatistics(matrices_or_networks_or_files)
        %   Creates an instance obj of the MetaStatistics class, where the
        %   argument matrices_or_networks_or_files can be a cell of
        %   bipartite matrices or a cell of Bipartite objects.
        %
        %   obj = MetaStatistics(matrices_or_networks_or_files,str_names)
        %   Creates an instance obj of the MetaStatistics class, where a
        %   cell of strings str_names is used for naming the bipartite
        %   networks.
        %
            assert(isa(matrices_or_networks_or_files,'cell') || ...
                isa(matrices_or_networks_or_files,'char'))
            
            if(isa(matrices_or_networks_or_files,'cell'))
                %If you are testing a single network why do group testing?
                obj.n_networks = length(matrices_or_networks_or_files);
                for i = 1:obj.n_networks
                    if(isa(matrices_or_networks_or_files{i},'Bipartite'))
                        obj.networks{i} = matrices_or_networks_or_files{i};
                        obj.matrices{i} = obj.networks{i}.matrix;
                    else
                        obj.networks{i} = Bipartite(MatrixFunctions.NON_ZERO_MATRIX(matrices_or_networks_or_files{i}));
                        %obj.networks{i} = Bipartite(matrices_or_networks_or_files{i});
                        obj.matrices{i} = obj.networks{i}.matrix;
                    end
                end
            end
            
            if(nargin == 2)
                assert(length(str_names) == obj.n_networks);
                obj.names = str_names;
                for i = 1:obj.n_networks
                    obj.networks{i}.name = obj.names{i};
                end
            else
                for i = 1:obj.n_networks
                    obj.networks{i}.name = sprintf('Network %i',i);
                    obj.names{i} = sprintf('Network %i',i);
                end
            end
            
            obj.plotter = MetaStatisticsPlotter(obj);
            
        end
        
        
        function obj = DoMetaAnalyisis(obj,ntrials,nullmodel)
        % DoMetaAnalyisis - Main method to perform an statistical analysis
        % in a set of matrices (networks). Remember to specify the flags about what
        % metrics you will test (see below).
        %
        %   obj = DoMetaAnalyisis(obj) - Perform the testing using default
        %   values for the number of replicates and null model.
        %
        %   obj = DoMetaAnalyisis(obj,replicates) - Perform the testing using
        %   replicates random matrices and the default null model.
        %
        %   obj = DoMetaAnalyisis(obj,replicates,nullmodel) - Perform the
        %   testing using the replicates random marices and the
        %   specified null model nullmodel.
        %
        % See also:
        %   MetaStatistics.do_modul, MetaStatistics.do_nodf, MetaStatistics.do_ntc
        
            if(nargin == 2)
                obj.replicates = ntrials;
            elseif(nargin == 3)
                obj.replicates = ntrials;
                obj.null_model = nullmodel;
            end
            
            n = obj.n_networks;
            net = obj.networks;
           
            for i = 1:n
                
                fprintf('Testing Matrix: %i . . .\n', i);
                
                %if(isempty(net{i}.webmatrix))
                %    continue;
                %end
                
                net{i}.statistics.DoNulls(obj.replicates,obj.null_model);
                net{i}.statistics.print_output = 0;
                net{i}.statistics.print_status = 0;
                
                if(obj.do_nodf == 1)
                    net{i}.statistics.TestNODF();
                    nest.value(i,1) = net{i}.statistics.nodf_values.value;
                    nest.mean(i,1) = net{i}.statistics.nodf_values.mean;
                    nest.std(i,1) = net{i}.statistics.nodf_values.std;
                    nest.zscore(i,1) = net{i}.statistics.nodf_values.zscore;
                    nest.percentile(i,1) = net{i}.statistics.nodf_values.percentile;
                    nest.random_values(i,:) = net{i}.statistics.nodf_values.random_values;
                end
                    
                if(obj.do_modul == 1)
                    net{i}.modules = obj.modularity_algorithm(net{i}.matrix);
                    net{i}.statistics.Modularity();
                    
                    qb.value(i,1) = net{i}.statistics.qb_values.value;
                    qb.mean(i,1) = net{i}.statistics.qb_values.mean;
                    qb.std(i,1) = net{i}.statistics.qb_values.std;
                    qb.zscore(i,1) = net{i}.statistics.qb_values.zscore;
                    qb.percentile(i,1) = net{i}.statistics.qb_values.percentile;
                    qb.random_values(i,:) = net{i}.statistics.qb_values.random_values;
                    
                    qr.value(i,1) = net{i}.statistics.qb_values.value;
                    qr.mean(i,1) = net{i}.statistics.qb_values.mean;
                    qr.std(i,1) = net{i}.statistics.qb_values.std;
                    qr.zscore(i,1) = net{i}.statistics.qb_values.zscore;
                    qr.percentile(i,1) = net{i}.statistics.qb_values.percentile;
                    qr.random_values(i,:) = net{i}.statistics.qb_values.random_values;
                    
                end
                
                if(obj.do_ntc == 1)
                    net{i}.statistics.TestNTC();
                    temp.value(i,1) = net{i}.statistics.ntc_values.value;
                    temp.mean(i,1) = net{i}.statistics.ntc_values.mean;
                    temp.std(i,1) = net{i}.statistics.ntc_values.std;
                    temp.zscore(i,1) = net{i}.statistics.ntc_values.zscore;
                    temp.percentile(i,1) = net{i}.statistics.ntc_values.percentile;
                    temp.random_values(i,:) = net{i}.statistics.ntc_values.random_values;
                end
                
                %Avoid memory overflow
                if(obj.clean_nulls == 1)
                    net{i}.statistics.nulls = {};
                end
            end
            
            if(obj.do_modul == 1)
                obj.qb_values = qb;
                obj.qb_values.algorithm = obj.modularity_algorithm;
                obj.qb_values.replicates = obj.replicates;
                obj.qb_values.null_model = obj.null_model;
                
                obj.qr_values = qr;
                obj.qr_values.algorithm = obj.modularity_algorithm;
                obj.qr_values.replicates = obj.replicates;
                obj.qr_values.null_model = obj.null_model;
            end
            
            if(obj.do_nodf == 1)
                obj.nodf_values = nest;
                obj.nodf_values.replicates = obj.replicates;
                obj.nodf_values.null_model = obj.null_model;
            end
            
            if(obj.do_ntc == 1)
                obj.ntc_values = temp;
                obj.ntc_values.replicates = obj.replicates;
                obj.ntc_values.null_model = obj.null_model;
            end
            
        end
        
        function q_values = TestCommunityStructureAlgorithms(obj)
            % TestCommunityStructureAlgorithms - Test all the modularity algorithms
            % with the idea to study what algorithm will be the more
            % appropiate for the corresponding group of matrices
            % (networks). 
            
            q_values = zeros(obj.n_networks,3);
            for i = 1:obj.n_networks
                
                fprintf('Testing Matrix: %i . . .\n', i);
                obj.networks{i}.modules = AdaptiveBrim(obj.matrices{i});
                obj.networks{i}.modules.Detect();
                q_values(i,1) = obj.networks{i}.modules.Qb;
                
                obj.networks{i}.modules = LPBrim(obj.matrices{i});
                obj.networks{i}.modules.Detect();
                q_values(i,2) = obj.networks{i}.modules.Qb;
                
                obj.networks{i}.modules = LeadingEigenvector(obj.matrices{i});
                obj.networks{i}.modules.Detect();
                q_values(i,3) = obj.networks{i}.modules.Qb;
                
            end
        end
        
        function Print(obj,filename)
        % Print - Print all analysed results
        %
        %   STR = Print(obj) Print the statistical results to screen and
        %   returns this information to the string STR
        %
        %   STR = Print(obj, FILE) Print the statistical resultss to screen and
        %   text file FILE and return this information to the string STR   
        %
        % See also: 
        %   Printer
        
            headers{1} = 'Network';
            columns = (1:obj.n_networks)';
            i = 2;
            if(~isempty(obj.qb_values))
                headers{i} = 'Qb';
                headers{i+1} = 'Qb mean';
                headers{i+2} = 'Qb z-score';
                headers{i+3} = 'Qb percent';
                columns = [columns obj.qb_values.value];
                columns = [columns obj.qb_values.mean];
                columns = [columns obj.qb_values.zscore];
                columns = [columns obj.qb_values.percentile];
                
                headers{i} = 'Qr';
                headers{i+1} = 'Qr mean';
                headers{i+2} = 'Qr z-score';
                headers{i+3} = 'Qr percent';
                columns = [columns obj.qr_values.value];
                columns = [columns obj.qr_values.mean];
                columns = [columns obj.qr_values.zscore];
                columns = [columns obj.qr_values.percentile];
                
                i = i+8;
            end
            
            if(~isempty(obj.nodf_values))
                headers{i} = 'NODF';
                headers{i+1} = 'NODF mean';
                headers{i+2} = 'NODF z-score';
                headers{i+3} = 'NODF percent';
                columns = [columns obj.nodf_values.value];
                columns = [columns obj.nodf_values.mean];
                columns = [columns obj.nodf_values.zscore];
                columns = [columns obj.nodf_values.percentile];
                i = i+4;
            end
            
            if(~isempty(obj.ntc_values))
                headers{i} = 'NTC';
                headers{i+1} = 'NTC mean';
                headers{i+2} = 'NTC z-score';
                headers{i+3} = 'NTC percent';
                        
                columns = [columns obj.ntc_values.value];
                columns = [columns obj.ntc_values.mean];
                columns = [columns obj.ntc_values.zscore];
                columns = [columns obj.ntc_values.percentile];
            end
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,',');
            
            fprintf(str)
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
    end
      
end
% 
%     function obj = DoTotalGroupTesting(obj,n_trials)
%             
%             n = obj.n_matrices;
%             
%             for i = 1:n
%                 bn = Bipartite(obj.matrices{i});
%                 bn.statistics.print_output = 0;
%                 bn.statistics.DoCompleteAnalysis(n_trials);
%                 obj.adaptive_modular.v(i) = bn.statistics.qb_values.Qb;
%                 obj.adaptive_modular.z(i) = bn.statistics.qb_values.zscore;
%                 obj.adaptive_modular.p(i) = bn.statistics.qb_values.percentile;
%                 
%                 obj.nodf_values.v(i) = bn.statistics.nodf_values.nodf;
%                 obj.nodf_values.z(i) = bn.statistics.nodf_values.zscore;
%                 obj.nodf_values.p(i) = bn.statistics.nodf_values.percentile;
%                 
%                 obj.ntc_nested.v(i) = bn.statistics.ntc_values.ntc;
%                 obj.ntc_nested.z(i) = bn.statistics.ntc_values.zscore;
%                 obj.ntc_nested.p(i) = bn.statistics.ntc_values.percentile;
%                 
%                 obj.eig_nested.v(i) = bn.statistics.eigvals.maxe;
%                 obj.eig_nested.z(i) = bn.statistics.eigvals.zscore;
%                 obj.eig_nested.p(i) = bn.statistics.eigvals.percentile;
%                 
%             end