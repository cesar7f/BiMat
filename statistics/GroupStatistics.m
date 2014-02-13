classdef GroupStatistics < handle
   
    properties
        matrices        = {}; % The matrices that will be tested
        names           = {}; % The names of the matrix (no required)
        qb_vals         = []; % Statistics for the modularity values
        nodf_vals        = []; % Statistics for the nodf values
        ntc_statistics        = []; % Statistics for the ntc values
        n_networks      = 0;  % Number of matrices (networks)
        do_modul        = 1;  % Flag to Perform the tests for modularity
        do_nest         = 1;  % Flag to Perform the tests for nodf
        do_temp         = 1;  % Flag to Perform the tests for temperature
        modul_class     = Options.MODULARITY_ALGORITHM; %Algorithm for modularity
        replicates      = Options.REPLICATES; %Number of replicates for the tests
        null_model      = Options.DEFAULT_NULL_MODEL; %Null model that will be used during the tests
        networks        = {}; % A cell that will contain the bipartite objects
        clean_nulls     = 1; % Clean the random matrices (nulls) after performing each test. Useful for not saturating memory
    end
    
    % Properties using during plotting
    properties
        p_value         = Options.P_VALUE;
        z_value         = Options.Z_VALUE;
        nest_mat_test   = 1;
        plotter   = {};
    end
    
    methods
        
        function obj = GroupStatistics(matrices_or_networks_or_files,str_names)
        % GroupStatistics - Main Constructor
        %   obj = GroupStatistics(matrices_or_networks_or_files)
        %   Creates an instance obj of the GroupStatistics class, where the
        %   argument matrices_or_networks_or_files can be a cell of
        %   matrices or a cell of Bipartite objects.
        %   obj = GroupStatistics(matrices_or_networks_or_files,str_names)
        %   Creates an instance obj of the GroupStatistics class, where the
        %   argument matrices_or_networks_or_files can be a cell of
        %   matrices or a cell of Bipartite objects. It also name each
        %   network using the corresponding cell of strings str_names
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
        
        
        function obj = DoGroupTesting(obj,ntrials,nullmodel)
        % DoGroupTesting - Main method to perform the testing in a group of
        % matrices (networks). Remember to specify the flags about what
        % metrics you will test (see below).
        %   obj = DoGroupTesting(obj) - Perform the testing using default
        %   values for the number of replicates and null model.
        %   obj = DoGroupTesting(obj,ntrials) - Perform the testing using
        %   ntrrials replicates and the default null model.
        %   obj = DoGroupTesting(obj,ntrials,nullmodel) - Perform the
        %   testing using the specified number of replicates ntrials and the
        %   specified null model nullmodel.
        % See also:
        %   GroupStatistics.do_modul, GroupStatistics.do_nest, GroupStatistics.do_temp
        
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
                
                if(obj.do_nest == 1)
                    net{i}.statistics.TestNODF();
                    nest.value(i,1) = net{i}.statistics.nodf_vals.value;
                    nest.mean(i,1) = net{i}.statistics.nodf_vals.mean;
                    nest.std(i,1) = net{i}.statistics.nodf_vals.std;
%                    nest.p(i,1) = net{i}.statistics.nodf_vals.p;
%                    nest.ci(i,:) = net{i}.statistics.nodf_vals.ci;
                    nest.zscore(i,1) = net{i}.statistics.nodf_vals.zscore;
                    nest.percent(i,1) = net{i}.statistics.nodf_vals.percent;
                    nest.random_values(i,:) = net{i}.statistics.nodf_vals.random_values;
                end
                    
                if(obj.do_modul == 1)
                    net{i}.modules = obj.modul_class(net{i}.matrix);
                    net{i}.statistics.Modularity();
                    qb.value(i,1) = net{i}.statistics.qb_vals.value;
                    qb.mean(i,1) = net{i}.statistics.qb_vals.mean;
                    qb.std(i,1) = net{i}.statistics.qb_vals.std;
%                    qb.p(i,1) = net{i}.statistics.qb_vals.p;
%                    qb.ci(i,:) = net{i}.statistics.qb_vals.ci;
                    qb.zscore(i,1) = net{i}.statistics.qb_vals.zscore;
                    qb.percent(i,1) = net{i}.statistics.qb_vals.percent;
                    qb.random_values(i,:) = net{i}.statistics.qb_vals.random_values;
                end
                
                if(obj.do_temp == 1)
                    net{i}.statistics.TestNTC();
                    temp.value(i,1) = net{i}.statistics.tempvals.value;
                    temp.mean(i,1) = net{i}.statistics.tempvals.mean;
                    temp.std(i,1) = net{i}.statistics.tempvals.std;
%                    temp.p(i,1) = net{i}.statistics.tempvals.p;
%                    temp.ci(i,:) = net{i}.statistics.tempvals.ci;
                    temp.zscore(i,1) = net{i}.statistics.tempvals.zscore;
                    temp.percent(i,1) = net{i}.statistics.tempvals.percent;
                    temp.random_values(i,:) = net{i}.statistics.tempvals.random_values;
                end
                
                %Avoid memory overflow
                if(obj.clean_nulls == 1)
                    net{i}.statistics.nulls = {};
                end
            end
            
            if(obj.do_modul == 1)
                obj.qb_vals = qb;
                obj.qb_vals.algorithm = obj.modul_class;
                obj.qb_vals.replicates = obj.replicates;
                obj.qb_vals.model = obj.null_model;
            end
            
            if(obj.do_nest == 1)
                obj.nodf_vals = nest;
                obj.nodf_vals.replicates = obj.replicates;
                obj.nodf_vals.model = obj.null_model;
            end
            
            if(obj.do_temp == 1)
                obj.tempvals = temp;
                obj.tempvals.replicates = obj.replicates;
                obj.tempvals.model = obj.null_model;
            end
            
        end
        
        function q_values = TestModularityAlgorithms(obj)
            % TestModularityAlgorithms - Test all the modularity algorithms
            % with the idea to study what algorithm will be the more
            % appropiate for the corresponding group of matrices
            % (networks). 
            %
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
        
        function PrintResults(obj)
           
            headers{1} = 'Network';
            columns = (1:obj.n_networks)';
            i = 2;
            if(~isempty(obj.qb_vals))
                headers{i} = 'Qb';
                headers{i+1} = 'Qb mean';
                headers{i+2} = 'Qb z-score';
                headers{i+3} = 'Qb percent';
                columns = [columns obj.qb_vals.value];
                columns = [columns obj.qb_vals.mean];
                columns = [columns obj.qb_vals.zscore];
                columns = [columns obj.qb_vals.percent];
                i = i+4;
            end
            
            if(~isempty(obj.nodf_vals))
                headers{i} = 'NODF';
                headers{i+1} = 'NODF mean';
                headers{i+2} = 'NODF z-score';
                headers{i+3} = 'NODF percent';
                columns = [columns obj.nodf_vals.value];
                columns = [columns obj.nodf_vals.mean];
                columns = [columns obj.nodf_vals.zscore];
                columns = [columns obj.nodf_vals.percent];
                i = i+4;
            end
            
            if(~isempty(obj.tempvals))
                headers{i} = 'NTC';
                headers{i+1} = 'NTC mean';
                headers{i+2} = 'NTC z-score';
                headers{i+3} = 'NTC percent';
                        
                columns = [columns obj.tempvals.value];
                columns = [columns obj.tempvals.mean];
                columns = [columns obj.tempvals.zscore];
                columns = [columns obj.tempvals.percent];
            end
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,',');
            
            fprintf(str)
            
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
%                 obj.adaptive_modular.v(i) = bn.statistics.qb_vals.Qb;
%                 obj.adaptive_modular.z(i) = bn.statistics.qb_vals.zscore;
%                 obj.adaptive_modular.p(i) = bn.statistics.qb_vals.percent;
%                 
%                 obj.nodf_vals.v(i) = bn.statistics.nodf_vals.nodf;
%                 obj.nodf_vals.z(i) = bn.statistics.nodf_vals.zscore;
%                 obj.nodf_vals.p(i) = bn.statistics.nodf_vals.percent;
%                 
%                 obj.ntc_nested.v(i) = bn.statistics.tempvals.ntc;
%                 obj.ntc_nested.z(i) = bn.statistics.tempvals.zscore;
%                 obj.ntc_nested.p(i) = bn.statistics.tempvals.percent;
%                 
%                 obj.eig_nested.v(i) = bn.statistics.eigvals.maxe;
%                 obj.eig_nested.z(i) = bn.statistics.eigvals.zscore;
%                 obj.eig_nested.p(i) = bn.statistics.eigvals.percent;
%                 
%             end