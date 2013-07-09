classdef GroupStatistics < handle
   
    properties
        matrices        = {};
        names           = {};
        qb_vals         = [];
        nestvals        = [];
        tempvals        = [];
        n_networks      = 0;
        do_modul        = 1;
        do_nest         = 1;
        do_temp         = 1;
        modul_class     = @AdaptiveBrim;
        replicates      = 100;
        null_model      = @NullModels.EQUIPROBABLE;
        networks        = {};
        clean_nulls     = 1;
    end
    
    methods
        
        function obj = GroupStatistics(matrices_or_networks_or_files,str_names)
            
            assert(isa(matrices_or_networks_or_files,'cell') || ...
                isa(matrices_or_networks_or_files,'char'))
            
            if(isa(matrices_or_networks_or_files,'cell'))
                %If you are testing a single network why do group testing?
                obj.n_networks = length(matrices_or_networks_or_files);
                for i = 1:obj.n_networks
                    if(isa(matrices_or_networks_or_files{i},'Bipartite'))
                        obj.networks{i} = matrices_or_networks_or_files{i};
                    else
                        obj.networks{i} = Bipartite(matrices_or_networks_or_files{i});
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
                    obj.networks{i}.name = ['Network %i', i];
                    obj.names{i} = ['Network %i', i];
                end
            end
            
        end
        
        
        function obj = DoGroupTesting(obj,ntrials,nullmodel)
            
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
                
                if(isempty(net{i}.webmatrix))
                    continue;
                end
                
                net{i}.statistics.DoNulls(obj.null_model,obj.replicates);
                net{i}.statistics.print_output = 0;
                
                if(obj.do_nest == 1)
                    net{i}.statistics.Nestedness();
                    nest.value(i,1) = net{i}.statistics.nestvals.value;
                    nest.mean(i,1) = net{i}.statistics.nestvals.mean;
                    nest.std(i,1) = net{i}.statistics.nestvals.std;
                    nest.p(i,1) = net{i}.statistics.nestvals.p;
                    nest.ci(i,:) = net{i}.statistics.nestvals.ci;
                    nest.zscore(i,1) = net{i}.statistics.nestvals.zscore;
                    nest.percent(i,1) = net{i}.statistics.nestvals.percent;
                    nest.random_values(i,:) = net{i}.statistics.nestvals.random_values;
                end
                    
                if(obj.do_modul == 1)
                    net{i}.modules = obj.modul_class(net{i}.matrix);
                    net{i}.statistics.Modularity();
                    qb.value(i,1) = net{i}.statistics.qb_vals.value;
                    qb.mean(i,1) = net{i}.statistics.qb_vals.mean;
                    qb.std(i,1) = net{i}.statistics.qb_vals.std;
                    qb.p(i,1) = net{i}.statistics.qb_vals.p;
                    qb.ci(i,:) = net{i}.statistics.qb_vals.ci;
                    qb.zscore(i,1) = net{i}.statistics.qb_vals.zscore;
                    qb.percent(i,1) = net{i}.statistics.qb_vals.percent;
                    qb.random_values(i,:) = net{i}.statistics.qb_vals.random_values;
                end
                
                if(obj.do_temp == 1)
                    net{i}.statistics.Temperature();
                    temp.value(i,1) = net{i}.statistics.tempvals.value;
                    temp.mean(i,1) = net{i}.statistics.tempvals.mean;
                    temp.std(i,1) = net{i}.statistics.tempvals.std;
                    temp.p(i,1) = net{i}.statistics.tempvals.p;
                    temp.ci(i,:) = net{i}.statistics.tempvals.ci;
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
                obj.nestvals = nest;
                obj.nestvals.replicates = obj.replicates;
                obj.nestvals.model = obj.null_model;
            end
            
            if(obj.do_temp == 1)
                obj.tempvals = temp;
                obj.tempvals.replicates = obj.replicates;
                obj.tempvals.model = obj.null_model;
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
            
            if(~isempty(obj.nestvals))
                headers{i} = 'NODF';
                headers{i+1} = 'NODF mean';
                headers{i+2} = 'NODF z-score';
                headers{i+3} = 'NODF percent';
                columns = [columns obj.nestvals.value];
                columns = [columns obj.nestvals.mean];
                columns = [columns obj.nestvals.zscore];
                columns = [columns obj.nestvals.percent];
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
    
    methods(Static)
       
%         function [matrices names] = GET_MATRICES_FROM_FILES(file_name)
%             
%             [pathstr, ~, ext] = fileparts(file_name);
%             files = dir(file_name);
%             
%             for i = 1:length(files)
%                
%                 fprintf('%s\n',files(i).name);
%                 matrices{i} = dlmread([pathstr,'/',files(i).name]);
%                 
%                 [~, name, ~] = fileparts(files(i).name);
%                 names{i} = name;
%                 
%             end
%             
%         end
        
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
%                 obj.nodf_nested.v(i) = bn.statistics.nestvals.nodf;
%                 obj.nodf_nested.z(i) = bn.statistics.nestvals.zscore;
%                 obj.nodf_nested.p(i) = bn.statistics.nestvals.percent;
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