classdef LPBrim < BipartiteModularity
% LPBrim - Main code class
% LPBrim algorithm that may work better than the adaptive brim version
% for large scale bipartite networks. To know how the algorithm works
% you can consult the following paper:
%
%    Liu, Xin and Murata, Tsuyoshi. Community detection in large-scale
%    bipartite networks. Web Intelligence and Intelligent Agent 
%    Technologies, 2009
%
% LPBrim Properties:
%    red_labels - Works as module identifiers for rows
%    blue_labels - Works as module identifiers for columns
%
% LPBrim Methods:
%    LPBrim - Main constructor    
%
% See also:
%    BipartiteModularity, AdaptiveBrim, and LeadingEigenvector
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        red_labels           = 0; %Works as module identifiers for rows
        blue_labels          = 0; %Works as module identifiers for columns
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURES ALGORITHM
    methods
        
        function obj = LPBrim(bipmatrix)
        % LPBrim - Main Constructor
        % 
        %   obj = LPBrim(MATRIX) Creates an LPBrim object obj
        %   using a bipartite adjacency matrix MATRIX that will be used to
        %   calculate modularity using the LPBrim Algorithm
        %
        % See also:
        %   LPBrim
            
            %Call the parent class    
            obj = obj@BipartiteModularity(bipmatrix);
           
        end
            
    end
    
    methods(Access = 'protected')
        
        function obj = DetectComponent(obj)
        % DetectComponent - Main method of the algorithm
        %
        %   obj = DetectComponent(obj) Detect the modularity in a specific
        %   component.
            
            %Apply the LP algorithm
            obj.LP();
            
            %Tunne the result by applying the BRIM algorithm in the found
            %configuration.
            [obj.rr_component obj.tt_component, obj.Qb] = BipartiteModularity.BRIM(obj.rr_component, obj.bb_component,obj.tt_component,sum(obj.matrix_component(:)));
            %obj.BRIM();
                        
            
        end
       
        function obj = LP(obj)
            
            %Find the initial maximal number of modules.
            cmax = max(obj.n_rows_component,obj.n_cols_component);
            
            %Create your module matrices with the appropiate size.
            obj.rr_component = zeros(obj.n_rows_component,cmax);
            obj.tt_component = zeros(obj.n_cols_component,cmax);
            
            %Assign each node in the module matrices to different module (column).
            obj.rr_component(1:obj.n_rows_component,1:obj.n_rows_component) = eye(obj.n_rows_component);
            obj.tt_component(1:obj.n_cols_component,1:obj.n_cols_component) = eye(obj.n_cols_component);
            
            nEdges = sum(sum(obj.matrix_component));
            
            %Calculate the initial modularity the previous configuration
            qmax = BipartiteModularity.CALCULATE_Qb_VALUE(obj.rr_component,obj.bb_component,obj.tt_component,nEdges);
            
            %Keep the current configuration as the best one so far
            qmaxglobal = qmax;
            ttMaxGlobal = obj.tt_component;
            rrMaxGlobal = obj.rr_component;
            
            %Repeat for obj.trials the algorithm and always keep the best
            %one
            for j = 1:obj.trials
                
                %Fin the maximal number of modules.
                cmax = max(obj.n_rows_component,obj.n_cols_component);
                
                %Assign the LP module indexes (each row and column node in
                %a independen module).
                obj.red_labels = 1:obj.n_rows_component;
                obj.blue_labels = 1:obj.n_cols_component;
                
                %Create your module matrices with the appropiate size and
                %assign each node to different module
                obj.rr_component = zeros(obj.n_rows_component,cmax);
                obj.tt_component = zeros(obj.n_cols_component,cmax);
                obj.rr_component(1:obj.n_rows_component,1:obj.n_rows_component) = eye(obj.n_rows_component);
                obj.tt_component(1:obj.n_cols_component,1:obj.n_cols_component) = eye(obj.n_cols_component);
                %qmax = obj.CalculateQValue();
                qmax = BipartiteModularity.CALCULATE_Qb_VALUE(obj.rr_component,obj.bb_component,obj.tt_component,nEdges);
                rrMax = obj.rr_component;
                ttMax = obj.tt_component;
                
                %LP BRIM main algorithm part
                while(1)

                    %Propagate red labels to blue labels
                    for i = 1:obj.n_cols_component
                        neighbors = obj.matrix_component(:,i);
                        labels = obj.red_labels(neighbors);
                        uniq = unique(labels);
                        nuniq = length(uniq);
                        rperm = randperm(nuniq);
                        uniq = uniq(rperm);
                        count = arrayfun(@(x) sum(x==labels), uniq);
                        [~, index] = max(count);
                        if(~isempty(uniq(index)))
                            obj.blue_labels(i) = uniq(index);
                        end
                    end

                    uniq = unique(obj.blue_labels);
                    nbluec = length(uniq);
                    newlab = arrayfun(@(x) find(uniq==x), obj.blue_labels);
                    obj.blue_labels = newlab;

                    %Propagate blue labels to red labels
                    for i = 1:obj.n_rows_component
                        neighbors = obj.matrix_component(i,:);
                        labels = obj.blue_labels(neighbors);
                        uniq = unique(labels);
                        nuniq = length(uniq);
                        rperm = randperm(nuniq);
                        uniq = uniq(rperm);
                        count = arrayfun(@(x) sum(x==labels), uniq);
                        [~, index] = max(count);
                        if(~isempty(uniq(index)))
                            obj.red_labels(i) = uniq(index);
                        end
                    end

                    uniq = unique(obj.red_labels);
                    nredc = length(uniq);
                    newlab = arrayfun(@(x) find(uniq==x), obj.red_labels);
                    obj.red_labels = newlab;

                    cmax = max(nbluec, nredc);
                    indxblue = sub2ind([obj.n_cols_component cmax], 1:obj.n_cols_component, obj.blue_labels);
                    indxred = sub2ind([obj.n_rows_component cmax], 1:obj.n_rows_component, obj.red_labels);

                    obj.rr_component = zeros(obj.n_rows_component,cmax);
                    obj.tt_component = zeros(obj.n_cols_component,cmax);
                    obj.rr_component(indxred) = 1;
                    obj.tt_component(indxblue) = 1;

                    q = BipartiteModularity.CALCULATE_Qb_VALUE(obj.rr_component,obj.bb_component,obj.tt_component,nEdges);
                    
                    %If the modularity configuration is better than the
                    %current one, replace it. Otherwise exit the LP loop.
                    if(q > qmax) 
                        qmax = q;
                        rrMax = obj.rr_component;
                        ttMax = obj.tt_component;
                    else
                        break;
                    end

                end
                
                %If the LP configuration found in trial i is better than
                %the global configuration replace the global with it.
                if(qmax > qmaxglobal)
                    qmaxglobal = qmax;
                    rrMaxGlobal = rrMax;
                    ttMaxGlobal = ttMax;
                end
            end
            
            %Assign the values found to the standart class variables.
            obj.Qb = qmaxglobal;
            obj.rr_component = rrMaxGlobal;
            obj.tt_component = ttMaxGlobal;
        end
        
        
    end
    
    
end