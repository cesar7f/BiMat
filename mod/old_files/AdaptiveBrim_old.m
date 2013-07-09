classdef AdaptiveBrim < BipartiteModularity
    % AdaptiveBrim Main code class. 
    % Adaptive Brim algorithm to evaluate modularity. However, this objects
    % is just a interface for the Barber adaptive brim algorithm code.
    % Thiss is a son class of the BipartiteModularity class. The adaptive
    % BRIM algorithm is described on the paper:
    %    
    %    Barber, Michael J. Modularity and community detection in bipartite
    %    networks. Physical Review E. 2007
    %
    % AdaptiveBrim Methods:
    %    AdaptiveBrim - Main constructor    
    %    DetectComponent - Detect the modularity in a single component of
    %    the network
    %
    % See BipartiteModularity.
    methods
        
        function obj = AdaptiveBrim(bipmatrix)
        % obj = AdaptiveBrim(bipmatrix) - Main constructor
            
            %Call the parent class
            obj = obj@BipartiteModularity(bipmatrix);
            
        end
        
        
        function obj = DetectComponent(obj)
        % obj = DetectComponent(obj) - Detect the modularity in a specific
        % component.
            
            %Get number of edges
            nEdges = sum(obj.matrix_component(:));
            
            qmax = -0.5; %Negative symbol avoids problems in full connected graphs
            
            %Repeat the algorithm obj.trials times and chose the best the
            %one with the highest modularity value.
            for i = 1:obj.trials
                
                %Call the specific Adaptive Brim algorithm coded by Barber
                %team.
                [R, T, Qhist1, ~, ~] = abrim(obj.bb_component, nEdges);
                
                %Keep always the best modularity configuration.
                if (Qhist1(end) > qmax)
                   
                    qmax = Qhist1(end);
                    RRmax = R;
                    TTmax = T;
                    
                end
                
            end
            
            %Assign the best configuration of modularity
            obj.rr_component = RRmax;
            obj.tt_component = TTmax;
            obj.Qb = qmax;

            %Clean empty columns in the community matrices.
            obj.rr_component = BipartiteModularity.CLEAR_MODULE_MATRIX(obj.rr_component);
            obj.tt_component = BipartiteModularity.CLEAR_MODULE_MATRIX(obj.tt_component);
            
            if(size(obj.rr_component,2) ~= size(obj.tt_component,2))
                display(3);
            end
            
            obj.N_component = size(obj.rr_component,2);
            
            %Assign each row and column in the matrix to its corresponding
            %modules.
            obj.row_modules_component = BipartiteModularity.ASSIGN_MODULES(obj.rr_component);
            obj.col_modules_component = BipartiteModularity.ASSIGN_MODULES(obj.tt_component);
            %obj.N_component = size(obj.row_modules_component,2); 
            %obj.AssignModules();
            
            %obj.CalculateQrValue();
            %obj.N = obj.CommunityQuantity;
            %obj.Qb = obj.Q;
            
        end
            
    end
    
    
end