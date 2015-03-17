classdef AdaptiveBrim < BipartiteModularity
	% AdaptiveBrim - Main code class
	% Adaptive Brim algorithm to evaluate modularity. However, this objects
	% is just a interface for the Barber adaptive brim algorithm code.
	% Thiss is a son class of the BipartiteModularity class. The adaptive
	% BRIM algorithm is described on the paper:
	%    
	%    Barber, Michael J. Modularity and community detection in bipartite
	%    networks. Physical Review E. 2007
	%
	% AdaptiveBrim Properties:
	%   prob_reassigment - probability of reassigment when new modules are created
	%   expansion_factor - factor of expansion for creating new modules (2.0=duplicate number of modules).
	%
	% AdaptiveBrim Methods:
	%    AdaptiveBrim - Main constructor 
	%    EXPAND_MODULE_QUANTITY - Expand the number of modules (columns) in a module matrix (rr or tt).
	%
	% See also:
	%    BipartiteModularity, LeadingEigenvector, and LPBrim
    
    properties
        prob_reassigment = 0.5; %probability of reassigment when new modules are created
        expansion_factor = 2.0; %factor of expansion for creating new modules (2.0=duplicate number of modules).
    end
    
    methods
        
        function obj = AdaptiveBrim(bipmatrix)
			% AdaptiveBrim - Main Constructor
			% 
			%   obj = AdaptiveBrim(MATRIX) Creates an AdaptiveBrim object obj
			%   using a bipartite adjacency matrix MATRIX that will be used to
			%   calculate modularity using the Adaptive Brim Algorithm
			%
			% See also:
			%   AdaptiveBrim
				
            %Call the parent class
            obj = obj@BipartiteModularity(bipmatrix);
            
        end
        
    end
    
    methods (Access=protected)
	
        function obj = DetectComponent(obj)
			% DetectComponent - Main method of the algorithm
			%
			%   obj = DetectComponent(obj) Detect the modularity in a specific
			%   component.
            
            %Get number of edges
            obj.n_edges_component = sum(obj.matrix_component(:));
            
            [n_rows n_cols] = size(obj.matrix_component);
            
            sol_inf.rr = ones(n_rows,1);
            sol_inf.tt = ones(n_cols,1);
            sol_inf.N = 1;
            sol_inf.Q = 0.0;
            
            sol_mid = sol_inf;
            
            n_max = max(n_rows,n_cols);
            sol_global.Q = -0.5;
            
            for i = 1:obj.trials
                while(sol_mid.N < n_max)
                    
                    n_end = min(obj.expansion_factor*sol_mid.N, n_max);
                    
                    sol_end = obj.FindExpandedSolution(sol_mid,n_end);
                    
                    if(sol_end.Q > sol_mid.Q)
                        sol_inf = sol_mid;
                        sol_mid = sol_end;
                        if(sol_end.N < n_end)
                            break;
                        end
                    else
                        break;
                    end
                end
                
                sol_opt = obj.RecursiveBisection(sol_inf,sol_mid,sol_end);
                
                if(sol_global.Q < sol_opt.Q)
                    sol_global = sol_opt;
                end
            end
            
            %Assign the best configuration of modularity
            obj.rr_component = sol_global.rr;
            obj.tt_component = sol_global.tt;
            obj.Qb = sol_global.Q;

            
            
        end     
        
        function sol_new = FindExpandedSolution(obj,prev_sol,N)
			% FindExpandedSolution - Expand the module matrices RR and TT to a
			% new number of modules
			%
			%   sol_new = FindExpandedSolution(obj,prev_sol,N) Expand the number of
			%   modules prev_sol.N to a new value N by first doing a simple
			%   expansion and later optimize the values of rr and tt using BRIM
			%   algorithm. The empty modules are deleted before returning the
			%   result in a structure sol_new.
        
            sol_new.N = N;
            sol_new.rr = AdaptiveBrim.EXPAND_MODULE_QUANTITY(prev_sol.rr, N, obj.prob_reassigment);
            sol_new.tt = AdaptiveBrim.EXPAND_MODULE_QUANTITY(prev_sol.tt, N, obj.prob_reassigment);
            %solt = sol_new;
            [sol_new.rr sol_new.tt sol_new.Q] = BipartiteModularity.BRIM(sol_new.rr,obj.bb_component,sol_new.tt,obj.n_edges_component);
            
            %[sol_new.rr sol_new.tt] = BipartiteModularity.CLEAN_MODULE_MATRICES(sol_new.rr,sol_new.tt);
            %sol_new.N = size(sol_new.rr,2);
            %sol_new.Q = BipartiteModularity.CALCULATE_Qb_VALUE(sol_new.rr,obj.bb_component,sol_new.tt,obj.n_edges_component);
            %if(solt.N ~=sol_new.N&&N==2)
            %    display(3);
            %end
        end
         
        
        function sol = RecursiveBisection(obj,sol_inf,sol_mid,sol_end)
			% RecursiveBisection - Recursively look for the best solution by
			% doing some kind of bisection method/ternary search.
			%
			%   sol = RecursiveBisection(obj,sol_inf,sol_mid,sol_end) Find for
			%   the best solution between sol_inf and sol_end by looking first
			%   between sol_mid and sol_end and later between sol_inf ans
			%   sol_mid and keeping the best value at sol_mid. This function
			%   assumes that the modularity value has only one maxima and is
			%   well behaved.
        
            n_try_right = floor((sol_end.N + sol_mid.N)/2);
            n_try_left = floor((sol_inf.N + sol_mid.N)/2);
            
            if(n_try_right ~= sol_mid.N)
                sol_try = obj.FindExpandedSolution(sol_mid,n_try_right);
                if(sol_try.Q > sol_mid.Q)
                    sol_inf = sol_mid;
                    sol_mid = sol_try;
                    sol = obj.RecursiveBisection(sol_inf,sol_mid,sol_end);
                else%/(sol_try.Q >= sol_end.Q)
                    sol_end = sol_try;
                    %if(sol_end.N == sol_mid.N)
                    %    if(sol_mid.Q > sol_end.Q); sol_mid = sol_end; else sol_end = sol_mid; end;
                    %end
                    sol = obj.RecursiveBisection(sol_inf,sol_mid,sol_end);
                end
            elseif(n_try_left ~= sol_inf.N)
                sol_try = obj.FindExpandedSolution(sol_inf,n_try_left);
                if(sol_try.Q > sol_mid.Q)
                    sol_end = sol_mid;
                    sol_mid = sol_try;
                    sol = obj.RecursiveBisection(sol_inf,sol_mid,sol_end);
                else%(sol_try.Q >= sol_end.Q)
                    sol_inf = sol_try;
                    sol = obj.RecursiveBisection(sol_inf,sol_mid,sol_end);
                end
            else
                sol = sol_mid;
            end
            
        end
    end
    
    methods(Static)


        
        function module_matrix = EXPAND_MODULE_QUANTITY(module_matrix,n_end_modules,prob)
			% EXPAND_MODULE_QUANTITY - Expand the number of modules (columns)
			% in a module matrix (rr or tt).
			%   module_matrix = EXPAND_MODULE_QUANTITY(module_matrix,n_end_modules,prob)
			%   Expand the number of modules of module_matrix to n_end_module
			%   using a probability of reassigment of prob. In other words the
			%   expected number of nodes that are reassigned to the new modules
			%   is equal to prob*size(module_matrix,1);
        
            [n_nodes n_modules] = size(module_matrix);
            
            module_matrix(:,n_modules+1:n_end_modules) = 0;
            new_modules_indexes = ceil( (-1+rand(n_nodes,1)) * n_modules + n_end_modules);

            switched_nodes = find(rand(n_nodes,1)<=prob);

            module_matrix(switched_nodes,:) = 0;
            module_matrix(sub2ind([n_nodes n_end_modules],...
                switched_nodes,new_modules_indexes(switched_nodes))) = 1;


        end
        
    end
    
    
end