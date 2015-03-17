% Diversity - Static class that share diversity functions used for calculating
% correlation between rows/columns labels and modularity.
%
% Diversity Methods:
%    SHANNON_INDEX - Shannon's diversity index
%    SIMPSON_INDEX - Simpson's diversity index
%
% See also:
%    InternalStatistics
classdef Diversity
   
    methods(Static)
            
        function H = SHANNON_INDEX(species,un)
            % SHANNON_INDEX - Shannon's diversity index
            %   H = SHANNON_INDEX(species) Calculate the Shannon's index H of
            %   vector species.
            %   H = SHANNON_INDEX(species,un) Calculate the shannon's index H of
            %   vector species using only specific indexes un
            %                 H = \sum_{i=1}^n - {p_i} \log p_i
            if(nargin==1)
                un = unique(species);
            end
            ni = length(un);
            H = 0;
            N = length(species);
            
            for i = 1:ni
               
                p = sum(species == un(i))/N;
                H = H - p * log(p);
                
            end
            
        end
        
        function D = SIMPSON_INDEX(species,un)
            % SIMPSON_INDEX - Simpson's diversity index
            %   H = SIMPSON_INDEX(species) Calculate the Simpson's index H of
            %   vector species.
            %   H = SIMPSON_INDEX(species,un) Calculate the Simpson's index H of
            %   vector species using only specific indexes un
            
            if(nargin==1)
                un = unique(species);
            end
            ni = length(un);
            D = 0;
            N = length(species);
            
            for i = 1:ni
                n = sum(species == un(i));
                D = D + (n*(n-1))/(N*(N-1));
                
            end
            D = 1 - D;
        end

        
    end
    
end