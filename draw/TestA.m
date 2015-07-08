classdef TestA
      
    %Algorithms
    properties(Constant)
        OPTIMIZE_COMPONENTS = 0; %Optimize by components or entire matrix in modularity
        % Modularity default algorithm. See also: AdaptiveBrim,
        % LPBrim, LeadingEigenvector
        MODULARITY_ALGORITHM = @AdaptiveBrim; % Modularity default algorithm.
        % Nestedness default algorithm. See also: NestednessNTC,
        % NestednessNODF
        NESTEDNESS_ALGORITHM = @NestednessNODF;
        TRIALS_MODULARITY = 20; %Number of random restarts for AdaptiveBrim and LPBrim algorithms. LeadingEigenvector is deterministic and therefore does not use this variable.
    end
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = Options

        end
    end 
end