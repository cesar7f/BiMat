classdef Options
% Options - Main default values for use in BiMat. Certain functions in
% BiMat include optional arguments. When those optional arguments are
% specified by the user, the values in this class are used instead.
    
    %Statistical TEST
    properties(Constant)
        %These values are used for testing for significance
        P_VALUE    = 0.05; % p_value for statistical test
        Z_VALUE    = 1.96; % z-value for statistical test
    end
    
    %Null Models
    properties(Constant)
        %Default null model for creating random networks. See also: NullModels
        DEFAULT_NULL_MODEL = @NullModels.EQUIPROBABLE; 
        ALLOW_ISOLATED_NODES = 1; %Create random networks that can or can not have empty rows or column nodes
        TRIALS_FOR_NON_EMPTY_NODES = 1000; %When ALLOW_ISOLATED_NODES = 0, this is the maximum number of random networks with empty rows/columns before giving up and relaxing that constraint
        INCLUDE_EMPTY_NODES = 1; %Keep  rows and columns from the matrix that are empty
        REPLICATES = 100;  %Number of random networks (or permutations )during statistical tests
    end
      
    %Algorithms
    properties(Constant)
        OPTIMIZE_COMPONENTS = 0; %Optimize by components or entire matrix in modularity
        % Modularity default algorithm. See also: AdaptiveBrim,
        % LPBrim, LeadingEigenvector
        MODULARITY_ALGORITHM = @AdaptiveBrim; % Modularity default algorithm.
        TRIALS_MODULARITY = 20; %Number of random restarts for AdaptiveBrim and LPBrim algorithms. LeadingEigenvector is deterministic and therefore does not use this variable.
    end
end