classdef Options
    %Statistical TEST
    properties(Constant)
        P_VALUE    = 0.05; %p_value for statistical test
        TWO_TAIL   = 1;    %Two or one tail for statistical test
        Z_VALUE    = 1.96; %z-value for statistical test
    end
    
    %Null Models
    properties(Constant)
        DEFAULT_NULL_MODEL = @NullModels.EQUIPROBABLE; %Default null model for creating random networks
        ALLOW_ISOLATED_NODES = 0; %Create random networks that can or can not have empty rows or column nodes
        TRIALS_FOR_NON_EMPTY_NODES = 1000; %When ALLOW_ISOLATED_NODES = 0, this is the maximum number of random networks with empty rows/columns before giving up
        INCLUDE_EMPTY_NODES = 0; %Delete rows and columns from the matrix (only for analysis) that are empty
        REPLICATES = 100;  %Number of random networks during tests
    end
    
    %Algorithms
    properties(Constant)
        OPTIMIZE_COMPONENTS = 0; %Optimize by components or entire matrix in modularity
        MODULARITY_ALGORITHM = @AdaptiveBrim;
    end
end