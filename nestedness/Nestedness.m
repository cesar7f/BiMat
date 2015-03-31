% Nestedness - Modularity parent class. 
% This class is the main part of all nested metrics.
% It contains all the shared methods of the nestedness algorithms.
%
% Nestedness Properties:
%    N - Nestedness value
%    matrix - Bipartite Adjacency Matrix
%    n_rows - Number of row nodes
%    n_cols - Number of column nodes
%    done - Flag to indicate if the algorith has ben performed
%    print_results - Flag to indicate if result output will be generated
%
% Nestedness Methods:
%    Detect - Main nestedness detection method
%    Print - Print nestedness information
%    nestedness = NTC(matrix) Calculate the nestedness using the NTC metric.    
%    nestedness = NODF(matrix) Calculate the nestedness using the NODF metric.    
%
% See also:
%    NestednessNODF, NestednessNTC
classdef Nestedness < handle
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        N                     = 0;    % Nestedness value
        matrix                = []    % Bipartite Adjacency Matrix
        n_rows                = 0;    % Number of row nodes
        n_cols                = 0;    % Number of column nodes
        %independent_rows_cols = false;% The nestedness algorithm can calculate nestedness for rows and columns independently
        done                  = 0;    % Flag to indicate if the algorith has been performed
    end
    
    properties
        print_results         = Options.PRINT_RESULTS % Flag to indicate if result output will be generated
    end
    
    methods(Abstract, Access = 'protected')
        
        % RunNestedAlgorithm - Abstract method to be implemented in all
        %    Nestedness son classes
        % See NestednessNODF, NestednessNTC
        obj = RunNestedAlgorithm(obj);

    end
    
    methods(Access= 'protected')

        function obj = Nestedness(bipmatrix)
        % Nestedness(bipmatrix) - Main constructor
        %   Can not be directly instantiated    
            obj.matrix = bipmatrix > 0; %Normalize the matrix
        end
        
    end
    
    methods(Access='public')
       
        function obj = Detect(obj)
            
            obj.RunNestedAlgorithm();
            
            if(obj.print_results)
                obj.Print();
            end
            
            obj.done = true;
            
        end
        
    end
    
    methods(Static)
       
        function nested = NTC(matrix)
        % NTC - Calculate the nestedness using the NTC metric
        %
        %   nestedness = NTC(matrix) Calculate the nestedness using the NTC metric.    
        
            nested = NestednessNTC(matrix);
            nested.print_results = false;
            nested.Detect();
            nested.Print();
            
        end

        function nested = NODF(matrix)
        % NODF - Calculate the nestedness using the NODF metric
        %
        %   nestedness = NODF(matrix) Calculate the nestedness using the NODF metric.    
        
            nested = NestednessNODF(matrix);
            nested.print_results = false;
            nested.Detect();
            nested.Print();
            
        end
        
    end
    
end

