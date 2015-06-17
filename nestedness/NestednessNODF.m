% NODF - Calculate normalized nodf value of a matrix. To know how this
% nestedness metric works you can consult the following paper:
%
%   Almeida-Neto, Mario and Guimaraes, Paulo and Guimaraes, Paulo R and 
%   Loyola, Rafael D and Ulrich, Werner. A consistent metric for nestedness
%   analysis in ecological systems: reconciling concept and measurement.
%   Oikos 2008
%
% NestednessNTC Properties:
%     N_rows - NODF value for rows
%     N_cols - NODF value for columns
%
% NestednessNTC Methods:
%     NestednessNODF - Main Constructor
%     Print - Print NODF nestedness information
%     NODF - Calculate the NODF nestedness o a matrix
%
% See also:
%    NestednessNTC, Nestedness
classdef NestednessNODF < Nestedness

    properties(GetAccess = 'public', SetAccess = 'protected')
        N_rows   = 0; % NODF value for rows
        N_cols   = 0; % NODF value for columns
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
    methods
        function obj = NestednessNODF(bipmatrix)
        % NestednessNODF - Main Constructor
        % 
        %   obj = NestednessNODF(MATRIX) Creates an NestednessNODF object obj
        %   using a bipartite adjacency matrix MATRIX that will be used to
        %   calculate nestedes using the NODF metric.
        %
        % See also:
        %    NestednessNTC
            
            obj = obj@Nestedness(bipmatrix);
            
            %obj.independent_rows_cols = true;
        end
       
        
        function obj = Print(obj,filename)
        % Print - Print NODF nestedness information
        %
        %   STR = Print(obj) Print the NODF information to screen and
        %   return this information to the string STR
        %
        %   STR = Print(obj, FILE) Print the NODF information to screen and
        %   text file FILE and return this information to the string STR   
        %
        % See also: 
        %   Printer     
            str = 'Nestedness NODF:\n';
            str = [str, '\tNODF (Nestedness value):    \t', sprintf('%20.4f',obj.N), '\n'];
            str = [str, '\tNODF (Rows value):          \t', sprintf('%20.4f',obj.N_rows), '\n'];
            str = [str, '\tNODF (Columns value):       \t', sprintf('%20.4f',obj.N_cols), '\n'];
            fprintf(str);  
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
           
        
        
    end

    methods(Access='protected')
       
                function obj = RunNestedAlgorithm(obj)
        % RunNestedAlgorithm - Main method for calculating NODF nestedness
        % values
        %
        %   obj = RunNestedAlgorithm(obj) Calculates the nestedness of the
        %   matrix. Use obj.N after calling this method to get the
        %   nestedness value. Aditionally, you can use obj.N_rows and
        %   obj.N_cols for the contribution of row and column to
        %   nestedness
        
            [obj.N,obj.N_rows,obj.N_cols] = NODF(obj.matrix);

        end
        
    end
    
    methods(Static)    
        
        function nodf_obj = NODF(matrix)
        % NODF - Calculate the NODF nestedness
        %
        %   nest = NODF(MATRIX) Calculate the NODF nestedness of MATRIX,
        %   print the basic information to
        %   screen and return an NestednessNODF object that contains such
        %   information in nest.
        
            nodf_obj = NestednessNODF(matrix);
            nodf_obj.Detect();
            nodf_obj.Print();
            
        end
        
    end
end
