%==========================================================================
% Name        : Nestedness.m
% Author      : Cesar Flores
% Created     : 22/Jun/2010
% Updated     : 22/Jun/2010
% Description : Represent the basic structure of a Nestedness algorithm
%==========================================================================
%
% Properties:
%   Network      -   A Network in wich the algorithm will work. Can be
%                    unipartite or bipartite
%
% Abstract Methods:
%
%   CalculateNestedness    -   Calculates the Modularity. It must be
%                              implemented in each subclass.
%
% Notes:
%
%   - THIS IS AN INTERFACE AND THEREFORE MUST NOT BE INSTANTIATED
%

classdef Nestedness < handle

    properties%(GetAccess = 'public', SetAccess = 'private')
        Matrix               = [];
        MatrixMinimal        = [];
        NetworkBipartite     = {};
        nRows                = 0;
        nCols                = 0;
        name                 = 'name'
        IndexRow             = []; %Register the row swaps.
        IndexCol             = []; %Register the columns swaps.
        SortingMethod        = 2; %1 for NTC, 2 for Sum Heuristic
        N                    = 0; %Nested Value
        Fill                 = 0;
        Name                 = '';
        MaxRandomStarts      = 5;                        %How many random starts we will have
        nOnes                = 0;
        nZeros               = 0;
        nRowSorts            = 0;
        nColSorts            = 0;
        DoSorting            = 1;
    end
    
    
    properties(GetAccess = 'private', SetAccess = 'private')
        %Debug                = 0;
    end
    
    methods(Abstract)
        CalculateNestedness(obj);
    end
    
    % SORTING AND MATRIX MANIPULATIONS
    methods
        
        function obj = SortMatrix(obj)
            %fprintf('I sort\n');
            switch obj.SortingMethod
                case 1
                    obj.SortMatrixByNTCHeuristic();
                otherwise
                    obj.SortMatrixBySUMHeuristic();
            end
        end
        
        function obj = SortMatrixByNTCHeuristic(obj)
                        
            %tmsRows = size(n);
            %tmsCols = size(n);
    
            %tmsRows(i) = sum(obj.GetRScore(1));
            %tmsCols(i) = sum(obj.GetRScore(2));

            if(obj.nCols >= obj.nRows) 
                obj.ReorderMatrixColumns();
                obj.ReorderMatrixRows();
            else
                obj.ReorderMatrixRows();
                obj.ReorderMatrixColumns();
            end
        end
        
        function obj = SortMatrixBySUMHeuristic(obj)
            
            sumCols = sum(obj.Matrix);
            sumRows = sum(obj.Matrix,2);
            %sumRows = sum(obj.Matrix');
            
            
            [~,tmpic]=sort(sumCols,'descend');
            [~,tmpir]=sort(sumRows,'descend');

            obj.ReorderCols(tmpic);
            obj.ReorderRows(tmpir);
            
        end
         
       
        function obj = ReorderMatrixColumns(obj)

            r = obj.GetRScore(2); %Get the vector r score in columns
            
            [~, indexr] = sort(r);
            
            obj.ReorderCols(indexr);
            
        end
        
        function obj = ReorderMatrixRows(obj)
            
            r = obj.GetRScore(1); %Get the vector r score in rows
            
            [~, indexr] = sort(r);
            
            obj.ReorderRows(indexr);
        end
        
        function r = GetRScore(obj, rowsOrCols)
            %Get the sequence or actual total r score of the rows or
            %columns.
            
            if(rowsOrCols == 1) %Case of Rows
            
                s = zeros(1,obj.nRows);
                for i = 1:obj.nRows
                    s(i) = 0;
                    for j = 1:obj.nCols
                        if(obj.Matrix(i,j) > 0)
                            s(i) = s(i) + (j*j);
                        end
                    end
                end

                t = zeros(1,obj.nRows);
                for i = 1:obj.nRows
                    t(i) = 0;
                    for j = 1:obj.nCols
                        if(obj.Matrix(i,j) == 0)
                            t(i) = t(i) + (obj.nCols-j+1)^2;
                        end
                    end
                end
                
            else %Case of Columns
                
                s = zeros(1,obj.nCols);
                for j = 1:obj.nCols
                    s(j) = 0;
                    for i = 1:obj.nRows
                        if(obj.Matrix(i,j) > 0)
                            s(j) = s(j) + (i*i);
                        end
                    end
                end

                t = zeros(1,obj.nCols);
                for j = 1:obj.nCols
                    t(j) = 0;
                    for i = 1:obj.nRows
                        if(obj.Matrix(i,j) == 0)
                            t(j) = t(j) + (obj.nRows-i+1)^2;
                        end
                    end
                end
            end
            r = 0.5*t - 0.5*s;
        end
        
        function obj = MirrorMatrix(obj)
            
            mirrorCols = fliplr(obj.IndexCol);
            mirrorRows = filplr(obj.IndexRow);
            
            obj.ReorderCols(mirrorCols);
            obj.ReorderRows(mirrorRows);
            
        end
        
        function obj = RandomizeMatrix(obj)
            %Creates a random permutation of the matrix
           
            colRandom = randperm(obj.nCols);
            rowRandom = randperm(obj.nRows);
            
            obj.ReorderCols(colRandom);
            obj.ReorderRows(rowRandom);
        end
      
        
        function obj = ReorderCols(obj,newIndexes)
             obj.Matrix  = obj.Matrix(:,newIndexes);
           %= temp;
            obj.IndexCol = obj.IndexCol(newIndexes);
        end
        
        function obj = ReorderRows(obj,newIndexes)
            temp = obj.Matrix(newIndexes,:);
            obj.Matrix = temp;
            obj.IndexRow = obj.IndexRow(newIndexes);
        end
        
      
        
        function obj = FixFinalSortIndexes(obj)
           
            for i = fliplr(2:obj.nRowSorts)
                
                a = obj.IndexRow{i-1};
                b = obj.IndexRow{i};
                obj.IndexRow{i-1} = a(b);
                %obj.IndexRow{i-1} = a(obj.IndexRow{i});
                
            end
            
            for i = fliplr(2:obj.nColSorts)
                
                a = obj.IndexCol{i-1};
                b = obj.IndexCol{i};
                obj.IndexCol{i-1} = a(b);
                
            end
            
        end
    end
    
end