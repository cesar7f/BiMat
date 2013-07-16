% NestednessBINMATNEST - NTC (Nestedness Temperature Calculator) algorithm
% This class calculates the nestedness of a matrix using the temperature.
% The value of nestedness if found by normalizing tempereature N = (100-T)
% / 100, and have values between 0 and 1, where 1 is perfectly nested and 0\
% perfectly anti-nested. For information about how this algorithm works you
% can consult the following papers:
%
%     Atmar, Wirt and Patterson, Bruce D. The measure of order and disorder
%     in the distribution of species in fragmented habitat. Oecologia 1993
%
%     Rodriguez-Girones, Miguel A and Santamaria, Luis. A new algorithm to
%     calculate the nestedness temperature of presence--absence matrices.
%     Journal of Biogeography 2006
%
% NestednessBINMATNEST Properties:
%     P - p parameter of the isocline function that is calculated using the fill of the matrix
%     PositionMatrixX - X coordinate of each matrix cell in a unit square
%     PositionMatrixY - Y coordinate of each matrix cell in a unit square
%     dMatrix - distances with the perfect nestedness line
%     DMatrix - Size of the diagonal that cross the matrix element
%     uMatrix - Unexpectedness Matrix
%     T - T emperature
%     CalculatedFill - Calculated fill using the integral of Fxp(P,X). Ideally must have the same value than the fill of the matrix.
%     X - Vector of X coordinate
%     Fxp - Vector of isocline values in y coordinate
%     UMin - Unexpectedness matrix
%     tunsorted - 
%     DoGeometry - Flag to indicate the calculus of geometry
%     done - Flag to indicate if the algorith has ben performed
%     PMax - Maximal P value for finding the isoclane.
%     PMin - Minimal p value for finding the isoclane.
%     UsedArea - Chose a value in 1,2,3
%     DeltaX - X Increment in order to get the vector of the Isoclane values (obj.Fxp)
%     DebugMessages - 1,0 Print Debug Messages
%     K  - Constant to normalize values of T in [0,100]
%     BreakRandom - How many initial random permutations in the matrix
%
% NestednessBINMATNEST Methods:
%    NestednessBINMATNEST - Main Constructor
%    SetMatrix - Change the matrix of the algorithm
%    CalculateNestedness - Main method calculating Nestedness
%    CalculateMatrixGeometry - Calculate all the geometry aspects
%    AssignMatrixPositions - Map the matrix elements to a unit
%    FindPValue - Find the parameter 'p' of the isocline function
%    GetFilledArea - Get the area above the isocline
%    CalculateDiagonalsAndDistances - Calculate diagonal and isocline distance size matrices
%    CalculateTemperature - Calculate the matrix temperature
%    CalculateUnexpectedness - Calculate the matrix unexpectedness
%    PERFECT_NESTED - Return a perfect nested matrix according to the NTC algorithm
%    FIND_UNEXPECTED_CELLS - Return a matrix that indicate what are the unexpected cells.
%    GET_ISOCLINE - Get the isocline function
%
% See also:
%    NODF
classdef NestednessBINMATNEST < Nestedness

    properties(GetAccess = 'public')%, SetAccess = 'private')
        P                  = 0;  % p parameter of the isocline function that is calculated using the fill of the matrix
        PositionMatrixX    = []; % X coordinate of each matrix cell in a unit square
        PositionMatrixY    = []; % Y coordinate of each matrix cell in a unit square
        dMatrix            = []; % distances with the perfect nestedness line
        DMatrix            = []; % Size of the diagonal that cross the matrix element
        uMatrix            = []; % Unexpectedness Matrix
        T                  = 0;  % T emperature
        CalculatedFill     = 0;  % Calculated fill using the integral of Fxp(P,X). Ideally must have the same value than the fill of the matrix.
        X                  = []; % Vector of X coordinate
        Fxp                = []; % Vector of isocline values in y coordinate
        UMin               = 0;  % Unexpectedness matrix
        tunsorted          = 0;  
        DoGeometry         = 1;  % Flag to indicate the calculus of geometry
        done               = 0;  % Flag to indicate if the algorith has ben performed.
    end
    %DEBUG Properties - Change to parametrize and Debug the algorithm;
    properties(GetAccess = 'public', SetAccess = 'private')
        PMax               = 99999;      % Maximal P value for finding the isoclane.
        PMin               = 0.0005;     % Minimal p value for finding the isoclane.
        UsedArea           = 2;          % Chose a value in 1,2,3
        DeltaX             = 0.001;      % X Increment in order to get the vector of the Isoclane values (obj.Fxp)
        DebugMessages      = 0;          % 1,0 Print Debug Messages
        K                  = 2.4125e+003 % 100 / 0.04145;  %Value found in the literature
        BreakRandom        = 20;         % How many initial random permutations in the matrix
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
    methods
        function obj = NestednessBINMATNEST(bipNetwork)
        % NestednessBINMATNEST - Main Constructor
        %   obj = NestednessBINMATNEST(bipNetwork) Create a NestednessBINMATNEST object
        %   called bp using a bipartite adjacency matrix or bipartite
        %   object named bipNetwork
            
            %Constructor method
            if(isa(bipNetwork, 'Bipartite'))
                obj.NetworkBipartite = bipNetwork;
                obj.Matrix = bipNetwork.webmatrix > 0; %Normalize the matrix
                obj.nRows = bipNetwork.n_rows;  %Number of Rows
                obj.nCols = bipNetwork.n_cols; %Number of Columns
                obj.Fill = bipNetwork.connectance;
            elseif(isa(bipNetwork, 'double')||isa(bipNetwork, 'logical'))
                matrix = bipNetwork;
                obj.Matrix = matrix > 0; %Normalize the matrix
                [obj.nRows obj.nCols] = size(matrix);
                obj.Fill = sum(sum(obj.Matrix))/(obj.nRows*obj.nCols);   
            end
            
            obj.IndexRow = 1:obj.nRows;
            obj.IndexCol = 1:obj.nCols;
            obj.MaxRandomStarts = 50;
            obj.UsedArea = 2;
            obj.BreakRandom     = 10;
            obj.SortingMethod        = 2;
        end
        
        function obj = SetMatrix(obj,matrix)
        % SetMatrix - Change the matrix of the algorithm
        %   obj = SetMatrix(obj,matrix) Change the matrix in which the
        %   algorithm will be performed. Useful only when the new matrix
        %   has the same size and similar connectance (fill) than a
        %   the old matrix. Usint this method, we do not have to perform
        %   the goemetrical pre-calculus another time (isocline, distances,
        %   diagonals, etc).
            obj.Matrix = matrix > 0; %Normalize the matrix
            [obj.nRows obj.nCols] = size(matrix);
            obj.Fill = sum(sum(obj.Matrix))/numel(obj.Matrix);
            obj.IndexRow = 1:obj.nRows;
            obj.IndexCol = 1:obj.nCols;
            obj.MaxRandomStarts = 50;
        end
        
        function obj = CalculateNestedness(obj)
        % CalculateNestedness - Main method calculating Nestedness
        % Temperature Calculator
        %   obj = CalculateNestedness(obj) Calculates the nestedness of the
        %   matrix. Use obj.N after calling this method to get the
        %   nestedness value, and obj.T for getting the temperature value.
        
            if(obj.nRows==1 || obj.nCols==1)
                obj.N = 0;
                return;
            end
            
            obj.nRowSorts = 1;
            obj.nColSorts = 1;
            
            % Perfrom the calculus of geometry (isocline, distances, etc)
            if(obj.DoGeometry)
                obj.CalculateMatrixGeometry();
            end

            %Calculate the temperature
            obj.CalculateTemperature();

            %Normalize the temperature to the nestedness value
            obj.N = (100-obj.T)/100;
            
            %If you want the calculus for a unsorted matrix you are
            %finished. Normally obj.DoSorting = 1, such that you want an optimal ordering. 
            if(obj.DoSorting == 0)
                return;
            end
                        
            % The next part focus in finding the ordering, such that you
            % will have the smaller possible value of temperature (and by
            % consequence the highest nestedness value)
            
            globalMinimalT = 500;
            matrixLocalMinima = [];
            indexRowLocalMinima = [];
            indexColLocalMinima = [];
            indexRowGlobalMinima = [];
            indexColGlobalMinima = [];
            
            
            failedtoincrease = 0; %Count if the next matrix randomization do an improvement
            % Do obj.MaxRandomStarts initial random permutations of the
            % matrix to be tested
            for i = 1:obj.MaxRandomStarts
                
                % If no increase is detected in obj.BreakRandom continuos
                % trials, no need for continue looking.
                if(failedtoincrease > obj.BreakRandom)
                    %fprintf('Break on i = %i\n',i);
                    break;
                end
                
                
                permutationMinimalT = 500; %temperature infinite
                obj.T = 500;
                
                obj.RandomizeMatrix();
%                i = 1;
                while(1)
                    %display(i);
                    %i = i+1;
                    obj.SortMatrix();
                    obj.CalculateTemperature();
                    
                    if(obj.DebugMessages == 1); fprintf('TLocal = %f T = %f\n', permutationMinimalT,obj.T); end;
                    
                    if(abs(permutationMinimalT - obj.T) <= 0.001 || obj.T > permutationMinimalT)
                        break;
                    end
                        
                    if(obj.T < permutationMinimalT)
                        permutationMinimalT = obj.T;
                        matrixLocalMinima = obj.Matrix;
                        indexRowLocalMinima = obj.IndexRow;
                        indexColLocalMinima = obj.IndexCol;
                    end
                    
                end
                if(obj.DebugMessages == 1); fprintf('finalizo ciclo\n'); end;
                
                %Save if permutation is smaller than the global minimal
                if(permutationMinimalT < globalMinimalT)
                    %fprintf('TMinimalGlob = %f\n', permutationMinimalT);
                    globalMinimalT = permutationMinimalT;
                    obj.MatrixMinimal = matrixLocalMinima;
                    indexRowGlobalMinima = indexRowLocalMinima;
                    indexColGlobalMinima = indexColLocalMinima;
                    failedtoincrease = 0;
                end
                
                failedtoincrease = failedtoincrease + 1;
            end
            
           %Keep the best sorting for NTC
           obj.Matrix = obj.MatrixMinimal;
           obj.IndexRow = indexRowGlobalMinima;
           obj.IndexCol = indexColGlobalMinima;
           obj.T = globalMinimalT;
           obj.N = (100-obj.T)/100;
            
           obj.done = 1;
           %obj.PrintOutput();
            
        end
           
    end

    % GEOMETRY DEFINITION SECTION
    methods
       
        function obj = CalculateMatrixGeometry(obj)
            % CalculateMatrixGeometry - Calculate all the geometry aspects
            % of the algorithm
            % obj = CalculateMatrixGeometry(obj)
            % This function calculate all the matrix geometry in the next
            % order:
            %   1.- Coordinate representation of of the matrix elements in a unit scuare.
            %   2.- Function of the isoclane f(x;p) based in the matrix density Fill
            %   3.- Main diagonal size for all the the matrix elements.
            %   4.- Distance along the main diagonal of all the matrix
            %   elements 
         
            %1.-Coordinate representation of of the matrix elements in a unit scuare.
            obj.AssignMatrixPositions();
            
            %2.- Function of the isoclane f(x;p) based in the matrix density Fill
            obj.X = (0.5/obj.nCols):obj.DeltaX:((obj.nCols-0.5)/obj.nCols); %Define the X Vector of the function
            obj.P = obj.FindPValue();
            obj.Fxp = 0.5/obj.nRows + ((obj.nRows-1)/obj.nRows) * (1-(1-(obj.nCols*(obj.X)-0.5)/(obj.nCols-1)).^(obj.P)).^(1/(obj.P));
            %3.-,4.-
            obj.CalculateDiagonalsAndDistances();
            
        end
        
        function obj = AssignMatrixPositions(obj)
            % AssignMatrixPositions - Map the matrix elements to a unit
            % square coordinate system.
            % obj = AssignMatrixPositiong(obj) - Map the matrix elements to a unit
            % square coordinate system.
            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    obj.PositionMatrixX(i,j) = (j-0.5)/obj.nCols;
                    obj.PositionMatrixY(i,j) = (obj.nRows-i+0.5)/obj.nRows;
                end 
            end
        end
        
        function p = FindPValue(obj)
            % FindPValue - Find the parameter 'p' of the isocline function
            % p = FindPValue(obj) - Get the parameter p of the isocline
            % function by doing a search in the p space and doing a
            % bisection method at the end, such that the area above the
            % isocline is the same than the connectance (fill) of the matrix.
            
            p = obj.PMin; %Starting with the minimal pre-defined value of p parameter
            filledarea = 0; %Area above the curve. The objective is to equalize to obj.Fill.
            while(p < obj.PMax) %After some predefined PMax the increase in p will not affect the form of the isocline
                filledarea = obj.GetFilledArea(p);
                if(obj.Fill > filledarea) 
                    break;
                end
                if(obj.DebugMessages); fprintf('area = %5.4f p = %5.4f\n', filledarea,p); end; 
                p = p*2;          
            end
            
            %if(obj.DebugMessages); fprintf('area = %10.9f p = %5.4f lastp = %5.4f\n', filledarea,upp,lowp); end;
            
            if(p < obj.PMax && p > obj.PMin) %If the parameter p is not an extreme case
                %BISECTION METHOD
                upp = p;
                lowp = p / 2;
                mid = 0;
                while( abs( obj.Fill - filledarea) > 0.001)
                    mid = (upp + lowp)/2;
                    filledarea = obj.GetFilledArea(mid);
                    if(filledarea < obj.Fill)
                        upp = mid;
                    else
                        lowp = mid;
                    end
                    if(obj.DebugMessages); fprintf('area = %10.9f p = %f\n', filledarea,mid); end;        
                end
                if(mid ~= 0)
                    p = mid;
                end;
            end
            
            obj.CalculatedFill = filledarea;
        end
        
        function Area = GetFilledArea(obj,p)
            % GetFilledArea - Get the area above the isocline
            %   Area = GetFilledArea(obj,p) - Ghet the area above the
            %   isocline with parameter p.
            
            
            %Isocline equation
            obj.Fxp = 0.5/obj.nRows + ((obj.nRows-1)/obj.nRows) * (1-(1-(obj.nCols*(obj.X)-0.5)/(obj.nCols-1)).^p).^(1/p);
            
            %Area below the isocline
            integral = trapz(obj.X,obj.Fxp);
            
            %Three ways of calculating the area (only important when the
            %matrix is small. Case 2 gives the best results.
            switch obj.UsedArea
                case 1
                    Area = 1 - real(integral);
                case 2
                    Area = 1 - real(integral) - (obj.nRows-0.5)*(0.5)/(obj.nRows*obj.nCols);
                otherwise
                    Area = (obj.nRows-0.5)/(obj.nRows-1) - real(integral) * obj.nRows * obj.nCols / ((obj.nCols-1)*(obj.nRows-1));
            end
        end
        
        function obj = CalculateDiagonalsAndDistances(obj)
            % CalculateDiagonalsAndDistances - Calculate diagonal and
            % isocline distance size matrices
            %   obj = CalculateDiagonalsAndDistances(obj) - Calculate diagonal and
            %   isocline distance size matrices
            obj.uMatrix = zeros(size(obj.Matrix));
            MaxDiag = sqrt(2);
            
            obj.DMatrix = zeros(size(obj.Matrix));
            obj.dMatrix = zeros(size(obj.Matrix));
            
            %For each row and column
            for i = 1:obj.nRows
                for j = 1:obj.nCols
                         
                    y1 = real(obj.PositionMatrixX(i,j) + obj.PositionMatrixY(i,j) - obj.X);
                    y2 = obj.Fxp;
                    
                    [~, index] = min(abs(y1-y2));
                    
                    %Intersection point between the diagonal and the
                    %iscoline
                    ycross = y1(index);
                    xcross = obj.X(index);

                    %Distance from the isocline to the matrix element
                    distance = sqrt( (obj.PositionMatrixX(i,j)-xcross)^2 + (obj.PositionMatrixY(i,j)-ycross)^2 );
                    obj.dMatrix(i,j) = distance;
                    obj.DMatrix(i,j) = (obj.PositionMatrixX(i,j) + obj.PositionMatrixY(i,j)) * sqrt(2);

                    if(obj.DMatrix(i,j) > MaxDiag)
                        obj.DMatrix(i,j) = abs(obj.PositionMatrixX(i,j) + obj.PositionMatrixY(i,j) - 2) * sqrt(2);
                    end
                    
                    % Change to negative elements below isocline, such that
                    % the sign will differentiate above vs below isocline
                    % elemnts.
                    if(obj.PositionMatrixY(i,j) < ycross)
                        obj.dMatrix(i,j) = -obj.dMatrix(i,j);
                    end
                end
            end       
        end 
    end
    
    
    
    % CALCULATE TEMPERATURE AND IMPORTANT VALUES
    methods
        
        function obj = CalculateTemperature(obj)
            % CalculateTemperature - Calculate the matrix temperature
            %   obj = CalculateTemperature(obj) - Calculate the temperature
            %   using the Atmar standard equation. The temperature is in
            %   the interval [0,100], while the nestedness in the interval
            %   [0,1]. High values of temperature corresponds to low values
            %   of nestedness.
            obj.UMin = obj.CalculateUnexpectedness();
            obj.T = obj.K*obj.UMin;
            obj.N = (100-obj.T)/100;
        end
        
        function unex = CalculateUnexpectedness(obj)
            % CalculateUnexpectedness - Calculate the matrix unexpectedness
            %   obj = CalculateUnexpectedness(obj) - Sum all temperature
            %   contributions from unexpected cells (absences below the
            %   matrix and presences above the matrix)
            obj.uMatrix = zeros(size(obj.Matrix));
            obj.uMatrix = ((obj.Matrix==0 & obj.dMatrix > 0) | (obj.Matrix ~=0 & obj.dMatrix < 0 )).*((obj.dMatrix./obj.DMatrix).^2);   
            unex = sum(sum(obj.uMatrix)) / (obj.nRows*obj.nCols);
        end
    end
    
    % PLOTTING AND OUTPUT
    methods
        
        function obj = PrintOutput(obj)
            %fprintf('NESTEDNESS ANALYSIS IN NETWORK : %s\n', obj.Name);
            %fprintf('T = %5.4f UMin = %5.4f UMax = %5.4f UNorm = %5.4f Fill = %5.4f P = %5.4f\n',...
            %        obj.T, obj.UMin, obj.UMax, obj.UNorm, obj.Fill, obj.P);
            
            fprintf('T = %5.4f Size = %i by %i\n', obj.T, obj.nRows, obj.nCols);
        end

    end
    
    methods(Static)
        
        function matrix = PERFECT_NESTED(nrows,ncols,fill)
        % PERFECT_NESTED - Return a perfect nested matrix according to the
        % NTC algorithm
        %   matrix = PERFECT_NESTED(nrows,ncols,fill) - Return a perfect
        %   nested matrix of size nrows by ncols and a connectance = fill.
        %   The perfect nested matrix follows the definition of the NTC
        %   algorithm (the isocline divide ones from zeros in the entire
        %   matrix).
            matrix = zeros(nrows,ncols);
        
            bnest = NetworkBipartite(matrix);
            
            nest = NestednessBINMATNEST(bnest);
            nest.Fill = fill;
            nest = nest.CalculateMatrixGeometry();
            nest = nest.CalculateDiagonalsAndDistances();
            
            for i = 1:nrows
                for j = 1:ncols
                    if(nest.dMatrix(i,j) > 0)
                        matrix(i,j) = 1;
                    else
                        matrix(i,j) = 0;
                    end
                end
            end
            
            matrix(nrows,1) = 1;
            matrix(1,ncols) = 1;
           
        end
        
        function matrix_unex = FIND_UNEXPECTED_CELLS(matrix)
        % FIND_UNEXPECTED_CELLS - Return a matrix that indicate what are the
        % unexpected cells.
        %   matrix = FIND_UNEXPECTED_CELLS(matrix) - For matrix 'matrix',
        %   calculate the geometry in order to return a matrix 'matrix_unex'
        %   with ones in the position of unexpected cells of the original
        %   matrix.
            nest = NestednessBINMATNEST(matrix);
            nest.CalculateMatrixGeometry();
            nest.CalculateUnexpectedness();
            
            matrix_unex = nest.uMatrix > 0;%.0005;
            
        end
        
        function [x y] = GET_ISOCLINE(n_rows,n_cols,p_value)
        % GET_ISOCLINE - Get the isocline function
        %   [x y] = GET_ISOCLINE(n_rows,n_cols,p_value) - Get the isocline
        %   function in x and y vectors for a matrix of size n_rows by
        %   n_cols and a connectance of p_value. Useful when the user is
        %   only interested in the isocline (e.g. plotting) and not the
        %   temperature value.
            if(nargin==1)
                matrix = n_rows>0;
            else
                matrix = zeros(n_rows,n_cols);
                len = n_rows*n_cols;
                matrix(1:round(len*p_value))=1;
            end
            
            nest = NestednessBINMATNEST(matrix);
            
            nest.CalculateMatrixGeometry();
            x = 0.5 + nest.nCols.*nest.X;
            y = 0.5 + nest.nRows.*nest.Fxp;
        end
        
    end
end
