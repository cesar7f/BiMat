% NestednessNTC - NTC (Nestedness Temperature Calculator) algorithm
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
% NestednessNTC Properties:
%     T - T emperature
%     do_geometry - Flag to indicate the calculus of geometry
%     index_rows - Register of the swaps in Rows.
%     index_cols - Register of the swaps in Cols.
%     connectance - Fill of the matrix
%     trials - Number of random initializations
%     do_sorting - Sort the matrix before calculating the temperature
%
% NestednessNTC Methods:
%     NestednessNTC - Main Constructor
%     SetMatrix - Change the adjacency matrix of the algorithm
%     Print - Print NTC nestedness information
%     NTC - Calculate the NTC nestedness
%     PERFECT_NESTED - Return a perfect nested matrix according to the NTC algorithm
%     FIND_UNEXPECTED_CELLS - Return a matrix that indicates what are the unexpected cells.
%     GET_ISOCLINE - Get the isocline function
%
% See also:
%    NestednessNODF, Nestedness
classdef NestednessNTC < Nestedness

    properties(GetAccess = 'public', SetAccess = 'protected')
        T                  = 0;    % T emperature
        do_geometry        = true; % Flag to indicate the calculus of geometry
        index_rows         = [];   % Register of the swaps in Rows.
        index_cols         = [];   % Register of the swaps in Cols.
        connectance        = 0;    % Fill of the matrix
        trials             = 5;    % Number of random initializations
        do_sorting         = true; % Sort the matrix before calculating the temperature
    end
    
    properties(Access = 'protected')
        matrix_minimal     = []; % The matrix with the smalles temperature (highest nestedness)
        P                  = 0;  % p parameter of the isocline function that is calculated using the fill of the matrix
        pos_x_matrix       = []; % X coordinate of each matrix cell in a unit square
        pos_y_matrix       = []; % Y coordinate of each matrix cell in a unit square
        d_matrix           = []; % distances with the perfect nestedness line
        diag_matrix        = []; % Size of the diagonal that cross the matrix element
        u_matrix           = []; % Unexpectedness Matrix
        calculated_fill    = 0;  % Calculated fill using the integral of fxp(P,X). Ideally must have the same value than the fill of the matrix.
        X                  = []; % Vector of X coordinate
        fxp                = []; % Vector of isocline values in y coordinate
        u_min              = 0;  % Unexpectedness matrix
        tunsorted          = 0;  
        sorting_method     = 2; %1 for NTC, 2 for Sum Heuristic
        n_row_sorts        = 0;
        n_col_sorts        = 0;
    end
    %DEBUG Properties - Change to parametrize and Debug the algorithm;
    properties(Access = 'protected')
        p_max               = 99999;      % Maximal P value for finding the isoclane.
        p_min               = 0.0005;     % Minimal p value for finding the isoclane.
        used_area           = 2;          % Chose a value in 1,2,3
        delta_x             = 0.001;      % X Increment in order to get the vector of the Isoclane values (obj.fxp)
        debug_messages      = 0;          % 1,0 Print Debug Messages
        K                   = 2.4125e+003 % 100 / 0.04145;  %Value found in the literature
        break_random        = 20;         % How many initial random permutations in the matrix
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
    methods

        
        function obj = NestednessNTC(bipmatrix)
        % NestednessNTC - Main Constructor
        % 
        %   obj = NestednessNTC(MATRIX) Creates an NestednessNTC object obj
        %   using a bipartite adjacency matrix MATRIX that will be used to
        %   calculate nestedes using the Nestedness Temperature Calculator
        %   (NTC).
        %
        % See also:
        %    NestednessNODF
            obj = obj@Nestedness(bipmatrix);
            
            [obj.n_rows, obj.n_cols] = size(obj.matrix);
            obj.connectance = sum(sum(obj.matrix))/(obj.n_rows*obj.n_cols);   
                        
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            
        end

        function obj = SetMatrix(obj,matrix)
        %   obj = SetMatrix(obj,matrix) Change the matrix in which the
        %   algorithm will be performed. Useful only when the new matrix
        %   has the same size and similar connectance (fill) than a
        %   the old matrix. Using this method, we do not have to perform
        %   the goemetrical pre-calculus another time (isocline, distances,
        %   diagonals, etc).
        %   
        %   Be cautious while changing the matrix. The change of the matrix
        %   will cause that NestednessNTC will not recalculate the values
        %   of geometry and distances. Therefore, the matrix must have the
        %   same size and similar connectance. Otherwise the results will
        %   be spurious.
            
            %Make sure that at least the size of the new matrix is the same
            assert(sum(size(matrix) == size(obj.matrix))==2);
            
            obj.matrix = matrix > 0; %Normalize the matrix
            [obj.n_rows, obj.n_cols] = size(matrix);
            obj.connectance = sum(sum(obj.matrix))/numel(obj.matrix);
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            obj.do_geometry = 0;

        end

        
        function str = Print(obj,filename)
        % Print - Print NTC nestedness information
        %
        %   STR = Print(obj) Print the NTC information to screen and
        %   return this information to the string STR
        %
        %   STR = Print(obj, FILE) Print the NTC information to screen and
        %   text file FILE and return this information to the string STR   
        %
        % See also: 
        %   Printer
            
            str = 'Nestedness NTC:\n';
            str = [str, '\tNTC (Nestedness value):     \t', sprintf('%20.4f',obj.N), '\n'];
            str = [str, '\tT (Temperature value):      \t', sprintf('%20.4f',obj.T), '\n'];
           
            fprintf(str);  
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
           
    end

    % GEOMETRY DEFINITION SECTION
    methods(Access = 'protected')
       
        function obj = RunNestedAlgorithm(obj)
        % RunNestedAlgorithm - Main method for calculating NTC nestedness
        % Temperature Calculator
        %
        %   obj = RunNestedAlgorithm(obj) Calculates the nestedness of the
        %   matrix. Use obj.N after calling this method to get the
        %   nestedness value, and obj.T for getting the temperature value.
        
            if(isempty(obj.matrix))
                obj.N= NaN;
                obj.T = NaN;
                return;
            end
        
            if(obj.n_rows==1 || obj.n_cols==1)
                obj.N= 0;
                return;
            end
            
            obj.n_row_sorts = 1;
            obj.n_col_sorts = 1;
            
            % Perfrom the calculus of geometry (isocline, distances, etc)
            if(obj.do_geometry)
                obj.CalculateMatrixGeometry();
            end

            %Calculate the temperature
            obj.CalculateTemperature();

            %Normalize the temperature to the nestedness value
            obj.N= (100-obj.T)/100;
            
            %If you want the calculus for a unsorted matrix you are
            %finished. Normally obj.do_sorting = 1, such that you want an optimal ordering. 
            if(obj.do_sorting == 0)
                return;
            end
                        
            % The next part focus in finding the ordering, such that you
            % will have the smaller possible value of temperature (and by
            % consequence the highest nestedness value)
            
            globalMinimalT = obj.T;
            obj.matrix_minimal = obj.matrix;
            indexRowGlobalMinima = obj.index_rows;
            indexColGlobalMinima = obj.index_cols;
            
            % Do obj.trials initial random permutations of the
            % matrix to be tested
            for i = 1:obj.trials
                
                [matrix_permuted,new_row_index,new_col_index] = MatrixFunctions.RANDOM_SORT(obj.matrix);
                obj.matrix = matrix_permuted;
                obj.index_rows(new_row_index);
                obj.index_cols(new_col_index);

                [sort_matrix, new_row_index, new_col_index] = MatrixFunctions.SORT_MATRIX(obj.matrix);
                obj.matrix = sort_matrix;
                obj.index_rows = obj.index_rows(new_row_index);
                obj.index_cols = obj.index_cols(new_col_index);

                obj.CalculateTemperature();

                if(obj.debug_messages == 1); fprintf('T = %f Tglobal = %f\n', obj.T, globalMinimalT); end;

                
                %Save if permutation is smaller than the global minimal
                if(obj.T < globalMinimalT)
                    %fprintf('TMinimalGlob = %f\n', permutationMinimalT);
                    globalMinimalT = obj.T;
                    obj.matrix_minimal = obj.matrix;
                    indexRowGlobalMinima = obj.index_rows;
                    indexColGlobalMinima = obj.index_cols;
                end

            end
            
           %Keep the best sorting for NTC
           obj.matrix = obj.matrix_minimal;
           obj.index_rows = indexRowGlobalMinima;
           obj.index_cols = indexColGlobalMinima;
           obj.T = globalMinimalT;
           obj.N= (100-obj.T)/100;
            
           obj.done = 1;
           %obj.PrintOutput();
            
        end

        
        function obj = CalculateMatrixGeometry(obj)
        % CalculateMatrixGeometry - Calculate all the geometry aspects
        % of the algorithm
        %
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
            obj.X = (0.5/obj.n_cols):obj.delta_x:((obj.n_cols-0.5)/obj.n_cols); %Define the X Vector of the function
            obj.P = obj.FindPValue();
            obj.fxp = 0.5/obj.n_rows + ((obj.n_rows-1)/obj.n_rows) * (1-(1-(obj.n_cols*(obj.X)-0.5)/(obj.n_cols-1)).^(obj.P)).^(1/(obj.P));
            %3.-,4.-
            obj.CalculateDiagonalsAndDistances();
            
        end
        
        function obj = AssignMatrixPositions(obj)
        % AssignMatrixPositions - Map the matrix elements to a unit
        % square coordinate system.
        %
        %   obj = AssignMatrixPositiong(obj) - Map the matrix elements to a unit
        %   square coordinate system.
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    obj.pos_x_matrix(i,j) = (j-0.5)/obj.n_cols;
                    obj.pos_y_matrix(i,j) = (obj.n_rows-i+0.5)/obj.n_rows;
                end 
            end
        end
        
        function p = FindPValue(obj)
            % FindPValue - Find the parameter 'p' of the isocline function
            %   p = FindPValue(obj) - Get the parameter p of the isocline
            %   function by doing a search in the p space and doing a
            %   bisection method at the end, such that the area above the
            %   isocline is the same than the connectance (fill) of the matrix.
            
            p = obj.p_min; %Starting with the minimal pre-defined value of p parameter
            filledarea = 0; %Area above the curve. The objective is to equalize to obj.connectance.
            while(p < obj.p_max) %After some predefined p_max the increase in p will not affect the form of the isocline
                filledarea = obj.GetFilledArea(p);
                if(obj.connectance > filledarea) 
                    break;
                end
                if(obj.debug_messages); fprintf('area = %5.4f p = %5.4f\n', filledarea,p); end; 
                p = p*2;          
            end
            
            %if(obj.debug_messages); fprintf('area = %10.9f p = %5.4f lastp = %5.4f\n', filledarea,upp,lowp); end;
            
            if(p < obj.p_max && p > obj.p_min) %If the parameter p is not an extreme case
                %BISECTION METHOD
                upp = p;
                lowp = p / 2;
                mid = 0;
                while( abs( obj.connectance - filledarea) > 0.001)
                    mid = (upp + lowp)/2;
                    filledarea = obj.GetFilledArea(mid);
                    if(filledarea < obj.connectance)
                        upp = mid;
                    else
                        lowp = mid;
                    end
                    if(obj.debug_messages); fprintf('area = %10.9f p = %f\n', filledarea,mid); end;        
                end
                if(mid ~= 0)
                    p = mid;
                end;
            end
            
            obj.calculated_fill = filledarea;
        end
        
        function Area = GetFilledArea(obj,p)
        % GetFilledArea - Get the area above the isocline
        %
        %   Area = GetFilledArea(obj,p) - Ghet the area above the
        %   isocline with parameter p.
            
            
            %Isocline equation
            obj.fxp = 0.5/obj.n_rows + ((obj.n_rows-1)/obj.n_rows) * (1-(1-(obj.n_cols*(obj.X)-0.5)/(obj.n_cols-1)).^p).^(1/p);
            
            %Area below the isocline
            integral = trapz(obj.X,obj.fxp);
            
            %Three ways of calculating the area (only important when the
            %matrix is small. Case 2 gives the best results.
            switch obj.used_area
                case 1
                    Area = 1 - real(integral);
                case 2
                    Area = 1 - real(integral) - (obj.n_rows-0.5)*(0.5)/(obj.n_rows*obj.n_cols);
                otherwise
                    Area = (obj.n_rows-0.5)/(obj.n_rows-1) - real(integral) * obj.n_rows * obj.n_cols / ((obj.n_cols-1)*(obj.n_rows-1));
            end
        end
        
        function obj = CalculateDiagonalsAndDistances(obj)
        % CalculateDiagonalsAndDistances - Calculate diagonal and
        % isocline distance size matrices
        %
        %   obj = CalculateDiagonalsAndDistances(obj) - Calculate diagonal and
        %   isocline distance size matrices
            obj.u_matrix = zeros(size(obj.matrix));
            MaxDiag = sqrt(2);
            
            obj.diag_matrix = zeros(size(obj.matrix));
            obj.d_matrix = zeros(size(obj.matrix));
            
            %For each row and column
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                         
                    y1 = real(obj.pos_x_matrix(i,j) + obj.pos_y_matrix(i,j) - obj.X);
                    y2 = obj.fxp;
                    
                    [~, index] = min(abs(y1-y2));
                    
                    %Intersection point between the diagonal and the
                    %iscoline
                    ycross = y1(index);
                    xcross = obj.X(index);

                    %Distance from the isocline to the matrix element
                    distance = sqrt( (obj.pos_x_matrix(i,j)-xcross)^2 + (obj.pos_y_matrix(i,j)-ycross)^2 );
                    obj.d_matrix(i,j) = distance;
                    obj.diag_matrix(i,j) = (obj.pos_x_matrix(i,j) + obj.pos_y_matrix(i,j)) * sqrt(2);

                    if(obj.diag_matrix(i,j) > MaxDiag)
                        obj.diag_matrix(i,j) = abs(obj.pos_x_matrix(i,j) + obj.pos_y_matrix(i,j) - 2) * sqrt(2);
                    end
                    
                    % Change to negative elements below isocline, such that
                    % the sign will differentiate above vs below isocline
                    % elemnts.
                    if(obj.pos_y_matrix(i,j) < ycross)
                        obj.d_matrix(i,j) = -obj.d_matrix(i,j);
                    end
                end
            end       
        end 
        
    end
    
    % CALCULATE TEMPERATURE AND IMPORTANT VALUES
    methods(Access = 'protected')
        
        function obj = CalculateTemperature(obj)
            % CalculateTemperature - Calculate the matrix temperature
            %   obj = CalculateTemperature(obj) - Calculate the temperature
            %   using the Atmar standard equation. The temperature is in
            %   the interval [0,100], while the nestedness in the interval
            %   [0,1]. High values of temperature corresponds to low values
            %   of nestedness.
            obj.u_min = obj.CalculateUnexpectedness();
            obj.T = obj.K*obj.u_min;
            obj.N= (100-obj.T)/100;
        end
        
        function unex = CalculateUnexpectedness(obj)
            % CalculateUnexpectedness - Calculate the matrix unexpectedness
            %   obj = CalculateUnexpectedness(obj) - Sum all temperature
            %   contributions from unexpected cells (absences below the
            %   matrix and presences above the matrix)
            obj.u_matrix = zeros(size(obj.matrix));
            obj.u_matrix = ((obj.matrix==0 & obj.d_matrix > 0) | (obj.matrix ~=0 & obj.d_matrix < 0 )).*((obj.d_matrix./obj.diag_matrix).^2);   
            unex = sum(sum(obj.u_matrix)) / (obj.n_rows*obj.n_cols);
        end
    end
    
    
    methods(Static)
        
        function ntc = NTC(matrix)
        % NTC - Calculate the NTC nestedness
        %
        %   nest = NTC(matrix) Calculate the NTC nestedness,
        %   print the basic information to
        %   screen and return an NestednessNTC object that contains such
        %   information in nest.
        
            ntc = NestednessNTC(matrix);
            ntc.Detect();
            ntc.Print();
        end
        
        function matrix = PERFECT_NESTED(nrows,ncols,fill)
        % PERFECT_NESTED - Return a perfect nested matrix according to the
        % NTC algorithm
        %
        %   matrix = PERFECT_NESTED(nrows,ncols,fill) - Return a perfect
        %   nested matrix of size nrows by ncols and a connectance = fill.
        %   The perfect nested matrix follows the definition of the NTC
        %   algorithm (the isocline divide ones from zeros in the entire
        %   matrix).
            matrix = zeros(nrows,ncols);
        
            bnest = NetworkBipartite(matrix);
            
            nest = NestednessNTC(bnest);
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
        %
        %   matrix = FIND_UNEXPECTED_CELLS(matrix) - For matrix 'matrix',
        %   calculate the geometry in order to return a matrix 'matrix_unex'
        %   with ones in the position of unexpected cells of the original
        %   matrix.
            nest = NestednessNTC(matrix);
            nest.CalculateMatrixGeometry();
            nest.CalculateUnexpectedness();
            
            matrix_unex = nest.uMatrix > 0;%.0005;
            
        end
        
        function [x,y] = GET_ISOCLINE(n_rows,n_cols,p_value)
        % GET_ISOCLINE - Get the isocline function
        %
        %   [x y] = GET_ISOCLINE(n_rows,n_cols,p_value) - Get the isocline
        %   function in x and y vectors for a matrix of size n_rows by
        %   n_cols and a connectance of p_value. Useful when the user is
        %   only interested in the isocline (e.g. plotting) and not the
        %   temperature value.
            if(nargin==1)
                matrix_loc = n_rows>0;
            else
                matrix_loc = zeros(n_rows,n_cols);
                len = n_rows*n_cols;
                matrix_loc(1:round(len*p_value))=1;
            end
            
            ntc = NestednessNTC(matrix_loc);
            
            ntc.CalculateMatrixGeometry();
            x = 0.5 + ntc.n_cols.*ntc.X;
            y = 0.5 + ntc.n_rows.*ntc.fxp;
        end
        
    end
    

end


% globalMinimalT = 500;
%             matrixLocalMinima = [];
%             indexRowLocalMinima = [];
%             indexColLocalMinima = [];
%             indexRowGlobalMinima = [];
%             indexColGlobalMinima = [];
%             
%             
%             failedtoincrease = 0; %Count if the next matrix randomization do an improvement
%             % Do obj.trials initial random permutations of the
%             % matrix to be tested
%             for i = 1:obj.trials
%                 
%                 % If no increase is detected in obj.break_random continuos
%                 % trials, no need for continue looking.
%                 if(failedtoincrease > obj.break_random)
%                     %fprintf('Break on i = %i\n',i);
%                     break;
%                 end
%                 
%                 
%                 permutationMinimalT = 500; %temperature infinite
%                 obj.T = 500;
%                 
%                 [matrix_permuted,new_row_index,new_col_index] = MatrixFunctions.RANDOM_SORT(obj.matrix);
%                 obj.matrix = matrix_permuted;
%                 obj.index_rows(new_row_index);
%                 obj.index_cols(new_col_index);
% %                i = 1;
%                 while(1)
%                     %display(i);
%                     %i = i+1;
%                     [sort_matrix, new_row_index, new_col_index] = MatrixFunctions.SORT_MATRIX(obj.matrix);
%                     obj.matrix = sort_matrix;
%                     obj.index_rows = obj.index_rows(new_row_index);
%                     obj.index_cols = obj.index_cols(new_col_index);
%                     
%                     obj.CalculateTemperature();
%                     
%                     if(obj.debug_messages == 1); fprintf('TLocal = %f T = %f\n', permutationMinimalT,obj.T); end;
%                     
%                     if(abs(permutationMinimalT - obj.T) <= 0.001 || obj.T > permutationMinimalT)
%                         break;
%                     end
%                         
%                     if(obj.T < permutationMinimalT)
%                         permutationMinimalT = obj.T;
%                         matrixLocalMinima = obj.matrix;
%                         indexRowLocalMinima = obj.index_rows;
%                         indexColLocalMinima = obj.index_cols;
%                     end
%                     
%                 end
%                 if(obj.debug_messages == 1); fprintf('finalizo ciclo\n'); end;
%                 
%                 %Save if permutation is smaller than the global minimal
%                 if(permutationMinimalT < globalMinimalT)
%                     %fprintf('TMinimalGlob = %f\n', permutationMinimalT);
%                     globalMinimalT = permutationMinimalT;
%                     obj.matrix_minimal = matrixLocalMinima;
%                     indexRowGlobalMinima = indexRowLocalMinima;
%                     indexColGlobalMinima = indexColLocalMinima;
%                     failedtoincrease = 0;
%                 end
%                 
%                 failedtoincrease = failedtoincrease + 1;
%             end
