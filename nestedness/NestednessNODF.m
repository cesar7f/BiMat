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
% 
% function [nodf nodf_rows nodf_cols] = NODF_WITHIN_MODULES(modules_or_bipweb_or_matrix)
%         %NOT A TESTED FUNCTION   
%             if(isa(modules_or_bipweb_or_matrix,'Bipartite'))
%                 modules = modules_or_bipweb_or_matrix.modules;
%             elseif(isa(modules_or_bipweb_or_matrix,'BipartiteModularity'))
%                 modules = modules_or_bipweb_or_matrix;
%             elseif(isa(modules_or_bipweb_or_matrix,'double'))
%                 modules = Options.MODULARITY_ALGORITHM(modules_or_bipweb_or_matrix);
%             end
%                
%             if(modules.done == 0); modules.Detect(); end;
%             
%             row_modules = modules.row_modules;
%             col_modules = modules.col_modules;
%             
%             
%             matrix = modules.matrix;
%                       
%             nodf_rows = 0; nodf_cols = 0;
%             deg_rows = sum(matrix,2);
%             deg_cols = sum(matrix,1);
%             for i = 1:modules.N
% 
%                 irow = find(row_modules==i);
%                 icol = find(col_modules==i);
%                 
%                 matrix_mod = matrix(irow,icol);
%                 
%                 for ii = 1:length(irow)
%                     for jj = ii+1:length(irow)
%                         if(deg_rows(ii)==0 || deg_rows(jj)==0 || deg_rows(jj) == deg_rows(ii)); continue; end;
%                         nodf_rows = nodf_rows + sum(matrix_mod(ii,:).*matrix_mod(jj,:))/min(deg_rows(ii),deg_rows(jj));
%                     end
%                 end
%                 
%                 for ii = 1:length(icol)
%                     for jj = ii+1:length(icol)
%                         if(deg_cols(ii)==0 || deg_cols(jj)==0 || deg_cols(jj) == deg_cols(ii)); continue; end;
%                         nodf_cols = nodf_cols + sum(matrix_mod(:,ii).*matrix_mod(:,jj))/min(deg_cols(ii),deg_cols(jj));
%                     end
%                 end
%                 
%             end
%             
%             [n_rows n_cols] = size(matrix);
%             denom = n_rows*(n_rows-1)/2 + n_cols*(n_cols-1)/2;
%             denom_rows = (n_rows*(n_rows-1)/2);
%             denom_cols = (n_cols*(n_cols-1)/2);
%             
%             nodf = (nodf_rows+nodf_cols)/denom;
%             nodf_rows = nodf_rows/denom_rows;
%             nodf_cols = nodf_cols/denom_cols;
%             
%         end
%         
%         function [nodf nodf_rows nodf_cols] = NODF_BETWEEN_MODULES(modules_or_bipweb_or_matrix)
%         %NOT A TESTED FUNCTION   
%             if(isa(modules_or_bipweb_or_matrix,'Bipartite'))
%                 modules = modules_or_bipweb_or_matrix.modules;
%             elseif(isa(modules_or_bipweb_or_matrix,'BipartiteModularity'))
%                 modules = modules_or_bipweb_or_matrix;
%             elseif(isa(modules_or_bipweb_or_matrix,'double'))
%                 modules = Options.MODULARITY_ALGORITHM(modules_or_bipweb_or_matrix);
%             end
%                
%             if(modules.done == 0); modules.Detect(); end;
%             
%             row_modules = modules.row_modules;
%             col_modules = modules.col_modules;
%             
%             
%             matrix = modules.matrix;
%                       
%             nodf = 0; nodf_rows = 0; nodf_cols = 0;
%                             
%             deg_rows = sum(matrix,2);
%             deg_cols = sum(matrix,1);
%             for i = 1:modules.N
% 
%                 irow = find(row_modules==i);
%                 icol = find(col_modules~=i);
%                 
%                 matrix_mod = matrix(irow,icol);
%                 
%                 for ii = 1:length(irow)
%                     for jj = ii+1:length(irow)
%                         if(deg_rows(ii)==0 || deg_rows(jj)==0); continue; end;
%                         nodf_rows = nodf_rows + sum(matrix_mod(ii,:).*matrix_mod(jj,:))/min(deg_rows(ii),deg_rows(jj));
%                     end
%                 end
%                 
%                 irow = find(row_modules~=i);
%                 icol = find(col_modules==i);
%                 
%                 matrix_mod = matrix(irow,icol);
%                 
%                 for ii = 1:length(icol)
%                     for jj = ii+1:length(icol)
%                         if(deg_cols(ii)==0 || deg_cols(jj)==0); continue; end;
%                         nodf_cols = nodf_cols + sum(matrix_mod(:,ii).*matrix_mod(:,jj))/min(deg_cols(ii),deg_cols(jj));
%                     end
%                 end
%                 
%             end
%             
%             [n_rows n_cols] = size(matrix);
%             denom = n_rows*(n_rows-1)/2 + n_cols*(n_cols-1)/2;
%             denom_rows = (n_rows*(n_rows-1)/2);
%             denom_cols = (n_cols*(n_cols-1)/2);
%             
%             nodf = (nodf_rows+nodf_cols)/denom;
%             nodf_rows = nodf_rows/denom_rows;
%             nodf_cols = nodf_cols/denom_cols;
%             
%         end
%    methods
       
                
%         function obj = CalculateNpairedExternalCommunities(obj,bipweb)
%             
%             obj.Matrix = bipweb.matrix;
%             
%             ssRows = bipweb.Modularity.SSRow;
%             ssCols = bipweb.Modularity.SSCol;
%             
%             sumrows = sum(obj.Matrix,2);
%             sumcols = sum(obj.Matrix,1);
%             
%             obj.Np = 0;
%             %nr = 0;
%             for j = 1:obj.nRows
%                 for k = j+1:obj.nRows
%                     if(ssRows(j) ~= ssRows(k))
%                         if(sumrows(j) > sumrows(k) && sumrows(k) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(j, obj.Matrix(k,:)==1))/sumrows(k);
%                         elseif(sumrows(j) < sumrows(k) && sumrows(j) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(k, obj.Matrix(j,:)==1))/sumrows(j);
%                         end
%                     
%                     end
%                    
%                 end
%             end
%             
%             for j = 1:obj.nCols
%                 for k = j+1:obj.nCols
%                     if(ssCols(j) ~= ssCols(k))
%                         if(sumcols(j) > sumcols(k) && sumcols(k) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,k)==1,j))/sumcols(k);
%                         elseif(sumcols(j) < sumcols(k) && sumcols(j) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,j)==1,k))/sumcols(j);
%                         end
%                     end
%                 end
%             end
%             
%         end
%         
%         function obj = CalculateNpairedInternalCommunities(obj,bipweb)
% 
%             bn = obj.NetworkBipartite;
%             obj.Matrix = bn.Modularity.SortedMatrix ~= 0;
%             
%             colComSize = bn.Modularity.ColComSize{bn.Modularity.FinalLevel};
%             rowComSize = bn.Modularity.RowComSize{bn.Modularity.FinalLevel};
%             ssRow = bn.Modularity.SSRow;
%             ssCol = bn.Modularity.SSCol;
%             ncom = bn.Modularity.CommunityQuantity(bn.Modularity.FinalLevel);
%             
%             sumrows = sum(obj.Matrix,2);
%             sumcols = sum(obj.Matrix,1);
%             
%             obj.Np = 0;
%             for i = 1:ncom
%                
%                 rowcomsize = rowComSize(i);
%                 colcomsize = colComSize(i);
% 
%                 irow = find(ssRow==i,1);
%                 icol = find(ssCol==i,1);
%                 
%                 for j = irow:(irow+rowcomsize-1)
%                     for k = j+1:(irow+rowcomsize-1)
%                         if(sumrows(j) > sumrows(k) && sumrows(k) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(j, obj.Matrix(k,:)==1))/sumrows(k);
%                         elseif(sumrows(j) < sumrows(k) && sumrows(j) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(k, obj.Matrix(j,:)==1))/sumrows(j);
%                         end
%                     end
%                 end
%                 
%                 for j = icol:(icol+colcomsize-1)
%                     for k = j+1:(icol+colcomsize-1)
%                         if(sumcols(j) > sumcols(k) && sumcols(k) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,k)==1,j))/sumcols(k);
%                         elseif(sumcols(j) < sumcols(k) && sumcols(j) > 0)
%                             obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,j)==1,k))/sumcols(j);
%                         end
%                     end
%                 end
%                 
%             end
%             
%         end
        
%    end
