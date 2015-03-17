% Reader - Static class for creating bipartite networks from different
% type of input files.
%
% Reader Methods:
%     READ_ADJACENCY_LIST - Create a Bipartite object from an adjacency list text file
%     READ_BIPARTITE_MATRIX - Create a Bipartite object from a bipartite adjacency matrix
%
% See also:
%    Printer
classdef Reader
   
    methods (Access = private)
    %private so that you can't instatiate.
        function out = Reader

        end
    end 
    
    methods(Static)
       
        function bp = READ_ADJACENCY_LIST(filename, delimiter)
        % READ_ADJACENCY_LIST - Create a Bipartite object from an adjacency
        % list text file
        %   BP = READ_ADJACENCY_LIST(FILE) Reads an adjacency list from
        %   text file FILE and use it for creating a bipartite object BP.
        %   The FILE must contain two or three columns separated by an
        %   space, such that the first colum will indicate node names that
        %   are linked to nodes in the last column. Values in an optional
        %   middle column will indicate the strenght of the link (only used
        %   for plotting).
        %   BP = READ_ADJACENCY_LIST(FILE, DELIM) Reads an adjacency list from
        %   text file FILE and use it for creating a bipartite object BP.
        %   The FILE must contain two or three columns separated by a
        %   single character defined in DELIM,
        %   such that the first colum will indicate node names that
        %   are linked to nodes in the last column. Values in an optional
        %   middle column will indicate the strenght of the link (only used
        %   for plotting).
        
            %Check how many columns
            fid = fopen(filename,'r');
            if(nargin == 1)
                delimiter = ' '; %or whatever
            end
            tLines = fgets(fid);
            numCols = numel(strfind(tLines,delimiter)) + 1;
            fclose(fid);
            
            fid = fopen(filename,'r');
            if(numCols == 2)
                adja_list = textscan(fid,'%s %s');
                row_names = adja_list{1};
                col_names = adja_list{2};
                use_weiths = 0;
            else
                adja_list = textscan(fid,'%s %d %s');
                row_names = adja_list{1};
                col_names = adja_list{3};
                weights = adja_list{2};
                use_weiths = 1;
            end
            fclose(fid);
            
            %A node can not be row and column node at the same time
            assert(isempty(intersect(row_names,col_names)));
            
            row_nodes = unique(row_names);
            col_nodes = unique(col_names);
            n_rows = length(row_nodes);
            n_cols = length(col_nodes);
            
            n_interactions = length(row_names);
            matrix = zeros(n_rows,n_cols);
            for i = 1:n_interactions
                if(use_weiths); inter = weights(i); else inter = 1; end;
                row_n = row_names{i};
                col_n = col_names{i};
                i_row = not(cellfun('isempty', strfind(row_nodes,row_n)));
                i_col = not(cellfun('isempty', strfind(col_nodes,col_n)));
                matrix(i_row,i_col) = inter;
            end
            
            bp = Bipartite(matrix);
            bp.row_labels = row_nodes;
            bp.col_labels = col_nodes;
            
        end
        
        function bp = READ_BIPARTITE_MATRIX(filename)
        % READ_BIPARTITE_MATRIX - Create a Bipartite object from a
        % bipartite adjacency matrix
        %   BP = READ_BIPARTITE_MATRIX(FILE) Reads a bipartite adjacency matrix from
        %   text file FILE and use it for creating a bipartite object BP.
        %   The values in the matrix must be separated either by space
        %   or comma character
        
            matrix = dlmread(filename);
            
            bp = Bipartite(matrix);
            
        end
        
    end
    
end

%  function bipartite_web = CREATE_FROM_MATRIX_WITH_LABELS(filename)
%             % CREATE_FROM_MATRIX_WITH_LABELS - Create a bipartite object
%             % using an input file
%             %   bipartite_web = CREATE_FROM_MATRIX_WITH_LABELS(filename)
%             %   Create a Bipartite object bipartite_web using filename
%             %   file. See rodents.web for an example of the format that
%             %   filename file must have.
%             
%             fid = fopen(filename); 
%             if(fid==-1)
%                 error('The file could not be open. Check the name of the file.');
%             end
%                 
%             row_labels = {}; col_labels = {};
%             lread = fgetl(fid);
%             idx = find(lread == '"');
%             
%             for i = 1:2:length(idx)-1
%                 col_labels{(i+1)/2} = lread(idx(i)+1:idx(i+1)-1);
%             end
%             
%             lread = fgetl(fid);
%             i = 1;
%             matrix = [];
%             while(lread~=-1)
%                 idx = find(lread == '"');
%                 row_labels{i} = lread(idx(1)+1:idx(2)-1);
%                 matrix(i,:) = str2num(lread(idx(2)+1:end));
%                 i = i + 1;
%                 lread = fgetl(fid);
%             end
%             
%             bipartite_web = Bipartite(matrix);
%             bipartite_web.row_labels = row_labels;
%             bipartite_web.col_labels = col_labels;
%             [paths namefile ext] = fileparts(filename);
%             bipartite_web.name = namefile;
%         end