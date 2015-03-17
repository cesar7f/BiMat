% Printer - Class used for printing output files on screen and files.
%
% Printer Properties:
%     delimiter - Delimiter of columns for information presented in tables
%     bipweb - Bipartite object for which the information will be printed
%
% Printer Methods:
%    Printer - Main Constructor
%    CreateCytoscapeData - Create input files for Cytoscape software
%    PrintAdjacencyList - Print the adjacency list
%    PrintGeneralProperties - Print the general properties 
%    PrintStructureValues - Print the sctructure values
%    PrintStructureStatistics - Print the statistics of structure values
%    PrintStructureStatisticsOfModules - Print the structure module
%    PrintRowModuleDiversity - Print the statistics of module
%    PrintColumnModuleDiversity - Print the statistics of module
%    PrintGeneralPropertiesOfModules - Print general properties of modules
%    PRINT_TO_FILE - Print a string to a text file
%    CREATE_FORMATED_STRING - Created a formated string
%
% See also:
%    Reader
classdef Printer < handle
    
    properties
        
        delimiter      = ','; % Delimiter of columns for information presented in tables
        bipweb         = {};  % Bipartite object for which the information will be printed
        
    end
    
    methods
        
        function obj = Printer(webbip)
        % Printer - Main Constructor
        %   pt = Printer(webbip) Create a Printer object pt wich will be
        %   related to Bipartite object webbip.
        
            obj.bipweb = webbip;
            
        end
                
        function obj = CreateCytoscapeData(obj,filename)
        % CreateCytoscapeData - Create input files for Cytoscape software
        %
        %   obj = CreateCytoscapeData(obj,filename) Create two input files
        %   filename_edges.csv and filename_nodes.csv for Cytoscape
        %   software. The first one is indispensable for Cytoscape and
        %   contains the list of edges of the bipartite network. The second
        %   contain additional information of the nodes such as module id,
        %   name label, etc. This additional information can be used in
        %   Cytoscape for creating layouts with more information.
            
            n_rows = obj.bipweb.n_rows;
            n_cols = obj.bipweb.n_cols;
            
            %inter = sum(sum(obj.bipweb.matrix));
            %adlist = zeros(inter,3);
            
            matrix = obj.bipweb.webmatrix;
            
            %nn = 1;

            fidedges = fopen([filename,'_edges.csv'],'w');
            fidnodes = fopen([filename,'_nodes.csv'],'w');
            
            for i = 1:n_rows; fprintf(fidedges,'%i\n',i); end;
            for j = 1:n_cols; fprintf(fidedges,'%i\n',j+n_cols); end;
            
            modules_done = obj.bipweb.community.done;
            
            row_modul = mod(find(obj.bipweb.community.rr'==1),obj.bipweb.community.N);
            col_modul = mod(find(obj.bipweb.community.tt'==1),obj.bipweb.community.N);
            row_modul = row_modul.*(row_modul>0)+(row_modul==0).*obj.bipweb.community.N;
            col_modul = col_modul.*(col_modul>0)+(col_modul==0).*obj.bipweb.community.N;
            
            for i = 1:n_rows
                for j = 1:n_cols
                    if(matrix(i,j) > 0)
                        if(~modules_done)
                            fprintf(fidedges,'%i %i %i\n',i,j+n_rows,matrix(i,j));
                        else
                            if(row_modul(i)==col_modul(j))
                                fprintf(fidedges,'%i %i %i %i\n',i,j+n_rows,matrix(i,j),row_modul(i));
                            else
                                fprintf(fidedges,'%i %i %i %i\n',i,j+n_rows,matrix(i,j),0);
                            end
                        end
                    end
                end
            end
            
            fclose(fidedges);
               
            fprintf(fidnodes,'ID,NAME,TYPE,MODULE\n');
            row_labels = obj.bipweb.row_labels;
            col_labels = obj.bipweb.col_labels;
            
            for i = 1:n_rows
                fprintf(fidnodes,'%i,%s,%i,%i\n',i,row_labels{i},1,row_modul(i));
            end
            for j = 1:n_cols
                fprintf(fidnodes,'%i,%s,%i,%i\n',j+n_rows,col_labels{j},2,col_modul(j));
            end
            
            fclose(fidnodes);

        end
        
    end
    
    methods
       
        function str = PrintAdjacencyList(obj,filename)
        % PrintAdjacencyList - Print the adjacency list
        %
        %   obj = PrintAdjacencyList(obj) - Print the adjacency list of a
        %   network to screen.
        %
        %   obj = PrintAdjacencyList(obj,FILE) - Print the adjacency list of a
        %   network to screen and text file FILE.
            str = '';
            for i = 1:obj.bipweb.n_rows
                for j = 1:obj.bipweb.n_cols
                    if(obj.bipweb.webmatrix(i,j)>0)
                        str = [str, sprintf('%s,%i,%s',obj.bipweb.row_labels{i},obj.bipweb.webmatrix(i,j),obj.bipweb.col_labels{j}), '\n'];
                    end
                end
            end
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
        function PrintGeneralProperties(obj,filename)
        % PrintGeneralProperties - Print the general properties 
        %
        %   obj = PrintGeneralProperties(obj) - Print the general 
        %   properties of a bipartite network to screen. 
        %
        %   obj = PrintGeneralProperties(obj,FILE) - Print the general 
        %   properties of a bipartite network to screen and text file FILE.
        
            str = 'General Properties\n';
            str = [str, '\t Number of species:       \t', sprintf('%6i',obj.bipweb.n_rows+obj.bipweb.n_cols), '\n'];
            str = [str, '\t Number of row species:   \t', sprintf('%6i',obj.bipweb.n_rows), '\n'];
            str = [str, '\t Number of column species:\t', sprintf('%6i',obj.bipweb.n_cols), '\n'];
            str = [str, '\t Number of Interactions:  \t', sprintf('%6i',obj.bipweb.n_edges), '\n'];
            str = [str, '\t Size:                    \t', sprintf('%6i',obj.bipweb.size_webmatrix), '\n'];
            str = [str, '\t Connectance or fill:     \t', sprintf('%6.3f',obj.bipweb.connectance), '\n'];
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
        function PrintStructureValues(obj,filename)
        % PrintStructureValues - Print the sctructure values
        %
        %   obj = PrintStructureValues(obj) - Print the structure values
        %   of nestedness and modularity to screen.
        %
        %   obj = PrintStructureValues(obj,FILE) - Print the structure values
        %   of nestedness and modularity to screen and text file FILE.
        
            if(obj.bipweb.community.done == 0)
                obj.bipweb.community.Detect();
            end
            
            if(obj.bipweb.nestedness.done == 0)
                obj.bipweb.nestedness.Detect();
            end
            
            if (nargin == 2);
                obj.bipweb.community.Print(filename);
                obj.bipweb.nestedness.Print(filename);
            else
                obj.bipweb.community.Print();
                obj.bipweb.nestedness.Print();
            end
            
        end
        
        function PrintStructureStatistics(obj,filename)
        % PrintStructureStatistics - Print the statistics of structure values
        %
        %   obj = PrintStructureStatistics(obj) - Print the statistics of
        %   nestedness and modularity to screen.
        %
        %   obj = PrintStructureStatistics(obj,FILE) - Print the statistics of
        %   nestedness and modularity to screen and text file FILE.
           
            if(obj.bipweb.statistics.community_done == 0)
                obj.bipweb.statistics.TestCommunityStructure();
            end
            
            if(obj.bipweb.statistics.nested_done == 0)
                obj.bipweb.statistics.TestNestedness();
            end
            
            if (nargin == 2);
                obj.bipweb.statistics.Print(filename);
            else
                obj.bipweb.statistics.Print();
            end
        end
        
        function PrintStructureStatisticsOfModules(obj,filename)
        % PrintStructureStatisticsOfModules - Print the structure module
        % statistics
        %
        %   obj = PrintStructureStatisticsOfModules(obj) - Print the statistics of
        %   nestedness and modularity for modules detected in the network
        %   to screen.
        %
        %   obj = PrintStructureStatisticsOfModules(obj,FILE) - Print the statistics of
        %   nestedness and modularity for modules detected in the network
        %   to screen and text file FILE.
               
            if(nargin == 2)
                obj.bipweb.internal_statistics.meta_statistics.Print(filename);
            else
                obj.bipweb.internal_statistics.meta_statistics.Print();
            end
            
        end
        
        function PrintRowModuleDiversity(obj, filename)
        % PrintRowModuleDiversity - Print the statistics of module
        % diversity for row nodes
        %
        %   obj = PrintRowModuleDiversity(obj) - Print the statistics of
        %   class node diversity for rows inside modules. This statistics can be
        %   used as a test to detect if correlation exist between the class
        %   of nodes and module configuration.
        %
        %   obj = PrintRowModuleDiversity(obj,FILE) - Print the statistics of
        %   class node diversity for rows inside modules to text file FILE. This statistics can be
        %   used as a test to detect if correlation exist between the class
        %   of nodes and module configuration.
        
            headers{1} = 'Module';
            headers{2} = 'index value';
            headers{3} = 'zscore';
            headers{4} = 'percent';
            
            columns = (1:length(obj.bipweb.internal_statistics.row_diversity.value))';
            columns = [columns obj.bipweb.internal_statistics.row_diversity.value'];
            columns = [columns obj.bipweb.internal_statistics.row_diversity.zscore'];
            columns = [columns obj.bipweb.internal_statistics.row_diversity.percentile'];
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,',');
            
            str = ['Random permutations:\t', sprintf('%25i',obj.bipweb.internal_statistics.row_diversity.n_permutations), '\n', str];
            str = ['Diversity index:    \t', sprintf('%25s',func2str(obj.bipweb.internal_statistics.row_diversity.diversity_index)), '\n', str];
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
        end
        
        function PrintColumnModuleDiversity(obj, filename)
        % PrintColumnModuleDiversity - Print the statistics of module
        % diversity for row nodes
        %
        %   obj = PrintColumnModuleDiversity(obj) - Print the statistics of
        %   class node diversity for columns inside modules. This statistics can be
        %   used as a test to detect if correlation exist between the class
        %   of nodes and module configuration.
        %
        %   obj = PrintColumnModuleDiversity(obj,FILE) - Print the statistics of
        %   class node diversity for columns inside modules to text file FILE. This statistics can be
        %   used as a test to detect if correlation exist between the class
        %   of nodes and module configuration.           
            
            headers{1} = 'Module';
            headers{2} = 'index value';
            headers{3} = 'zscore';
            headers{4} = 'percent';
            
            columns = (1:length(obj.bipweb.internal_statistics.col_diversity.value))';
            columns = [columns obj.bipweb.internal_statistics.col_diversity.value'];
            columns = [columns obj.bipweb.internal_statistics.col_diversity.zscore'];
            columns = [columns obj.bipweb.internal_statistics.col_diversity.percentile'];
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,',');
            
            str = ['Random permutations:\t', sprintf('%25i',obj.bipweb.internal_statistics.col_diversity.n_permutations), '\n', str];
            str = ['Diversity index:    \t', sprintf('%25s',func2str(obj.bipweb.internal_statistics.col_diversity.diversity_index)), '\n', str];
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
        function PrintGeneralPropertiesOfModules(obj, filename)
        % PrintGeneralPropertiesOfModules - Print general properties of modules
        %
        %   obj = PrintGeneralPropertiesOfModules(obj) - Print the general properties
        %   of modules to screen.
        %
        %   obj = PrintGeneralPropertiesOfModules(obj,FILE) - Print the general
        %   properties of modules to screen and text file FILE. 
        
            fprintf('No. \t H \t P \t S \t I \t M \t C \t Lh \t Lp\n');
            
            headers{1} = 'Module';
            headers{2} = 'm';
            headers{3} = 'n';
            headers{4} = 'S';
            headers{5} = 'I';
            headers{6} = 'M';
            headers{7} = 'C';
            headers{8} = 'Lr';
            headers{9} = 'Lc';
            
            if(obj.bipweb.community.done == 0)
                obj.bipweb.community.Detect();
            end
            
            obj.bipweb.internal_statistics.module_networks = obj.bipweb.community.ExtractCommunityModules();
            
            columns = (1:length(obj.bipweb.internal_statistics.module_networks))';
            
            columns = [columns cellfun(@(x) x.n_rows, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_cols, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_cols+x.n_rows, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_edges, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.size_webmatrix, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.connectance, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_edges/x.n_rows, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_edges/x.n_cols, obj.bipweb.internal_statistics.module_networks)];
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,obj.delimiter);
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
    end
    
    methods(Static)
        
        function PRINT_TO_FILE(str,filename)
        % PRINT_TO_FILE - Print a string to a text file
        %
        %   PRINT_TO_FILE(STR, FILE) Print a string STR to text file FILE
            fid = fopen(filename,'w');
            fprintf(fid,str);            
            fclose(fid);
        end
        
        function formated_string = CREATE_FORMATED_STRING(headers,columns,delimiter)
        % CREATE_FORMATED_STRING - Created a formated string
        %
        %   STR = CREATE_FORMATED_STRING(HEADERS, COLUMNS, DELIM) Create a
        %   formated single string STR that will contain a table with
        %   headers HEADERS using DELIM as column delimiter. The
        %   information of columns will be extracted from COLUMS.

        
            assert((isa(columns,'cell') && length(headers) == length(columns)) ...
                || size(columns,2)==length(headers));
            
            n_cols = length(headers);
            
            if(nargin==2)
                delimiter = ',';
            end
            
            
            
            if(isa(columns,'double'))
                columns_str = arrayfun(@num2str, columns, 'unif', 0);
                [n_r n_c] = size(columns_str);
                spacing = zeros(n_cols,1);
                for i = 1:n_cols
                    spacing(i) = max(max(cellfun('length',columns_str(:,i))),length(headers{i}));
                end
                
                
                row_size = (sum(spacing)+n_cols+1)*(n_r+1);
                %allocate memory for char
                formated_string(row_size) = char(0);
                i_char = 1;
                for i = 1:n_cols
                    format_data = ['%',num2str(spacing(i)),'s'];
                    formated_string(i_char:i_char+spacing(i)-1) = sprintf(format_data, headers{i});
                    if(i ~= n_cols)
                        formated_string(i_char+spacing(i)) = delimiter;
                    end
                    i_char = i_char + spacing(i)+1;
                end
                formated_string(i_char-1:i_char) = '\n';
                i_char = i_char+1;
                for j = 1:n_r
                    for i = 1:n_c
                        format_data = ['%',num2str(spacing(i)),'s'];
                        formated_string(i_char:i_char+spacing(i)-1) = ...
                            sprintf(format_data, columns_str{j,i});
                        if(i ~= n_cols)
                            formated_string(i_char+spacing(i)) = delimiter;
                        end
                        i_char = i_char + spacing(i)+1;
                    end
                    formated_string(i_char-1:i_char) = '\n';
                    i_char = i_char+1;
                end
            end
            
            
        end
        
    end
    
end