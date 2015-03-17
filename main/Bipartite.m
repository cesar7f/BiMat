classdef Bipartite < handle
    % Bipartite Main code class. 
    % This class is the main part of the code and
    % is used as bridge between all the algorithms and funtionalities of
    % this software.
    %
    % Bipartite Properties:
    %    webmatrix - Interaction matrix (not neccesearly a binary matrix).
    %    matrix - Adjacency matrix (binary matrix)
    %    community - Instance of the modularity algorithm
    %    nestedness - Instance of a nestedness algorithm
    %    n_edges - Number of edges
    %    connectance - Fill of the webmatrix
    %    n_rows - Number of rows
    %    n_cols - Number of columns
    %    size_webmatrix - Size of the matrix (n_rows times n_cols);
    %    row_degrees - Degrees of each row node.
    %    col_degrees - Degrees of each col node.
    %    name - Name of the network
    %    statistics - Object of the statistical analysis class (StatisticalTest).
    %    row_labels - Labels for row nodes.
    %    col_labels - Labels for col nodes.
    %    printer - Object of the class Print.m for output.
    %    plotter - Instance of the class PlotWebs.m for plotting matrices and graphs.
    %    row_class - Id's of the row classes
    %    col_class - Id's of the col classes.
    %    internal_statistics - Instance of the class InternalStatistics.m for multi-scale analysis.
    %    print_results - Flag to indicate if result output will be generated
    %
    % Bipartite Methods:
    %    Bipartite - Main constructor
    properties
        webmatrix             = []; % Interaction matrix (not neccesearly a binary matrix).
        matrix                = []; % Adjacency matrix (binary matrix)
        % Instance of the modularity algorithm. 
        % See also AdaptiveBrim, LeadingEigenvector, and LPBrim
        community               = {}; 
        % Instance of the nestedness algorithm. 
        % See also NestednessNTC, NestednessNODF
        nestedness            = {}; 
        n_edges               = 0;  % Number of edges
        connectance           = 0;  % Fill of the webmatrix
        n_rows                = 0;  % Number of rows
        n_cols                = 0;  % Number of columns
        size_webmatrix        = 0;  % Size of the matrix (n_rows times n_cols);
        row_degrees           = []; % Degrees of each row node.
        col_degrees           = []; % Degrees of each col node.
        name                  = {}; % Name of the network
        statistics            = {}; % Object of the statistical analysis class (StatisticalTest).
        row_labels            = {}; % Labels for row nodes.
        col_labels            = {}; % Labels for col nodes.
        printer               = {}; % Object of the class Print.m for output.
        plotter               = {}; % Instance of the class PlotWebs.m for plotting matrices and graphs.
        row_class             = []; % Id's of the row nodes
        col_class             = []; % Id's of the col nodes.
        internal_statistics   = {}; % Instance of the class InternalStatistics.m for multi-scale analysis.
        print_results         = Options.PRINT_RESULTS % Flag to indicate if result output will be generated
    end
    
    methods
        
        function obj = Bipartite(web,namebip)
            % Bipartite - Main Constructor
            %   bp = BIPARTITE(web) Create a bipartite object using a
            %   quatitative matrix web
            %   bp = BIPARTITE(web,namebip) Create a bipartite instance using a
            %   quatitative matrix web and name the object with namebip.
            
            if(nargin == 0)
               error('You need to specify a double matrix or a txt file in matrix format');
            end
            
            if(nargin == 1 && (isa(web,'double')||isa(web,'logical'))); namebip = 'No name'; end;
            
            if(nargin == 1 && isa(web,'char'))
                [~, namefile, ~] = fileparts(web);
                namebip = namefile;
                web = dlmread(web,' ');
            end
            web = 1.0*web;
            
            obj.name = namebip;
            obj.webmatrix = web;
            
            %Delete empty rows/columns according to options
            if(Options.INCLUDE_EMPTY_NODES == 0)
                obj.webmatrix = MatrixFunctions.NON_ZERO_MATRIX(web);
            end
            
            %General Properties
            [obj.n_rows, obj.n_cols] = size(web);
            %if(obj.n_rows == 0 || obj.n_cols == 0)
            %    return;
            %end;
            obj.size_webmatrix = obj.n_rows * obj.n_cols;
            obj.matrix = 1.0*(obj.webmatrix > 0);
            obj.n_edges = sum(obj.matrix(:));
            obj.connectance = sum(obj.matrix(:))/numel(obj.matrix);
            
            %Degree
            obj.row_degrees = sum(obj.matrix,2);
            obj.col_degrees = sum(obj.matrix,1)';
            
            %Nestedness
            obj.nestedness = Options.NESTEDNESS_ALGORITHM(obj.webmatrix);
            
            %Modularity
            obj.community = Options.MODULARITY_ALGORITHM(obj.webmatrix);
            
            %Statistical tests
            obj.statistics = StatisticalTest(obj);
            
            %Output
            obj.printer = Printer(obj);
            
            %Matrix and graph layouts
            %obj.plotter = PlotWebs(obj);
            
            %Multi-scale analysis
            obj.internal_statistics = InternalStatistics(obj);
            
            %Default labels
            for i = 1:obj.n_rows; obj.row_labels{i} = sprintf('row%03i',i); end;
            for i = 1:obj.n_cols; obj.col_labels{i} = sprintf('col%03i',i); end;
            
            %Default classification
            for i = 1:obj.n_rows; obj.row_class(i) = 1; end;
            for i = 1:obj.n_cols; obj.col_class(i) = 1; end;
            
            
        end
        
        
        
    end
       
    % SET/GET methods
    methods
        function value = get.plotter(obj)
            %The plotter is created the first time than is called
            if(isempty(obj.plotter))
                obj.plotter = PlotWebs(obj);
            end
            value = obj.plotter;
        end
        
%         function obj = set.print_results(obj,value)
%             The plotter is created the first time than is called
%             if(~isempty(obj.community))
%                 obj.community.print_results = value;
%             end
%             value = obj.plotter;
%         end
    end
    
end

