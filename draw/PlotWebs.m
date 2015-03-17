% PlotWebs - Matrix and graph plotting class. 
% This class allows the plotting of matrix and graph layouts that represent
% the bipartite network.
%
% PlotWebs Properties:
%   GENERAL
%     matrix - Adjacency matrix in boolean version
%     webmatrix - Adjacency matrix with edge weight information
%     biweb - Bipartite object that an object of the class will reference to
%     n_rows - Number of row nodes
%     n_cols - Number of column nodes
%     row_labels - Text labels for row nodes
%     col_labels - Text labels for column nodes
%     use_labels - Flag that indicate the use of text labels
%     font_size - Font size for text labels. Change according to network size
%     colors - Colors used for draw modules
%     interpreter - none or latex for text labels
%     ax - Axis of the plot
%   MATRIX LAYOUT
%     cell_color - Cell color
%     use_empty_cell - Flag to indicate if we want to plot empty cells
%     cell_empty_color - Color of the empty cell 
%     back_color - Back color
%     margin - Margin between cells
%     isocline_color - Color of the isocline
%     division_color - Color used for the division of the modules
%     line_width - Line width used in the isocline
%     use_isocline - Flag that indicated the plotting of isocline in a nested graph
%     color_interactions - Color of interactions. Used only if use_type_interaction = true
%     use_module_format - Flag to give appropiate color format to modules
%     use_type_interaction - Flag to color cells according to weight (only if discrete)
%     use_isocline_module - Flag to indicate if isocline will be plotted inside module sorting
%   GRAPH LAYOUT
%     radius - Radius of the nodes for graph layouts.
%     vertical_margin - Vertical margin between nodes for graph layouts.
%     horizontal_proportion - Horizontal margin (proportional to the y total size) between nodes for graph layouts. 
%     bead_color_rows - Color of the row nodes.
%     bead_color_columns - Color of the column nodes.
%     link_color - Color of the links.
%     link_width - Edge width
%
% PlotWebs Methods:
%     PlotWebs - Main Constructor
%     PlotGraph - Plot a original sorting graph layout
%     PlotNestedGraph - Plot a nested sorting graph layout
%     PlotModularGraph - Plot a modular sorting graph layout
%     PlotMatrix - Plot a original sorting matrix layout
%     PlotNestedMatrix - Plot a nested sorting matrix layout
%     PlotModularMatrix - Plot a modular sorting matrix layout
%     PLOT_TO_PDF - Plot to PDF file the current Figure
%     PLOT_TO_EPS - Plot to EPS file the current Figure
%     PLOT_MATRIX - Plot a matrix layout of a bipartite network
%     PLOT_NESTED_MATRIX - Plot a nested matrix layout of a bipartite network
%     PLOT_MODULAR_MATRIX - Plot a modular matrix layout of a bipartite network
%     PLOT_GRAPH - Plot a graph layout of a bipartite network
%     PLOT_NESTED_GRAPH - Plot a nested graph layout of a bipartite network
%     PLOT_MODULAR_GRAPH - Plot a modular graph layout of a bipartite network
%     
%
% See also:
%    Printer.CreateCytoscapeData
classdef PlotWebs < handle;
    
    % GENERAL PROPERTIES
    properties
        
        matrix               = [];       % Adjacency matrix in boolean version
        webmatrix            = [];       % Adjacency matrix with edge weight information
        biweb                = {};       % Bipartite object that an object of the class will reference to
        n_rows               = [];       % Number of row nodes
        n_cols               = [];       % Number of column nodes
        row_labels           = {};       % Text labels for row nodes
        col_labels           = {};       % Text labels for column nodes
        use_labels           = true;     % Flag that indicate the use of text labels
        font_size            = 5;        % Font size for text labels. Change according to network size.
        colors               = [];       % Colors used for draw modules.
        interpreter          = 'none';   % none or latex for text labels.
        ax                   = 'image';  % Axis of the plot.
        
    end
    
    properties(Access = 'protected')
        index_rows           = [];       % Index position of rows after modular or nested sorting
        index_cols           = [];       % Index position of columns after modular or nested sorting.
    end
    
    % MATRIX LAYOUT PROPERTIES
    properties
       
        cell_color           = [0 0 0];  % Cell color
        use_empty_cell       = false;    % Flag to indicate if we want to plot empty cells
        cell_empty_color     = 'white';  % Color of the empty cell        
        back_color           = [1 1 1];  % Back color
        margin               = 0.12;     % Margin between cells
        isocline_color       = 'red';    % Color of the isocline
        division_color       = 'blue';   % Color used for the division of the modules
        line_width           = 1.5;      % Line width used in the isocline
        use_isocline         = true;     % Flag that indicated the plotting of isocline in a nested graph
        color_interactions   = [];       % Color of interactions. Used only if use_type_interaction = true
        use_type_interaction = false;    % Flag to color cells according to weight (only if discrete)
        use_module_format    = true;     % Flag to give appropiate color format to modules
        use_isocline_module  = false;    % Flag to indicate if isocline will be plotted inside module sorting
        use_module_division = true;     % Plot Divisions for modules

    end
    
    properties(Access = 'protected')
        pos_x                = [];       % X coordinate of each cell in a unit square
        pos_y                = [];       % Y coordinate of each cell in a unit square
        border_color         = 'none';   % Border color
        border_line_width    = 0.0005;   % Border line width
    end
    
    
    % GRAPH LAYOUT PROPERTIES
    properties
       
        radius                = 0.5;     % Radius of the nodes for graph layouts.
        vertical_margin       = 0.12;    % Vertical margin between nodes for graph layouts.
        horizontal_proportion = 0.5;     % Horizontal margin (proportional to the y total size) between nodes for graph layouts. 
        bead_color_rows       = [1 0 0]; % Color of the row nodes.
        bead_color_columns    = [0 0 1]; % Color of the column nodes.
        link_color            = [0 0 0]; % Color of the links.
        link_width            = 0.75;     % Edge width
        
    end
    
    properties(Access = 'protected')
        row_pos               = [];      % Positions of the row nodes.
        col_pos               = [];      % Positions of the column nodes.
    end
    
    % PUBLIC METHODS
    methods
       
        function obj = PlotWebs(bipmatrix_or_biweb,plot_format)
        % PlotWebs - Main Constructor
        %
        %   bp = PlotWebs(bipmatrix_or_biweb) Create a PlotWebs object
        %   called bp using a bipartite adjacency matrix or bipartite
        %   object named bipmatrix_or_biweb
        
            if(isa(bipmatrix_or_biweb,'Bipartite'))
                obj.biweb = bipmatrix_or_biweb;
                obj.matrix = obj.biweb.matrix;
                obj.webmatrix = obj.biweb.webmatrix;
            else
                obj.matrix = bipmatrix_or_biweb;
                obj.webmatrix = bipmatrix_or_biweb;
            end
            
            if(nargin == 2)
                obj.SetPlotFormat(plot_format)
            end
            
            [obj.n_rows obj.n_cols] = size(obj.matrix);
            
            if(length(unique(obj.webmatrix))>2)
                obj.use_type_interaction = 1;
            end
            
                       
            
            obj.StartUp();
            
        end  
        
        function obj = SetPlotFormat(obj, plot_format)
           
            if(~isa(plot_format,'PlotFormat'))
               throw(MException('PlotWebs:SetPlotFormat', 'You must send a PlotFormat object')); 
            end
            
            % GENERAL LAYOUT PROPERTIES
            obj.use_labels           = plot_format.use_labels;
            obj.font_size            = plot_format.font_size;
            obj.colors               = plot_format.colors;
            obj.interpreter          = plot_format.interpreter;
            obj.ax                   = plot_format.ax;
    
            % MATRIX LAYOUT PROPERTIES
            obj.cell_color           = plot_format.cell_color;
            obj.use_empty_cell       = plot_format.use_empty_cell;
            obj.cell_empty_color     = plot_format.cell_empty_color;
            obj.back_color           = plot_format.back_color;
            obj.margin               = plot_format.margin;
            obj.isocline_color       = plot_format.isocline_color;
            obj.division_color       = plot_format.division_color;
            obj.line_width           = plot_format.line_width;
            obj.use_isocline         = plot_format.use_isocline;
            obj.color_interactions   = plot_format.color_interactions;
            obj.use_type_interaction = plot_format.use_type_interaction;
            obj.use_module_format    = plot_format.use_module_format;
    
            % GRAPH LAYOUT PROPERTIES
            obj.radius                = plot_format.radius;
            obj.vertical_margin       = plot_format.vertical_margin;
            obj.horizontal_proportion = plot_format.horizontal_proportion;
            obj.bead_color_rows       = plot_format.bead_color_rows;
            obj.bead_color_columns    = plot_format.bead_color_columns;
            obj.link_color            = plot_format.link_color;
            obj.link_width            = plot_format.link_width;
            
            if(isempty(obj.colors))
                obj.colors = colormap(jet(64));%colormap('jet');
                obj.colors = obj.colors([23,2,13,42,57,20,15,11,9,16,3,28,26,24,46,59,41,18,56,40,17,48,27,53,6,62,5,60,14,32,64,19,36,58,39,21,4,8,35,30,50,63,25,51,55,34,61,37,47,44,54,43,38,12,52,33,31,1,22,29,10,45,49,7],:);
            end
            
        end
        
        function PlotGraph(obj)
        % PlotGraph - Plot a original sorting graph layout
        %
        %   obj = PlotGraph(obj) - Plot a graph layout using the original sorting of nodes.
            cla;
            maxd = max(obj.n_rows,obj.n_cols);
            x1 = 1; x2 = 1+maxd*obj.horizontal_proportion;
            
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;            
            
            local_matrix = obj.webmatrix;
            
            hold on;
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j))
                        if(obj.use_type_interaction)
                            %plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.link_color,'LineWidth',obj.link_width);
                            plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.color_interactions(obj.webmatrix(i,j),:),'LineWidth',obj.link_width);
                        else
                            plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.link_color,'LineWidth',obj.link_width);
                        end
                        
                    end
                end
            end
            hold off;
            
            for i = 1:obj.n_rows
                obj.DrawCircle(x1,obj.row_pos(i),obj.bead_color_rows);
            end
            
            for j = 1:obj.n_cols
                obj.DrawCircle(x2,obj.col_pos(j),obj.bead_color_columns);
            end
            
            
            
            obj.ApplyBasicBFormat();
            
        end
        
        function PlotNestedGraph(obj)
        % PlotNestedGraph - Plot a nested sorting graph layout
        %
        %   obj = PlotNestedGraph(obj) - Plot a graph layout where nodes
        %   are sorted such that the nestedness pattern is more visible
            
            cla;
            maxd = max(obj.n_rows,obj.n_cols);
            x1 = 1; x2 = 1+maxd*obj.horizontal_proportion;
            
            [~, obj.index_rows] = sort(sum(obj.matrix,2),'descend');
            [~, obj.index_cols] = sort(sum(obj.matrix,1),'descend'); 
            
            local_matrix = obj.matrix(obj.index_rows,obj.index_cols);
            
            hold on;
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j)==1)
                        plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.link_color,'LineWidth',obj.link_width);
                    end
                end
            end
            hold off;
            
            for i = 1:obj.n_rows
                obj.DrawCircle(x1,obj.row_pos(i),obj.bead_color_rows);
            end
            
            for j = 1:obj.n_cols
                obj.DrawCircle(x2,obj.col_pos(j),obj.bead_color_columns);
            end
            
            obj.ApplyBasicBFormat();
        end
        
        function PlotModularGraph(obj)
        % PlotModularGraph - Plot a modular sorting graph layout
        %
        %   obj = PlotModularGraph(obj) - Plot a graph layout where nodes
        %   are sorted such that the modular pattern is more visible            
            cla;
            maxd = max(obj.n_rows,obj.n_cols);
            x1 = 1; x2 = 1+maxd*obj.horizontal_proportion;
            
            if(isempty(obj.biweb))
                obj.biweb = Bipartite(obj.matrix);
            end
            if(obj.biweb.community.done == 0)
                tmp = obj.biweb.community.print_results;
                obj.biweb.community.print_results = 0;
                obj.biweb.community.Detect();
                obj.biweb.community.print_results = tmp;
            end
            
            obj.index_rows = obj.biweb.community.index_rows;
            obj.index_cols = flipud(obj.biweb.community.index_cols);
            
            row_mod = obj.biweb.community.row_modules(obj.index_rows);
            col_mod = obj.biweb.community.col_modules(obj.index_cols);
            
            local_matrix = obj.matrix(obj.index_rows,obj.index_cols);
            n_col = length(obj.colors);
                
            hold on;
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j)==1 && row_mod(i)~=col_mod(j))
                        plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.link_color,'LineWidth',obj.link_width);
                    end
                    
                    if(local_matrix(i,j)==1 && row_mod(i)==col_mod(j))
                        idx_col = mod(row_mod(i),n_col); if(idx_col==0); idx_col = n_col; end;
                        plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.colors(idx_col,:),'LineWidth',obj.link_width);
                    end
                    
                end
            end
            hold off;
            
            for i = 1:obj.n_rows
                idx_col = mod(row_mod(i),n_col); if(idx_col==0); idx_col = n_col; end;
                obj.DrawCircle(x1,obj.row_pos(i),obj.colors(idx_col,:));
            end
            
            for j = 1:obj.n_cols
                idx_col = mod(col_mod(j),n_col); if(idx_col==0); idx_col = n_col; end;
                obj.DrawCircle(x2,obj.col_pos(j),obj.colors(idx_col,:));
            end
            
            obj.ApplyBasicBFormat();
        end    
        
        
        
        function obj = PlotMatrix(obj)
        % PlotMatrix - Plot a original sorting matrix layout
        %
        %   obj = PlotMatrix(obj) - Plot a matrix layout using the original sorting of nodes.    
            cla;
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(obj.matrix(i,j) > 0)
                        if(obj.use_type_interaction)
                            obj.DrawCell(i,j,obj.color_interactions(obj.webmatrix(i,j),:));
                        else
                            obj.DrawCell(i,j,obj.cell_color);
                        end
                    else
                        if(obj.use_empty_cell == true)
                            obj.DrawCell(i,j,obj.cell_empty_color);
                        end
                    end
                end
            end
          
            obj.ApplyBasicFormat();
            
        end
        
        function obj = PlotNestedMatrix(obj)
        % PlotNestedMatrix - Plot a nested sorting matrix layout
        %
        %   obj = PlotNestedMatrix(obj) - Plot a matrix layout using the nested sorting of nodes.    
        
            cla;
            [~, obj.index_rows] = sort(sum(obj.matrix,2),'descend');
            [~, obj.index_cols] = sort(sum(obj.matrix,1),'descend');
                        
            local_matrix = obj.webmatrix(obj.index_rows,obj.index_cols);
                       
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j) > 0)
                        if(obj.use_type_interaction)
                            obj.DrawCell(i,j,obj.color_interactions(local_matrix(i,j),:));
                        else
                            obj.DrawCell(i,j,obj.cell_color);
                        end
                    else
                        if(obj.use_empty_cell == true)
                            obj.DrawCell(i,j,obj.cell_empty_color);
                        end
                    end
                end
            end
            
            if(obj.use_isocline)
                hold on;
                [x,y] = NestednessNTC.GET_ISOCLINE(obj.matrix);
                p = plot(x,y);
                set(p,'Color',obj.isocline_color,'LineWidth',obj.line_width);
                hold off;
            end
            
            obj.ApplyBasicFormat();
            
        end
        
        function obj = PlotModularMatrix(obj)
        % PlotModularMatrix - Plot a modular sorting matrix layout
        %
        %   obj = PlotModularMatrix(obj) - Plot a matrix layout using the modular sorting of nodes.    
        
            cla;
            if(isempty(obj.biweb))
                obj.biweb = Bipartite(obj.matrix);
            end
            if(obj.biweb.community.done == 0)
                tmp = obj.biweb.community.print_results;
                obj.biweb.community.print_results = 0;
                obj.biweb.community.Detect();
                obj.biweb.community.print_results = tmp;
            end
            
            obj.index_rows = obj.biweb.community.index_rows;
            obj.index_cols = obj.biweb.community.index_cols;
            
            row_mod = obj.biweb.community.row_modules(obj.index_rows);
            col_mod = obj.biweb.community.col_modules(obj.index_cols);
            
            local_matrix = obj.webmatrix(obj.index_rows,obj.index_cols);
            n_col = length(obj.colors);
            
            start_y = 1;
            start_x = obj.n_cols;
            if(obj.use_module_division == true)
                for i = 1:obj.biweb.community.N

                    dx = length(find(col_mod==i));
                    dy = length(find(row_mod==i));
                    if(dx>0 && dy>0)
                        if(obj.use_module_format)
                            idx_col = mod(i,n_col); if(idx_col==0); idx_col = n_col; end;
                            obj.MakeDivision(start_y+dy-1,start_x-dx+1,dy,dx,obj.colors(idx_col,:),obj.back_color);
                        else
                            obj.MakeDivision(start_y+dy-1,start_x-dx+1,dy,dx,obj.division_color,obj.back_color);
                        end
                    end
                    start_x = start_x-dx;
                    start_y = start_y+dy;

                end
            end
            
            inter_ones = 0;
            inter_total = 0;
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j)>0)
                        if(obj.use_type_interaction)
                            obj.DrawCell(i,j,obj.color_interactions(local_matrix(i,j),:));
                            if(row_mod(i)~=col_mod(j) && local_matrix(i,j)==1)
                                inter_ones = inter_ones + 1;
                            elseif(row_mod(i)~=col_mod(j) && local_matrix(i,j)==2)
                                inter_total = inter_total + 1;
                            end
                        else
                            if(row_mod(i)==col_mod(j))
                                idx_col = mod(row_mod(i),n_col); if(idx_col==0); idx_col = n_col; end;
                                obj.DrawCell(i,j,obj.colors(idx_col,:));
                            else
                                obj.DrawCell(i,j,obj.cell_color);
                            end
                        end
                        
                    else
                        if(obj.use_empty_cell == true)
                            obj.DrawCell(i,j,obj.cell_empty_color);
                        end
                    end
                end

            end

            if(obj.use_isocline_module)
                matrices = obj.biweb.community.ExtractCommunityMatrices();
                start_y = obj.n_rows;
                start_x = obj.n_cols;
                hold on;

                for i = 1:obj.biweb.community.N

                    dx = length(find(col_mod==i));
                    dy = length(find(row_mod==i));
                    
                    if(isempty(matrices{i}) || sum(sum(matrices{i})) == numel(matrices{i}) || sum(sum(matrices{i})) == 0)
                        start_x = start_x-dx;
                        start_y = start_y-dy;
                        continue;
                    end
                    
                    [x,y] = NestednessNTC.GET_ISOCLINE(matrices{i});
                    p = plot(start_x-dx+x,start_y-dy+y);
                    %p = plot(x,y);
                    if(obj.use_module_format)
                        idx_col = mod(i,n_col); if(idx_col==0); idx_col = n_col; end;
                        set(p,'Color',obj.colors(idx_col,:),'LineWidth',obj.line_width);
                    else
                        set(p,'Color',obj.isocline_color,'LineWidth',obj.line_width);
                    end

                    %obj.MakeDivision(start_y+dy-1,start_x-dx+1,dy,dx,'black',obj.back_color);

                    start_x = start_x-dx;
                    start_y = start_y-dy;

                end
                hold off;
            end
            
            obj.ApplyBasicFormat();
        end
        
    end
    
    % PROTECTED METHODS
    methods(Access = 'protected')
        function obj = StartUp(obj)
            obj.colors = colormap(jet(64));%colormap('jet');
            obj.colors = obj.colors([23,2,13,42,57,20,15,11,9,16,3,28,26,24,46,59,41,18,56,40,17,48,27,53,6,62,5,60,14,32,64,19,36,58,39,21,4,8,35,30,50,63,25,51,55,34,61,37,47,44,54,43,38,12,52,33,31,1,22,29,10,45,49,7],:);
            obj.color_interactions = obj.colors;
            obj.color_interactions(1,:) = [1 0 0];
            obj.color_interactions(2,:) = [0 0 1];
            
            
           
            %Delete This section thereafter--------------------------
%             obj.color_interactions(1,:) = [0 255 255]/255;
%             obj.color_interactions(2,:) = [0 160 160]/255;
%             obj.color_interactions(3,:) = [0 80 80]/255;
%             obj.colors(2,:) = [0    0.7500    1.0000];
%             obj.colors(3,:) = [1.0000    0.8750         0];
%             obj.colors(4,:) = [186,85,211]/255;
%             obj.colors(5,:) = [255,165,0]/255;
%             obj.colors(6,:) = [34 139 34]/255;
            %--------------------------------------------------------
            
            obj.FillPositions();
        end
        
        function obj = FillPositions(obj)
        % FillPositions - Fill positions
        % obj = FillPositions(obj) - Calculate the geometrical position of matrix
        % cells for rows and columns.
            obj.pos_x = zeros(obj.n_rows,obj.n_cols);
            obj.pos_y = zeros(obj.n_rows,obj.n_cols);

            obj.pos_x = repmat(1:obj.n_cols, obj.n_rows, 1);
            obj.pos_y = repmat(((obj.n_rows+1)-(1:obj.n_rows))',1,obj.n_cols);
            
            maxd = max(obj.n_rows,obj.n_cols);
            obj.row_pos = linspace(maxd,1,obj.n_rows);
            obj.col_pos = linspace(maxd,1,obj.n_cols);
            
        end
        
        function DrawLine(start_cord,end_cord,color)
        % DrawLine - Draw and edge for graph layouts
        % obj = DrawLine(start_cord,end_cord,color) - Draw an edge between
        % start_cord and end_cord geometrical coordinates (positions) using
        % color as color. Function used for graph layout functions
        
            plot([start_cord(1) end_cord(1)],[start_cord(2) end_cord(2)],color,'LineWidth',1.0);
            
        end
        
        function obj = DrawCircle(obj,x,y,color)
        % DrawCircle - Draw a node for graph layouts
        % obj = DrawCircle(obj,x,y,color) - Draw a node in coordinate (x,y)
        % using color as color.
            r = obj.radius;
            marg = obj.vertical_margin;
            rec = rectangle('Position',[x-r+marg,y-r+marg,2*(r-marg),2*(r-marg)],'Curvature',[1 1]);
            set(rec,'FaceColor',color);
            set(rec,'EdgeColor',color);
        end
        
        function obj = ApplyBasicFormat(obj)
        % ApplyBasicFormat - Apply final format to the matrix layout
        % obj = ApplyBasicFormat(obj) - Print the text labels and give a
        % final format to the matrix layout plot.
        
        %obj.FillLabels();
            
 %           idx_rows = obj.index_rows;
%            idx_type = obj.biweb.row_class;
 %           idx_type = idx_type(idx_rows);
            
            if(obj.use_labels)
                if(~isempty(obj.biweb))
                    obj.col_labels = obj.biweb.col_labels;
                    obj.row_labels = obj.biweb.row_labels;
                end
                
                obj.FillLabels();
                
                for i = 1:obj.n_rows
                    %if(idx_type(i) == 2)
                        text(0,obj.pos_y(i,1),obj.row_labels{i},'HorizontalAlignment','right','FontSize',obj.font_size,'interpreter',obj.interpreter);
                    %else
                    %    text(0,obj.pos_y(i,1),obj.row_labels{i},'HorizontalAlignment','right','FontSize',obj.font_size,'interpreter',obj.interpreter,'Color',[0 128 0]/255,'FontAngle','italic');
                    %end
                end
                for j = 1:obj.n_cols
                    text(obj.pos_x(1,j),0,obj.col_labels{j},'HorizontalAlignment','right','Rotation',90,'FontSize',obj.font_size,'interpreter',obj.interpreter);
                end
            else
                set(gca,'xticklabel',[]);
                set(gca,'yticklabel',[]);
                set(gca,'YTick',[]);
                set(gca,'XTick',[]);
            end
            
            axis(obj.ax);
            set(gca,'Color',obj.back_color);
            set(gcf,'Color','white');
            set(gcf, 'InvertHardCopy', 'off'); %Do not plot in white the background
            xlim([0.5-1.1*obj.margin,0.5+obj.n_cols+obj.margin]);
            ylim([0.5-obj.margin,0.5+obj.n_rows+obj.margin]);
            box on;
            
        end
        
        function obj = ApplyBasicBFormat(obj)
        % ApplyBasicBFormat - Apply final format to the graph layout
        % obj = ApplyBasicBFormat(obj) - Print the text labels and give a
        % final format to the graph layout plot.            
            
            axis(obj.ax);
            maxd = max(obj.n_rows,obj.n_cols);
            x1 = 1; x2 = 1+maxd*obj.horizontal_proportion;
            xlim([x1-obj.radius x2+obj.radius]);
            ylim([1-obj.radius maxd+obj.radius]);
            
            
            
            if(obj.use_labels)
               
                obj.FillLabels();
                
                for i = 1:obj.n_rows
                    text(x1-obj.radius,obj.row_pos(i),obj.row_labels{i},'HorizontalAlignment','right','FontSize',obj.font_size,'interpreter',obj.interpreter);
                end
                for j = 1:obj.n_cols
                    text(x2+obj.radius,obj.col_pos(j),obj.col_labels{j},'HorizontalAlignment','left','FontSize',obj.font_size,'interpreter',obj.interpreter);
                end
                
            end
            
            set(gca,'Color','white');
            set(gcf,'Color','white');
            set(gca,'xcolor','white');
            set(gca,'ycolor','white');
            set(gcf, 'InvertHardCopy', 'off'); %Do not plot in white the background
        end
        
        function obj = FillLabels(obj)
        % FillLabels - Create text labels
        % obj = FillLabels(obj) - Create text labels according to the
        % labels of the property biweb. If the property was not defined
        % (e.g. the PlotWebs object obj was created from a matrix and not
        % from a Bipartite object), the labels are specified using default
        % names.
        
            c = findall(gca,'Type','text');
            delete(c);
            
            if(isempty(obj.biweb))
                obj.row_labels = cell(obj.n_rows,1);
                obj.col_labels = cell(obj.n_cols,1);
                for i = 1:obj.n_rows; obj.row_labels{i} = sprintf('row%03i',i); end;
                for j = 1:obj.n_cols; obj.col_labels{j} = sprintf('col%03i',j); end;
            else
                obj.row_labels = obj.biweb.row_labels;
                obj.col_labels = obj.biweb.col_labels;
            end
            
            obj.row_labels = obj.row_labels(obj.index_rows);
            obj.col_labels = obj.col_labels(obj.index_cols);
            
            set(gca,'xticklabel',[]);
            set(gca,'yticklabel',[]);
            set(gca,'YTick',[]);
            set(gca,'XTick',[]);
            
        end
        
        function obj = DrawBack(obj)
        % DrawBack - Chose the color for the background
        % obj = DrawBack(obj) - Changes the color for the background of the
        % plot.
        
            set(gca,'Color',obj.back_color);
            
        end
        
        function obj = DrawCell(obj, i, j, color)
        % DrawCell - Draw a cell for matrix layouts
        % obj = DrawCell(obj, i, j, color) - Draw a cell in position i,j
        % (in integer units) using color as color. This function is used
        % only for matrix layouts.
            
            %x1 = obj.pos_x(i,j)-0.5+obj.margin;
            %x2 = x1 + 1-2*obj.margin;
            %y1 = obj.pos_y(i,j)-0.5+obj.margin;
            %y2 = y1 + 1-2*obj.margin;
            rec1 = rectangle('Position',[obj.pos_x(i,j)-0.5+obj.margin,obj.pos_y(i,j)-0.5+obj.margin,1-2*obj.margin,1-2*obj.margin]);
            %rec1 = patch([x1,x2,x2,x1],[y1,y1,y2,y2],'r');
            set(rec1,'FaceColor',color);
            set(rec1,'EdgeColor','none');
%            set(rec1,'facealpha',1.0);
        end
        
        function obj = MakeDivision(obj,i,j,nrow,ncol,bordercolor,backcolor)
        % MakeDivision - Make a module division.
        % obj = MakeDivision(obj,i,j,nrow,ncol,bordercolor,backcolor) -
        % Make a module division starting in cell positioned at i,j of size
        % nrow by ncol cells. bordercolor and backcolor are used for border
        % and back color, respectivally.
        
        %TODO: alpha_value. This will be use for doing transparencies.
        %However it doest not work currently
            
             rec1 = rectangle('Position',[obj.pos_x(i,j)-0.5-obj.margin,obj.pos_y(i,j)-0.5-obj.margin, ncol+2*obj.margin, nrow+2*obj.margin],'LineWidth',1);
             %x1 = obj.pos_x(i,j)-0.5-obj.margin;
             %x2 = x1+ncol+2*obj.margin;
             %y1 = obj.pos_y(i,j)-0.5-obj.margin;
             %y2 = y1+nrow+2*obj.margin;
             %rec1 = patch([x1,x2,x2,x1],[y1,y1,y2,y2],'r','linewidth',1);
             set(rec1,'FaceColor',bordercolor);
             set(rec1,'EdgeColor','none');
             
             %set(rec,'LineStyle','--');
             %x1 = obj.pos_x(i,j)-0.5+obj.margin;
             %x2 = x1 + ncol-2*obj.margin;
             %y1 = obj.pos_y(i,j)-0.5+obj.margin;
             %y2 = y1 + nrow-2*obj.margin;
             rec2 = rectangle('Position',[obj.pos_x(i,j)-0.5+obj.margin,obj.pos_y(i,j)-0.5+obj.margin, ncol-2*obj.margin, nrow-2*obj.margin]);
             %rec2 = patch([x1,x2,x2,x1],[y1,y1,y2,y2],'r');
             set(rec2,'FaceColor',backcolor);
             set(rec2,'EdgeColor','none');
             
%             if(nargin == 8)
%                set(rec1,'facealpha',alpha_value);
%                set(rec2,'facealpha',alpha_value);
%             end
     
             %rec = rectangle('Position',[obj.PosX(i,j)-0.5,obj.PosY(i,j)-0.5, ncol, nrow]);
             %set(rec,'FaceColor',obj.CellColor);
             %set(rec,'EdgeColor','black');

        end
    end
    
    % STATIC METHODS
    methods(Static)

        function PLOT_TO_PDF(filenampdf)
        % PLOT_TO_PDF - Plot to PDF file the current Figure
        %
        %   PLOT_TO_PDF(filenampdf) Print to PDF file filenampdf what is
        %   plotted in the current figure
            set(gcf,'PaperPositionMode','auto');
            print(filenampdf,'-dpdf');           
        end
        
        function PLOT_TO_EPS(filenameps)
        % PLOT_TO_EPS - Plot to EPS file the current Figure
        %
        %   PLOT_TO_EPS(filenampdf) Print to EPS file filenameps what is
        %   plotted in the current figure
            set(gcf,'PaperPositionMode','auto');
            print(filenameps,'-depsc2');           
        end
        
        function PLOT_TO_JPG(filejpg)
        % PLOT_TO_JPG - Plot to JPG file the current Figure
        %
        %   PLOT_TO_JPG(filenampdf) Print to JPG file filenameps what is
        %   plotted in the current figure            
            set(gcf,'PaperPositionMode','auto');
            print(filejpg,'-djpeg2'); 
        end
        
        function plotter = PLOT_MATRIX(matrix, plot_format)
        % PLOT_MATRIX - Plot a matrix layout of a bipartite network
        %
        %   PW = PLOT_MATRIX(MATRIX) Plot a bipartite network with bipartite
        %   adjacency matrix MATRIX in matrix layout using the original
        %   sorting and returns a PlotWebs object PW
            plotter = PlotWebs(matrix);
            if(nargin==2)
                plotter.SetPlotFormat(plot_format);
            end
            plotter.PlotMatrix();
        end
        
        function plotter = PLOT_NESTED_MATRIX(matrix, plot_format)
        % PLOT_NESTED_MATRIX - Plot a nested matrix layout of a bipartite network
        %
        %   PW = PLOT_NESTED_MATRIX(MATRIX) Plot a bipartite network with bipartite
        %   adjacency matrix MATRIX in matrix layout using the nested
        %   sorting and returns a PlotWebs object PW
            plotter = PlotWebs(matrix);
            if(nargin==2)
                plotter.SetPlotFormat(plot_format);
            end
            plotter.PlotNestedMatrix();
        end
        
        function plotter = PLOT_MODULAR_MATRIX(matrix, plot_format)
        % PLOT_MODULAR_MATRIX - Plot a modular matrix layout of a bipartite network
        %
        %   PW = PLOT_MODULAR_MATRIX(MATRIX) Plot a bipartite network with bipartite
        %   adjacency matrix MATRIX in matrix layout using the modular
        %   sorting and returns a PlotWebs object PW            
            plotter = PlotWebs(matrix);
            if(nargin==2)
                plotter.SetPlotFormat(plot_format);
            end
            plotter.PlotModularMatrix();
        end
        
        function plotter = PLOT_GRAPH(matrix, plot_format)
        % PLOT_GRAPH - Plot a graph layout of a bipartite network
        %
        %   PW = PLOT_GRAPH(MATRIX) Plot a bipartite network with bipartite
        %   adjacency matrix MATRIX in graph layout using the original
        %   sorting and returns a PlotWebs object PW             
            plotter = PlotWebs(matrix);
            if(nargin==2)
                plotter.SetPlotFormat(plot_format);
            end
            plotter.PlotGraph();
        end
        
        function plotter = PLOT_NESTED_GRAPH(matrix, plot_format)
        % PLOT_NESTED_GRAPH - Plot a nested graph layout of a bipartite network
        %
        %   PW = PLOT_NESTED_GRAPH(MATRIX) Plot a bipartite network with bipartite
        %   adjacency matrix MATRIX in graph layout using the nested
        %   sorting and returns a PlotWebs object PW              
            plotter = PlotWebs(matrix);
            if(nargin==2)
                plotter.SetPlotFormat(plot_format);
            end
            plotter.PlotNestedGraph();
        end
        
        function plotter = PLOT_MODULAR_GRAPH(matrix, plot_format)
        % PLOT_MODULAR_GRAPH - Plot a modular graph layout of a bipartite network
        %
        %   PW = PLOT_MODULAR_GRAPH(MATRIX) Plot a bipartite network with bipartite
        %   adjacency matrix MATRIX in graph layout using the modular
        %   sorting and returns a PlotWebs object PW              
            plotter = PlotWebs(matrix);
            if(nargin==2)
                plotter.SetPlotFormat(plot_format);
            end
            plotter.PlotModularGraph();
        end      
        
    end
    
end