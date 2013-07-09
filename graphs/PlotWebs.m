% PlotWebs Matrix and graph plotting class. 
% This class allows the plotting of matrix and graph layouts that represent
% the bipartite network.
%
% PlotWebs Properties:
%     pos_x - X coordinate of each cell in a unit square
%     pos_y - Y coordinate of each cell in a unit square
%     cell_color - Cell color
%     back_color - Back color
%     margin - Margin between cels
%     ax - Axis of the plot.
%     border_color -
%     border_line_width -
%     isocline_color - Color of the isocline
%     division_color - Color used for the division of the modules
%     line_width - Line width used in the isocline
%     matrix - Adjacency matrix in boolean version
%     webmatrix - Adjacency matrix with edge weight information
%     use_type_interaction - Flag to color cells according to weight (only if discrete)
%     use_row_sections - Flag to color sections
%     use_module_format - Flag to give appropiate color format to modules
%     biweb - Bipartite object that an object of the class will reference to
%     row_labels - Text labels for row nodes
%     col_labels - Text labels fof column nodes
%     use_labels - Flag that indicate the use of text labels
%     FontSize - Font size for text labels. Change according to number of rows and columns.
%     index_rows - Index position of rows after modular or nested sorting
%     index_cols - Index position of columns after modular or nested sorting.
%     colors - Colors used for draw modules.
%     color_interactions - Color of interactions. Used only if use_type_interaction = 1
%     color_section - Color of the section. Used only if use_row_sections = 1
%     n_rows - Number of rows in the matrix
%     n_cols - Number of columns in the matrix
%     interpreter - none or latex for text labels.
%     plot_iso_modules - Flag that indicate the plotting of the isocline inside modules
%     use_isocline - Flag that indicated the plotting of isocline in a nested graph    
%     radius - Radius of the nodes for graph layouts.
%     vertical_margin - Vertical margin between nodes for graph layouts.
%     horizontal_proportion - Horizontal margin (proportional to the y total size) between nodes for graph layouts. 
%     row_pos - Positions of the row nodes.
%     col_pos - Positions of the column nodes.
%     bead_color - Color of the nodes.
%     link_color - Color of the links.        
%
% PlotWebs Methods:
%     PlotWebs - Main Constructor
%     FillPositions - Fill positions
%     PlotBMatrix - Plot a original sorting graph layout
%     PlotBNestedMatrix - Plot a nested sorting graph layout
%     PlotBModularMatrix - Plot a modular sorting graph layout
%     DrawLine - Draw and edge for graph layouts
%     DrawCircle - Draw a node for graph layouts
%     PlotMatrix - Plot a original sorting matrix layout
%     PlotNestedMatrix - Plot a nested sorting matrix layout
%     PlotModularMatrix - Plot a modular sorting matrix layout
%     ApplyBasicFormat - Apply final format to the matrix layout
%     ApplyBasicBFormat - Apply final format to the graph layout
%     FillLabels - Create text labels
%     DrawBack - Chose the color for the background
%     DrawCell - Draw a cell for matrix layouts
%     MakeDivision - Make a module division.
%
% See also:
%    Printer.CreateCytoscapeData
classdef PlotWebs < handle;
   
    properties
       
        pos_x                = [];       % X coordinate of each cell in a unit square
        pos_y                = [];       % Y coordinate of each cell in a unit square
        cell_color           = [0 0 0];  % [1 1 1]; %Cell color
        back_color           = [1 1 1];  % [0 0 128]/255; %Back color
        margin               = 0.12;     % Margin between cels
        ax                   = 'image';  % Axis of the plot.
        border_color         = 'none'; 
        border_line_width    = 0.0005;
        isocline_color       = 'red';    % Color of the isocline
        division_color       = 'blue';   % Color used for the division of the modules
        line_width           = 1.5;      % Line width used in the isocline
        
    end
    
    properties
        
        matrix               = [];       % Adjacency matrix in boolean version
        webmatrix            = [];       % Adjacency matrix with edge weight information
        use_type_interaction = 0;        % Flag to color cells according to weight (only if discrete)
        use_row_sections     = 0;        % Flag to color sections
        use_module_format    = 1;        % Flag to give appropiate color format to modules
        biweb                = {};       % Bipartite object that an object of the class will reference to
        row_labels           = {};       % Text labels for row nodes
        col_labels           = {};       % Text labels fof column nodes
        use_labels           = 1;        % Flag that indicate the use of text labels
        FontSize             = 5;        % Font size for text labels. Change according to number of rows and columns.
        index_rows           = [];       % Index position of rows after modular or nested sorting
        index_cols           = [];       % Index position of columns after modular or nested sorting.
        colors               = [];       % Colors used for draw modules.
        color_interactions   = [];       % Color of interactions. Used only if use_type_interaction = 1
        color_section        = [];       % Color of the section. Used only if use_row_sections = 1
        n_rows               = [];       % Number of rows in the matrix
        n_cols               = [];       % Number of columns in the matrix
        interpreter          = 'none';   % none or latex for text labels.
        plot_iso_modules     = 1;        % Flag that indicate the plotting of the isocline inside modules
        use_isocline         = 1;        % Flag that indicated the plotting of isocline in a nested graph
    end
    
    properties
       
        radius                = 0.5;     % Radius of the nodes for graph layouts.
        vertical_margin       = 0.12;    % Vertical margin between nodes for graph layouts.
        horizontal_proportion = 0.5;     % Horizontal margin (proportional to the y total size) between nodes for graph layouts. 
        row_pos               = [];      % Positions of the row nodes.
        col_pos               = [];      % Positions of the column nodes.
        bead_color            = [0 0 0]; % Color of the nodes.
        link_color            = [0 0 0]; % Color of the links.
        
    end
    
    methods
       
        function obj = PlotWebs(bipmatrix_or_biweb)
        % PlotWebs - Main Constructor
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
            
            [obj.n_rows obj.n_cols] = size(obj.matrix);
            
            obj.colors = colormap('jet');
            obj.colors = obj.colors([23,2,13,42,57,20,15,11,9,16,3,28,26,24,46,59,41,18,56,40,17,48,27,53,6,62,5,60,14,32,64,19,36,58,39,21,4,8,35,30,50,63,25,51,55,34,61,37,47,44,54,43,38,12,52,33,31,1,22,29,10,45,49,7],:);
            obj.color_interactions(1,:) = [1 0 0];
            obj.color_interactions(2,:) = [0 0 1];
            obj.color_section = [0.75 0.75 0.75];
            
            if(length(unique(obj.webmatrix))>2)
                obj.use_type_interaction = 1;
            end
            
            if(length(unique(obj.biweb.row_class))>1)
                obj.use_row_sections = 1;
            end
            
            %Delete This section thereafter--------------------------
            obj.color_interactions(1,:) = [255 0 0]/255;
            obj.color_interactions(2,:) = [0 0 255]/255;
            obj.colors(2,:) = [0    0.7500    1.0000];
            obj.colors(3,:) = [1.0000    0.8750         0];
            obj.colors(4,:) = [186,85,211]/255;
            obj.colors(5,:) = [255,165,0]/255;
            obj.colors(6,:) = [34 139 34]/255;
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
        
        function PlotBMatrix(obj)
        % PlotBMatrix - Plot a original sorting graph layout
        % obj = PlotBMatrix(obj) - Plot a graph layout using the original sorting of nodes.
            cla;
            maxd = max(obj.n_rows,obj.n_cols);
            x1 = 1; x2 = 1+maxd*obj.horizontal_proportion;
            
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;            
            
            local_matrix = obj.matrix;
            
            hold on;
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j)==1)
                        plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'black');
                    end
                end
            end
            hold off;
            
            for i = 1:obj.n_rows
                obj.DrawCircle(x1,obj.row_pos(i),'black');
            end
            
            for j = 1:obj.n_cols
                obj.DrawCircle(x2,obj.col_pos(j),'black');
            end
            
            
            
            obj.ApplyBasicBFormat();
            
        end
        
        function PlotBNestedMatrix(obj)
        % PlotBNestedMatrix - Plot a nested sorting graph layout
        % obj = PlotBNestedMatrix(obj) - Plot a graph layout where nodes
        % are sorted such that the nestedness pattern is more visible
            
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
                        plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.link_color);
                    end
                end
            end
            hold off;
            
            for i = 1:obj.n_rows
                obj.DrawCircle(x1,obj.row_pos(i),obj.bead_color);
            end
            
            for j = 1:obj.n_cols
                obj.DrawCircle(x2,obj.col_pos(j),obj.bead_color);
            end
            
            obj.ApplyBasicBFormat();
        end
        
        function PlotBModularMatrix(obj)
        % PlotBModularMatrix - Plot a modular sorting graph layout
        % obj = PlotBModularMatrix(obj) - Plot a graph layout where nodes
        % are sorted such that the modular pattern is more visible            
            cla;
            maxd = max(obj.n_rows,obj.n_cols);
            x1 = 1; x2 = 1+maxd*obj.horizontal_proportion;
            
            if(isempty(obj.biweb))
                obj.biweb = Bipartite(obj.matrix);
            end
            if(obj.biweb.modules.done == 0)
                obj.biweb.modules.Detect();
            end
            
            obj.index_rows = obj.biweb.modules.index_rows;
            obj.index_cols = flipud(obj.biweb.modules.index_cols);
            
            row_mod = obj.biweb.modules.row_modules;
            col_mod = flipud(obj.biweb.modules.col_modules);
            
            local_matrix = obj.matrix(obj.index_rows,obj.index_cols);
            n_col = length(obj.colors);
                
            hold on;
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j)==1 && row_mod(i)~=col_mod(j))
                        plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color','black');
                    end
                end
            end
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(local_matrix(i,j)==1 && row_mod(i)==col_mod(j))
                        plot([1 x2],[obj.row_pos(i) obj.col_pos(j)],'Color',obj.colors(mod(row_mod(i),n_col),:));
                    end
                end
            end
            hold off;
            
            for i = 1:obj.n_rows
                obj.DrawCircle(x1,obj.row_pos(i),obj.colors(mod(row_mod(i),n_col),:));
            end
            
            for j = 1:obj.n_cols
                obj.DrawCircle(x2,obj.col_pos(j),obj.colors(mod(col_mod(j),n_col),:));
            end
            
            obj.ApplyBasicBFormat();
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
        
        function obj = PlotMatrix(obj)
        % PlotMatrix - Plot a original sorting matrix layout
        % obj = PlotMatrix(obj) - Plot a matrix layout using the original sorting of nodes.    
            cla;
            obj.index_rows = 1:obj.n_rows;
            obj.index_cols = 1:obj.n_cols;
            
            if(obj.use_row_sections)
                row_section = obj.biweb.rows_idx;
                row_section_unique = unique(row_section);
                for i = 1:length(row_section_unique)
                    if(mod(i,2)==1)
                        x = find(row_section==row_section_unique(i));
                        %MakeDivision(obj,i,j,nrow,ncol,bordercolor,backcolor)
                        obj.MakeDivision(x(end),1,length(x),obj.n_cols,obj.color_section,obj.color_section);
                    end
                end
            end
            
            for i = 1:obj.n_rows
                for j = 1:obj.n_cols
                    if(obj.matrix(i,j) > 0)
                        if(obj.use_type_interaction)
                            obj.DrawCell(i,j,obj.color_interactions(obj.webmatrix(i,j),:));
                        else
                            obj.DrawCell(i,j,obj.cell_color);
                        end
                    end
                end
            end
            
            obj.ApplyBasicFormat();
            
        end
        
        function obj = PlotNestedMatrix(obj)
        % PlotNestedMatrix - Plot a nested sorting matrix layout
        % obj = PlotNestedMatrix(obj) - Plot a matrix layout using the nested sorting of nodes.    
        
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
                    end
                end
            end
            
            hold on;
            [x,y] = NestednessBINMATNEST.GET_ISOCLINE(obj.matrix);
            p = plot(x,y);
            set(p,'Color',obj.isocline_color,'LineWidth',obj.line_width);
            hold off;
            
            obj.ApplyBasicFormat();
            
        end
        
        function obj = PlotModularMatrix(obj)
        % PlotModularMatrix - Plot a modular sorting matrix layout
        % obj = PlotModularMatrix(obj) - Plot a matrix layout using the modular sorting of nodes.    
        
            cla;
            if(isempty(obj.biweb))
                obj.biweb = Bipartite(obj.matrix);
            end
            if(obj.biweb.modules.done == 0)
                obj.biweb.modules.Detect();
            end
            
            obj.index_rows = obj.biweb.modules.index_rows;
            obj.index_cols = obj.biweb.modules.index_cols;
            
            row_mod = obj.biweb.modules.row_modules;
            col_mod = obj.biweb.modules.col_modules;
            
            local_matrix = obj.webmatrix(obj.index_rows,obj.index_cols);
            n_col = length(obj.colors);
            
            start_y = 1;
            start_x = obj.n_cols;
            
            for i = 1:obj.biweb.modules.N
               
                dx = length(find(col_mod==i));
                dy = length(find(row_mod==i));
                if(dx>0 && dy>0)
                    if(obj.use_module_format)
                        obj.MakeDivision(start_y+dy-1,start_x-dx+1,dy,dx,obj.colors(mod(i,n_col)+1,:),obj.back_color);
                    else
                        obj.MakeDivision(start_y+dy-1,start_x-dx+1,dy,dx,obj.division_color,obj.back_color);
                    end
                end
                start_x = start_x-dx;
                start_y = start_y+dy;
                
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
                                obj.DrawCell(i,j,obj.colors(mod(row_mod(i),n_col)+1,:));
                            else
                                obj.DrawCell(i,j,obj.cell_color);
                            end
                        end
                        
                    end
                end

            end

            if(obj.plot_iso_modules)
                matrices = obj.biweb.modules.ExtractCommunityMatrices();
                start_y = obj.n_rows;
                start_x = obj.n_cols;
                hold on;

                for i = 1:obj.biweb.modules.N

                    dx = length(find(col_mod==i));
                    dy = length(find(row_mod==i));
                    
                    if(isempty(matrices{i}) || sum(sum(matrices{i})) == numel(matrices{i}) || sum(sum(matrices{i})) == 0)
                        start_x = start_x-dx;
                        start_y = start_y-dy;
                        continue;
                    end
                    
                    [x,y] = NestednessBINMATNEST.GET_ISOCLINE(matrices{i});
                    p = plot(start_x-dx+x,start_y-dy+y);
                    %p = plot(x,y);
                    if(obj.use_module_format)
                        set(p,'Color',obj.colors(mod(i,n_col)+1,:),'LineWidth',obj.line_width);
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
        
        
        
        function obj = ApplyBasicFormat(obj)
        % ApplyBasicFormat - Apply final format to the matrix layout
        % obj = ApplyBasicFormat(obj) - Print the text labels and give a
        % final format to the matrix layout plot.
        
        %obj.FillLabels();
            
 %           idx_rows = obj.index_rows;
%            idx_type = obj.biweb.rows_idx;
 %           idx_type = idx_type(idx_rows);
            
            if(obj.use_labels)
                if(~isempty(obj.biweb))
                    obj.col_labels = obj.biweb.col_labels;
                    obj.row_labels = obj.biweb.row_labels;
                end
                
                obj.FillLabels();
                
                for i = 1:obj.n_rows
                    %if(idx_type(i) == 2)
                        text(0,obj.pos_y(i,1),obj.row_labels{i},'HorizontalAlignment','right','FontSize',obj.FontSize,'interpreter',obj.interpreter);
                    %else
                    %    text(0,obj.pos_y(i,1),obj.row_labels{i},'HorizontalAlignment','right','FontSize',obj.FontSize,'interpreter',obj.interpreter,'Color',[0 128 0]/255,'FontAngle','italic');
                    %end
                end
                for j = 1:obj.n_cols
                    text(obj.pos_x(1,j),0,obj.col_labels{j},'HorizontalAlignment','right','Rotation',90,'FontSize',obj.FontSize,'interpreter',obj.interpreter);
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
                    text(x1-obj.radius,obj.row_pos(i),obj.row_labels{i},'HorizontalAlignment','right','FontSize',obj.FontSize,'interpreter',obj.interpreter);
                end
                for j = 1:obj.n_cols
                    text(x2+obj.radius,obj.col_pos(j),obj.col_labels{j},'HorizontalAlignment','left','FontSize',obj.FontSize,'interpreter',obj.interpreter);
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
        
            rec = rectangle('Position',[obj.pos_x(i,j)-0.5+obj.margin,obj.pos_y(i,j)-0.5+obj.margin,1-2*obj.margin,1-2*obj.margin]);
            set(rec,'FaceColor',color);
            set(rec,'EdgeColor','none');

        end
        
        function obj = MakeDivision(obj,i,j,nrow,ncol,bordercolor,backcolor)
        % MakeDivision - Make a module division.
        % obj = MakeDivision(obj,i,j,nrow,ncol,bordercolor,backcolor) -
        % Make a module division starting in cell positioned at i,j of size
        % nrow by ncol cells. bordercolor and backcolor are used for border
        % and back color, respectivally.
            
             rec = rectangle('Position',[obj.pos_x(i,j)-0.5-obj.margin,obj.pos_y(i,j)-0.5-obj.margin, ncol+2*obj.margin, nrow+2*obj.margin],'LineWidth',1);
             set(rec,'FaceColor',bordercolor);
             set(rec,'EdgeColor','none');
             %set(rec,'LineStyle','--');
             rec = rectangle('Position',[obj.pos_x(i,j)-0.5+obj.margin,obj.pos_y(i,j)-0.5+obj.margin, ncol-2*obj.margin, nrow-2*obj.margin]);
             set(rec,'FaceColor',backcolor);
             set(rec,'EdgeColor','none');
     
             %rec = rectangle('Position',[obj.PosX(i,j)-0.5,obj.PosY(i,j)-0.5, ncol, nrow]);
             %set(rec,'FaceColor',obj.CellColor);
             %set(rec,'EdgeColor','black');

        end
        
    end
    
end