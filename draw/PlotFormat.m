% PlotWebs - Matrix and graph plotting class. 
% This class allows the plotting of matrix and graph layouts that represent
% the bipartite network.
%
% PlotFormat Properties:
%   GENERAL
%     use_labels - Flag that indicate the use of text labels
%     font_size - Font size for text labels. Change according to network size
%     colors - Colors used for draw modules
%     interpreter - none or latex for text labels
%     ax - Axis of the plot
%   MATRIX LAYOUT
%     cell_color - Cell color
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
% See also:
%    PlotWebs
classdef PlotFormat < handle;
    
    % GENERAL PROPERTIES
    properties
        
        use_labels           = 0;       % Flag that indicate the use of text labels
        font_size            = 5;        % Font size for text labels. Change according to network size.
        colors               = [];       % Colors used for draw modules.
        interpreter          = 'none';   % none or latex for text labels.
        ax                   = 'image';  % Axis of the plot.
        
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
    
    methods(Access = 'private')
        function valid = isValidColor(~,value)
            valid_str_values = {'yellow', 'magenta', 'cyan', 'red', 'green', 'blue', 'white', 'black', 'y', 'm', 'c', 'r', 'g','b','w', 'k'};
            valid = true;
            if(~(isa(value,'char') || isnumeric(value)))
                valid = false;
            end
            
            if(isa(value,'char')  && sum(ismember(valid_str_values,lower(value))) == 0)
                valid = false;
            elseif(isnumeric(value))
                % Test for a vector of 1x3 or 3x1
                if((size(value,1)+size(value,2) ~= 4) || (mod(size(value,1),2)==0) )
                    valid = false;
                % Test for all values in the range (0,1)
                elseif(sum((value>=0) + (value <= 1)) ~= 6)
                    valid = false;
                end
            end
        end
    end
    
    methods

        function set.use_labels(obj,value)
            if(value ~= 0 && value ~= 1); throw(MException('PlotFormat:use_labels', 'Value must be a boolean value')); end;
            obj.use_labels = value;
        end
       
        function set.font_size(obj,value)
            if(isnumeric(value) == false || value <= 0); throw(MException('PlotFormat:font_size', 'Font Size must be a positive number')); end;
            obj.font_size = value;
        end
        
        function set.colors(obj, value)
            obj.colors = value;
        end
        
        function set.interpreter(obj,value)
            value = lower(value);
            if(strcmp(value,'tex') == 0 && strcmp(value,'latex') == 0 && strcmp(value,'none') == 0)
                throw(MException('PlotFormat:interpreter', 'interpreter must be ''none'', ''latex'', or ''tex'' '))
            end
            obj.interpreter = value;
        end
        
        function set.ax(obj,value)
            valid_values = {'auto', 'manual', 'tight', 'fill', 'ij', 'xy', 'equal', 'image', 'square', 'vis3d', 'normal', 'off', 'on'};
            value = lower(value);
            
            if(sum(ismember(valid_values,value)) == 0)
                throw(MException('PlotFormat:ax', 'AX must be a valid Matlab axis'))
            end
            obj.ax = value;
        end
        
      
         
        
        %cell_color;%           = [0 0 0];  % Cell color
        %back_color;%           = [1 1 1];  % Back color
        %margin;%               = 0.12;     % Margin between cells
        %isocline_color;%       = 'red';    % Color of the isocline
        %division_color;%       = 'blue';   % Color used for the division of the modules
       
        function set.cell_color(obj,value)
            if(obj.isValidColor(value) == false)
                throw(MException('PlotFormat:cell_color', 'Cell Color must be a valid Matlab color'));
            end
            
            obj.cell_color = value;
        end       
        
        function set.back_color(obj,value)
            if(obj.isValidColor(value) == false)
                throw(MException('PlotFormat:back_color', 'Back Color must be a valid Matlab color'));
            end
            
            obj.back_color = value;
        end    
        
        function set.margin(obj,value)
            if(isnumeric(value) == false || value > 0.5 || value < 0); throw(MException('PlotFormat:margin', 'Margin between matrix cells must be in the range (0,0.5)')); end;
            
            obj.margin = value;
        end       
        
        function set.isocline_color(obj,value)
            if(obj.isValidColor(value) == false)
                throw(MException('PlotFormat:isocline_color', 'Isocline color must be a valid Matlab color'));
            end
            
            obj.isocline_color = value;
        end   
        
        function set.division_color(obj,value)
            if(obj.isValidColor(value) == false)
                throw(MException('PlotFormat:division_color', 'Division color must be a valid Matlab color'));
            end
            
            obj.division_color = value;
        end   
        
%         line_width;%           = 1.5;      % Line width used in the isocline
%         use_isocline;%         = true;     % Flag that indicated the plotting of isocline in a nested graph
%         color_interactions;%   = [];       % Color of interactions. Used only if use_type_interaction = true
%         use_type_interaction;% = false;    % Flag to color cells according to weight (only if discrete)
%         use_module_format;%    = true;     % Flag to give appropiate color format to modules
%         
         function set.line_width(obj,value)
            if(isnumeric(value) == false || value <= 0); throw(MException('PlotFormat:line_width', 'Line width must be a positive number')); end;
            obj.line_width = value;
        end  
        
        function set.use_isocline(obj,value)
            if(value ~= 0 && value ~= 1); throw(MException('PlotFormat:use_isocline', 'use_isocline value must be a boolean value')); end;
            obj.use_isocline = value;
        end
        
        function set.color_interactions(obj,value)
            obj.color_interactions = value;
        end
        
        function set.use_type_interaction(obj,value)
            if(value ~= 0 && value ~= 1); throw(MException('PlotFormat:use_type_interaction', 'Use of type interactions must be a boolean value')); end;
            obj.use_type_interaction = value;
        end
        
        function set.use_module_format(obj,value)
            if(value ~= 0 && value ~= 1); throw(MException('PlotFormat:use_module_format', 'Use Module Format value must be a boolean value')); end;
            obj.use_module_format = value;
        end
        
        function set.use_isocline_module(obj,value)
            if(value ~= 0 && value ~= 1); throw(MException('PlotFormat:use_isocline_module', 'Use isocline module value must be a boolean value')); end;
            obj.use_isocline_module = value;
        end
        %     use_isocline_module - Flag to indicate if isocline will be plotted inside module sorting
        
%         radius;%                = 0.5;     % Radius of the nodes for graph layouts.
%         vertical_margin;%       = 0.12;    % Vertical margin between nodes for graph layouts.
%         horizontal_proportion;% = 0.5;     % Horizontal margin (proportional to the y total size) between nodes for graph layouts. 
%         bead_color_rows;%       = [1 0 0]; % Color of the row nodes.
%         bead_color_columns;%    = [0 0 1]; % Color of the column nodes.
%         link_color;%            = [0 0 0]; % Color of the links.
%         link_width;%            = 0.75;     % Edge width
       
        function set.radius(obj,value)
            if(isnumeric(value) == false || value <= 0); throw(MException('PlotFormat:radius', 'Node Radius must be a positive number')); end;
            obj.radius = value;
        end  
        
        function set.vertical_margin(obj,value)
            if(isnumeric(value) == false || value <= 0); throw(MException('PlotFormat:vertical_margin', 'Vertical vargin must be a positive number')); end;
            obj.vertical_margin = value;
        end  
        
        function set.horizontal_proportion(obj,value)
            if(isnumeric(value) == false || value <= 0); throw(MException('PlotFormat:horizontal_proportion', 'Horizontal proportion must be a positive number')); end;
            obj.horizontal_proportion = value;
        end  

        function set.bead_color_rows(obj,value)
            if(obj.isValidColor(value) == false)
                throw(MException('PlotFormat:bead_color_rows', 'Row nodse color must be a valid Matlab color'));
            end
            
            obj.bead_color_rows = value;
        end   
        
        function set.bead_color_columns(obj,value)
            if(obj.isValidColor(value) == false)
                throw(MException('PlotFormat:bead_color_columns', 'Column nodes color must be a valid Matlab color'));
            end
            
            obj.bead_color_columns = value;
        end           
        
        function set.link_color(obj,value)
            if(obj.isValidColor(value) == false)
                throw(MException('PlotFormat:link_color', 'Link (edge) color must be a valid Matlab color'));
            end
            
            obj.link_color = value;
        end           
        
        function set.link_width(obj,value)
            if(isnumeric(value) == false || value <= 0); throw(MException('PlotFormat:link_width', 'Link (edge) witdh must be a positive number')); end;
            obj.link_width = value;
        end  
        
    end

    
end