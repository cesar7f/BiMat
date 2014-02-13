% MetaStatisticsPlotter - Meta-analysis plotter. 
% This class allows the plotting of matrix and graph layouts when many
% networks are analysed and compared at the same time. Further, it includes
% some extra functions for a better visualization of the statistics
% results.
%
% MetaStatisticsPlotter Properties:
%   GENERAL
%     meta_statistics - MetaStatistics object that an object of this class will reference to
%     p_value - p-value for plotting random value bars and significance
%     z_value - z-value for plotting random value bars and significance
%     use_labels - Flag that indicate the use of text labels
%     font_size - Font size for text labels. Change according to number of rows and columns.
%     nest_test - 1 - for NTC, 2 - for NODF during graph and plot layout.
%     do_test_in_plots - Color labels according to statistical significance.
%     idx_to_be_ploted - Indexes of the networks that will be plotted
%   MATRIX LAYOUT
%     cell_color - Cell color
%     back_color - Back color
%     line_width - Line width used in the isocline
%     use_isocline - Flag that indicated the plotting of isocline in a nested graph
%   GRAPH LAYOUT
%     radius - Radius of the nodes for graph layouts.
%     vertical_margin - Vertical margin between nodes for graph layouts.
%     horizontal_proportion - Horizontal margin (proportional to the y total size) between nodes for graph layouts. 
%     bead_color_rows - Color of the row nodes.
%     bead_color_columns - Color of the column nodes.
%     link_color - Color of the links.
%     link_width - Edge width for grah payout plots
%
% PlotWebs Methods:
%     MetaStatisticsPlotter - Main Constructor
%     PlotNTCValues - Plot NTC nestedness values
%     PlotNODFValues - Plot NODF nestedness values
%     PlotModularValues - Plot modularity values
%     PlotMatrices - Plot all the networks in matrix layout in original sorting
%     PlotNestedMatrices - Plot all the networks in matrix layout in nested sorting
%     PlotModularMatrices - Plot all the networks in matrix layout in modular sorting
%     PlotGraphs - Plot all the networks in graph layout in original sorting
%     PlotNestedGraphs - Plot all the networks in graph layout in nested sorting
%     PlotModularGraphs - Plot all the networks in graph layout in modular sorting
%     
% See also:
%    PlotWebs, GroupStatistics
classdef MetaStatisticsPlotter;
   
    % GENERAL
    properties
        meta_statistics       = {}; % GroupStatistics object that an object of this class will reference to
        p_value               = Options.P_VALUE; % p-value for plotting random value bars and significance
        z_value               = Options.Z_VALUE; % z-value for plotting random value bars and significance
        use_labels            = false;    % Flag that indicate the use of text labels
        font_size             = 10;       % Font size for text labels. Change according to number of rows and columns.
        nest_test             = 1;        % 1 - for NTC, 2 - for NODF during graph and plot layout.
        do_test_in_plots      = true;     % Color labels according to statistical significance.
        idx_to_be_ploted      = [];       % Indexes of the networks that will be plotted
    end
    
    % MATRIX LAYOUT
    properties
        cell_color            = [0 0 0];  % Cell color
        back_color            = [1 1 1];  % Back color
        line_width            = 1.5;      % Line width used in the isocline
        use_isocline          = true;     % Flag that indicated the plotting of isocline in a nested graph
    end
    
    % GRAPH LAYOUT
    properties
        radius                = 0.5;      % Radius of the nodes for graph layouts.
        vertical_margin       = 0.12;     % Vertical margin between nodes for graph layouts.
        horizontal_proportion = 0.5;      % Horizontal margin (proportional to the y total size) between nodes for graph layouts. 
        bead_color_rows       = [1 0 0];  % Color of the row nodes.
        bead_color_columns    = [0 0 1];  % Color of the column nodes.
        link_color            = [0 0 0];  % Color of the links.
        link_width            = 0.5;      % Edge width for grah payout plots
    end
    
    methods
       
        function obj = MetaStatisticsPlotter(group_statistics)
        % MetaStatisticsPlotter - Main Constructor
        %
        %   MP = MetaStatisticsPlotter(META_STAT) Create a MetaStatisticsPlotter object
        %   called MP using a MetaStatistics object called META_STAT that
        %   can be used to perform diverse kinds of meta-analysis plots.
        
            obj.meta_statistics = group_statistics;
            obj.idx_to_be_ploted = 1:obj.meta_statistics.n_networks;
            
        end
        
        function obj = PlotNTCValues(obj,pvalue)
        % PlotNTCValues - Plot NTC nestedness values
        %
        %   obj = PlotNTCValues(obj) Plot NTC and random
        %   expectation values for all the networks that are being
        %   analyzed using a double tail p-value defined on obj.p_value
        %
        %   obj = PlotNTCValues(obj,P_VAL) Plot NTC and random
        %   expectation values for all the networks that are being
        %   analyzed using a double tail p-value defined on P_VAL.
        
            if(nargin == 1)
                obj.PlotValues(obj.meta_statistics.ntc_values.value, obj.meta_statistics.ntc_values.random_values);
            else
                obj.PlotValues(obj.meta_statistics.ntc_values.value, obj.meta_statistics.ntc_values.random_values,pvalue);
            end
            %Labels in title, y-axis and legends
            legend('Measured NTC','Random expectation',1,'Location','NorthWest');
            legend('boxoff')
            ylabel('Nestedness (NTC)','fontsize',16);

        end
        
        function obj = PlotNODFValues(obj,pvalue)
        % PlotNODFValues - Plot NODF nestedness values
        %
        %   obj = PlotNODFValues(obj) Plot NODF and random
        %   expectation values for all the networks that are being
        %   analyzed using a double tail p-value defined on obj.p_value
        %
        %   obj = PlotNODFValues(obj,P_VAL) Plot NODF and random
        %   expectation values for all the networks that are being
        %   analyzed using a double tail p-value defined on P_VAL.     
        
            if(nargin==1)
                obj.PlotValues(obj.meta_statistics.nodf_values.value, obj.meta_statistics.nodf_values.random_values);
            else
                obj.PlotValues(obj.meta_statistics.nodf_values.value, obj.meta_statistics.nodf_values.random_values,pvalue);
            end
            %Labels in title, y-axis and legends
            legend('Measured NODF','Random expectation',1,'Location','NorthWest');
            legend('boxoff')
            ylabel('Nestedness (NODF)','fontsize',16);
        end
        
        function obj = PlotModularValues(obj,pvalue)
        % PlotModularValues - Plot modularity values
        %
        %   obj = PlotModularValues(obj) Plot Qb and random
        %   expectation values for all the networks that are being
        %   analyzed using a double tail p-value defined on obj.p_value
        %
        %   obj = PlotModularValues(obj,P_VAL) Plot Qb and random
        %   expectation values for all the networks that are being
        %   analyzed using a double tail p-value defined on P_VAL.     
        
            if(nargin==1)
                obj.PlotValues(obj.meta_statistics.qb_values.value, obj.meta_statistics.qb_values.random_values);
            else
                obj.PlotValues(obj.meta_statistics.qb_values.value, obj.meta_statistics.qb_values.random_values,pvalue);
            end
            
            %Labels in title, y-axis and legends
            legend('Measured Modularity','Random expectation',1,'Location','NorthWest');
            legend('boxoff');
            ylabel('Modularity (Q)','fontsize',16);
        end
        
        
        
        function PlotMatrices(obj,n_rows,n_cols)
        % PlotMatrices - Plot all the networks in matrix layout in original
        % sorting
        %
        %   obj = PlotMatrices(obj) Perform a plot of all the bipartite
        %   network in matrix layout using original sorting. The number of
        %   subplots across rows and columns is estimated to equal or
        %   approximatly equal.
        %
        %   obj = PlotMatrices(obj, N_ROWS, N_COLS) Perform a plot of all the bipartite
        %   network in matrix layout using original sorting. The number of
        %   subplots across rows and columns is N_ROWS x N_COLS.
            
            if(nargin == 1)
                n_rows = floor(sqrt(obj.meta_statistics.n_networks));
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
            elseif(nargin == 2)
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
            end
            
            clf;
            obj.FormatPlotters();
            
            for i = 1:obj.meta_statistics.n_networks
                subplot(n_rows, n_cols, i);
                obj.meta_statistics.networks{i}.plotter.use_labels = 0;
                obj.meta_statistics.networks{i}.plotter.PlotMatrix();
                title(obj.meta_statistics.names{i}, 'FontSize',10);
            end
            
            set(gcf,'Position', [148         213        1142         746]);

        end
        
        function PlotNestedMatrices(obj,n_rows,n_cols,pvalue)
        % PlotNestedMatrices - Plot all the networks in matrix layout in
        % nested sorting
        %
        %   obj = PlotNestedMatrices(obj) Perform a plot of all the bipartite
        %   network in matrix layout using nested sorting. The number of
        %   subplots across rows and columns is estimated to equal or
        %   approximatly equal.
        %
        %   obj = PlotNestedMatrices(obj, N_ROWS, N_COLS) Perform a plot of all the bipartite
        %   network in matrix layout using nested sorting. The number of
        %   subplots across rows and columns is N_ROWS x N_COLS.
        
            if(nargin == 1)
                n_rows = floor(sqrt(obj.meta_statistics.n_networks));
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;
            elseif(nargin == 2)
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;   
            elseif(nargin==3)
                pvalue = obj.p_value;   
            end
            
            if(obj.nest_test == 1)
                sig_indices = obj.meta_statistics.ntc_values.percentile >= 100*(1-pvalue/2.0);
                no_sig_indices = obj.meta_statistics.ntc_values.percentile <= 50*pvalue;
            elseif(obj.nest_test == 2)
                sig_indices = obj.meta_statistics.nodf_values.percentile >= 100*(1-pvalue/2.0);
                no_sig_indices = obj.meta_statistics.nodf_values.percentile <= 50*pvalue;
            end
            
            clf;
            obj.FormatPlotters();
            
            for i = 1:obj.meta_statistics.n_networks
                subplot(n_rows, n_cols, i);
                obj.meta_statistics.networks{i}.plotter.use_labels = 0; %Do not show row/column labels
                obj.meta_statistics.networks{i}.plotter.use_isocline = 1;%No isocline inside modules
                obj.meta_statistics.networks{i}.plotter.PlotNestedMatrix();
                col = 'black'; % Color for not significance
                
                if(obj.do_test_in_plots == 1)
                    if(sig_indices(i) == 1) % Color for significant modularity
                        col = 'red';
                    elseif(no_sig_indices(i) == 1)%Color for significant antimodularity 
                        col = 'blue';
                    end
                end
                
                title(obj.meta_statistics.names{i},'Color',col, 'FontSize',10);
            end
            set(gcf,'Position', [148         213        1142         746]);
        end
        
        function PlotModularMatrices(obj,n_rows,n_cols,pvalue)
        % PlotModularMatrices - Plot all the networks in matrix layout in
        % modular sorting
        %
        %   obj = PlotModularMatrices(obj) Perform a plot of all the bipartite
        %   network in matrix layout using modular sorting. The number of
        %   subplots across rows and columns is estimated to equal or
        %   approximatly equal.
        %
        %   obj = PlotModularMatrices(obj, N_ROWS, N_COLS) Perform a plot of all the bipartite
        %   network in matrix layout using modular sorting. The number of
        %   subplots across rows and columns is N_ROWS x N_COLS.
                    
            
            if(nargin == 1)
                n_rows = floor(sqrt(obj.meta_statistics.n_networks));
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;
            elseif(nargin == 2)
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;   
            elseif(nargin==3)
                pvalue = obj.p_value;   
            end
            
            sig_indices = obj.meta_statistics.qb_values.percentile >= 100*(1-pvalue/2.0);
            no_sig_indices = obj.meta_statistics.qb_values.percentile <= 50*pvalue;
            
            clf;
            obj.FormatPlotters();
            
            for i = 1:obj.meta_statistics.n_networks
                subplot(n_rows, n_cols, i);
                obj.meta_statistics.networks{i}.plotter.use_labels = 0; %Do not show row/column labels
                obj.meta_statistics.networks{i}.plotter.plot_iso_modules = 0;%No isocline inside modules
                obj.meta_statistics.networks{i}.plotter.PlotModularMatrix();
                col = 'black'; % Color for not significance
                
                if(obj.do_test_in_plots == 1)
                    if(sig_indices(i) == 1) % Color for significant modularity
                        col = 'red';
                    elseif(no_sig_indices(i) == 1)%Color for significant antimodularity 
                        col = 'blue';
                    end
                end
                
                title(obj.meta_statistics.names{i},'Color',col, 'FontSize',10);
            end
            set(gcf,'Position', [148         213        1142         746]);

        end
        
        
        function PlotGraphs(obj,n_rows,n_cols)
        % PlotGraphs - Plot all the networks in graph layout in original
        % sorting
        %
        %   obj = PlotGraphs(obj) Perform a plot of all the bipartite
        %   network in graph layout using original sorting. The number of
        %   subplots across rows and columns is estimated to equal or
        %   approximatly equal.
        %
        %   obj = PlotGraphs(obj, N_ROWS, N_COLS) Perform a plot of all the bipartite
        %   network in graph layout using original sorting. The number of
        %   subplots across rows and columns is N_ROWS x N_COLS.
            
            if(nargin == 1)
                n_rows = floor(sqrt(obj.meta_statistics.n_networks));
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
            elseif(nargin == 2)
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
            end
            
            clf;
            obj.FormatPlotters();
            
            for i = 1:obj.meta_statistics.n_networks
                subplot(n_rows, n_cols, i);
                obj.meta_statistics.networks{i}.plotter.use_labels = 0;
                obj.meta_statistics.networks{i}.plotter.PlotGraph();
                title(obj.meta_statistics.names{i}, 'FontSize',10);
            end
            
            set(gcf,'Position', [148         213        1142         746]);

        end
        
        
        
        function PlotNestedGraphs(obj,n_rows,n_cols,pvalue)
        % PlotNestedGraphs - Plot all the networks in graph layout in
        % nested sorting
        %
        %   obj = PlotNestedGraphs(obj) Perform a plot of all the bipartite
        %   network in graph layout using nested sorting. The number of
        %   subplots across rows and columns is estimated to equal or
        %   approximatly equal.
        %
        %   obj = PlotNestedGraphs(obj, N_ROWS, N_COLS) Perform a plot of all the bipartite
        %   network in graph layout using nested sorting. The number of
        %   subplots across rows and columns is N_ROWS x N_COLS.
            
            if(nargin == 1)
                n_rows = floor(sqrt(obj.meta_statistics.n_networks));
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;
            elseif(nargin == 2)
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;   
            elseif(nargin==3)
                pvalue = obj.p_value;   
            end
            
            if(obj.nest_test == 1)
                sig_indices = obj.meta_statistics.ntc_values.percentile >= 100*(1-pvalue/2.0);
                no_sig_indices = obj.meta_statistics.ntc_values.percentile <= 50*pvalue;
            elseif(obj.nest_test == 2)
                sig_indices = obj.meta_statistics.nodf_values.percentile >= 100*(1-pvalue/2.0);
                no_sig_indices = obj.meta_statistics.nodf_values.percentile <= 50*pvalue;
            end
                
            clf;
            obj.FormatPlotters();
            
            for i = 1:obj.meta_statistics.n_networks
                subplot(n_rows, n_cols, i);
                obj.meta_statistics.networks{i}.plotter.use_labels = 0; %Do not show row/column labels
                obj.meta_statistics.networks{i}.plotter.PlotNestedGraph();
                col = 'black'; % Color for not significance
                
                if(obj.do_test_in_plots == 1)
                    if(sig_indices(i) == 1) % Color for significant modularity
                        col = 'red';
                    elseif(no_sig_indices(i) == 1)%Color for significant antimodularity 
                        col = 'blue';
                    end
                end
                
                title(obj.meta_statistics.names{i},'Color',col, 'FontSize',10);
            end
            set(gcf,'Position', [148         213        1142         746]);
        end
        
        
        
        function PlotModularGraphs(obj,n_rows,n_cols,pvalue)
        % PlotModularGraphs - Plot all the networks in graph layout in
        % modular sorting
        %
        %   obj = PlotModularGraphs(obj) Perform a plot of all the bipartite
        %   network in graph layout using modular sorting. The number of
        %   subplots across rows and columns is estimated to equal or
        %   approximatly equal.
        %
        %   obj = PlotModularGraphs(obj, N_ROWS, N_COLS) Perform a plot of all the bipartite
        %   network in graph layout using modular sorting. The number of
        %   subplots across rows and columns is N_ROWS x N_COLS.
            
            if(nargin == 1)
                n_rows = floor(sqrt(obj.meta_statistics.n_networks));
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;
            elseif(nargin == 2)
                n_cols = ceil(obj.meta_statistics.n_networks/n_rows);
                pvalue = obj.p_value;   
            elseif(nargin==3)
                pvalue = obj.p_value;   
            end
            
            sig_indices = obj.meta_statistics.qb_values.percentile >= 100*(1-pvalue/2.0);
            no_sig_indices = obj.meta_statistics.qb_values.percentile <= 50*pvalue;
            clf;
            obj.FormatPlotters();
            for i = obj.idx_to_be_ploted 
                subplot(n_rows, n_cols, i);
                obj.meta_statistics.networks{i}.plotter.PlotModularGraph();
                col = 'black'; % Color for not significance
                
                if(obj.do_test_in_plots == 1)
                    if(sig_indices(i) == 1) % Color for significant modularity
                        col = 'red';
                    elseif(no_sig_indices(i) == 1)%Color for significant antimodularity 
                        col = 'blue';
                    end
                end
                
                title(obj.meta_statistics.names{i},'Color',col, 'FontSize',10);
            end
            set(gcf,'Position', [148         213        1142         746]);
        end
        



    end
    
    methods(Access = 'protected')
        
        function obj = FormatPlotters(obj)
        % FormatPlotters - Format all the PlotWeb objects with the format
        % of the current MetaStatisticsPlotter object.
        %
        %   obj = FormatPlotters(obj) Format all the PlotWeb objects with the format
        %   of the current MetaStatisticsPlotter object.
        

            for i = 1:obj.meta_statistics.n_networks
               
                obj.meta_statistics.networks{i}.plotter.cell_color = obj.cell_color;
                obj.meta_statistics.networks{i}.plotter.back_color = obj.back_color;
                obj.meta_statistics.networks{i}.plotter.line_width = obj.line_width;
                obj.meta_statistics.networks{i}.plotter.use_labels = obj.use_labels;
                obj.meta_statistics.networks{i}.plotter.font_size = obj.font_size;
                obj.meta_statistics.networks{i}.plotter.use_isocline = obj.use_isocline;
                obj.meta_statistics.networks{i}.plotter.vertical_margin = obj.vertical_margin;
                obj.meta_statistics.networks{i}.plotter.horizontal_proportion = obj.horizontal_proportion;
                obj.meta_statistics.networks{i}.plotter.bead_color_rows = obj.bead_color_rows;
                
                obj.meta_statistics.networks{i}.plotter.bead_color_columns = obj.bead_color_columns;
                obj.meta_statistics.networks{i}.plotter.link_color = obj.link_color;
                obj.meta_statistics.networks{i}.plotter.link_width = obj.link_width;
                
            end
            
        end
        
        function obj = PlotValues(obj, values, random_values, pvalue)
        % PlotValues - Plot values with random expectations error bars
        %
        %   obj = PlotValues(obj,VALS, RANDOM_VALS) Perfom a plot, in which
        %   VALS are the points to be plotted, and RANDOM_VALS are values
        %   which corresponds to a null model. RANDOM_VALS are used for
        %   drawing error bars in which a double two-tail p-value =
        %   obj.p_value is used.
        %
        %   obj = PlotValues(obj,VALS, RANDOM_VALS,P_VAL) Perfom a plot, in which
        %   VALS are the points to be plotted, and RANDOM_VALS are values
        %   which corresponds to a null model. RANDOM_VALS are used for
        %   drawing error bars in which a double two-tail p-value =
        %   P_VAL is used.
   
            
            if(nargin == 3)
                pvalue = obj.p_value/2.0;
            end
            
            %pvalue can be at most 1.0
            assert(pvalue <= 1.0);
            
            sup_limit = 1 - (pvalue/2.0);
            inf_limit = pvalue/2.0;
            
            [~,sorted_indexes] = sort(values); % I will plot in increasing NTC value

            %Get random values and sort according to sorted_indexes
            values = values(sorted_indexes);
            mean_random_vals = mean(random_values,2);
            mean_random_vals = mean_random_vals(sorted_indexes);
            random_values = random_values(sorted_indexes,:); %sort in rows
            names_matrix = obj.meta_statistics.names(sorted_indexes);

            %Find the limits of the error bars using two tail p-value=0.05
            sup_bound = random_values(:,round(obj.meta_statistics.replicates * sup_limit));
            low_bound = random_values(:,round(obj.meta_statistics.replicates * inf_limit));

            %Plot the data of the real matrices
            clf;
            plot(1:obj.meta_statistics.n_networks, values,'o','MarkerFaceColor','red','MarkerEdgeColor','red');
            hold on;
            %Plot the data of the random values
            errorbar(1:obj.meta_statistics.n_networks, mean_random_vals, ... 
                mean_random_vals - low_bound, sup_bound - mean_random_vals, ...
                'o','MarkerFaceColor','white','MarkerEdgeColor','black');
            hold off;

            %Write the labels
            set(gca,'xticklabel',[]);
            for i=1:obj.meta_statistics.n_networks
                tmph=text(i,-0.01,names_matrix(i));
                set(tmph,'HorizontalAlignment','right');
                set(tmph,'rotation',90);
                set(tmph,'fontsize',obj.font_size);
            end


            %Give format to the matrix
            xlim([0 1+obj.meta_statistics.n_networks]);
            ylim([0 1]);

            %Give appropiate size to the figure window
            set(gcf,'Position',[91   135   859   505]);
            set(gca,'Units','pixels');
            set(gcf,'Position', [91   135   859   505+150])
            apos = get(gca,'position');
            apos(2) = apos(2) + 82;
            set(gca,'position',apos);
            set(gcf,'position',[91   135   859   596]);
        end
        
    end
    
end