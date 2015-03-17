%This file includes the code for creating the figures that are shown on the
%paper ....
classdef FiguresPaper
    
    
    methods(Static)
       
       
        function CreateFigure1()
        % Create the panels of Figure 1.
        
            bip = Reader.READ_BIPARTITE_MATRIX('memmott_matrix.txt');
            
            % Map visit frequency to log scale. We do this in order to
            % decrease the number of colors in the interaction type
            matrix = bip.webmatrix;
            matrix = log2(matrix+2);matrix(~isfinite(matrix))=0;matrix=ceil(matrix);
            matrix = matrix - 1; %Avoid having 1'in 0's
            bip.webmatrix = matrix;
            % Specify the color of interactions
            bip.plotter.color_interactions = colormap(jet(max(matrix(:))));
            bip.plotter.use_labels = false;
            
            % Format plotter
            bip.plotter.line_width = 0.5; %isocline line width
            bip.plotter.link_width = 0.25; %Link width in graph layout
            bip.plotter.margin = 0.07; %Margin between cells in matrix. Normalized units to cell size.
            
            %Make graph layout size proportional to matrix layout
            bip.plotter.horizontal_proportion = bip.n_cols/bip.n_rows; 
            
            figure(1);
            clf;
            set(gcf,'position',[13 71 564 926]);
            
            %Font sizes
            titles_fsize = 16;
            label_fsize = 10;
            layout_fsize = 23;
            move_in_y = 0.07; % units are in normalized matlab figure units
            
            %MATRICES
            subplot(2,3,1);
            bip.plotter.PlotMatrix();
            set(gca,'position',get(gca,'position')-[0 move_in_y 0 0]);
            title('Original','FontSize',titles_fsize);
            text(-0.1,0.5,'Matrix Layout','units','normalized','rotation',90,'HorizontalAlignment',...,
                'center','VerticalAlignment','Bottom','FontSize',layout_fsize,'FontName','Times');
            
            subplot(2,3,2);
            bip.plotter.PlotNestedMatrix();
            set(gca,'position',get(gca,'position')-[0 move_in_y 0 0]);
            title('Nested','FontSize',titles_fsize);
            xlabel('Plants','FontSize',label_fsize);
            ylabel('Pollinators','FontSize',label_fsize);
            subplot(2,3,3);
            bip.plotter.PlotModularMatrix();
            set(gca,'position',get(gca,'position')-[0 move_in_y 0 0]);
            title('Modular','FontSize',titles_fsize);
            
            %GRAPHS
            subplot(2,3,4);
            bip.plotter.PlotGraph();
            text(-0.1,0.5,'Graph Layout','units','normalized','rotation',90,'HorizontalAlignment',...
                'center','VerticalAlignment','Bottom','FontSize',layout_fsize,'FontName','Times');
            
            subplot(2,3,5);
            bip.plotter.PlotNestedGraph();
            text(0,0.5,'Pollinators','units','normalized','rotation',90,'HorizontalAlignment',...
                'center','VerticalAlignment','Bottom','FontSize',label_fsize);
            text(1.0,0.5,'Plants','units','normalized','rotation',-90,'HorizontalAlignment',...
                'center','VerticalAlignment','Bottom','FontSize',label_fsize);
            xlabel('Plants','FontSize',label_fsize);
            ylabel('Pollinators','FontSize',label_fsize);
            subplot(2,3,6);
            bip.plotter.PlotModularGraph();
            
            
            
        end
       
        function mstat = CreateFigure3(mstat)
        % Create the two panels that are show in Figure 3. Due to the
        % randomness in the null model, the figures may not look as they
        % appear in the paper.
        
            close all;
            
            if nargin == 0
                
                load group_testing_data.mat;
                mstat = MetaStatistics(grouptesting.matrices); 
                mstat.replicates = 100; %How many random networks we want for each matrix
                mstat.null_model = @NullModels.EQUIPROBABLE; %Our Null model
                mstat.modularity_algorithm = @AdaptiveBrim; %Algorithm for modularity.
                mstat.nestedness_algorithm = @NestednessNTC; %Algorithm for nestedness.
                mstat.do_community = 1; % Perform Modularity analysis (default)
                mstat.do_nestedness = 1; % Perform Nestedness analysis (default)
                mstat.names = grouptesting.name;
                mstat.DoMetaAnalyisis(); % Perform the analysis.
                
            end
            
            mstat.plotter.p_value = 0.05; %p-value for error bars
            mstat.plotter.font_size = 10; %Size for x-labels.
            figure(1);
            mstat.plotter.PlotModularValues();
            figure(2);
            mstat.plotter.PlotNestednessValues();
            
        end
        
        function mstat = CreateFigure4(mstat)
        % Create Figure 4. Due to the
        % randomness in the null model, the figures may not look as they
        % appear in the paper.
        
            close all;
            
            if nargin == 0
                
                load group_testing_data.mat;
                mstat = MetaStatistics(grouptesting.matrices); 
                mstat.replicates = 100; %How many random networks we want for each matrix
                mstat.null_model = @NullModels.EQUIPROBABLE; %Our Null model
                mstat.modularity_algorithm = @AdaptiveBrim; %Algorithm for modularity.
                mstat.nestedness_algorithm = @NestednessNTC; %Algorithm for nestedness.
                mstat.do_community = 1; % Perform Modularity analysis (default)
                mstat.do_nestedness = 1; % Perform Nestedness analysis (default)
                mstat.names = grouptesting.name;
                mstat.DoMetaAnalyisis(); % Perform the analysis.
                
            end
            
            mstat.plotter.p_value = 0.05; %p-value for color labeling
            mstat.plotter.link_width = 0.5; %Edge width
            mstat.plotter.use_isocline = false; %Do no show isocline inside modules
            mstat.plotter.PlotModularMatrices(5,8); %Use a grid of 5 x 8


        end
        
        function bp = CreateFigure5(bp)
        % Create all plots show in Figure 5. Due to the randomness in the
        % algorithms, the figures may not look exactly as in the paper.
        %
        % 
        
            %Meta Analisis Figure
            close all;
            
            % Do initial analysis if no bipartite object was passed
            if nargin == 0
                load moebus_use_case;
                bp = Bipartite(moebus.weight_matrix);
                bp.nestedness = NestednessNTC(bp.webmatrix); %Use NTC instead of default (NODF)
                bp.community = AdaptiveBrim(bp.webmatrix); %Default algorithm
                bp.community.Detect();
                %focus only
                bp.internal_statistics.idx_to_focus_on = 1:15;
                bp.internal_statistics.TestInternalModules();
            end
            bp.plotter.font_size = 2.0;
            
            figure(1);
            bp.plotter.PlotModularMatrix();
            set(gcf,'position',[9    71   902   926]);
            
            figure(2);
            bp.internal_statistics.meta_statistics.plotter.font_size = 16;
            bp.internal_statistics.meta_statistics.plotter.PlotNestednessValues();
            set(get(gca,'ylabel'),'fontsize',24)
            set(gca,'fontsize',16);
            
            
            figure(3);
            %Plot cell color according to interaction type (strong/weak).
            bp.internal_statistics.meta_statistics.plotter.use_type_interaction = true;
            bp.internal_statistics.meta_statistics.plotter.font_size = 16;
            bp.internal_statistics.meta_statistics.plotter.use_specific_colors = true;
            bp.internal_statistics.meta_statistics.plotter.colors = bp.plotter.colors;
            bp.internal_statistics.meta_statistics.plotter.isocline_color = 'black';
            bp.internal_statistics.meta_statistics.plotter.PlotNestedMatrices();
            set(gcf,'position',[880   116   951   563]);

        end
        
         function bp = CreateFigureMainPaper(bp)
        % Create all plots show in Figure 5. Due to the randomness in the
        % algorithms, the figures may not look exactly as in the paper.
        %
        % 
        
            %Meta Analisis Figure
            close all;
            
            % Do initial analysis if no bipartite object was passed
            if nargin == 0
                load moebus_use_case;
                bp = Bipartite(moebus.weight_matrix);
                bp.nestedness = NestednessNTC(bp.webmatrix); %Use NTC instead of default (NODF)
                bp.community = AdaptiveBrim(bp.webmatrix); %Default algorithm
                bp.community.Detect();
                %focus only
                bp.internal_statistics.idx_to_focus_on = 1:15;
                bp.internal_statistics.TestInternalModules();
            end
            bp.plotter.font_size = 2.0;
            
            figure(1);
            bp.plotter.use_module_format = false;
            bp.plotter.isocline_color = 'black';
            bp.plotter.division_color = 'black';
            bp.plotter.link_width = 0.53;
            bp.plotter.use_labels = false;
            bp.plotter.PlotModularMatrix();
            set(gcf,'position',[9    71   902   926]);
            
            figure(2);
            bp.internal_statistics.meta_statistics.plotter.font_size = 16;
            bp.internal_statistics.meta_statistics.plotter.PlotNestednessValues();
            set(get(gca,'ylabel'),'fontsize',24)
            set(gca,'fontsize',16);
            
            
            figure(3);
            %Plot cell color according to interaction type (strong/weak).
            bp.internal_statistics.meta_statistics.plotter.use_type_interaction = true;
            bp.internal_statistics.meta_statistics.plotter.font_size = 16;
            bp.internal_statistics.meta_statistics.plotter.use_specific_colors = false;
            bp.internal_statistics.meta_statistics.plotter.colors = bp.plotter.colors;
            bp.internal_statistics.meta_statistics.plotter.isocline_color = 'black';
            bp.internal_statistics.meta_statistics.plotter.PlotNestedMatrices();
            set(gcf,'position',[880   116   951   563]);

        end
        
    end
    
end