%% BiMat Use case using Moebus cross-infection matrix data
% This use case will introduce the user to the most basic functionalities
% of the BiMat Software. In order to do that we will calculate some of the
% results presented on the Flores et al 2012 paper (Multi-scale structure
% and geographic drivers of cross-infection within marine bacteria and
% phages). We will show how to perform the next things:
%

%% Add the source to the matlab path
%Assuming that you run this script from examples directory
g = genpath('../'); addpath(g);
close all;
%%%
% We need also to load the data from which we will be working on
load moebus_use_case.mat;

%% Creating the Bipartite network object
bp = Bipartite(moebus.weight_matrix);    % Create the main object
bp.row_labels = moebus.bacteria_labels;  % Updating node labels
bp.col_labels = moebus.phage_labels;     
bp.row_class = moebus.bacteria_stations;   % Updating node ids
bp.col_class = moebus.phage_stations;
%%%
% We can print the general properties of the network with:
bp.printer.PrintGeneralProperties();
%% Calculating Modularity
% The modularity algorithm is encoded in the property modules of the
% bipartite object (|bp.modules|). Tree algorithms are available:
%
% # Adaptive BRIM (|AdaptiveBrim.m|)
% # LP&BRIM (|LPBrim.m|)
% # Newman Algorithm (|NewmanAlgorithm.m|)
%
% Adaptive BRIM algorithm is assigned by default during the creation of the
% Bipartite object. However, we can assign another algorithm dinamically.
% For example to change to Newman Algorithm:
bp.modules = NewmanModularity(bp.matrix);
bp.modules.DoKernighanLinTunning = true; % This flag is exclusive of the Newman Algorithm.
bp.modules = AdaptiveBrim(bp.matrix); % Return to the default algorithm
%%%
% We need to calculate the modularity by calling:
bp.modules.Detect();
%%
% If we are interested only in node community indexes we can use
% |bp.modules.row_modules| and |bp.modules.col_modules|. However for
% modularity values we need to call |bp.modules.Qb| or |bp.modules.Qr| as
% is:
% the next lines of code:
fprintf('The modularity value Qb is %f\n', bp.modules.Qb);
fprintf('The fraction of interactions inside modules Qr is %f\n',bp.modules.Qr);
%%
% The value 0<=Qb<=1 is calculated using the standard bipartite modularity
% function (introduced by Barber):
% $$Q_b = \frac{1}{m} Trace(R^T (B-P) T)$
%
% while the value 0<=Qr<=1 represents the fraction of interactions that fall
% inside modules:
% $$Q_r = \frac{1}{m} \sum_{i=1}^m \sum_{j=1}^n B_{ij}\delta(g_i,g_j)$

%% Calculating Nestedness
% Two algorithms exist for calculating nestedness. Contrary to the case
% of modularity, in this case there is no need to change the algorithm
% because all the algorithms have an independent property in the Bipartite
% object. These algorithms are:
%
% * NODF (Nestedness metric based on overlap and decreasing filling). With
% value in the interval [0,1].
% * NTC (Nestedness Temperature Calculator) With value in the interval [0
% 1].
%
% The first algorithm is runned during the creation of the Bipartite
% object, but because the NTC algorithm is slow, you need to run the
% algorithm explicitally:
bp.ntc.CalculateNestedness();
%%
% Finally to show acces the values of the two algorithms you need to call:
fprintf('The nestedness NODF value is %f\n', bp.nodf.nodf); %The same value will be printed all the times
fprintf('The nestedness NTC value is %f\n', bp.ntc.N); %Because the value depends in the sorting of rows and columns, it may variate from trial to trial.
%%%
% We can print all the values of structure values by just calling:
bp.printer.PrintStructureValues
%% Plotting in Matrix Layout
% You can print the layout of the original, nestedness and modular sorting.
% If you matrix is weighted in a categorical way using integers
% (0,1,2...) you can visualize a different color for each
% interaction, where 0 is no interaction. For using this functionality you
% need to assign a color for each interaction and specifically indicate
% that you want a color for each interaction before calling the plot function:
figure(1);
set(gcf,'Position',[19    72   932   922]); % Matlab command to change the figure window;
bp.plotter.FontSize = 2.5; %Change the font size of the rows and labels
bp.plotter.use_type_interaction = true; % I want to use different color for each kind of interaction
bp.plotter.color_interactions(1,:) = [1 0 0]; %Red color for clear lysis
bp.plotter.color_interactions(2,:) = [0 0 1]; %Blue color for turbid spots
bp.plotter.back_color = 'white';
bp.plotter.PlotMatrix(); %After changing all the format we finally can call the plotting function.
%%
% For plotting the nestedness matrix you may decide to use or not an
% iscoline. The nestedness pattern is just the matrix sorted in decreasing
% degree for row and column nodes.
figure(2);
set(gcf,'Position',[19+200    72   932   922]); % Matlab command to change the figure window;
bp.plotter.use_isocline = true; % The isocline is used in the NTC algorithm.
bp.plotter.isocline_color = 'red'; %Decide the color of the isocline.
bp.plotter.PlotNestedMatrix();
%%
% For plotting the modularity sort, the plot function will calculate the
% modularity (call |bp.modules.Detect()|) if you have not previouslly call
% it. 
figure(3);
set(gcf,'Position',[19+400    72   932   922]); % Matlab command to change the figure window;
bp.plotter.plot_iso_modules = true; %This will plot isoclines inside modules.
bp.plotter.PlotModularMatrix();

%% Statistical analysis in the entire network
% We can perform an statistical analysis in the entire network for
% nestedness and modularity. In order to make an statistical analysis of
% the structure values we need to decide how many replicates we will need
% and what null model is more convenient for what we need. File
% |NullModels.m| contain all the availaible null models, while file
% |StatisticalTest.m| contains all the functions required for performing this
% analysis. We can perfom the analysis in all the structure values by the
% use of a single function call the required functions for nestedness and
% modularity analysis:
bp.statistics.print_output = 0; % Do not show outputs
bp.statistics.DoCompleteAnalysis(10, @NullModels.EQUIPROBABLE); %Bernoulli
bp.statistics.DoCompleteAnalysis(10, @NullModels.AVERAGE); %Probability degree
bp.printer.PrintStructureStatistics(); %Print the statistical values
%%
% All structure statistics calculate the next numbers:
%
% # value: Structure value of the tested bipartite network
% # p: p-value of a t-test
% # ci: Confidence interval of the mean of the random bipartite networks
% using a t-test
% # percent: Ratio of how many random matrices have an structure value
% smallen than the tested network
% # z-score: Z-score of the value in the tested network using the values of
% the random networks. Be careful when using this value because the random
% distribution may not be normal.
% # replicates: The quantity of random networks used in the analysis.
%
%% Statistical Analysis of the internal modules
% In addition to be able to perform structure analysis in the entire
% network, we may be able (depending in the size and module configuration
% of the tested matrix) to perform structure analysis in the internal
% modules. We will show next (i) how to do an analysis of modularity and
% nestedness in the internal modules and (ii) how to test for a possible
% correlation between node labeling and module configuration.
% All the functions for performing this kind of analysis is encoded in file
% |InternalStatistics.m|. For calculating the the statistical structure of the
% internal modules we just need to call:
bp.internal_statistics.TestInternalModules(20,@NullModels.EQUIPROBABLE); %100 replicates and Bernoulli as null model
bp.printer.PrintStructureStatisticsOfModules(); % Print the results
%%
% We can also study if a correlation exist between the the row labeling and
% the module configuration. For performing this analysis we always will
% need a labaling for rows and/or columns that group them in different
% sets. In this case we have as labeling the station number from which the
% bacteria and phages were extracted. Therefore what we will study is if
% there exist a correlation between the station location (geography) and
% the module configuration. We will use the same method that was used in
% Flores et al 2012. Given the labeling this method calculates the
% diversity index of the labeling inside each module and compare it with
% random permutations of the labeling across the matrix.
bp.internal_statistics.TestDiversityRows(1000); %Using the labeling of bp and 1000 random permutations
bp.internal_statistics.TestDiversityColumns(1000,moebus.phage_stations,@Diversity.SHANNON_INDEX); % Using specific labeling and Shannon index
bp.printer.PrintColumnModuleDiversity(); %Print the information of column diversity
