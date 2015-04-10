%% BiMat - Group Testing Use case
% This use case will introduce the user to the functionalities about
% how to perform an statistical analysis of a group of bipartite networks
% (matrices). For doing that we will use the data from Flores et Al 2011.
% This data consist of 38 bipartite adjacency matrices of different sizes.
% Each matrix is named according to the first author paper
% from which it was extracted.
% We will perform an analysis of modularity and nestedness in the entire set.
%
%% Add the source to the matlab path
%Assuming that you run this script from examples directory
g = genpath('../'); addpath(g);
close all; %close all open figures
%%%
% We need also to load the data from which we will be working on
load phage_bacteria_matrices.mat;

%% Creating a Group Testing object
% If the number of random matrices and the null model are not assigned, 100 and
% AVERAGE are used as default. Here we will use 100 random matrices with
% the EQUIPROBABLE null model
mstat = MetaStatistics(phage_bacteria_matrices.matrices);    % Create the main object

%% Perform an statististical analysis in the set of matrices
% Suppose that we are interested in calculating the modularity and
% nestedness using the NTC algorithm as Flores et Al 2011 did. In addition,
% following the approach of Flores et Al 2011, we want to use the
% equiprobable model as null model in our random networks. The way to
% perform this analysis is by running the next lines:
mstat.replicates = 100; %How many random networks we want for each matrix
mstat.null_model = @NullModels.EQUIPROBABLE; %Our Null model
mstat.modularity_algorithm = @AdaptiveBrim; %Algorithm for modularity.
mstat.nestedness_algorithm = @NestednessNTC; %Algorithm for nestedness.
mstat.do_community = 1; % Perform Modularity analysis (default)
mstat.do_nestedness = 1; % Perform Nestedness analysis (default)
mstat.names = phage_bacteria_matrices.name;
mstat.DoMetaAnalyisis(); % Perform the analysis.
%%
% Notice that |DoMetaAnalyisis| method prints informatino about the current
% networks that is being analyzed, such that the user will know at every
% moment the current status of the analysis.
% After the analysis is finished a simple statistical measure to say that a
% matrix is nested and/or modular is to chose a two tail p-value = 0.05 as
% Flores et al 2011 did. Therefore, the next lines of code will show how
% many matrices are found nested and/or modular
fprintf('Number of nested matrices: %i\n',sum(mstat.N_values.percentile >= 97.5));
fprintf('Number of modular matrices: %i\n',sum(mstat.Qb_values.percentile >= 97.5));
%%
% Because we only did 100 random matrices you may get different results. For
% a more accurate result you may try 1.000 or even 10,000.
% We can also show the entire set of results by calling:
mstat.Print();

%% Ploting results
%The user can visualize the results of the last output in a graphical way. For example for visualizing
%the results of modularity and NTC nestedness value, the user can type:
mstat.plotter.p_value = 0.05; %p-value for error bars
mstat.plotter.font_size = 10; %Size for x-labels.
figure(1);
mstat.plotter.PlotModularValues();
figure(2);
mstat.plotter.PlotNestednessValues();
%%
% In addition the user can also plot the data in either graph or matrix
% layout. Here we show for graph nested layout and modular matrix layout.
% As in the case of a single networks, is possible to specify some of the
% most fundamental format properties.
mstat.plotter.p_value = 0.05; %p-value for color labeling
% Plot of nested graphs
mstat.plotter.bead_color_rows = 'blue'; %Color of row nodes
mstat.plotter.bead_color_columns = 'red'; %Color of column nodes
mstat.plotter.link_width = 0.5; %Edge width
mstat.plotter.use_isocline = false; %Do no show isocline inside modules
figure(3);
mstat.plotter.PlotModularMatrices(5,8); %Use a grid of 5 x 8
%Plot of modular matrices
figure(4);
mstat.plotter.PlotNestedGraphs(5,8);


