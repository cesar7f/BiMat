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
load group_testing_data.mat;

%% Creating a Group Testing object
% If the number of random matrices and the null model are not assigned, 100 and
% AVERAGE are used as default. Here we will use 100 random matrices with
% the EQUIPROBABLE null model
gp = GroupStatistics(grouptesting.matrices);    % Create the main object

%% Perform an statististical analysis in the set of matrices
% Suppose that we are interested in calculating the modularity and
% nestedness using the NTC algorithm as Flores et Al 2011 did. In addition,
% following the approach of Flores et Al 2011, we want to use the
% equiprobable model as null model in our random networks. The way to
% perform this analysis is by running the next lines:
gp.replicates = 100; %How many random networks we want for each matrix
gp.null_model = @NullModels.EQUIPROBABLE; %Our Null model
gp.modul_class = @AdaptiveBrim; %Algorithm for modularity.
gp.do_temp = 1; % Perform NTC analysis (default)
gp.do_modul = 1; % Perform Modularity analysis (default)
gp.do_nest = 0; % Do not perform NODF analysis
gp.DoGroupTesting(); % Perform the analysis.
gp.names = grouptesting.name;
%%
% Notice that |DoGroupTesting| method prints informatino about the current
% networks that is being analyzed, such that the user will know at every
% moment the current status of the analysis.
% After the analysis is finished a simple statistical measure to say that a
% matrix is nested and/or modular is to chose a two tail p-value = 0.05 as
% Flores et al 2011 did. Therefore, the next lines of code will show how
% many matrices are found nested and/or modular
fprintf('Number of nested matrices: %i\n',sum(gp.tempvals.percent >= 97.5));
fprintf('Number of modular matrices: %i\n',sum(gp.qb_vals.percent >= 97.5));
%%
% We can also show the entire set of results by calling:
gp.PrintResults();

%% Using a GroupStatistics object to create your own plots
% We can use a GroupStatistics object (gp in this case) to create specific
% plots. Suppose we are interested in plotting all the matrices in modular
% sorting, such that the lables in the matrix are in red if the matrix is
% modular and in blue if it is antimodular. A simple script for performing
% this task will be:
n_rows = 5;
n_cols = 8;
modular_indices = gp.qb_vals.percent >= 97.5;
no_modular_indices = gp.qb_vals.percent <= 2.5;
figure(1);
for i = 1:gp.n_networks
    subplot(n_rows, n_cols, i);
    gp.networks{i}.plotter.use_labels = 0; %Do not show row/column labels
    gp.networks{i}.plotter.use_isocline = 0; %No isocline inside modules
    gp.networks{i}.plotter.PlotModularMatrix();
    col = 'black'; % Color for not significance
    if(modular_indices(i) == 1) % Color for significant modularity
        col = 'red';
    elseif(no_modular_indices(i) == 1) % Color for significant antimodularity
        col = 'blue';
    end
    title(gp.names{i},'Color',col, 'FontSize',10);
end
set(gcf,'Position', [148         213        1142         746]);
%%
% We may want to create a plot that compare the values of the networks with
% the random values of the null model. The next lines will show how to
% create such plot for the case of the NTC results
ntc_vals = gp.tempvals.value;
[~,sorted_indexes] = sort(ntc_vals); % I will plot in increasing NTC value

%Get random values and sort according to sorted_indexes
ntc_vals = ntc_vals(sorted_indexes);
mean_random_vals = gp.tempvals.mean(sorted_indexes);
random_values = gp.tempvals.random_values; %variable already sorted in rows
random_values = random_values(sorted_indexes,:); %sort in rows
names = gp.names(sorted_indexes);

%Find the limits of the error bars using two tail p-value=0.05
sup_bound = random_values(:,round(gp.replicates * 0.975));
low_bound = random_values(:,round(gp.replicates * 0.025));

%Plot the data of the real matrices
figure(2);
plot(1:gp.n_networks, ntc_vals,'o','MarkerFaceColor','red','MarkerEdgeColor','red');
hold on;
%Plot the data of the random values
errorbar(1:gp.n_networks, mean_random_vals, mean_random_vals - low_bound, ...,
    sup_bound - mean_random_vals, 'o','MarkerFaceColor','white','MarkerEdgeColor','black');
hold off;

%Write the labels
set(gca,'xticklabel',[]);
for i=1:gp.n_networks
    tmph=text(i,-0.01,names(i));
    set(tmph,'HorizontalAlignment','right');
    set(tmph,'rotation',90);
    set(tmph,'fontsize',10);
end

%Labels in title, y-axis and legends
tmplh = legend('Measured modularity','Random expectation',1,'Location','NorthWest');
legend('boxoff')
title('Nestedness in Bacteria-Phage Networks','fontsize',20);
ylabel('Nestedness (NTC)','fontsize',16);

%Give format to the matrix
xlim([1 1+gp.n_networks]);
ylim([0 1]);

%Give appropiate size to the figure window
set(gcf,'Position',[91   135   859   505]);
set(gca,'Units','pixels');
set(gcf,'Position', [91   135   859   505+150])
apos = get(gca,'position');
apos(2) = apos(2) + 82;
set(gca,'position',apos);
set(gcf,'position',[91   135   859   596]);


