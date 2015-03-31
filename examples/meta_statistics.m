%% Meta Statistics
% Some times the user may be in need of not only evaluating nestedness and/or
% modularity in a single bipartite network, but in a set of related
% networks. There exists many examples in the literature of how these kinds 
% of analysis are performed 
% (i.e. see <http://www.pnas.org/content/108/28/E288.abstract Flores et al
% 2011> for Phage-bacteria and
% <http://www.pnas.org/content/100/16/9383.full Bascompte et al 2003> for
% plant-pollinator networks). 
% |BiMat| comes with a couple of classes that the user can use in order to
% perform this kind of analysis. The name of these classes are
% |MetaStatistics| and |MetaStatisticsPlotter|. The main functionality of
% the first one is to take as input an entire set of bipartite networks and
% perform modularity and nested analysis among them. Finally, the second
% class is used in order to have a better visual understanding of the
% results.


%% Example: Seed dispersal networks
% This example will show how to perform this kind of analysis. In order to
% do that, we will use seed-dispersal bipartite networks. These networks
% were collected for Bascompte's lab website <http://www.web-of-life.es/
% The Web of Life>. We will explain the process of this kind of analysis
% from the point in which the user have already a set of data sets, which
% in these case correspond to bipartite matrices encoded in text files.
% The data used in this example can be located inside the |examples/data/|
% directory. 
%We will need first to read the data and create all the required matrices.
%This part of the code need to be executed inside examples directory
matrix_files = dir('data/seed_dispersal_matrices/M_SD_*');
n_files = length(matrix_files);
bip_networks = cell(n_files,1);
for i = 1:n_files
    bip_networks{i} = Reader.READ_BIPARTITE_MATRIX(['data/seed_dispersal_matrices/',matrix_files(i).name]);
end

%Now we can create our MetaStatistics instance:
meta_stat = MetaStatistics(bip_networks);
%%
% Once we have the |MetaStatistics| instance, we can configure it to
% perform modularity and/or nestedness tests. The default is to perform
% both, but the user can change the behavior by updating the properties
% |do_community| or |do_nestedness| to 0 (false). Equally possible is to
% chose the number of replicates by matrix, the null model, and the nested
% and modularity algorithms. If the properties are not updated, the |BiMat|
% will as usual use the default values specified on the file
% |base/Options.m|. Continuing with the experiment, we will just focus on
% evaluating Nestedness using NODF algorithm, with 100 replicates and the
% most simple null model:
%Configure the meta-analysis experiment
meta_stat.do_community = false;
meta_stat.do_nestedness = true; %default value
meta_stat.nestedness_algorithm = @NestednessNODF;
meta_stat.replicates = 100;
%%
% After configuring the experiment, we can just launch it by typing:
meta_stat.DoMetaAnalyisis();
%%
% After performing the experiment we can print or plot the results:
%print results
meta_stat.Print;
%%
%The next value will be used for the error bars of the plot
meta_stat.plotter.p_value = 0.05; %Default value
%plot the results
figure(1);meta_stat.plotter.PlotNestednessValues();
%%
% Finally, we can get a better idea of the results by plotting the nested
% matrices.
meta_stat.plotter.back_color = [255 229 204]/255;
meta_stat.plotter.cell_color = [96 96 96]/255;
meta_stat.plotter.isocline_color = [200 100 0]/255;
figure(2); meta_stat.plotter.PlotNestedMatrices;
%%
% As we can see, most of the matrices are nested (the ones with (+) on the
% label name. They are nested under the |p_value = 0.05| that we chose
% before using a two side test. This means that the nested matrices (+) are
% have a value larger than 97.5 % of their corresponding random networks.
% The user must know that we just repeated part of the results of a very
% famous paper on the field. If he has the interest, we would like to
% recommend him to read the original paper 
% <http://www.pnas.org/content/100/16/9383.full Bascompte et al 2003>.
%%
% Finally, we must say that the user can continue the analysis for
% modularity. The next lines could complete the work:
%%
% 
%   meta_stat.do_community = true;
%   meta_stat.modularity_algorithm = @AdaptiveBrim;
%   meta_stat.DoMetaAnalyisis();
%   figure(3); meta_stat.plotter.PlotModularValues();
%   figure(4); meta_stat.plotter.PlotModularMatrices();
% 

