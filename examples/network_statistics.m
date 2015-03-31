%% Network Statistics
% As described in the null models section, modularity and nestedness values
% lack statistical significance if they are not compared with an
% appropiate null model. Further, an appropiate ensamble of those
% random networks need to be created in order to analyze those networks.
% Altought the user can use |BiMat| null models for designing his own
% statistical tests, |BiMat| comes already with an |StatisticalTest| that
% helps the user to perform such tests.
% |BiMat| comes with statistical test for both modularity and nestedness.
% These tests can be executed directly on an instance of the |Bipartite|
% class (where the property |statistics| is an instance of the
% |StatisticalTest| class, or directly in the adjacency matrix.


%% Example: Working directly on a |Bipartite| instance
% This example will show the usual flow to evaluate the significance of
% modularity and nestedness on an specific matrix. We will use one of the
% phage bacteria matrices 
% (<http://www.sciencedirect.com/science/article/pii/S0168160509006576 Zinno 2010>) collected on 
% <http://www.pnas.org/content/108/28/E288.abstract Flores et al 2012>.
%loading and creating data
load phage_bacteria_matrices.mat;
bp = Bipartite(phage_bacteria_matrices.matrices{38}); %Matrices 38 corresponds to Zinno 2010
%A quick visual of how the data looks like:
pformat = PlotFormat();
pformat.use_labels = false;
pformat.back_color = [110,110,20]/255;
pformat.cell_color = 'white';
pformat.use_isocline_module = false;
bp.plotter.SetPlotFormat(pformat);

font_size = 16;

figure(1); set(gcf,'Position', [38,64,1290,435]);
subplot(1,3,1); bp.plotter.PlotMatrix; title('Original','FontSize',font_size);
subplot(1,3,2); bp.plotter.PlotNestedMatrix; title('Nested','FontSize',font_size);xlabel('Zinno 2010','FontSize',font_size+6);
subplot(1,3,3); bp.plotter.PlotModularMatrix; title('Modular','FontSize',font_size);

%% 
% By looking at the last plots, we can infer that this matrix have a
% strong community structure pattern, while nestedness is not apparent. In
% order to confirm these observations, we can finally perform nestedness
% and modularity tests on this matrix (using default values):
bp.statistics.TestCommunityStructure();
%%
bp.statistics.TestNestedness();
%%
% As we can see, |BiMat| shows the status of the current
% evaluation, which in this case was performed using 100 random matrices
% (the default value). In order to print the results, the user
% only need to call the |Print| method:
bp.statistics.Print;
%%
% The output just shows the configuration of each evaluation, which in this
% case was performed using the default values, which may not be a strong
% statistical test. Let's repeat the experiment using a larger number of
% random matrices and a different null model:
%The null model will be the AVERAGE one with 1000 replicates:
bp.statistics.DoNulls(1000, @NullModels.AVERAGE);
bp.statistics.TestCommunityStructure();
bp.statistics.TestNestedness();
%%
% And finally printing the results:
bp.statistics.Print();

%% Example: Working using the functional approach
% Sometimes the user just want to do a quick analysis, without the extra
% functionalities that working with a |Bipartite| instance provides. For
% those cases, the user can just execute the functional calls of the
% |StatisticalTest| class on a specific matrix. Using the previous matrix
% as an example:
matrix = phage_bacteria_matrices.matrices{38};
stest_modul = StatisticalTest.TEST_COMMUNITY_STRUCTURE(matrix,200,@NullModels.AVERAGE,@LeadingEigenvector);
%%
stest_nest = StatisticalTest.TEST_NESTEDNESS(matrix,200,@NullModels.EQUIPROBABLE,@NestednessNODF);
