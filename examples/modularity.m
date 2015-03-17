%% Modularity

%% Description
% Modularity indicates the presence of dense clusters of related nodes embedded 
% within the network. In many systems, we can 
% find a partition of nodes into specific communities or modules.  
% This modularity metric an be expressed for a bipartite network as:
% <latex>
% \begin{equation}
% Q = \frac{1}{E} \sum_{ij} \left( B_{ij} - \frac{k_i d_j}{E} \right) \delta(g_i,h_j),
% \end{equation}
% </latex>
% where $B_ij$ is the element in the bipartite matrix representing a link
% (1)
% or no link (0) between nodes $i$ and $j$, $g_i$ and $h_i$ are the module indexes of nodes $i$ (that belongs to set
% $R$) and $j$ (that belongs to set $C$), $k_i$ is the degree of node $i$, $d_j$ is
% the degree of node $j$ and $E$ is the number of links in the network.  
% |BiMat| contains three algorithms to maximize the last equation, all of
% them containing different heuristics. 
%
% * |AdaptiveBrim|: The algorithm uses the BRIM algorithm
%   with an heuristic for looking the optimal value of modules. This heuristic consists
%   in multiply the number of modules $N$ by a factor of two at each interaction as long
%   as the modularity value $Q_b$ continues to increase, at which time the heuristic 
%   perform a bisection method between the last values $N$ that were visited.
% * |LPBrim|: The algorithm is a combination between the BRIM and
%   LP algorithms. The heuristic of this algorithms consist in searching for the
%   best module configuration by first using the LP algorithm. The BRIM algorithm
%   is used at the end to refine the results.
% * |LeadingEigenvector|: Altough this algorithm was defined for unipartite
%   networks, for certain bipartite networks, this algorithm can give the
%   best optimization.
%   equation, when defined in the unipartite version for two modules, can be expressed
%   as a linear combination of the normalized eigenvectors of the modularity matrix. Therefore,
%   this method divides the network in two modules by choosing the eigenvector with the
%   largest eigenvalue, such that negative and positive components of this eigenvector 
%   specifies the division. To implement this algorithm,
%   |BiMat| first converts the bipartite matrix into its unipartite version.


%% Example 1: Calculating Modularity
% The next example shows how to detect modularity using the default algorithm:
%
matrix = MatrixFunctions.BLOCK_MATRIX(2,10);
bp = Bipartite(matrix);
bp.plotter.PlotMatrix();
bp.community.Detect();
%%
% We can also change the default algorithm before detecting the community
% structure, choosing among one of the three options described before:
bp.community = LPBrim(matrix);
bp.community.Detect();
%%
% Further, there is no need to work directly with a |bipartite| instance. The user can
% also chose to work with a |BipartiteModularity| instance instead:
%By creating an instance and then calculating modularity
com = LeadingEigenvector(matrix);
com.Detect();
%Or by calling a static method:
com2 = BipartiteModularity.LEADING_EIGENVECTOR(matrix);

%% Example 2: Accesing detailed results
% Altough by calculating the modularity we can already see what are the
% modularity results, sometimes we may need to know detailed values. By
% just typing the name of the |BipartiteModularity| instance we have access
% to these values:
%A list of all properties
com2
%%
%Number of modules
com2.N
%%
%Module indices for row nodes
com2.row_modules'
%%
%Moudlarity value
com2.Qb

%% Example 3: Optimizing at the component level
% The default behavior of |BiMat| is to optimize modularity using the entire
% adjacency matrix. However, some times the network may have
% isolated components (sub-graphs with no connections between them). For
% those instances, the user may want to optimize each sub-graph independent
% of the others, and then calculate the modularity:
%Loading the data
load data/moebus_data.mat;
%default optimization
bp = Bipartite(moebus.weight_matrix > 0);
bp2.community.optimize_by_component = false;
%sub-graph optimization
bp2 = Bipartite(moebus.weight_matrix > 0);
bp2.community.optimize_by_component = true;
%visualizing the results:
plotFormat = PlotFormat();
plotFormat.back_color = [41,38,78]/255;
plotFormat.cell_color = 'white';
plotFormat.use_labels = false;
subplot(1,2,1);
bp.plotter.SetPlotFormat(plotFormat);
bp.plotter.PlotModularMatrix();
title(['$Q = ', num2str(bp.community.Qb),'$, $N = ',num2str(bp.community.N),'$'],...
    'FontSize',16, 'interpreter','latex');
xlabel('Default Optimization','FontSize', 20);
subplot(1,2,2);
bp2.plotter.SetPlotFormat(plotFormat);
bp2.plotter.PlotModularMatrix();
title(['$Q = ', num2str(bp2.community.Qb),'$, $N = ',num2str(bp2.community.N),'$'],...
    'FontSize',16, 'interpreter','latex');
xlabel('Optimizaing by isolated components','FontSize', 20);
set(gcf,'position', [70, 223, 1224, 673]);
%%
% We can see that the default optimization gives a higher modularity value
% (we are optimizing at the entire matrix) but lower number of nodes than
% the component optimization. By optimizing at the component level, we do
% not have any more the constraint of the entire network optimization and
% therefore it is possible to break into more modules. However, the
% reported value is still the value of the entire network partition.

