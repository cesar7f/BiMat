%% Modularity

%% Description
% Modularity indicates the presence of dense clusters of related nodes embedded 
% within the network. In many systems, we can 
% find a partition of nodes into specific communities or modules.  
% This modularity metric can be expressed for a bipartite network as:
% <latex>
% \begin{equation}
% Q_b = \frac{1}{E} \sum_{ij} \left( B_{ij} - \frac{k_i d_j}{E} \right) \delta(g_i,h_j),
% \end{equation}
% </latex>
% where @B_ij@ is the element in the bipartite matrix representing a link
% (1)
% or no link (0) between nodes @i@ and @j@, @g_i@ and @h_i@ are the module indexes of nodes @i@ (that belongs to set
% @R@) and @j@ (that belongs to set @C@), @k_i@ is the degree of node @i@, @d_j@ is
% the degree of node @j@ and @E@ is the number of links in the network.  
% |BiMat| contains three algorithms to maximize the last equation, all of
% them containing different heuristics. 
%
% * |AdaptiveBrim|: The algorithm uses the BRIM algorithm
%   with an heuristic to look for the optimal value of modules. This heuristic consists
%   in multiply the number of modules @N@ by a factor of two at each interaction as long
%   as the modularity value @Q_b@ continues to increase, at which time the heuristic 
%   perform a bisection method between the last values @N@ that were visited.
% * |LPBrim|: The algorithm is a combination between the BRIM and
%   LP algorithms. The heuristic of this algorithm consists in searching for the
%   best module configuration by first using the LP algorithm. The BRIM algorithm
%   is used at the end to refine the results.
% * |LeadingEigenvector|: Although this algorithm was defined for unipartite
%   networks, for certain bipartite networks, this algorithm can give the
%   highest modularity value.
%   The bipartite modularity equation, when defined in the unipartite version for two modules, can be expressed
%   as a linear combination of the normalized eigenvectors of the modularity matrix. Therefore,
%   this method divides the network in two modules by choosing the eigenvector with the
%   largest eigenvalue, such that negative and positive components of this eigenvector 
%   specifies the division. To implement this algorithm,
%   |BiMat| first converts the bipartite matrix into its unipartite version
%   and the initial step.
%
% In addition to optimize the standard modularity @Q_b@ \bimat also evaluates
% (after optimizing @Q_b@) an a posteriori measure of modularity @Q_r@
% introduced in <http://f1000research.com/articles/2-130/v3 Poisot 2013> and defined as:
% <latex>
% \begin{equation}
% Q_r = 2\times\frac{W}{E}-1 %% This version returns values in 0;1
% \end{equation}
% </latex>
% where @W = \sum_{ij} B_{ij} \delta(g_i,h_j)@ is the number of edges that
% are inside modules. Alternatively, @Q_r \equiv \frac{W-T}{W+T}@ where
% @T@ is the number of edges that are between modules.  In other words, 
% this quantity maps the relative difference of edges that are
% within modules to those between modules on a scale 
% from 1 (all edges are within modules) to @-1@ (all edges are between modules). 
% This measure allows to compare the output of different algorithms.

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
% The description of these properties is showed in the next table:
%
% <html>
% <table class="tftable" border="1">
% <tr><th>Property</th><th>Algorithm</th><th>Description</th></tr>
% <tr> <td><tt>Qb</tt></td> <td>All</td> <td>standard modularity value between 0 and 1</td> </tr>
% <tr> <td><tt>Qr</tt></td> <td>All</td> <td>a-posteriori modularity value between 0 and 1</td> </tr>
% <tr> <td><tt>N</tt></td> <td>All</td> <td>number of modules</td> </tr>
% <tr> <td><tt>matrix</tt></td> <td>All</td> <td>bipartite adjacency matrix</td> </tr>
% <tr> <td><tt>row_modules</tt></td> <td>All</td> <td>module index for each row</td> </tr>
% <tr> <td><tt>col_modules</tt></td> <td>All</td> <td>module index for each column</td> </tr>
% <tr> <td><tt>bb</tt></td> <td>All</td> <td>modularity matrix</td> </tr>
% <tr> <td><tt>n_rows</tt></td> <td>All</td> <td>number of rows</td> </tr>
% <tr> <td><tt>n_cols</tt></td> <td>All</td> <td>number of columns</td> </tr>
% <tr> <td><tt>n_edges</tt></td> <td>All</td> <td>number of edges</td> </tr>
% <tr> <td><tt>index_rows</tt></td> <td>All</td> <td>row indexes after sorting by modules</td> </tr>
% <tr> <td><tt>index_cols</tt></td> <td>All</td> <td>column indexes after sorting by modules</td> </tr>
% <tr> <td><tt>done</tt></td> <td>All</td> <td>if modularity has been already detected</td> </tr>
% <tr> <td><tt>optimize_by_component</tt></td> <td>All</td> <td>should modularity be optimized by isolated component first? (default <tt>false</tt>)</td> </tr>
% <tr> <td><tt>print_results</tt></td> <td>All</td> <td>if results will be printed after detection (default <tt>true</tt>)</td> </tr>
% <tr> <td><tt>trials</tt></td> <td><tt>AdaptiveBrim/LPBrim</tt></td> <td>number of initialization trials (default see <tt>Options.TRIALS_MODULARITY)</tt></td> </tr>
% <tr> <td><tt>prob_reassigment</tt></td> <td><tt>AdaptiveBrim</tt></td> <td>probability of reassigment when new modules are created (default 0.5)</td> </tr>
% <tr> <td><tt>expansion_factor</tt></td> <td><tt>AdaptiveBrim</tt></td> <td>factor of expansion for creating new modules (default 2.0)</td> </tr>
% <tr> <td><tt>red_labels</tt></td> <td><tt>LPBrim</tt></td> <td>Similar than <tt>row_modules</tt></td> </tr>
% <tr> <td><tt>col_labels</tt></td> <td><tt>LPBrim</tt></td> <td>Similar than <tt>col_modules</tt></td> </tr>
% <tr> <td><tt>DoKernighanLinTunning</tt></td> <td><tt>LeadingEigenvector</tt></td> <td>execute Kernighan tunning after modularity detection? (default <tt>true</tt>)</tt></td> </tr>
% </table>
% </html>
%
% Now, we can access a property by just typing it:
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
xlabel('Optimizing by isolated components','FontSize', 20);
set(gcf,'position', [70, 223, 1224, 673]);
%%
% We can see that the default optimization gives a higher modularity value
% (we are optimizing at the entire matrix) but lower number of modules than
% the component optimization. By optimizing at the component level, we do
% not have any more the constraint of the entire network optimization and
% therefore it is possible to break into more modules. However, the
% reported value is still the value of the entire network partition.

