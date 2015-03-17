% StatisticalTest - Statistical analysis class for a bipartite object in terms of
% modularity and nestedness
%
% StatisticalTest Properties:
%     bipweb - A bipartite network object in which the analysis will be done
%     nulls - A set of random matrices used as null model
%     N_values - Results for the Nestedness value
%     Qb_values - Results for the standard modularity value
%     Qr_values - Results for the ratio of interactions inside modules
%     null_model - Null Model that will be used for creating the random networks nulls
%     replicates - Number of random networks used for the null model
%     community_done - The analysis in modularity was performed
%     nested_done      = 0;   % The analysis in nestedness was performed
%     print_status - Print status of the statistical test (print the random matrix that is tested)
%
% StatisticalTest Methods:
%     StatisticalTest - Main Constructor
%     DoNulls - Create random matrices for the statistical analysis
%     CleanNulls - Function to delete the nulls (random matrices)
%     DoCompleteAnalysis - Perform the entire modularity and nestedness analysis
%     TestNestedness - Perform the Nestedness Statistical Analysis
%     TestCommunityStructure - Perform the Modularity Statistical Analysis
%     GET_STATISTICAL_VALUES - Get the statistics of a value with respect to some random replicates
%     TEST_COMMUNITY_STRUCTURE - Perform a test of the modularity in a        % bipartite matrix
%     TEST_NESTEDNESS - Perform a test of the nestedness in a bipartite matrix
%
% See also:
%     MetaStatistics, NullModels, BipartiteModularity, Nestedness
classdef StatisticalTest < handle

    properties(GetAccess = 'public', SetAccess = 'protected')
        bipweb           = {};     % A bipartite network object in which the analysis will be done
        nulls            = {};     % A set of random matrices used as null model
        N_values         = [];     % Results for the Nestedness value
        Qb_values        = [];     % Results for the standard modularity value
        Qr_values        = [];     % Results for the ratio of interactions inside communities
        null_model       = Options.DEFAULT_NULL_MODEL; % Null Model that will be used for creating the random networks nulls
        replicates       = Options.REPLICATES;    % Number of random networks used for the null model
        community_done   = 0;      % The analysis in modularity was performed
        nested_done      = 0;   % The analysis in nestedness was performed
    end
    
    properties(Access = 'public')
        print_status     = 1;      % Print status of the statistical test (print the random matrix that is tested) 
    end
    
    methods

        function obj = StatisticalTest(webbip)
        % StatisticalTest - Main Constructor
        %   obj = StatisticalTest(webbip) Create a StatisticalTest object that makes
        %   reference to the Bipartite object webbip
        
            obj.bipweb = webbip;
        end
        

        function obj = DoNulls(obj,replic,nullmodel)
        % DoNulls - Create random matrices for the statistical analysis
        %   obj = DoNulls(obj) Create Options.REPLICATES random matrices using the
        %   default null model (Options.DEFAULT_NULL_MODEL).
        %
        %   obj = DoNulls(obj,replic) Create replic random matrices using
        %   the default null model (Options.DEFAULT_NULL_MODEL).
        %
        %   obj = DoNulls(obj,replic,nullmodel) Create replic random
        %   matrices using the null model indicated in the variable
        %   nullmodel
        %
        % See also:
        %   NullModels, Options.REPLICATES, Options.DEFAULT_NULL_MODEL
        
            %obj.community_done = 0;
            %obj.nested_done = 0;
            
            if(nargin == 1)
                obj.null_model = Options.DEFAULT_NULL_MODEL;
                obj.replicates = Options.REPLICATES;
            elseif(nargin == 2)
                obj.null_model = Options.DEFAULT_NULL_MODEL;
                obj.replicates = replic;
            else
                obj.null_model = nullmodel;
                obj.replicates = replic;
            end
            
            obj.nulls = NullModels.NULL_MODEL(obj.bipweb.matrix,obj.null_model,obj.replicates);
            
        end
        
        function obj = CleanNulls(obj)
        % CleanNulls(obj) - Function to delete the nulls (random matrices)
        % Useful when the user do not need to keep record of the created
        % random matrices and hence, avoid saturing memory.
            obj.nulls = {};
        end
        
        
        function obj = DoCompleteAnalysis(obj, replic, nullmodel)
        % DoCompleteAnalysis - Perform the entire modularity and nestedness analysis
        %
        %   obj = DoCompleteAnalysis(obj) Perform the entire analysis for
        %   nestedness and modularity using default values for null model
        %   and number of replicates (Options.REPLICATES,
        %   Options.DEFAULT_NULL_MODEL).
        %
        %   obj = DoCompleteAnalysis(obj,replic) Perform the entire analysis for
        %   nestedness and modularity using Options.DEFAULT_NULL_MODEL for
        %   creating replic random matrices
        %
        %   obj = DoCompleteAnalysis(obj,replic,nullmodel) Perform the entire analysis for
        %   nestedness and modularity using the the specified Null model
        %   and a total of replic random matrices
        %
        % See also:
        %   NullModels, Options.REPLICATES, Options.DEFAULT_NULL_MODEL
        
            
            if(nargin == 1)
                nullmodel = Options.DEFAULT_NULL_MODEL;
                replic = Options.REPLICATES;
            elseif(nargin == 2)
                nullmodel = Options.DEFAULT_NULL_MODEL;
            end
            
            fprintf('Creating %i null random matrices...\n', replic);
            obj.DoNulls(replic,nullmodel);
            fprintf('Performing Nestedness statistical analysis...\n');
            obj.TestNestedness();
            fprintf('Performing Modularity statistical analysis...\n');
            obj.TestCommunityStructure();
            
            
        end
                
        function obj = TestNestedness(obj)
        % TestNestedness - Perform the Nestedness Statistical Analysis
        %
        %   obj = TestNestedness(obj) Perform the Nestedness Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise Options.REPLICATES and Options.DEFAULT_NULL_MODEL
        %   will be used as number of rando matrices and null model.
        %
        % See also:
        %   StatisticalTest.DoNulls
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            
            n = length(obj.nulls);
            
            nestedness = obj.bipweb.nestedness;
            
            if(nestedness.done == 0)
                tmp = nestedness.print_results;
                nestedness.print_results = false;
                nestedness.Detect();
                nestedness.print_results = tmp;
                N_val = obj.bipweb.nestedness.N;
            else
                N_val = obj.bipweb.nestedness.N;
            end
            
            expect = zeros(n,1);
            nested_class = str2func(class(obj.bipweb.nestedness));
            
            print_stat = floor(linspace(1,n+1,11));
            jp = 1;
            
            for i = 1:n
                if(obj.print_status && print_stat(jp)==i)
                    fprintf('Evaluating Nestedness in random matrices: %i - %i\n', print_stat(jp),print_stat(jp+1)-1);
                    jp = jp+1;
                end
                nested = feval(nested_class,obj.nulls{i});
                nested.print_results = false;
                nested.Detect();
                expect(i) = nested.N;
            end
                       
            obj.N_values = StatisticalTest.GET_STATISTICAL_VALUES(N_val,expect);
            obj.N_values.algorithm = class(obj.bipweb.nestedness);
            obj.N_values.replicates = obj.replicates;
            obj.N_values.null_model = obj.null_model;
                        
            obj.nested_done = 1;
            
        end
        
        
        function obj = TestCommunityStructure(obj)
        % TestCommunityStructure - Perform the Modularity Statistical Analysis
        %
        %   obj = TestCommunityStructure(obj) Perform the Modularity Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise Options.REPLICATES and Options.DEFAULT_NULL_MODEL
        %   will be used as number of rando matrices and null model.
        %
        % See also:
        %   StatisticalTest.DoNulls
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            
            %Calculate the modularity of the Bipartite object
            if(obj.bipweb.community.done == 0)
                tmp = obj.bipweb.community.print_results;
                obj.bipweb.community.print_results = false;
                obj.bipweb.community.Detect(100);
                obj.bipweb.community.print_results = tmp;
            end
            
            wQr = obj.bipweb.community.Qr;
            wQb = obj.bipweb.community.Qb;
            n = obj.replicates;
            
            Qb_random = zeros(n,1);
            Qr_random = zeros(n,1);
            
            if(isprop(obj.bipweb.community,'DoKernighanLinTunning'))
                exist_kernig = 1;
                do_kernig = obj.bipweb.community.DoKernighanLinTunning;
            else
                exist_kernig = 0;
                do_kernig = 0;
            end
            
            modul_class = str2func(class(obj.bipweb.community));
            n_trials = obj.bipweb.community.trials;
            
            print_stat = floor(linspace(1,n+1,11));
            jp = 1;
            
            for i = 1:n
                if(obj.print_status && print_stat(jp)==i)
                    fprintf('Evaluating Modularity in random matrices: %i - %i\n', print_stat(jp),print_stat(jp+1)-1);
                    jp = jp+1;
                end
                modularity = feval(modul_class,obj.nulls{i});
                modularity.trials = n_trials;
                if(exist_kernig == 1)
                    modularity.DoKernighanLinTunning = do_kernig;
                end
                modularity.print_results = false;
                modularity.Detect(10);
                Qb_random(i) = modularity.Qb;
                Qr_random(i) = modularity.Qr; 
            end
            
            obj.Qb_values = StatisticalTest.GET_STATISTICAL_VALUES(wQb,Qb_random);
            obj.Qr_values = StatisticalTest.GET_STATISTICAL_VALUES(wQr,Qr_random);
            
            obj.Qb_values.algorithm = class(obj.bipweb.community);
            obj.Qb_values.replicates = obj.replicates;
            obj.Qb_values.null_model = obj.null_model;
            
            obj.Qr_values.algorithm = class(obj.bipweb.community);
            obj.Qr_values.replicates = obj.replicates;
            obj.Qr_values.null_model = obj.null_model;
            
            obj.community_done = 1;
        end
        
        
        
        function obj = Print(obj,filename)
        % Print - Print all previouslly statistically tested values
        %
        %   STR = Print(obj) Print the statistical results to screen and
        %   returns this information to the string STR
        %
        %   STR = Print(obj, FILE) Print the statistical resultss to screen and
        %   text file FILE and return this information to the string STR   
        %
        % See also: 
        %   Printer    
            str = '';
            if(obj.community_done == 1)
                str = 'Modularity\n';
                str = [str, '\t Used algorithm:\t', sprintf('%30s',obj.Qb_values.algorithm), '\n'];
                str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.Qb_values.null_model)), '\n'];
                str = [str, '\t Replicates:    \t', sprintf('%30i',obj.Qb_values.replicates), '\n'];
                str = [str, '\t Qb value:      \t', sprintf('%30.4f',obj.Qb_values.value), '\n'];
                str = [str, '\t     mean:      \t', sprintf('%30.4f',obj.Qb_values.mean), '\n'];
                str = [str, '\t     std:       \t', sprintf('%30.4f',obj.Qb_values.std), '\n'];
                str = [str, '\t     z-score:   \t', sprintf('%30.4f',obj.Qb_values.zscore), '\n'];
                str = [str, '\t     t-score:   \t', sprintf('%30.4f',obj.Qb_values.tscore), '\n'];
                str = [str, '\t     percentil: \t', sprintf('%30.4f',obj.Qb_values.percentile), '\n'];
                
                str = [str, '\t Qr value:      \t', sprintf('%30.4f',obj.Qr_values.value), '\n'];
                str = [str, '\t     mean:      \t', sprintf('%30.4f',obj.Qr_values.mean), '\n'];
                str = [str, '\t     std:       \t', sprintf('%30.4f',obj.Qr_values.std), '\n'];
                str = [str, '\t     z-score:   \t', sprintf('%30.4f',obj.Qr_values.zscore), '\n'];
                str = [str, '\t     t-score:   \t', sprintf('%30.4f',obj.Qr_values.tscore), '\n'];
                str = [str, '\t     percentil: \t', sprintf('%30.4f',obj.Qr_values.percentile), '\n'];
            end
            
            if(obj.nested_done == 1)
                str = [str, 'Nestedness\n'];
                str = [str, '\t Used algorithm:\t', sprintf('%30s',obj.N_values.algorithm), '\n'];
                str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.N_values.null_model)), '\n'];
                str = [str, '\t Replicates:    \t', sprintf('%30i',obj.N_values.replicates), '\n'];
                str = [str, '\t Nestedness value:\t', sprintf('%30.4f',obj.N_values.value), '\n'];
                str = [str, '\t     mean:      \t', sprintf('%30.4f',obj.N_values.mean), '\n'];
                str = [str, '\t     std:       \t', sprintf('%30.4f',obj.N_values.std), '\n'];
                str = [str, '\t     z-score:   \t', sprintf('%30.4f',obj.N_values.zscore), '\n'];
                str = [str, '\t     t-score:   \t', sprintf('%30.4f',obj.N_values.tscore), '\n'];
                str = [str, '\t     percentil: \t', sprintf('%30.4f',obj.N_values.percentile), '\n'];
            end
            

            fprintf(str);  
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
    end
    


    methods(Static)
        
        function out = GET_STATISTICAL_VALUES(real_value, random_values)
        % GET_STATISTICAL_VALUES - Get the statistics of a value with
        % respect to some random replicates
        %
        %   out = GET_STATISTICAL_VALUES(real_value, random_values) Get
        %   some statistics for real_value using the vector random_values. It returns
        %   an structure out with the folloging variables:
        %
        %     value: Just real_value
        %     mean: mean(random_values)
        %     std: std(random_values)
        %     zscore: the z-score of real_value in the distribution of
        %             random_values
        %     percentile: the percentile of real_value in the distribution
        %                 of random_values
        %     random_values: the sorted version of random_values
        %     
            n = length(random_values);
            me = mean(random_values);
            stad = std(random_values);
            zscore = (real_value - me)/stad;
            tscore = (real_value - me)/(stad/sqrt(n));
            
            percent = 100.0*sum(real_value>random_values+0.00001)/n;
            
            out.value = real_value; out.mean = me; out.std = stad;
            out.zscore = zscore; out.percentile = percent;
            out.tscore = tscore;
            out.random_values = sort(random_values);
            
        end
        
        function stest = TEST_COMMUNITY_STRUCTURE(matrix,replicates,null_model,modul_algorithm)
        % TEST_COMMUNITY_STRUCTURE - Perform a test of the modularity in a
        % bipartite matrix
        %
        %   stest = TEST_COMMUNITY_STRUCTURE(MATRIX) Perform a statistical
        %   analysis of modularity in MATRIX using default values for null
        %   model, modularity algorithm, and replicates. It returns a
        %   StatisticalTest object.
        %
        %   obj = TEST_COMMUNITY_STRUCTURE(MATRIX,REPLICATES) Perform a statistical
        %   analysis of modularity in MATRIX using REPLICATES random matrices
        %
        %   obj = TEST_COMMUNITY_STRUCTURE(MATRIX,REPLICATES,NULL_MODEL)
        %   Perform an statistical analysis of modularity in MATRIX unsing
        %   NULL_MODEL as null model.
        %   
        %   stest = TEST_COMMUNITY_STRUCTURE(MATRIX,REPLICATES,NULL_MODEL,MODUL_ALGORITHM)
        %   Perform an statistical analysis of modularity in MATRIX unsing
        %   MODUL_ALGORITHM as modularity algorithm.
        %
        % See also:
        %   NullModels, Options.REPLICATES, Options.DEFAULT_NULL_MODEL,
        %   Options.MODULARITY_ALGORITHM
        
            
            if(nargin==1)
                replicates = Options.REPLICATES;
                null_model = Options.DEFAULT_NULL_MODEL;
                modul_algorithm = Options.MODULARITY_ALGORITHM;
            elseif(nargin==2)
                null_model = Options.DEFAULT_NULL_MODEL;
                modul_algorithm = Options.MODULARITY_ALGORITHM;
            elseif(nargin==3)
                modul_algorithm = Options.MODULARITY_ALGORITHM;
            end
            
            bp = Bipartite(matrix);
            bp.community = modul_algorithm(bp.matrix);
            stest = StatisticalTest(bp);
            stest.DoNulls(replicates,null_model);
            stest.TestCommunityStructure();
            
            stest.Print();
        end
        
        function stest = TEST_NESTEDNESS(matrix,replicates,null_model,nested_algorithm)
        % TEST_NESTEDNESS - Perform a test of the nestedness in a
        % bipartite matrix
        %
        %   stest = TEST_NESTEDNESS(MATRIX) Perform a statistical
        %   analysis of Nestedness in MATRIX using default values for null
        %   model, nested algorithm, and replicates. It returns a
        %   StatisticalTest object.
        %
        %   stest = TEST_NESTEDNESS(MATRIX,REPLICATES) Perform a statistical
        %   analysis of Nestedness in MATRIX using REPLICATES random matrices
        %
        %   stest = TEST_NESTEDNESS(MATRIX,REPLICATES,NULL_MODEL)
        %   Perform an statistical analysis of Nestedness NODF in MATRIX unsing
        %   NULL_MODEL as null model.
        %
        %   stest = TEST_NESTEDNESS(MATRIX,REPLICATES,NULL_MODEL,NESTED_ALGORITHM)
        %   Perform an statistical analysis of modularity in MATRIX unsing
        %   MODUL_ALGORITHM as modularity algorithm.
        %
        % See also:
        %   NullModels, Options.REPLICATES, Options.DEFAULT_NULL_MODEL    
            if(nargin==1)
                replicates = Options.REPLICATES;
                null_model = Options.DEFAULT_NULL_MODEL;
                nested_algorithm = Options.NESTEDNESS_ALGORITHM;
            elseif(nargin==2)
                null_model = Options.DEFAULT_NULL_MODEL;
                nested_algorithm = Options.NESTEDNESS_ALGORITHM;
            elseif(nargin==3)
                nested_algorithm = Options.NESTEDNESS_ALGORITHM;
            end
            
            bp = Bipartite(matrix);
            bp.nestedness = nested_algorithm(matrix);
            stest = StatisticalTest(bp);
            stest.DoNulls(replicates,null_model);
            stest.TestNestedness();
            
            stest.Print();
            
        end
        
    end

end

%   NEXT COMMENTED CODE INCLUDES FUNCTIONS TO CALCULATE THE NODF
%   CONTRIBUTION OF BOTH ROWS AND COLUMNS, AND THE MAXIUM EIGENVALUE NEW
%   APPROACH TO EVALUATE NESTEDNESS IN A BIPARTITE MATRIX. THEY WERE
%   COMMENTED BECAUSE NO EXAHUSTIVE TEST WAS PERFORMED IN THESE FUNCTIONS.
%     properties(GetAccess = 'private', SetAccess = 'protected')
%     %Untested algorithms
%         nodf_row_contrib = []; % Row NODF contributions
%         nodf_col_contrib = []; % Column Nodf Contributions
%         eigvals          = [];    % Results for the espectral Radius algorithm value
%         nest_contrib_done= 0; % The analysis on nestedness contribution was performed
%         eig_done         = 0;      % The analysis en expectral radius was performed
%     end
%     
%     methods(Access = 'private')
%         
%         function obj = MaxEigenvalue(obj)
%         %NOT TESTED!!!!!
%         % MaxEigenvalue - Perform the spectral radius Statistical Analysis
%         %
%         %   obj = MaxEigenvalue(obj) Perform the NTC Statistical analsysis. Be
%         %   sure to create the random matrices before calling this
%         %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
%         %   be used for the analysis
%             
%             if(isempty(obj.nulls))
%                 obj.DoNulls();
%             end
%             
%             [obj.eigvals] = StatisticalTest.GET_DEV_EIG(obj.bipweb,obj.nulls);
%             
%             obj.eig_done = 1;
%         end
%         
%     end
%     

% 
%     methods
%         
%                
%         
%        function obj = NestednessContributions(obj)
%         % NOT TESTED!!!!!!
%         % NestednessContributions - Perform the nestedness contribution Statistical Analysis
%         %
%         %   obj = NestednessContributions(obj) Perform the nestedness contribution Statistical analsysis. Be
%         %   sure to create the random matrices before calling this
%         %   function. Otherwise default number of replicates and null model
%         %   will be used (see main/Options.m)
%         
%             if(isempty(obj.nulls))
%                 obj.DoNulls();
%             end
%             [obj.nodf_row_contrib,obj.nodf_col_contrib] = StatisticalTest.GET_NEST_CONTRIBUTIONS(obj.bipweb.matrix,obj.nulls);
%             
%             obj.nodf_row_contrib.replicates = obj.replicates;
%             obj.nodf_row_contrib.null_model = obj.null_model;
%             
%             obj.nodf_col_contrib.replicates = obj.replicates;
%             obj.nodf_col_contrib.null_model = obj.null_model;
%             
%             obj.nest_contrib_done = 1;
%         end
%         
%         function [out] = GET_DEV_EIG(webbip,rmatrices)
%         % NOT TESTED!!!!!!
%         % GET_DEV_EIG - Perform a spectral radius Statistical Analysis
%         %
%         %   [out out_row out_col] =  GET_DEV_EIG(webbip,rmatrices) Perform t-test and
%         %   z-test in the NTC value of the bipartite object webbip using
%         %   the ser of random matrices rmatrices. Return a
%         %   structure for NTC statistical values in the entire matrix(out)
%         %   with the next elements:
%         %      value   - The value in the empirical matrix
%         %      p       - p-value of the performed t-test
%         %      ci      - Confidence interval of the performet t-test
%         %      percent - The percent of random networks for which  the
%         %                empirical value is bigger than the value of the random
%         %                networks
%         %      z       - z-score of the empirical value      
%             n = length(rmatrices);
%             expect = zeros(n,1);
%             
%             maxe = MatrixNull.GetBiggestEigenvalue(webbip.matrix);
%             
%             for i = 1:n
%                 expect(i) = MatrixNull.GetBiggestEigenvalue(rmatrices{i});
%             end
%             
%             out = StatisticalTest.GET_STATISTICAL_VALUES(maxe,expect);
%             
%         end
%         
%         function [c_rows, c_cols] = GET_NEST_CONTRIBUTIONS(matrix, rmatrices)
%         % NOT TESTED!!!!!!
%         % GET_NEST_CONTRIBUTIONS - Get NODF nestedness contributions of a
%         % matrix
%         %
%         %   [c_rows c_cols] = GET_NEST_CONTRIBUTIONS(matrix, rmatrices) Get
%         %   the NODF contributions for both rows and columns. The
%         %   contribution is calculated as the z-score of each individually
%         %   randomly permuted row and columns
%             nodf = NODF(matrix);
%             N = nodf.nodf;
%             [n_rows, n_cols] = size(matrix);
%             nn = length(rmatrices);
%             row_contrib = zeros(nn,n_rows);
%             col_contrib = zeros(nn,n_cols);
%             %orig_matrix = matrix;
%             for i = 1:n_rows
%                 temp_row = matrix(i,:);
%                 for j = 1:nn
%                     matrix(i,:) = rmatrices{j}(i,:);
%                     nodf = NODF(matrix,0);
%                     row_contrib(j,i) = nodf.nodf;
%                     if(isnan(nodf.nodf))
%                         continue;
%                     end
%                 end
%                 matrix(i,:) = temp_row;
%             end
%             
%             for i = 1:n_cols
%                 temp_col = matrix(:,i);
%                 for j = 1:nn
%                     matrix(:,i) = rmatrices{j}(:,i);
%                     nodf = NODF(matrix,0);
%                     col_contrib(j,i) = nodf.nodf;
%                 end
%                 matrix(:,i) = temp_col;
%             end
%             std_rows = std(row_contrib);
%             std_cols = std(col_contrib);
%             mean_rows = mean(row_contrib);
%             mean_cols = mean(col_contrib);
%             
%             c_rows = (N - mean_rows)./std_rows;
%             c_cols = (N - mean_cols)./std_cols;
%             
%         end
%         
%         
%     end