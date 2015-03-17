classdef LeadingEigenvector < BipartiteModularity
% LeadingEigenvector - Main code class
% Newman's leading eigenvector algorithm to calculate modularity in unipartite
% networks. In this case, the bipartite network is converted to
% unipartite and then, the Newman algorithm is applied in the
% unipartite adjacency matrix of the network. The algorithm is
% explained in the paper:
%
%    Newman, Mark EJ. Modularity and community structure in networks.
%    PNAS 2006
%
% LeadingEigenvector Properties:
%    DoKernighanLinTunning - Do a final tuning in the modularity configuration to improve the results.
%
% LeadingEigenvector Methods:
%    LeadingEigenvector - Main constructor    
%
% See also:
%    BipartiteModularity, AdaptiveBrim, and LPBrim
    
    properties
        DoKernighanLinTunning= true;   %Do a final tuning in the modularity configuration to improve the results.
    end
    
    properties(Access = 'protected')
        adjacency            = [];  %Unipartite version of the adjacency matrix
        n_nodes              = 0;   %n_cols + n_rows
        pp                   = [];  %Null model matrix
        kk                   = 0;   %Node Degrees
        ss                   = 0;   %Community indices of the unipartite version
        Q                    = 1;   %Modularity value in uniparte version.
    end
    
    methods
       
        function obj = LeadingEigenvector(bipmatrix)
        % LeadingEigenvector - Main Constructor
        % 
        %   obj = LeadingEigenvector(MATRIX) Creates an LeadingEigenvector object obj
        %   using a bipartite adjacency matrix MATRIX that will be used to
        %   calculate modularity using the Newman's Leading Eigenvector
        %   algorithm
        %
        % See also:
        %   LeadingEigenvector
            
            %Call the parent class
            obj = obj@BipartiteModularity(bipmatrix);
            
        end
        
    end
        
    methods(Access = 'protected')
        
        function obj = DetectComponent(obj)
        % DetectComponent - Main method of the algorithm
        %
        %   obj = DetectComponent(obj) Detect the modularity in a specific
        %   component.
        
            %Convert the bipartite matrix component to unipartite version.
            obj.adjacency = [zeros(obj.n_rows_component) obj.matrix_component; obj.matrix_component' zeros(obj.n_cols_component)];
            %obj.adjacency = bipmatrix > 0; %Uncomment if you want to work directly with unipartite networks
            
            %Get the number of total nodes, degrees and number of edges
            obj.n_nodes = size(obj.adjacency,1);
            obj.kk = sum(obj.adjacency,2);
            obj.n_edges_component =  sum(obj.kk)/2;
            
            %If there is no interaction, the modular matrix is the same
            %than the adjacency matrix
            if all(all(obj.matrix_component == 0))
                obj.bb_component = obj.adjacency;
            %otherwise calculate the modular matrix using the unipartite
            %version equation
            else
                %obj.bb_component = obj.adjacency - (obj.kk*obj.kk')/(2*obj.n_edges_component);
                if(obj.optimize_by_component == false)
                    obj.bb_component = obj.adjacency - (obj.kk*obj.kk')/(2*obj.n_edges);
                else
                    obj.bb_component = obj.adjacency - (obj.kk*obj.kk')/(2*obj.n_edges_component);
                end
            end
            
            %Run the newman algorithm.
            obj.NewmanAlgorithm();
            ncom = length(unique(obj.ss));
            
            obj.rr_component = zeros(obj.n_rows_component,ncom);
            obj.tt_component = zeros(obj.n_cols_component,ncom);
            
            for i = 1:obj.n_nodes
                if(i <= obj.n_rows_component); obj.rr_component(i,obj.ss(i)) = 1;
                else obj.tt_component(i-obj.n_rows_component,obj.ss(i)) = 1; end;
            end
            
        end
        
        function obj = NewmanAlgorithm(obj)
        % obj = NewmanAlgorithm(obj) - Detect the modularity using Newman
        % algorithm.
            
            %Start with a single component.
            obj.N_component = 1;
            obj.Q = 0;
            
            %Community indexes of all nodes correspond to module 1.
            obj.ss = ones(obj.n_nodes,1);
            
            %Try to divide all the existing modules while increase in Q is
            %detected.
            qinc = 1;
            while(qinc == 1)
                   
                qinc = 0;
                nn = obj.N_component;
                
                %For each of the current modules, divide them only if
                %increase in Q is detected.
                for i = 1:nn
                    
                    %Get from the modular matrix the nodes that correspond
                    %to module i and find a new local modular matrix.
                    nindices = find(obj.ss==i);
                    nsize = length(nindices);
                    blocal = obj.bb_component(nindices,nindices);
                    blocal(1:nsize+1:end) = blocal(1:nsize+1:end)' - sum(blocal,2);
                    
                    %Use the biggest eigenvalue vector of the local matrix
                    %to divide the module in two.
                    [evec eval] = eig(blocal);
                    eval = diag(eval);
                    [~, maxindex] = max(eval);
                    evecmax = evec(:,maxindex);
                    
                    sslocal = ones(nsize,1);
                    sslocal(evecmax<=0) = -1;
                    
                    %Do a final tuning to improve the modularity
                    %configuration.
                    if(obj.DoKernighanLinTunning)
                        sslocal = obj.KernighanLin(sslocal,blocal);
                    end
                    
                    %Calculate the improvement in Q if module i is divided
                    %in two.
                    deltaQ = (sslocal' * blocal * sslocal) / (4*obj.n_edges_component);
                    
                    %If improvement is detected, divide the module in two,
                    %otherwise reject the division.
                    if(deltaQ > 0.000001)
                        qinc = 1;
                        newind = sslocal==1;
                        obj.ss(nindices(newind)) = obj.N_component+1;
                        obj.N_component = obj.N_component+1;
                        obj.Q = deltaQ + obj.Q;
                    end
                    
                end
                
            end
            
        end
        
        function ssnew = KernighanLin(obj,ss,blocal)
        % obj = KernighanLin(obj,ss,blocal) - Improve the modularity
        % configuration by applying KernighanLin algorithm. The idea of
        % this algorithm is just swapping always one node from one module
        % to the other such that it will be the best increase in
        % modularity. At each step the module which swapping gives the
        % biggest increase in modularity is swapped and the process
        % continue as long as modularity increase.
        
            %Calculate the modularity of a two module local network.
            nn = length(ss);
            ss_states = cell(nn+1,1);
            num_edges = obj.n_edges;
            %qglobal = obj.ModularityTwo(ss,blocal);
            
            qglobal = modul_two(ss,blocal,num_edges);
            
            qinc = 1;
            
            %Repeat while increase in modularity Q is detected.
            while(qinc == 1)
                
                qinc = 0;
                ss_states{1} = ss;
                q_values = zeros(nn,1);
                moved = zeros(nn,1);
                
                %Move all the nodes at least one time
                for i = 1:nn   

                    sslocal = ss_states{i};
                    ssbb = sslocal'*blocal;
                    qloc = -10;
                    %Move the node that has not been moved and has the biggest Q increase
                    for j = 1:nn

                        if(moved(j) == 1); continue; end;

                        sslocal(j)  = sslocal(j) * -1;
                        %newq = (sslocal'*blocal*sslocal) / (4*num_edges);
                        newq = (ssbb + 2*sslocal(j)*blocal(j,:))*sslocal / (4*num_edges);

                        if(newq > qloc)
                            qindex = j;
                            qloc = newq;
                        end
                        sslocal(j)  = sslocal(j) * -1;

                    end
                    sslocal(qindex) = sslocal(qindex) * -1;
                    ss_states{i+1} = sslocal;
                    moved(qindex) = 1;
                    q_values(i) = qloc; 
                end
                
                %Find the intermidiate state with the biggest Q
                [max_q max_state] = max(q_values);
                
                %Check if increase in global Q is detected
                if(max_q > qglobal)
                    qglobal = max_q;
                    ss = ss_states{max_state+1};
                    qinc = 1;
                end
            end
            
            ssnew = ss;
            
            function qq = modul_two(ss,bb,num_edges)
            %Nested function (faster than callin an outsider function)
                qq = (ss'*bb*ss) / (4*num_edges);
            end
        end
        
        
        function qq = Modularity(obj,ss,bb)
        %qq = Modularity(obj,ss,bb) Get the modularity of a unipartite
        %network using ss as module indesing and bb as modular matrix.
            ncom = length(unique(ss));
            nn = length(bb);
            S = zeros(nn,ncom);
            for i = 1:nn
                S(i,ss(i)) = 1;
            end
            %qq = trace(S'*bb*S) / (2*obj.n_edges_component);
            qq = trace(S'*bb*S) / (2*obj.n_edges);
        end
        
        
    end
    
    
end