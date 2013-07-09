classdef NewmanModularity < handle

    properties
        matrix               = [];  %Bipartite adjacency matrix
        adjacency            = [];  %Unipartite version of the adjacency matrix
        %rr                   = [];  %Red Nodes (rows) Communities matrix. Size = n_rows*CommunityQuantity
        %tt                   = [];  %Blue Nodes (columns) Communities. Size = n_cols*CommunityQuantity
        n_rows               = 0;   %Number of rows
        n_cols               = 0;   %Number of columns
        n_edges              = 0;
        n_nodes              = 0;   %n_cols + n_rows
        pp                   = [];  %Null model matrix
        bb                   = [];  %Original - Null.
        index_rows           = [];  %Register of the swaps in Rows.
        index_cols           = [];  %Register of the swaps in Cols.
        red_labels           = 0;
        blue_labels          = 0;
        trials               = 10;
        Q                    = 0;
        Qb                   = 0;
        Qr                   = 0;
        N                    = 0;
        row_modules          = [];
        col_modules          = [];
        kk                   = 0;   %Node Degrees
        done                 = 0;
        ss                   = 0;   %Community indices of the unipartite version
        rr                   = 0;
        tt                   = 0;
        DoKernighanLinTunning= 1;
    end
    
    methods
       
        function obj = NewmanModularity(bipmatrix)
            
            obj.matrix = bipmatrix > 0;
            [obj.n_rows obj.n_cols] = size(obj.matrix);
            
            obj.adjacency = [zeros(obj.n_rows) obj.matrix; obj.matrix' zeros(obj.n_cols)];
            %obj.adjacency = bipmatrix > 0; %Uncomment if you want to work directly with unipartite networks
            obj.n_nodes = size(obj.adjacency,1);
            
            obj.trials = 10;
            
            obj.kk = sum(obj.adjacency,2);
            obj.n_edges = sum(obj.kk)/2;
            
            if all(all(obj.matrix == 0))
                obj.bb = obj.adjacency;
            else
                obj.bb = obj.adjacency - (obj.kk*obj.kk')/(2*obj.n_edges);
            end
            
        end
        
        function obj = Detect(obj)
            
            obj.NewmanAlgorithm();
            ncom = length(unique(obj.ss));
            obj.rr = zeros(obj.n_rows,ncom);
            obj.tt = zeros(obj.n_cols,ncom);
            
            for i = 1:obj.n_nodes
                if(i <= obj.n_rows); obj.rr(i,obj.ss(i)) = 1;
                else obj.tt(i-obj.n_rows,obj.ss(i)) = 1; end;
            end
            coldeg = sum(obj.matrix, 1);
            rowdeg = sum(obj.matrix, 2);
            nedges = sum(rowdeg);
            bbbip = obj.matrix - (1/nedges) * rowdeg * coldeg;
            obj.Qb = trace(obj.rr' * bbbip * obj.tt) / nedges;
            obj.CalculateQrValue(nedges);
        end
        
        function obj = NewmanAlgorithm(obj)
            
            obj.N = 1;
            obj.Q = 0;
            
            obj.ss = ones(obj.n_nodes,1);
            
            qinc = 1;
            while(qinc == 1)
                   
                qinc = 0;
                nn = obj.N;
                
                for i = 1:nn
                    
                    nindices = find(obj.ss==i);
                    nsize = length(nindices);
                    blocal = obj.bb(nindices,nindices);
                    blocal(1:nsize+1:end) = blocal(1:nsize+1:end)' - sum(blocal,2);
                    
                    [evec eval] = eig(blocal);
                    eval = diag(eval);
                    [meval maxindex] = max(eval);
                    evecmax = evec(:,maxindex);
                    
                    sslocal = ones(nsize,1);
                    sslocal(evecmax<=0) = -1;
                    
                    if(obj.DoKernighanLinTunning)
                        sslocal = obj.KernighanLin(sslocal,blocal);
                    end
                    
                    deltaQ = (sslocal' * blocal * sslocal) / (4*obj.n_edges);
                    
                    if(deltaQ > 0.000001)
                        qinc = 1;
                        newind = sslocal==1;
                        obj.ss(nindices(newind)) = obj.N+1;
                        obj.N = obj.N+1;
                        obj.Q = deltaQ + obj.Q;
                    end
                    
                end
                
            end
            
        end
        
        function ssnew = KernighanLin(obj,ss,blocal)
            
            nn = length(ss);
            ss_states = cell(nn+1,1);
            qglobal = obj.ModularityTwo(ss,blocal);
            
            qinc = 1;
            
            while(qinc == 1)
                
                qinc = 0;
                ss_states{1} = ss;
                q_values = zeros(nn,1);
                moved = zeros(nn,1);
                for i = 1:nn   

                    sslocal = ss_states{i};
                    qloc = -10;
                    for j = 1:nn

                        if(moved(j) == 1); continue; end;

                        sslocal(j)  = sslocal(j) * -1;
                        newq = obj.ModularityTwo(sslocal,blocal);

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
                
                [max_q max_state] = max(q_values);
                
                if(max_q > qglobal)
                    qglobal = max_q;
                    ss = ss_states{max_state+1};
                    qinc = 1;
                end
            end
            
            ssnew = ss;
            
        end
        
        function qq = ModularityTwo(obj,ss,bb)
            qq = (ss'*bb*ss) / (4*obj.n_edges);
        end
        
        function qq = Modularity(obj,ss,bb)
            ncom = length(unique(ss));
            nn = length(bb);
            S = zeros(nn,ncom);
            for i = 1:nn
                S(i,ss(i)) = 1;
            end
            qq = trace(S'*bb*S) / (2*obj.n_edges);
        end

        function obj = CalculateQrValue(obj,nedges)
           
            obj.Qr = 0;
            
            for i = 1:obj.N
                row_index = find(obj.rr(:,i));
                col_index = find(obj.tt(:,i));
                nr = length(row_index);
                nc = length(col_index);
                
                for j = 1:nr
                    for k = 1:nc
                        if(obj.matrix(row_index(j),col_index(k)) > 0)
                            obj.Qr = obj.Qr + 1;
                        end
                    end
                end
            end
            
            obj.Qr = obj.Qr / nedges;
        end
        
    end
    
end