%==========================================================================
% Name        : ModularityBRIMBarberVersion.m
% Author      : C�sar Flores
% Created     : 14/Apr/2011
% Updated     : 14/Apr/2011
% Description : It is just an interface for using the Barber version algorithm
%==========================================================================
%              

classdef AdaptiveBrim < Modularity

    properties%(GetAccess = 'private', SetAccess = 'private')
        RR                   = [];  %Red Nodes (rows) Communities matrix. Size = nRows*CommunityQuantity
        TT                   = [];  %Blue Nodes (columns) Communities. Size = nCols*CommunityQuantity
        nRows                = 0;   %Number of rows
        nCols                = 0;   %Number of columns
        nEdges               = 0;
        PP                   = [];  %Null model matrix
        BB                   = [];  %Original - Null.
        RowComSize           = [];
        ColComSize           = [];
        IndexRow             = [];  %Register of the swaps in Rows.
        IndexCol             = [];  %Register of the swaps in Cols.
        SortedMatrix         = [];
        Sorted               = 0;
        Qr                   = 0;
        Qb                   = 0;
        Name = 'name';
        N                    = 0;
    end
    
    %DEBUG Properties - Change to parametrize and Debug the algorithm;
    properties
        Trials               = 100;  %How many trials by iteraction.
        QHistTrials          = [];
        ComQuantityHist      = [];
        RRHistoryBef         = {};
        RRHistoryAf          = {};
        TTHistoryBef         = {};
        TTHistoryAf          = {};
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURES ALGORITHM
    methods
        
        function obj = AdaptiveBrim(bipweb)
            %ModularityBRIM(bipartiteNetwork): Creates an object of the
            %class ModularityBRIM in order to use the BRIM algorithm in the
            %object netBipartite. The last one must be an instance of
            %the class NetworkBipartite
            %
            % Input arguments -
            %   netBipartite  = This an object of the class
            %                   NetworkBipartite
            % Output          -
            %   obj           = An object of the class ModularityBrim
            % 
            % Example:
            %   A = stats.phage_host;
            %   netbip = NetworkBipartite(A);
            %   brim = ModularityBRIM(netbip);
            
           
            
            %obj.Network = netBipartite;
            
            obj.Matrix = bipweb > 0;
            
            obj.Trials = 100;
            
        end
        
        
        function obj = Detect(obj)
           
            obj.nEdges = sum(obj.Matrix(:));
            [obj.nRows, obj.nCols] = size(obj.Matrix);
            obj.BB = bpmodmat(obj.Matrix);
            
            qmax = -0.5; %Negative symbol avoids problems in full connected graphs
            for i = 1:obj.Trials
          
                [R, T, Qhist1, Ahist1, Nhist1] = abrim(obj.BB, obj.nEdges);
                
                if (Qhist1(end) > qmax)
                   
                    qmax = Qhist1(end);
                    RRmax = R;
                    TTmax = T;
                    
                end
                
            end

            obj.RR = RRmax;
            obj.TT = TTmax;
            obj.Q = qmax;

            obj.CleanCommunities();
            
            obj.CalculateQrValue();
            obj.N = obj.CommunityQuantity;
            obj.Qb = obj.Q;
            
        end
            
    end
    
    methods
        
        function obj = CalculateQrValue(obj)
           
            obj.Qr = 0;
            rr = obj.RR;
            tt = obj.TT;
            nedges = sum(sum(obj.Matrix>0));
            
            for i = 1:obj.CommunityQuantity
                row_index = find(rr(:,i));
                col_index = find(tt(:,i));
                nr = length(row_index);
                nc = length(col_index);
                
                for j = 1:nr
                    for k = 1:nc
                        if(obj.Matrix(row_index(j),col_index(k)) > 0)
                            obj.Qr = obj.Qr + 1;
                        end
                    end
                end
            end
            
            obj.Qr = obj.Qr / nedges;
        end
       
        function obj = CleanCommunities(obj)
           
            %CLEAN THE COMMUNITIES
            com1 = find(any(obj.TT));
            com2 = find(any(obj.RR));
            
            if(length(com1) >= length(com2)) %Not sure if this validation is necesarry.
                obj.RR = obj.RR(:,com1);
                obj.TT = obj.TT(:,com1);
            else
                obj.RR = obj.RR(:,com2);
                obj.TT = obj.TT(:,com2);
            end
            
            obj.CommunityQuantity = size(obj.RR,2);
            
        end
        
        function obj = SortCommunities(obj,sortbyroworcol)
            %Not a part of the main Barber code.
            %This is a function that returns an useful matrix for plotting
            %the modularity pattern.
            
            
            if(nargin == 1)
               sortbyroworcol = 1;      
            end
            
            obj.Sorted = 1;
            
            obj.IndexRow = 1:obj.nRows;
            obj.IndexCol = 1:obj.nCols;
            
            rrTemp = obj.RR;
            ttTemp = obj.TT;
            
            %Sizes of the communities.
            obj.RowComSize = sum(rrTemp);
            obj.ColComSize = sum(ttTemp);
                                    
            %Variables that will store the community index of each row,
            %column.
            ssRow = 1:obj.nRows;
            ssCol = 1:obj.nCols;            
            
            %We sort by decrease order in community size in rows or cols
            if(sortbyroworcol == 1)
                [sorted index] = sort(obj.RowComSize,'descend');
                obj.RowComSize = sorted;
                obj.ColComSize = obj.ColComSize(index);
            else
                [sorted index] = sort(obj.ColComSize,'descend');
                obj.ColComSize = sorted;
                obj.RowComSize = obj.RowComSize(index);
            end
                
            obj.SortedMatrix = double(obj.Matrix);
            
            %First sort the community matrices
            rrTemp = rrTemp(:,index);
            ttTemp = ttTemp(:,index);
            
            %Sort Rows
            [rrTemp ind] = sortrows(-rrTemp);
            rrTemp = (rrTemp < 0);
            %rrTemp = flipud(rrTemp);
            %ind = flipud(ind);
            obj.SortedMatrix = double(obj.Matrix(ind,:));
            obj.IndexRow(:) = ind;
            [i,j] = ind2sub(size(rrTemp),find(rrTemp>0));
            ssRow(:) = j;

            %Sort Columns
            [ttTemp ind] = sortrows(-ttTemp);
            ttTemp = (ttTemp < 0);
            %ttTemp = flipud(ttTemp);
            %ind = flipud(ind);
            obj.SortedMatrix = obj.SortedMatrix(:,ind);
            obj.IndexCol(:) = ind;
            [i,j] = ind2sub(size(ttTemp),find(ttTemp>0));
            ssCol(:) = j;
            
            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    if(obj.SortedMatrix(i,j) == 1)
                       
                        if(ssRow(i) == ssCol(j))
                            obj.SortedMatrix(i,j) = ssRow(i);
                        else
                            obj.SortedMatrix(i,j) = -1;
                        end
                        
                    end
                    
                end
            end

        end
        
    end
    
end