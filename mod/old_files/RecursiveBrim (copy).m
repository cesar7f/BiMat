%==========================================================================
% Name        : ModularityBRIMRecursive.m
% Author      : Cesar
% Created     : 30/May/2011
% Updated     : 30/May/2011
% Description : Represents a recursive version of the BRIM algorithm
%==========================================================================
%0.9167

classdef RecursiveBrim < Modularity

    properties%(GetAccess = 'private', SetAccess = 'private')
        RR                   = [];  %Red Nodes (rows) Communities matrix. Size = nRows*CommunityQuantity
        TT                   = [];  %Blue Nodes (columns) Communities. Size = nCols*CommunityQuantity
        nRows                = 0;   %Number of rows
        nCols                = 0;   %Number of columns
        PP                   = [];  %Null model matrix
        BB                   = [];  %Original - Null.
        RowComSize           = [];
        ColComSize           = [];
        IndexRow             = [];  %Register of the swaps in Rows.
        IndexCol             = [];  %Register of the swaps in Cols.
        SortedMatrix         = [];
        SortByRows           = 1;
        Sorted               = 0;
        nLevels              = 10000;
        FinalLevel           = 0;
        Parents              = {};
        Heuristic            = 1;
        M                    = 0;
        Mij                  = [];
        I                    = [];
        Iij                  = [];
        FinalLevelTop        = 0;
        Rho                  = 0;
        Probs                = [];
        ProbAcum             = [];
        Name = 'name';
        SSRow                = [];
        SSCol                = [];
        Qr                   = 0;
        Qb                   = 0;
        N                    = 0;
    end
    
    %DEBUG Properties - Change to parametrize and Debug the algorithm;
    properties
        PortionNewComm       = 2;   %New size of communities based in the previous size. Default = 2
        PortionRandom        = 0.5; %Pourcentage of the nodes that are assigned to the communities created recently. Default = 0.5
        PortionCommAssigment = 0.5; %Portion of the communities that will be used to reassign the nodes. 
                                    %This pourcentage starts in the last communities. 
                                    %For example with 10 communities, a
                                    %value of 0.5 will use the communities
                                    %from 6 to 10. Using the default values
                                    %will create a situation where the new
                                    %created communities will be the
                                    %choseen comunities.
                                    %Default = 0.5
        Trials               = 100;  %How many trials by iteraction.
        tolerance            = 0.0001
        DebugMessages        = 0;  %1 for printing messages, 0 otherwhise
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURES ALGORITHM
    methods
        function obj = RecursiveBrim(bipmatrix)
                      
            obj.Matrix = bipmatrix > 0;
            [obj.nRows obj.nCols] = size(obj.Matrix);
            %obj.nLevels = 2;
            obj.Trials = 200;
            %obj.Network = netBipartite;
            
        end
        
        function obj = Detect(obj)
            obj.Heuristic = HierarchicalHeuristic.OPTIMAL_GLOBAL;
            %obj.HierarchicalBrimOptimization();
            obj.HierarchicalBRIM();
            %obj.DivideInIndependentComponents();
            %obj.FinalLevel = 1;
            %obj.Q = q;
            %display(obj.nRows);
            obj.CalculateQrValue();
            obj.N = obj.CommunityQuantity(end);
            obj.Qb = obj.Q(end);
        end
        
        function obj = DivideInIndependentComponents(obj)
           
            mm = [zeros(obj.nRows, obj.nRows) obj.Matrix; obj.Matrix' zeros(obj.nCols, obj.nCols)];
            [S, C] = graphconncomp(sparse(mm),'Weak', true);
            
            idrows = 1:obj.nRows;
            idcols = (obj.nRows+1):(obj.nRows+obj.nCols);
            
            rowsisolated = setdiff(C(idrows),C(idcols));
            colsisolated = setdiff(C(idcols),C(idrows));
            
            rr = zeros(obj.nRows, S);
            tt = zeros(obj.nCols, S);
            
            
            rr(sub2ind(size(rr),idrows,C(idrows))) = 1;
            tt(sub2ind(size(tt),1:obj.nCols,C(idcols))) = 1;
            
            %Clean communities with isolated nodes
            rr(:,[rowsisolated colsisolated]) = [];
            tt(:,[rowsisolated colsisolated]) = [];
            
            rr(:,end+1) = 0;
            tt(:,end+1) = 0;
            
            rowsisolated = sum(rr,2)==0;
            colsisolated = sum(tt,2)==0;
            
            rr(rowsisolated,end) = 1;
            tt(colsisolated,end) = 1;
            
            [rr tt cq]  = obj.CleanCommunities(rr,tt);

            obj.RR{1} = rr;
            obj.TT{1} = tt;
            %obj.FinalLevel = 1;
            obj.RowComSize{1} = sum(obj.RR{1});
            obj.ColComSize{1} = sum(obj.TT{1});
            obj.CommunityQuantity(1) = cq;
            
        end
        
        function obj = DivideInOneCommunity(obj)
           
            obj.RR{1} = ones(obj.nRows,1);
            obj.TT{1} = ones(obj.nCols,1);
            
            obj.RowComSize{1} = sum(obj.RR{1});
            obj.ColComSize{1} = sum(obj.TT{1});
            obj.CommunityQuantity(1) = 1;
            
        end
        
        function [rmax tmax cquant] = CleanCommunities(obj,rmax,tmax)
           
            %CLEAN THE COMMUNITIES
            com1 = find(any(tmax));
            com2 = find(any(rmax));
            
            if(length(com1) >= length(com2)) %Not sure if this validation is necesarry.
                rmax = rmax(:,com1);
                tmax = tmax(:,com1);
            else
                rmax = rmax(:,com2);
                tmax = tmax(:,com2);
            end
            
            if(size(rmax,2)~=size(tmax,2)); fprintf('eureka'); end;
            
            cquant = size(rmax,2);
            
            if(length(com1) ~= length(com2)) %This avoids having a empty columns in one of the matrices.
                cquant = -1;
            end;
            
        end


        
        function obj = CalculateQrValue(obj)
           
            obj.Qr = 0;
            rr = obj.RR{end};
            tt = obj.TT{end};
            nedges = sum(sum(obj.Matrix>0));
            
            for i = 1:obj.CommunityQuantity(end)
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

        
        function initialQ = GetInitialQ(obj, localbb)
                    
            [nrows ncols] = size(localbb);
            
            rr = assignsamemodule(nrows,1);
            tt = assignsamemodule(ncols,1);
            
            initialQ = bpmodularity(rr, localbb, tt, sum(obj.Matrix(:)));
                
        end
        
        function obj = HierarchicalBRIM(obj)
           
            obj.BB = bpmodmat(obj.Matrix); %Global ~B
            obj.DivideInIndependentComponents();
            %obj.DivideInOneCommunity();
            
            %obj.RR{1} = assignsamemodule(obj.nRows,1);
            %obj.TT{1} = assignsamemodule(obj.nCols,1);
            obj.Q(1) = bpmodularity(obj.RR{1}, obj.BB, obj.TT{1}, sum(obj.Matrix(:)));% = -0.5; %Avoid complete connected graph problems;
            %obj.RowComSize{1} = sum(obj.RR{1});
            %obj.ColComSize{1} = sum(obj.TT{1});
            %obj.CommunityQuantity(1) = 1;
            numedges = sum(obj.Matrix(:));
            %obj.nLevels
            succes = 1;
            obj.Parents{1} = 0;
            for i = 1:obj.nLevels

                obj.RR{i+1} = [];%zeros(obj.nRows,1);
                obj.TT{i+1} = [];%zeros(obj.nCols,1);
                rrnext = zeros(obj.nRows,1);
                ttnext = zeros(obj.nCols,1);
                comquantitynext = 0;
                obj.CommunityQuantity(i+1) = 0;
                rr = obj.RR{i};
                tt = obj.TT{i};
                parents = [];
                
                
                succes = 0;
                for j = 1:obj.CommunityQuantity(i);
                
                    indrow = find(rr(:,j)==1);
                    indcol = find(tt(:,j)==1);

                    localmatrix = obj.Matrix(indrow,indcol);
                    [~, ncols] = size(localmatrix);
                    
                    bb = [];
                    nedges = 0;
                    if(obj.Heuristic == HierarchicalHeuristic.SIMPLE_GLOBAL || obj.Heuristic == HierarchicalHeuristic.SIMPLE_LOCAL)
                        bb =  bpmodmat(localmatrix);
                        nedges = sum(sum(localmatrix));
                    elseif (obj.Heuristic == HierarchicalHeuristic.OPTIMAL_GLOBAL || obj.Heuristic == HierarchicalHeuristic.OPTIMAL_LOCAL)
                        bb = obj.BB(indrow,indcol);
                        nedges = numedges;
                    end
                    
                          
                    qmax = -0.5;
                    rlocal = [];
                    tlocal = [];
                    for k = 1:obj.Trials
                        
                        tt0 = assignrandmodule(ncols,2);
                        try
                            [rlocaltrial tlocaltrial qbrim] = brim(bb,tt0,nedges);
                        catch err
                            display('error');
                        end;
                        
                        if(qbrim(end) > qmax)
                            
                            qmax = qbrim(end);
                            rlocal = rlocaltrial;
                            tlocal = tlocaltrial;
                            
                        end
                        
                    end
                    
                    deltaq = 0;
                    if(obj.Heuristic == HierarchicalHeuristic.SIMPLE_LOCAL || obj.Heuristic == HierarchicalHeuristic.OPTIMAL_LOCAL)

                        if(obj.Heuristic == HierarchicalHeuristic.SIMPLE_LOCAL)
                            qmax =  bpmodularity(rlocal, obj.BB(indrow,indcol), tlocal, numedges);
                        end

                        deltaq = qmax - obj.GetInitialQ(obj.BB(indrow,indcol));

                    else
                        deltaq = 10; %Any positive quantity will do the trick
                    end
                    %qmax
                    %deltaq
                    [rlocal tlocal cquantity]  = obj.CleanCommunities(rlocal,tlocal);
                    if(cquantity == 2 && deltaq >= 0) %Divide in two only if there is an icrease in the modularity.
                        
                        index = [];
                        if(obj.SortByRows == 1)
                            [~, index] = sort(sum(rlocal),'descend');
                        else
                            [~, index] = sort(sum(tlocal),'descend');
                        end
                        
                        rlocal = rlocal(:,index);
                        tlocal = tlocal(:,index);
                        
                        rrnext(indrow,comquantitynext+1:comquantitynext+2) = rlocal;
                        ttnext(indcol,comquantitynext+1:comquantitynext+2) = tlocal;
                        succes = 1;
                        comquantitynext = cquantity+comquantitynext;
                        
                        parents = [parents j j];
                    else %Do not divide if deltaq is negative or there is not division
                        rlocal = rr(indrow,j); %Retreive rlocal in case there is a division
                        tlocal = tt(indcol,j); %''
                        rrnext(indrow,comquantitynext+1) = rlocal;
                        ttnext(indcol,comquantitynext+1) = tlocal;
                        comquantitynext = 1+comquantitynext;
                        parents = [parents j];
                    end
                    
                    
                end
                obj.RR{i+1} = rrnext;
                obj.TT{i+1} = ttnext;
                obj.CommunityQuantity(i+1) = comquantitynext;
                obj.Q(i+1) = bpmodularity(rrnext, obj.BB, ttnext, sum(obj.Matrix(:)));
                obj.Parents{i+1} = parents;
                
                obj.RowComSize{i+1} = sum(rrnext);
                obj.ColComSize{i+1} = sum(ttnext);
                %qmax
                if( (obj.Heuristic == HierarchicalHeuristic.OPTIMAL_LOCAL || obj.Heuristic == HierarchicalHeuristic.SIMPLE_LOCAL) && succes == 0)
                    %fprintf('salgo en el primer if');
                    break;
                elseif( (obj.Heuristic == HierarchicalHeuristic.OPTIMAL_GLOBAL || obj.Heuristic == HierarchicalHeuristic.SIMPLE_GLOBAL) && (obj.Q(i+1) < obj.Q(i) || succes == 0))
                    %fprintf('salgo en el segundo if');
                    %obj.Q(i+1)
                    %obj.Q(i)
                    
                    break;
                end
            end
            
            obj.RR(end) = [];
            obj.Parents(end) = [];
            obj.TT(end) = [];
            obj.Q(end) = [];
            %obj.ColComSize = [];
            %obj.RowComSize = [];
            obj.CommunityQuantity(end) = [];
            obj.FinalLevel = length(obj.Q);
            
            
        end

    end
end
% OLD CODE
% function obj = HierarchicalBrimRandomDivision(obj)
%             
%             obj.BB = bpmodmat(obj.Matrix); %Global ~B
%             
%             obj.Q(1) = 0;% = -0.5; %Avoid complete connected graph problems;
%             obj.RR{1} = assignsamemodule(obj.nRows,1);
%             obj.TT{1} = assignsamemodule(obj.nCols,1);
%             obj.CommunityQuantity(1) = 1;
%             %numedges = sum(obj.Matrix(:));
%             %obj.nLevels
%             succes = 1;
%             for i = 1:obj.nLevels
%                
%                 obj.RR{i+1} = [];%zeros(obj.nRows,1);
%                 obj.TT{i+1} = [];%zeros(obj.nCols,1);
%                 rrnext = zeros(obj.nRows,1);
%                 ttnext = zeros(obj.nCols,1);
%                 comquantitynext = 0;
%                 obj.CommunityQuantity(i+1) = 0;
%                 rr = obj.RR{i};
%                 tt = obj.TT{i};
%                 
%                 if(succes == 0); break; end;
%                 
%                 succes = 0;
%                 for j = 1:obj.CommunityQuantity(i);
%                 
%                     indrow = find(rr(:,j)==1);
%                     indcol = find(tt(:,j)==1);
% 
%                     localmatrix = obj.Matrix(indrow,indcol);
%                     localbb = bpmodmat(localmatrix);
%                     localnumedges = sum(sum(localmatrix));
%                     
%                     [~, ncols] = size(localmatrix);
%                     
%                     qmax = -0.5;
%                     rlocal = [];
%                     tlocal = [];
%                     for k = 1:obj.Trials
%                         
%                         tt0 = assignrandmodule(ncols,2);
%                         [rlocaltrial tlocaltrial qbrim] = brim(localbb,tt0,localnumedges);
%                         
%                         if(qbrim(end) > qmax)
%                             
%                             qmax = qbrim(end);
%                             rlocal = rlocaltrial;
%                             tlocal = tlocaltrial;
%                             
%                         end
%                         
%                     end
%                     
%                     qrandom = NullModelGenerator.GetModularitySample(localmatrix,100,NullModel.BERNOULLI_CONSTRAINED);
%                     
%                     qsorted = sort(qrandom);
%                     
%                     [rlocal tlocal cquantity]  = obj.CleanCommunities(rlocal,tlocal);
%                     
%                     if(cquantity == 2 && qmax > qsorted(3)) %Divide in two only if there is an icrease in the modularity.
%                         rrnext(indrow,comquantitynext+1:comquantitynext+2) = rlocal;
%                         ttnext(indcol,comquantitynext+1:comquantitynext+2) = tlocal;
%                         succes = 1;
%                         comquantitynext = cquantity+comquantitynext;
%                     else %Do not divide if deltaq is negative or there is not division
%                         rlocal = rr(indrow,j); %Retreive rlocal in case there is a division
%                         tlocal = tt(indcol,j); %''
%                         rrnext(indrow,comquantitynext+1) = rlocal;
%                         ttnext(indcol,comquantitynext+1) = tlocal;
%                         comquantitynext = 1+comquantitynext;
%                     end
%                     
%                     
%                 end
%                 obj.RR{i+1} = rrnext;
%                 obj.TT{i+1} = ttnext;
%                 obj.CommunityQuantity(i+1) = comquantitynext;
%                 obj.Q(i+1) = bpmodularity(rrnext, obj.BB, ttnext, sum(obj.Matrix(:)));
%             end
%             
%             obj.RR(end) = [];
%             obj.Parents(end) = [];
%             obj.TT(end) = [];
%             obj.Q(end) = [];
%             obj.CommunityQuantity(end) = [];
% end
%         
%         function obj = HierarchicalBrimSimpleDivision(obj)
%             
%             obj.BB = bpmodmat(obj.Matrix); %Global ~B
%             
%             obj.Q(1) = 0;% = -0.5; %Avoid complete connected graph problems;
%             obj.RR{1} = assignsamemodule(obj.nRows,1);
%             obj.TT{1} = assignsamemodule(obj.nCols,1);
%             obj.RowComSize{1} = sum(obj.RR{1});
%             obj.ColComSize{1} = sum(obj.TT{1});
%             obj.CommunityQuantity(1) = 1;
%             %numedges = sum(obj.Matrix(:));
%             %obj.nLevels
%             succes = 1;
%             obj.Parents{1} = 0;
%             for i = 1:obj.nLevels
%                
%                 if(succes == 0); break; end;
%                 
%                 obj.RR{i+1} = [];%zeros(obj.nRows,1);
%                 obj.TT{i+1} = [];%zeros(obj.nCols,1);
%                 rrnext = zeros(obj.nRows,1);
%                 ttnext = zeros(obj.nCols,1);
%                 comquantitynext = 0;
%                 obj.CommunityQuantity(i+1) = 0;
%                 rr = obj.RR{i};
%                 tt = obj.TT{i};
%                 parents = [];
%                 
%                 
%                 succes = 0;
%                 for j = 1:obj.CommunityQuantity(i);
%                 
%                     indrow = find(rr(:,j)==1);
%                     indcol = find(tt(:,j)==1);
% 
%                     localmatrix = obj.Matrix(indrow,indcol);
%                     localbb = bpmodmat(localmatrix);
%                     localnumedges = sum(sum(localmatrix));
%                     
%                     [~, ncols] = size(localmatrix);
%                     
%                     qmax = -0.5;
%                     rlocal = [];
%                     tlocal = [];
%                     for k = 1:obj.Trials
%                         
%                         tt0 = assignrandmodule(ncols,2);
%                         [rlocaltrial tlocaltrial qbrim] = brim(localbb,tt0,localnumedges);
%                         
%                         if(qbrim(end) > qmax)
%                             
%                             qmax = qbrim(end);
%                             rlocal = rlocaltrial;
%                             tlocal = tlocaltrial;
%                             
%                         end
%                         
%                     end
%                     
%                     %deltaq = qmax - obj.GetInitialQ(localmatrix,obj.BB(indrow,indcol));
%                     
%                     %rtemp = rlocal;
%                     %ttemp = tlocal;
%                     [rlocal tlocal cquantity]  = obj.CleanCommunities(rlocal,tlocal);
%                     
%                     if(cquantity == 2) %Divide in two only if there is an icrease in the modularity.
%                         
%                         index = [];
%                         if(obj.SortByRows == 1)
%                             [~, index] = sort(sum(rlocal),'descend');
%                         else
%                             [~, index] = sort(sum(tlocal),'descend');
%                         end
%                         
%                         rlocal = rlocal(:,index);
%                         tlocal = tlocal(:,index);
%                         
%                         rrnext(indrow,comquantitynext+1:comquantitynext+2) = rlocal;
%                         ttnext(indcol,comquantitynext+1:comquantitynext+2) = tlocal;
%                         succes = 1;
%                         comquantitynext = cquantity+comquantitynext;
%                         
%                         parents = [parents j j];
%                     else %Do not divide if deltaq is negative or there is not division
%                         rlocal = rr(indrow,j); %Retreive rlocal in case there is a division
%                         tlocal = tt(indcol,j); %''
%                         rrnext(indrow,comquantitynext+1) = rlocal;
%                         ttnext(indcol,comquantitynext+1) = tlocal;
%                         comquantitynext = 1+comquantitynext;
%                         parents = [parents j];
%                     end
%                     
%                     
%                 end
%                 obj.RR{i+1} = rrnext;
%                 obj.TT{i+1} = ttnext;
%                 obj.CommunityQuantity(i+1) = comquantitynext;
%                 obj.Q(i+1) = bpmodularity(rrnext, obj.BB, ttnext, sum(obj.Matrix(:)));
%                 obj.Parents{i+1} = parents;
%                 
%                 obj.RowComSize{i+1} = sum(rrnext);
%                 obj.ColComSize{i+1} = sum(ttnext);
%             end
%             
%             obj.RR(end) = [];
%             obj.Parents(end) = [];
%             obj.TT(end) = [];
%             obj.Q(end) = [];
%             obj.CommunityQuantity(end) = [];
%             obj.FinalLevel = length(obj.Q);
%         end
%         
%         function obj = HierarchicalBrimOptimization(obj)
%             
%             obj.BB = bpmodmat(obj.Matrix); %Global ~B
%             
%             obj.Q(1) = 0;% = -0.5; %Avoid complete connected graph problems;
%             obj.RR{1} = assignsamemodule(obj.nRows,1);
%             obj.TT{1} = assignsamemodule(obj.nCols,1);
%             obj.CommunityQuantity(1) = 1;
%             numedges = sum(obj.Matrix(:));
%             %obj.nLevels
%             succes = 1;
%             for i = 1:obj.nLevels
%                
%                 obj.RR{i+1} = [];%zeros(obj.nRows,1);
%                 obj.TT{i+1} = [];%zeros(obj.nCols,1);
%                 rrnext = zeros(obj.nRows,1);
%                 ttnext = zeros(obj.nCols,1);
%                 comquantitynext = 0;
%                 obj.CommunityQuantity(i+1) = 0;
%                 rr = obj.RR{i};
%                 tt = obj.TT{i};
%                 
%                 if(succes == 0); break; end;
%                 
%                 succes = 0;
%                 for j = 1:obj.CommunityQuantity(i);
%                 
%                     indrow = find(rr(:,j)==1);
%                     indcol = find(tt(:,j)==1);
% 
%                     localmatrix = obj.Matrix(indrow,indcol);
% 
%                     [nrows ncols] = size(localmatrix);
%                           
%                     imagesc(localmatrix);
%                     qmax = -0.5;
%                     rlocal = [];
%                     tlocal = [];
%                     for k = 1:obj.Trials
%                         
%                         tt0 = assignrandmodule(ncols,2);
%                         [rlocaltrial tlocaltrial qbrim] = brim(obj.BB(indrow,indcol),tt0,numedges);
%                         
%                         if(qbrim(end) > qmax)
%                             
%                             qmax = qbrim(end);
%                             rlocal = rlocaltrial;
%                             tlocal = tlocaltrial;
%                             
%                         end
%                         
%                     end
%                     
%                     deltaq = qmax - obj.GetInitialQ(localmatrix,obj.BB(indrow,indcol));
%                     
%                     rtemp = rlocal;
%                     ttemp = tlocal;
%                     [rlocal tlocal cquantity]  = obj.CleanCommunities(rlocal,tlocal);
%                     
%                     if(cquantity == 2 && deltaq >= 0) %Divide in two only if there is an icrease in the modularity.
%                         rrnext(indrow,comquantitynext+1:comquantitynext+2) = rlocal;
%                         ttnext(indcol,comquantitynext+1:comquantitynext+2) = tlocal;
%                         succes = 1;
%                         comquantitynext = cquantity+comquantitynext;
%                     else %Do not divide if deltaq is negative or there is not division
%                         rlocal = rr(indrow,j); %Retreive rlocal in case there is a division
%                         tlocal = tt(indcol,j); %''
%                         rrnext(indrow,comquantitynext+1) = rlocal;
%                         ttnext(indcol,comquantitynext+1) = tlocal;
%                         comquantitynext = 1+comquantitynext;
%                     end
%                     
%                     
%                 end
%                 obj.RR{i+1} = rrnext;
%                 obj.TT{i+1} = ttnext;
%                 obj.CommunityQuantity(i+1) = comquantitynext;
%                 obj.Q(i+1) = bpmodularity(rrnext, obj.BB, ttnext, sum(obj.Matrix(:)));
%             end
%             
%         end
%         
%         function [rrec trec cquantity q] = RecursiveBRIM(obj,matrix,level)
% 
%             bodModmat = bpmodmat(matrix);
%             [nrows ncols] = size(matrix);
%             numedges = sum(matrix(:));
%             
%             Qmax = -0.5; %Avoid complete connected graphs problems
%             Tmax = [];
%             Rmax = [];
%             for i = 1:obj.Trials
%                 T0 = assignrandmodule(ncols, 2);
%                 [R, T, Qbrim] = brim(bodModmat, T0, numedges);
%                 
%                 if(Qbrim(end) > Qmax)
%                     Qmax = Qbrim(end);
%                     Tmax = T;
%                     Rmax = R;
%                 end
%             end            
%             
%             obj.Q = Qbrim(end);
%             q = obj.Q;
%             
%             [Rmax Tmax cquantity] = obj.CleanCommunities(Rmax,Tmax);
%             
%             trec = Tmax;
%             rrec = Rmax;
%             
%             qbef = bpmodularity(Rmax, bodModmat, Tmax, numedges);
%             
%             if(qbef < 0.3)
%                 return;
%             end
% 
%             if(cquantity > 1 && level < obj.nLevels)
%             
%                 a = {};
%                 b = {};
%                 rmaxloc = {};%zeros(size(Rmax));
%                 tmaxloc = {};%zeros(size(Tmax));
%                 for i = 1:2
%                 
%                     matrixnew = matrix(Rmax(:,i)==1,Tmax(:,i)==1); 
%                 
%                     [rloc tloc cq qq] = obj.RecursiveBRIM(matrixnew,level+1);
%                     a{i} = find(Rmax(:,i)==1);
%                     b{i} = find(Tmax(:,i)==1);
%                     
%                     rmaxloc{i} = rloc;
%                     tmaxloc{i} = tloc;
%         
%                 end
%                 
%                 Rmax = zeros(size(rmaxloc{1},1)+size(rmaxloc{2},1),size(rmaxloc{1},2)+size(rmaxloc{2},2));
%                 Tmax = zeros(size(tmaxloc{1},1)+size(tmaxloc{2},1),size(tmaxloc{1},2)+size(tmaxloc{2},2));
%                 
%                 Rmax(a{1},1:size(rmaxloc{1},2)) = rmaxloc{1};
%                 Tmax(b{1},1:size(tmaxloc{1},2)) = tmaxloc{1};
%                 Rmax(a{2},size(rmaxloc{1},2)+1:size(rmaxloc{1},2)+size(rmaxloc{2},2)) = rmaxloc{2};
%                 Tmax(b{2},size(tmaxloc{1},2)+1:size(tmaxloc{1},2)+size(tmaxloc{2},2)) = tmaxloc{2};
%                 
%             end
%             
%             rrec = Rmax;
%             trec = Tmax;
%             q = bpmodularity(Rmax, bodModmat, Tmax, numedges);
%             
%         end
%                
%        
%         
%         function obj = SortCommunitiesTemp(obj,level)
%             %Not a part of the main Barber code.
%             %This is a function that returns an useful matrix for plotting
%             %the modularity pattern.
%             
%             
%             sortbyroworcol = 1;      
%             
%             obj.Sorted = 0;
%             
%             obj.IndexRow = 1:obj.nRows;
%             obj.IndexCol = 1:obj.nCols;
%             
%             rrTemp = obj.RR{level};
%             ttTemp = obj.TT{level};
%             
%             %Sizes of the communities.
%             obj.RowComSize = sum(rrTemp);
%             obj.ColComSize = sum(ttTemp);
%             obj.RowComSize                        
%             %Variables that will store the community index of each row,
%             %column.
%             ssRow = 1:obj.nRows;
%             ssCol = 1:obj.nCols;            
%             
%             %We sort by decrease order in community size in rows or cols
%             if(sortbyroworcol == 1)
%                 %obj.ColComSize
%                 %rrTemp
%                 %obj.RowComSize
%                 [sorted index] = sort(obj.RowComSize,'descend');
%                 obj.RowComSize = sorted;
%                 sorted
%                 %index
%                 
%                 obj.ColComSize = obj.ColComSize(index);
%             else
%                 [sorted index] = sort(obj.ColComSize,'descend');
%                 obj.ColComSize = sorted;
%                 obj.RowComSize = obj.RowComSize(index);
%             end
%                 
%             obj.SortedMatrix = double(obj.Matrix);
%             
%             %First sort the community matrices
%             rrTemp = rrTemp(:,index);
%             ttTemp = ttTemp(:,index);
%             
%             %Sort Rows
%             [rrTemp ind] = sortrows(-rrTemp);
%             rrTemp = (rrTemp < 0);
%             %rrTemp = flipud(rrTemp);
%             %ind = flipud(ind);
%             %ind
%             %obj.IndexRow
%             %ind
%             obj.SortedMatrix = double(obj.Matrix(ind,:));
%             obj.IndexRow(:) = ind;
%             [i,j] = ind2sub(size(rrTemp),find(rrTemp>0));
%             ssRow(:) = j;
% 
%             %Sort Columns
%             [ttTemp ind] = sortrows(-ttTemp);
%             ttTemp = (ttTemp < 0);
%             %ttTemp = flipud(ttTemp);
%             %ind = flipud(ind);
%             obj.SortedMatrix = obj.SortedMatrix(:,ind);
%             obj.IndexCol(:) = ind;
%             [i,j] = ind2sub(size(ttTemp),find(ttTemp>0));
%             ssCol(:) = j;
%             
% %             start = 1;
% %             for i = 1:length(obj.RowComSize)
% %                
% %                 smallmatrix = obj.SortedMatrix(start:sum(obj.RowComSize(1:i)),:);
% %                 %smallmatrix
% %                 [values sorted] = sort(sum(smallmatrix,2),'descend');
% %                 
% %                 smallmatrix = smallmatrix(sorted,:);
% %                 
% %                 obj.SortedMatrix(start:sum(obj.RowComSize(1:i)),:) = smallmatrix;
% %                 
% %                 start = start + obj.RowComSize(i);
% %                 start
% %             end
% %             
% %             start = 1;
% %             for i = 1:length(obj.ColComSize)
% %                
% %                 %start:sum(obj.ColComSize(i))
% %                 smallmatrix = obj.SortedMatrix(:,start:sum(obj.ColComSize(1:i)));
% %                 %smallmatrix
% %                 [values sorted] = sort(sum(smallmatrix),'descend');
% %                 
% %                 smallmatrix = smallmatrix(:,sorted);
% %                 
% %                 obj.SortedMatrix(:,start:sum(obj.ColComSize(1:i))) = smallmatrix;
% %                 
% %                 start = start + obj.ColComSize(i);
% %                 start
% %             end
% 
% %             startrow = 1;
% %             startcol = 1;
% %             for i = 1:length(obj.RowComSize)
% %                
% %                 rowcomsize = obj.RowComSize(i);
% %                 colcomsize = obj.ColComSize(i);
% %                 
% %                 smallmatrix = obj.SortedMatrix(startrow:startrow+rowcomsize-1,startcol:startcol+colcomsize-1);
% %                 %smallmatrix
% %                 [values sortedrow] = sort(sum(smallmatrix,2),'descend');
% %                 [values sortedcol] = sort(sum(smallmatrix,1),'descend');
% %                 
% %                 smallmatrix = smallmatrix(sortedrow,sortedcol);
% %                 
% %                 obj.SortedMatrix(startrow:startrow+rowcomsize-1,startcol:startcol+colcomsize-1) = smallmatrix;
% %                 
% %                 startrow = startrow + obj.RowComSize(i);
% %                 startcol = startcol + obj.ColComSize(i);
% %            end
%             
%             for i = 1:obj.nRows
%                 for j = 1:obj.nCols
%                     if(obj.SortedMatrix(i,j) ~= 0)
%                        
%                         if(ssRow(i) == ssCol(j))
%                             obj.SortedMatrix(i,j) = ssRow(i);
%                         else
%                             obj.SortedMatrix(i,j) = -1;
%                         end
%                         
%                     end
%                     
%                 end
%             end
% 
%         end
        