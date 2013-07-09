classdef Descriptions
   
    methods(Static)
        
        function adj = ADJACENCY(web)
            adj = web > 0;
        end
        
        function inter = INTERACTIONS(web)
           
            adj = web > 0;
            inter = adj > 0;
            
        end
        
        function web_size = WEB_SIZE(web)
            
            [n_rows n_cols] = size(web);
            web_size = n_rows*n_cols;
            
        end
        
        function connectance = CONNECTANCE(web)
            
            adj = web > 0;
            web_size = Descriptions.WEB_SIZE(web);
            connectance = sum(adj(:)) / web_size;
            
        end
        
        function generality = GENERALITY(web)
            
            generality = sum(web > 0, 2);
            
        end
        
        function vulnerability = VULNERALITY(web)
           
            vulnerability = sum(web > 0, 1);
            
        end
        

        
    end
    
end