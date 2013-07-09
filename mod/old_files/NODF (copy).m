function nodf = NODF(matrix, nodf_strict)
%TO DO, implement the decreasing strict (nodf_sctrict) section
    
    [n_rows n_cols] = size(matrix);

    C1 = (1.0*matrix)*matrix';
    D1 = repmat(diag(C1),1,n_rows);
    nr = (D1~=D1');%
    nr = C1.*nr;
    nr = nr./(min(D1,D1'));

    C2 = (1.0*matrix')*matrix;
    D2 = repmat(diag(C2),1,n_cols);
    nc = C2.*(D2~=D2');
    nc = nc./(min(D2,D2'));

    nr = 100*nansum(nr(:))/2;
    nc = 100*nansum(nc(:))/2;
    nvalue = nr+nc;
    
    denom = n_rows*(n_rows-1)/2 + n_cols*(n_cols-1)/2;
    nr = nr/(n_rows*(n_rows-1)/2);
    nc = nc/(n_cols*(n_cols-1)/2);   
    nvalue = nvalue/denom;
    
    nodf.nodf_rows = nr;
    nodf.nodf_cols = nc;
    nodf.nodf = nvalue;
    
end