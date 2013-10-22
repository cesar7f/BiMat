function commod(obj,eventdata,tanal)
%Does modularity according to Newman 2006 PNAS
global assocm
rnumind=length(assocm(1,:));
td=assocm-diag(diag(assocm));
ki=sum(td);m=sum(ki)/2;
switch tanal
    case 'Associations'
        typmod=questdlg('Modularity using?','Type of modularity','Modularity 1 (from gregariousness)','Modularity 2 (from 
permutations)','Cancel','Modularity 1 (from gregariousness)');
    case 'Interactions'
        typmod='Modularity 1 (from gregariousness)';
end
if ~strcmp(typmod,'Cancel')
    if strcmp(typmod,'Modularity 1 (from gregariousness)')
        stqd=ones(rnumind,1)*ki;
        stqd=stqd.*stqd'/(2*m);
        B=td-stqd;
    else
        pmatr=pertest(0,0,2);
        pmatr=pmatr-diag(diag(pmatr));
        B=full(td-sum(sum(td))*pmatr/sum(sum(pmatr)));
        B=B-diag(diag(B));%after suggestion from Newman to Lusseau
    end
    bq=B;
    numclus=1;
    clusass=ones(rnumind,1);
    clusfin(numclus)=1;
    prnam=labelmake;
    while 1
        uq=find(clusfin);
        if isempty(uq);break;end
        useind=find(clusass==uq(1));
        bq=B(useind,useind);
        kg=sum(td(useind,useind))-ki(useind)*sum(ki(useind))/(2*m);
        bq=bq-diag(kg);
        [V,D]=eig(bq);
        [dum,bige]=max(diag(D));
        cln=V(:,bige);%cluster assignments
        clna=(cln>0);nw=length(clna);
        while 1
            qq=clna*ones(1,nw);modu=sum(sum(bq.*(qq==qq')));
            sw=1;
            for i=1:length(clna)
                clnu=clna;clnu(i)=~clna(i);
                qqu=clnu*ones(1,nw);moduu=sum(sum(bq.*(qqu==qqu')));
                if moduu>modu
                    sw=0;
                    clna=clnu;modu=moduu;
                end
            end
            if sw;break;end
        end
        if sum(clna)*sum(~clna)
            numclus=numclus+1;
            clusass(useind(find(clna)))=numclus;
            clusval(useind)=cln;
            clusfin(numclus)=1;
        else
            clusfin(uq(1))=0;%finished with that cluster
        end
    end
    [dum,ql]=sort(clusass);
    disp('  ')
    disp('Assignment of individuals to clusters (eigs near zero indicate uncertainty)')
    switch tanal
        case 'Associations'
            disp(['  using eigenvector method of Newman (2006) and ' typmod])
        case 'Interactions'
            disp('  using eigenvector method of Newman (2006')
    end
    disp('  Individual   eig  Cluster')
    for j=1:rnumind
        disp(sprintf('%12s %7.4f %2.0f',prnam{ql(j)},clusval(ql(j)),clusass(ql(j))));
    end
    disp('  ')
    qq=clusass*ones(1,rnumind);moduf=B.*(qq==qq');moduf=sum(sum(moduf-diag(diag(moduf))))/(2*m);%moduf=sum(sum(B.*(qq==qq')))/(2*m);
    disp(sprintf('Modularity of this arrangement is %5.3f',moduf))

    global afields idat ridat atypfield rnam nam asf1
    butt=questdlg('Save Clusters as Supplemental Fields');
    if strcmp(butt,'Yes')
        if isempty(find(strcmp(afields,'Cluster')))
            defstr={'Cluster'};
        else
            defstr={''};
        end
        while 1
            namclus=inputdlg('Name of cluster field','',1,defstr);
            if isempty(find(strcmp(afields,namclus{1})))
                break
            else
                warndlg({[namclus{1} ' is already a supplemental field name--try again']});
                uiwait;
            end
        end
        asf1=asf1+1;
        afields{asf1}=namclus{1};
        atypfield(asf1)=1;
        ridat{asf1}=clusass;
        idat{asf1}=NaN*zeros(length(nam),1);
        for i=1:length(rnam)
            matid=find(strcmp(rnam{i},nam));
            idat{asf1}(matid)=clusass(i);
        end
        set(findobj('tag','uclass'),'string',namclus{1});
    end
end