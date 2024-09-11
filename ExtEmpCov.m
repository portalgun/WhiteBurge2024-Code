classdef ExtEmpCov < handle
properties
    moude
    bAll
    loc
    bCorr
    bBar

    AColor
    MColor
    ANColor
    MNColor
    bNorm
    bA
    bM
    bAN
    bMN

    cv
    edges
    bins
    counts
    locs
    ylocs
    ymax

    N
    lineWidth
end
properties(Hidden)
    parent
end
methods(Static)
    function P=getP()
        P={...
            'lineWidth',2;
            'bAll',true;
            'loc','median';
            'bCorr',false;
            'bBar',true;
            'moude','all';
        };

    end
end
methods
    function obj=ExtEmpCov(parent,varargin)
        obj.parent=parent;
        obj.parse(varargin{:});
    end
    function parse(obj,varargin)
        P=obj.getP();
        Args.simple(obj,P,varargin{:});
        Args.empty2Def(P,obj);

    end
    function set.moude(obj,moude)
        switch moude
        case 'all'
            obj.bA=true;
            obj.bM=false;

            obj.bAN=true;
            obj.bMN=false;

            obj.AColor='k';
            obj.ANColor={'b','r'};
            obj.bNorm=false;
        case 'comb'
            obj.bA=false;
            obj.bM=true;

            obj.bAN=false;
            obj.bMN=true;

            obj.MColor='k';
            obj.MNColor={'b','r'};
            obj.bNorm=false;
        case 'cmp'
            obj.bA=true;
            obj.bM=true;

            obj.bAN=false;
            obj.bMN=false;

            obj.AColor='k';
            obj.MColor='m';
            obj.bNorm=true;
        case {'cmpDiff','CmpDiff'}
            obj.bA=false;
            obj.bM=false;

            obj.bAN=true;
            obj.bMN=true;

            obj.ANColor={'b','r'}; % XXX
            obj.MNColor={'g','m'}; % XXX
            obj.bNorm=true;
        otherwise
            error('invalid mode')
        end
        obj.moude=moude;
    end
    function cv=get.cv(obj)
        cv=obj.parent.M.covLB;
        if obj.bCorr
            varB=obj.parent.M.varB;
            varLD=obj.parent.M.varLD;
            nInd=(varB <= 0);
            varB(nInd)=0;
            cv=cv./sqrt(varB .* varLD);
            cv(nInd)=nan;
        end
    end
    function edges=get.edges(obj)
        if obj.bCorr
            edges={linspace(-1,1,10)};
        else
            edges={};

        end
        cl=obj.parent.tmpSet('bAvg',false);
        SInds=obj.parent.sInds;
        KInds=obj.parent.kInds;
        [~,edges]=histcounts(obj.cv(:,:,SInds,KInds),edges{:});
        d=edges-1;
        if obj.bCorr
            d=d(find(abs(d)==min(abs(d)),1,'first'));
            edges=edges-d;
            edges=Hist.padEdges(edges,2);
        else
            edges=Hist.padEdges(edges,1);
        end
        obj.bins=Hist.edges2ctrs(edges)';
    end
    function N=get.N(obj)
        if obj.bAll
            N=1;
        else
            N=obj.parent.nStd;
        end
    end
    function get_counts(obj);
        cv=obj.cv;
        edges=obj.edges;

        cl=obj.parent.tmpSet('bAvg',false);
        SInds=obj.parent.sInds;
        KInds=obj.parent.kInds;
        n    =obj.parent.nStd;
        N=obj.N;
        nB=numel(edges)-1;

        switch obj.loc
        case 'mean'
            locFun=@(x,y) nanmean(x);
        case 'median'
            locFun=@(x,y) nanmedian(x);
        case 'mode'
            locFun=@(x,y) Hist.mode(x,y);
        end

        z=zeros(nB,N);
        c0=z;
        c1=z;
        c2=z;

        cm0=z;
        cm1=z;
        cm2=z;

        z=zeros(N,1);
        m0=z;
        m1=z;
        m2=z;

        mc0=z;
        mc1=z;
        mc2=z;

        for i = 1:N
            if obj.bAll
                inds=1:n;
            else
                inds=i;
            end

            cvm0=cv(inds,:,SInds,1);
            cvm1=cv(inds,1,SInds,1);
            cvm2=cv(inds,2,SInds,1);

            cv0=cv(inds,:,SInds,KInds);
            cv1=cv(inds,1,SInds,KInds);
            cv2=cv(inds,2,SInds,KInds);

            % COUNTS
            [c0(:,i)]=histcounts(cv0(:),edges);
            [c1(:,i)]=histcounts(cv1(:),edges);
            [c2(:,i)]=histcounts(cv2(:),edges);

            [cm0(:,i)]=histcounts(cvm0(:),edges);
            [cm1(:,i)]=histcounts(cvm1(:),edges);
            [cm2(:,i)]=histcounts(cvm2(:),edges);
            if obj.bNorm
                c0(:,i)=c0(:,i)./sum(c0(:,i));
                c1(:,i)=c1(:,i)./sum(c1(:,i));
                c2(:,i)=c2(:,i)./sum(c2(:,i));

                cm0(:,i) =cm0(:,i)./sum(cm0(:,i));
                cm1(:,i) =cm1(:,i)./sum(cm1(:,i));
                cm2(:,i) =cm2(:,i)./sum(cm2(:,i));
            end

            % LOC PARAM
            m0(i)=locFun(cv0(:),edges);
            m1(i)=locFun(cv1(:),edges);
            m2(i)=locFun(cv2(:),edges);

            mc0(i)=locFun(cvm0(:),edges);
            mc1(i)=locFun(cvm1(:),edges);
            mc2(i)=locFun(cvm2(:),edges);
        end

        C=struct();
        C.c0=c0;
        C.c1=c1;
        C.c2=c2;

        C.cm0=cm0;
        C.cm1=cm1;
        C.cm2=cm2;

        M=struct();
        M.m0=m0;
        M.m1=m1;
        M.m2=m2;

        M.mc0=mc0;
        M.mc1=mc1;
        M.mc2=mc2;

        obj.counts=C;
        obj.locs=M;

        ymx=obj.ymax;
        if obj.bNorm
            if obj.bA && obj.bM
                c0mx=max( c0,[],1);
                cmmx=max(cm0,[],1);
                c0 = c0./c0mx.*ymx;
                cm0=cm0./cmmx.*ymx;
            elseif obj.bAN && obj.bAN
                cm1mx=max(cm1,[],1);
                cm2mx=max(cm2,[],1);
                c1mx=max(c1,[],1);
                c2mx=max(c2,[],1);

                cm1=cm1./cm1mx.*ymx;
                cm2=cm2./cm2mx.*ymx;
                c1=c1./c1mx.*ymx;
                c2=c2./c2mx.*ymx;
            end
        end
        C.c0=c0;
        C.c1=c1;
        C.c2=c2;

        C.cm0=cm0;
        C.cm1=cm1;
        C.cm2=cm2;
    end
    function out=get.ylocs(obj)
        N=obj.N;
        z=zeros(N,1);
        bins=obj.bins;
        C=obj.counts;
        M=obj.locs;

        if obj.bBar
            intMeth='nearest';
        else
            intMeth='linear';
        end

        y0=z;
        y1=z;
        y2=z;

        ym0=z;
        ym1=z;
        ym2=z;
        % Y INTERP FOR MEAN
        for i=1:N
            y0(i) =interp1(bins(:),C.c0(:,i),  M.m0(i),intMeth);
            y1(i) =interp1(bins(:),C.c1(:,i),  M.m1(i),intMeth);
            y2(i) =interp1(bins(:),C.c2(:,i),  M.m2(i),intMeth);

            ym0(i)=interp1(bins(:),C.cm0(:,i), M.mc0(i),intMeth);
            ym1(i)=interp1(bins(:),C.cm1(:,i), M.mc1(i),intMeth);
            ym2(i)=interp1(bins(:),C.cm2(:,i), M.mc2(i),intMeth);
        end
        y=struct();

        y.y0=y0;
        y.y1=y1;
        y.y2=y2;
        y.ym0=ym0;
        y.ym1=ym1;
        y.ym2=ym2;

        out=y;
    end
    function out=get.ymax(obj)
        C=obj.counts;
        if obj.bA && obj.bM
            ymx=max([C.c0(:); C.cm0(:)]);
        elseif obj.bAN && obj.bMN
            ymx=max([C.cm1(:); C.cm2(:); C.c1(:); C.c2(:)]);
        elseif obj.bA
            ymx=C.c0(:);
        elseif obj.bM
            ymx=C.cm(:);
        elseif obj.bAN
            ymx=[C.c1(:); C.c2(:)];
        elseif obj.bMN
            ymx=[C.cm1(:); C.cm2(:)];
        end
        out=max(ymx(:))*1.05;
    end
    function plot(obj)
        if obj.bAll
            obj.plot2D();
        else
            obj.plot3D();
        end

    end
    function newFig(obj)
        moude=[upper(obj.moude(1)) obj.moude(2:end)];
        moude
        obj.parent.newFig(['histCov' moude num2str(obj.bAll) ]);
        hold off;
    end
    function plot2D(obj)
        obj.bAll=true;
        obj.get_counts();

        edges=obj.edges;
        bins=obj.bins;
        if obj.bBar
            plotfun=@(x,l,c) Plot.bar([],edges,x,'LineWidth',obj.lineWidth,'LineStyle',l,'Color',c);
        else
            plotfun=@(x,l,c) plot(bins,x,'LineWidth',obj.lineWidth,'LineStyle',l,'Color',c);
        end
        locfun=@(x,y, l,c) plot([x x],[0 y],'LineStyle',l,'Color',c);
        C=obj.counts;
        x=obj.locs;
        y=obj.ylocs;

        obj.newFig();

        % HIST
        lStyle='-';
        if obj.bA
            plotfun(C.c0,lStyle,obj.AColor);  hold on
        end
        if obj.bAN
            plotfun(C.c1,lStyle,obj.ANColor{1}); hold on
            plotfun(C.c2,lStyle,obj.ANColor{2});
        end

        if obj.bM
            plotfun(C.cm,lStyle,obj.MColor); hold on
        end
        if obj.bMN
            plotfun(C.cm1,lStyle,obj.ANColor{1}); hold on
            plotfun(C.cm2,lStyle,obj.ANColor{2});
        end

        % LINES
        lStyle='--';
        if obj.bA
            locfun(x.m0,y.y0,lStyle,obj.AColor);
        end
        if obj.bAN
            locfun(x.m1,y.y1,lStyle,obj.ANColor{1});
            locfun(x.m2,y.y2,lStyle,obj.ANColor{2});
        end

        if obj.bM
            locfun(x.mc0,y.ym0,lStyle,obj.MColor);
        end
        if obj.bMN
            locfun(x.mc1,y.ym1,lStyle,obj.MNColor{1});
            locfun(x.mc2,y.ym2,lStyle,obj.MNColor{2});
        end
        obj.format();
    end
    function plot3D(obj)
        obj.bAll=false;

        obj.get_counts();
        edges=obj.edges;
        bins=obj.bins;
        N=obj.N;
        nB=numel(edges)-1;
        C=obj.counts;
        M=obj.locs;
        Y=obj.ylocs;

        X=permute(obj.parent.E.X,[3 1 2]);
        X=repmat(X,nB,1,1);

        if obj.bBar
            plotfun=@(x,z,c) obj.bar_fun(x,edges,z,c);
        else
            plotfun=@(x,z,c) plot3(x,bins,z,'LineWidth',2,'LineStyle','-','Color',c);
        end

        obj.newFig();

        % Counts
        for i = 1:N
            x=X(:,i);
            if obj.bA
                plotfun(x,C.c0(:,i),obj.AColor); hold on
            end
            if obj.bAN
                plotfun(x,C.c1(:,i),obj.ANColor{1}); hold on
                plotfun(x,C.c2(:,i),obj.ANColor{2});
            end

            if obj.bM
                plotfun(x,C.cm0(:,i),obj.MColor); hold on
            end
            if obj.bMN
                plotfun(x,C.cm1(:,i),obj.MNColor{1}); hold on
                plotfun(x,C.cm2(:,i),obj.MNColor{2});
            end
        end

        args={'LineStyle','--','LineWidth',2,'Color'};

        % Locs
        x=X(1,:);
        if obj.bA
            plot3(x,M.m0,Y.y0,args{:},obj.AColor);
        end
        if obj.bAN
            plot3(x,M.m1,Y.y1,args{:},obj.ANColor{1});
            plot3(x,M.m2,Y.y2,args{:},obj.ANColor{2});
        end

        if obj.bM
            plot3(x,M.mc,Y.ym0,args{:},obj.MColor);
        end
        if obj.bMN
            plot3(x,M.mc1,Y.ym1,args{:},obj.MNColor{1});
            plot3(x,M.mc2,Y.ym2,args{:},obj.MNColor{2});
        end
        obj.format();
    end
    function bar_fun(obj,x,edges,z,c)
        z=Vec.row(z);

        X=[x; x];
        Y=[edges(1:end-1); edges(2:end)];
        Z=[z; z];

        Y=Y(:);
        X=X(:);
        Z=Z(:);
        plot3(X,Y,Z,'LineWidth',2,'LineStyle','-','Color',c);
    end
    function label(obj)
        stdlbl='Disparity (arcmin)';
        cntlbl='count';
        if obj.bCorr
            covlbl='Corr(L,B)';
        else
            covlbl='Cov(L,B)';
        end

        if obj.bAll
            xlabel(covlbl);
            ylabel(cntlbl);
        else
            xlabel(stdlbl);
            ylabel(covlbl);
            zlabel(cntlbl);
        end
    end
    function format(obj)
        obj.lim();
        if obj.bAll
            Axis.format();
            axis square;
        else
            set(gca,'YDir','Reverse');
        end
        obj.label();
    end
    function lim(obj)
        if obj.bAll
            ylim([0 obj.ymax]);
        end
    end
end
end
