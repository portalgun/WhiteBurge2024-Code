classdef ExtP < handle
properties
    dims
    Subjs
    model
    passes
    bNew
    fitType

    eModel
    M
    VarL
    nCov
    lCov

    covLB

    E
    O
    P
    D
    I
    Err

    binLabels
    binColors

    figNums
    FRC
    sp

    flds
    fldsX

    bLog
    oMargin
    iMargin
    oUnits
    iUnits

    Fit
    Boot
end
properties(Hidden)
    Fit_

    nM
    nFlds
    nFigs

    N
    sz0
    ndim0
    ndim
    sz

    Scale

    nSubj
    nBin

    ES
    clabels
    syms

    is
    ib
    im

    s

    sub0
    sub

    % FITTING
end
methods(Static)
    function P=getP()
        P={...
           'dims',[];
           'Subjs',[];
           'model',[];
           'passes',[];
           'bNew',[];
           'fitType', 0;
           ...
           'eModel', 'LmCR';
           'bLog',    0;
           'binNums', [];
           'flds',[];
           ...
           'oMargin',[2 1.5 3 4];
           'iMargin',[1 0 2 1];
           'oUnits','char';
           'iUnits','char';
           ...
           'nCov',1000;
           'lCov',[-4 4];
           'VarL',1;
            'M',[];
        };
    end
end
methods
    function obj=ExtP(dims,SUBJS,passes,model,bNew,varargin)
        EPFigs.parser(nargin,4);
        if isempty(SUBJS)
            SUBJS=Figs.getSUBJS(dims(1));
        end

        varargin=['dims',dims,'Subjs',{SUBJS},'passes',passes,'model',model,'bNew',bNew,varargin];
        Args.simple(obj,ExtP.getP(),varargin{:});

        % XXX
        obj.bNew=false;
        obj.fitType=0;
        obj.bLog=false;
        obj.binColors={'k','r'};
        obj.binLabels={'1','13'};

        if isempty(obj.M)
            switch obj.eModel
            case 'LmCR'
                obj.M=[.25 .5 1 1.5 2 2.5 3 4];
            case 'LMCR'
                obj.M=[.25 .5 1 1.5 2 2.5 3 4];
            end
        end
        obj.get_ext_data();
        if isempty(obj.flds)
            obj.get_default_flds();
        end

        % XXX
        obj.get_pred();

    end
%- GET
    function set.lCov(obj,val)
        obj.lCov=val;
        if ~isempty(obj.nCov)
            obj.covLB=linspace(obj.lCov(1),obj.lCov(2),obj.nCov);
        end
    end
    function set.nCov(obj,val)
        obj.nCov=val;
        if ~isempty(obj.lCov)
            obj.covLB=linspace(obj.lCov(1),obj.lCov(2),obj.nCov);
        end
    end

    function out=get.Fit(obj)
        if isempty(obj.Fit_)
            obj.Fit_=ExtPFit(obj);
        end
        out=obj.Fit_;
    end
%- SET
    function out=get.sp(obj)
        out=obj.FRC.sp;
    end
    function set.bLog(obj,in)
        obj.bLog=in;
        if obj.bLog
            obj.Scale='log';
        else
            obj.Scale='linear';
        end

    end
    function set.M(obj,in)
        obj.M=in;
        obj.nM=numel(obj.M);
    end
    function set.figNums(obj,in)
        obj.figNums=Vec.row(in);
        obj.nFigs=numel(obj.figNums);
    end
    function set.Subjs(obj,in)
        obj.Subjs=in;
        obj.nSubj=numel(obj.Subjs);


        E=ESubjs();
        ES=cell(obj.nSubj,1);
        syms=cell(obj.nSubj,1);
        for s =1:obj.nSubj
            ES{s}=E.(obj.Subjs{s});
            syms{s}=ES{s}.Marker;
        end
        obj.ES=ES;
        obj.syms=syms;
    end
    function set.flds(obj,in)
        obj.flds=in;
        obj.nFlds=numel(obj.flds);
        if isempty(obj.flds)
            return
        end
        obj.fldsX=regexprep(obj.flds,'[a-zA-Z]+','X');
    end
    function get_default_flds(obj)
        switch obj.eModel
        case 'LmCR'
            flds={'Rho23','Rho33'};
        case 'LMCR'
            flds={'Rho23'};
        end
        obj.flds=flds;
    end
%- DATA
    function get_ext_data(obj)
        E=EPs.ext_data(obj.dims,obj.Subjs,obj.passes,obj.model,obj.bNew,obj.fitType);

        obj.nBin= size(E.varLD,2);

        O=struct();
        O.Rho22=E.rho22;
        O.Rho33=E.rho33;
        O.Rho23=E.rho23;
        O.T3=E.T3;
        O.T2=E.T2;


        obj.sz0    =size(O.Rho22);

        obj.N      =numel(O.Rho22);
        obj.ndim0  =numel(obj.sz0);
        obj.ib     =obj.ndim0-1;
        obj.is     =obj.ndim0;
        obj.im     =obj.ndim0+1;

        obj.E=E;
        obj.O=O;

    end
    function get_pred(obj)
        obj.sz  =[obj.sz0 obj.nM];
        obj.ndim=numel(obj.sz);
        obj.sub0=substruct('()',repmat({':'},1,obj.ndim));

        %% INIT D & P
        z=zeros(obj.sz);

        P=struct();
        E=struct();
        D=struct();
        I=struct();
        for i = 1:obj.nFlds
            P.(obj.flds{i})=z;
            D.(obj.flds{i})=z;
            E.(obj.flds{i})=z;
            I.(obj.flds{i})=z;
        end
        for i = 1:obj.nFlds
            P.(obj.fldsX{i})=z;
        end

        O=obj.O;

        covLB=linspace(obj.lCov(1),obj.lCov(2),obj.nCov);

        sub=obj.sub0;
        clabels=cell(obj.nM,1);

        switch obj.eModel
        case 'LmCR'
            b23_33=true;
            b33_23=true;
        case 'LMCR'
            b23_33=false;
            b33_23=true;

        end

        for i = 1:obj.N
        for j = 1:obj.nM
            [sub.subs{1:obj.ndim0}]=ind2sub(obj.sz0,i);
            sub.subs{obj.im}=j;

            switch obj.eModel
            case {'LmCR','LMCR'}
                obji=Ext.LmCR(obj.VarL,obj.M(j),covLB,O.Rho22(i));
            end
            if i == 1
                clabels{j}=obji.getMTitle();
            end

            if b23_33
                [e33,x33,rho33p,rho33i]=obji.rho23_rho33_error(O.Rho23(i),0,true);
                P.X33  =subsasgn(P.X33,  sub, x33);

                P.Rho33=subsasgn(P.Rho33,sub, rho33p);
                I.Rho33=subsasgn(I.Rho33,sub, rho33i);

                E.Rho33=subsasgn(E.Rho33,sub, e33);

                D.Rho33=subsasgn(D.Rho33,sub, O.Rho33(i)-rho33p);
            end
            if b33_23
                [e23,x23,rho23p,rho23i]=obji.rho33_rho23_error(O.Rho33(i),0,true);
                P.X23  =subsasgn(P.X23,  sub, x23);

                P.Rho23=subsasgn(P.Rho23,sub, rho23p);
                I.Rho23=subsasgn(I.Rho23,sub, rho23i);

                E.Rho23=subsasgn(E.Rho23,sub, e23);

                D.Rho23=subsasgn(D.Rho23,sub, O.Rho23(i)-rho23p);
            end

        end
        end
        obj.clabels=clabels;
        obj.D=D;
        obj.P=P;
        obj.Err=E;
        obj.I=I;
    end
    function errors(obj)
        %% Ideal - Predicted
        inds=1:obj.nSubj;

        flds={'Rho23','Rho33'};
        E_IP=[];
        E_IO=[];

        nF=length(flds);
        E=cell(nF);
        for i = 1:nF
            O=obj.O.(flds{i})(:,:,inds);

            P=obj.P.(flds{i})(:,:,inds,:);

            I=obj.I.(flds{i})(:,:,inds,:);

            e_IP=squeeze(sum(abs(bsxfun(@minus,I,P)),[1,3]));
            e_IO=squeeze(sum(abs(bsxfun(@minus,P,O)),[1,3]));

            E_IP=cat(3,E_IP,e_IP);
            E_IO=cat(3,E_IO,e_IO);
        end
        E_IP=sum(E_IP,3);
        E_IO=sum(E_IO,3);

        dispV(efun(E_IP),'m_IP');
        dispV(efun(E_IO),'m_IO');

        function m=efun(E)
            for i = 1:size(E,1)
                mn(i)=min(E(i,:));
                ind(i)=find(E==mn(i));

                m(1,i)=obj.M(ind(i));
            end
        end


    end
%- FITTING
%- FIGURES
    function out=getSPArgs(obj,N,M)
        out={'iMargin',obj.iMargin, ...
                   'oMargin',obj.oMargin, ...
                   'iUnits',obj.iUnits, ...
                   'oUnits',obj.oUnits ...
        };
    end
    function init_SP(obj)
        args=obj.getSPArgs();
        RC=[obj.nFlds,obj.nM];

        for f = 1:numel(obj.figNums);
            names{f}=['ExtPredict ' num2str(obj.figNums(f)) ' ' num2str(obj.fitType) num2str(obj.VarL)];
        end
        obj.FRC=FRCs(obj,names,RC,args);
    end
%- SPEC SETTINGS
    function get_setting(obj,K)
        if obj.bLog
            lims=[.05 1];
        else
            lims=[-1 1];
        end
        s=struct();
        switch K
        case 0
            % SCATTER OBSERSERVED | PREDICTED
            s.indexMode='sb';
            s.plotMode='scatter';
            s.vPlotZLine=2;

            s.xlbl='observed';
            s.ylbl='predicted';
            s.xFld='i';
            s.yFld='p';

            s.xl=lims;
            s.yl=lims;
            s.YScale=obj.Scale;
            s.XScale=obj.Scale;
        case 1
            % SCATTER OBSERSERVED | PREDICTED
            s.indexMode='sb';
            s.plotMode='scatter';
            s.vPlotZLine=2;

            s.xlbl='observed';
            s.ylbl='predicted';
            s.xFld='o';
            s.yFld='p';

            s.xl=lims;
            s.yl=lims;
            s.YScale=obj.Scale;
            s.XScale=obj.Scale;
        case 2
            % SCATTER COV(CORR) |  PREDICTED
            s.indexMode='sb';
            s.plotMode='scatter';
            s.vPlotZLine=2;
            s.xFld='x';
            s.yFld='p';

            s.xlbl='corr(L,B)';
            s.ylbl='\rho';
            s.xl=[-1 1];
            s.yl=lims;

            s.YScale=obj.Scale;
            s.XScale='linear';
        case 3
            % HIST RHO Observed - PREDICTED
            s.indexMode='';
            s.plotMode='hist';
            s.vPlotZLine=-1;
            s.edges=linspace(-1,1,10);
            s.xFld='d';
            s.yFld='';

            s.xlbl='Observed - Predicted \rho';
            s.ylbl='count';
            s.xl=[-1 1];
            s.yl=[];

            s.YScale=obj.Scale;
            s.XScale='linear';
        case 4
            % HIST COV(CORR)
            s.indexMode='';
            s.plotMode='hist';
            s.vPlotZLine=-1;
            s.edges=linspace(-1,1,10);
            s.xFld='x';
            s.yFld='';

            s.xlbl='corr(L,B)';
            s.ylbl='count';
            s.xl=[-1 1];
            s.yl=[];

            s.YScale=obj.Scale;
            s.XScale='linear';
        end
        obj.s=s;
    end
%- PLOT
    function plotMo(obj);
        MO=obj.E.varB./obj.E.varLD;
        sA=spArgs;
        sA{1}=[1 2];
        SP(end+1)=SubPlot(sA{:});

        x1=MO(:,1,1:nSubj);
        x2=MO(:,2,1:nSubj);
        SP(end).sel([1 1]);
        hist(x1(:));
        SP(end).sel([1 2]);
        hist(x2(:));
        dispV([mean(x1(:)) median(x1(:)) exp(mean(log(x1(:)))) ],'Rho22_1')
        dispV([mean(x2(:)) median(x2(:)) exp(mean(log(x2(:)))) ],'Rho22_2')
    end
    function plot(obj,figNums)
        if nargin < 1 || isempty(figNums)
            figNums=1:4;
        end
        obj.figNums=figNums;

        obj.init_SP();
        if isempty(obj.flds)
            obj.get_default_flds();
        end

        for f=1:obj.nFigs
            obj.sub=obj.sub0;

            obj.get_setting(obj.figNums(f));
            plotMeth=['plot_' obj.s.indexMode];

            for j=1:obj.nM
                obj.sub.subs{obj.im}=j;
                for i=1:obj.nFlds

                    obj.FRC.sel(f,i,j);
                    hold off;

                    obj.s.fldx=obj.fldsX{i};
                    obj.s.fld =obj.flds{i};

                    obj.(plotMeth)();

                    if isempty(obj.s.yl)
                        yl=ylim();
                    else
                        yl=obj.s.yl;
                    end
                    switch obj.s.vPlotZLine
                    case -1
                        plot([0 0],   yl,'k:');
                    case 2
                        plot(obj.s.xl,yl,'k:');
                    end

                    obj.format(f);
                end
            end
            obj.fformat(f);
        end
    end
%- SPEC PLOT
    function varargout=get_v(obj,sub)
        varargout=cell(1,2);
        flds={obj.s.xFld, obj.s.yFld};
        S=obj.s;
        for i = 1:2
            switch flds{i}
            case 'd'
                varargout{i}=subsref(obj.D.(S.fld), sub);
            case 'e'
                varargout{i}=subsref(obj.Err.(S.fld), sub);
            case 'i'
                varargout{i}=subsref(obj.I.(S.fld), sub);
            case 'p'
                varargout{i}=subsref(obj.P.(S.fld ),sub);

            case 'x'
                varargout{i}=subsref(obj.P.(S.fldx),sub);
            case 'o'
                %obj.O.(fld)(:,b,u);
                s=sub;
                s.subs=s.subs(1:obj.ndim0);
                varargout{i}=subsref(obj.O.(S.fld), s);
            case ''
                varargout{i}=S.edges;
            end
        end
    end
    function plot_(obj,X,Y,color,marker)
        if nargin < 2
            [X,Y]=obj.get_v(obj.sub);
            if nargin < 4
                color=[0 0 0];
                if nargin < 5
                    marker='o';
                end
            end
        end
        switch obj.s.plotMode;
        case 'hist'
            obj.hist_(X,Y,color);
        case 'scatter'
            obj.scatter_(X,Y,color,marker);
        end
        hold on;
    end
    function plot_b(obj)
        S  =obj.s;
        sub=obj.sub;
        for b = 1:obj.nBin
            sub.subs{obj.ib}=b;

            [X,Y]=obj.get_v(sub);
            marker='o';
            color=obj.binColors{b};

            obj.plot_(X,Y,color,marker);
        end
    end
    function plot_sb(obj)
        S  =obj.s;
        sub=obj.sub;
        for s = 1:obj.nSubj
        for b = 1:obj.nBin
            sub.subs{obj.is}=s;
            sub.subs{obj.ib}=b;

            [X,Y]=obj.get_v(sub);
            marker=obj.syms{s};
            color=obj.binColors{b};

            obj.plot_(X,Y,color,marker);

        end
        end
    end
    function scatter_(obj,X,Y,color,marker);
        plot(X(:), Y(:), marker,   ...
            'MarkerFaceColor',color,...
            'MarkerEdgeColor','none'...
        ); hold on
    end
    function edges=hist_(obj,data,edges,color)
        if nargin < 3
            color='k';
        end
        if isempty(edges)
            [counts,edges]=histcounts(data);
        else
            [counts,edges]=histcounts(data,edges);
        end
        [ctrs,widths] =Hist.edges2ctrs(edges);
        plot(ctrs,counts,'-','Color',color); hold on
    end
%- FORMAT
    function format(obj,f)
        S=obj.s;
        axis square;
        if ~isempty(S.xl)
            xlim(S.xl);
        end
        if ~isempty(S.yl)
            Axis.ylim(S.yl);
        end
        obj.sp.xlabel(S.xlbl,true,[]);
        obj.sp.ylabel(S.ylbl,true,[],[]);
        set(gca,'YScale',S.YScale,'XScale',S.XScale);
        obj.sp.xticks();
        obj.sp.yticks();
        Axis.format();
    end
    function fformat(obj,f)
        obj.sp.rlabel(obj.flds);
        obj.sp.clabel(obj.clabels);
    end
end
end
