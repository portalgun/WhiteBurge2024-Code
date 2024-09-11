classdef ExtPFit < handle
properties
    VarL
    nCov
    lCov
    covLB

    nSubj
    sz0
    nBin
    N

    O
    E
    P

    x1
    x2
    e1
    e2
    r1p
    r2p
    r1i
    r2i
    subjInds

    typ
    All
    Subj
    Each
    Obs

    typs={'F0','C0X','C0Y','PO','IO','DX','DY'};
    expn=2

    syms
    binColors
end
properties(Access=protected)
    z
    parent
end
methods
%- CON
    function obj=ExtPFit(parent)
        obj.parent=parent;
    end
%- GET
    function out=get.sz0(obj)
        out=obj.parent.sz0;
    end
    function out=get.N(obj)
        out=obj.parent.N;
    end
    function out=get.nBin(obj)
        out=obj.parent.nBin;
    end
    function out=get.nSubj(obj)
        out=obj.parent.nSubj();
    end
    function out=get.covLB(obj)
        out=obj.parent.covLB;
    end
    function out=get.VarL(obj)
        out=obj.parent.VarL;
    end
    function out=get.O(obj)
        out=obj.parent.O;
    end
    function out=get.P(obj)
        out=obj.parent.P;
    end
    function out=get.E(obj)
        out=obj.parent.E;
    end
    function out=get.syms(obj)
        out=obj.parent.syms;
    end
    function out=get.binColors(obj)
        out=obj.parent.binColors;
    end
    function varargout=hist_(obj,varargin)
        varargout=cell(1,nargout);
        [varargout{:}]=obj.parent.hist_(varargin{:});
    end
    function scatter_(obj,varargin)
        obj.parent.scatter_(varargin{:});
    end
%- UTIL
    function C=m2C(obj,m,num,bCorr)
        if nargin < 4
            bCorr=false;
        end
        C=zeros(size(m));

        r22=obj.O.Rho22;
        r23=obj.O.Rho23;
        r33=obj.O.Rho33;


        for i = 1:numel(m)
            if m(i) <= 0
                C(i)=nan;
                continue
            end
            obji=Ext.LmCR(obj.E.varLD(i),m(i),obj.covLB,r22(i));
            if num==23
                c=obji.rho33_rho23(r33(i),bCorr,false);
            elseif num==33
                c=obji.rho23_rho33(r23(i),bCorr,false);
            end
            if ~isreal(c)
                c=nan;
            end
            C(i)=c;
        end
    end
    function R=m2Rho(obj,m,num)
        R=zeros(size(m));
        bCorr=true;

        r22=obj.O.Rho22;
        r23=obj.O.Rho23;
        r33=obj.O.Rho33;

        for i = 1:numel(m)
            if m(i) <= 0
                R(i)=nan;
                continue
            end
            obji=Ext.LmCR(obj.VarL,m(i),obj.covLB,r22(i));
            if num==23
                [~,r]=obji.rho33_rho23(r33(i),bCorr,true);
            elseif num==33
                [~,r]=obji.rho23_rho33(r23(i),bCorr,true);
            end
            if ~isreal(r)
                r=nan;
            end
            R(i)=r;
        end
    end
    function out=index(obj,in,bin,subj)
        if nargin < 3 || isempty(bin)
            bin=1:obj.nBin;
        end
        if nargin < 4 || isempty(subj)
            subj=1:obj.nSubj;
        end
        in=in(:,bin,subj);
        out=in(:);
    end
%- FIT HELPERS
    function X=fit(obj,fun,con,typ)
        if nargin < 3
            con=[];
        end

        A.UB=[];
        A.LB=[];
        A.A=[];
        A.B=[];
        A.Aeq=[];
        A.Beq=[];
        switch typ
        case 'LMCR'
            A.UB=10;
            A.LB=1/10;
        case 'LmCR'
            A.UB=10;
            A.LB=1/A.UB;
        end

        %dis='iter';
        dis='off';
        Opts=optimoptions('fmincon','Display',dis);
        f0=(rand(1)*(A.UB-A.LB))+A.LB;
        X = fmincon(fun,f0, A.A, A.B, A.Aeq, A.Beq, A.LB, A.UB,con,Opts);
    end
    function init(obj,name,Typ,typs)
        if nargin < 4 || isempty(typs)
            typs=obj.typs;
        end
        obj.z=zeros(obj.sz0);
        flds={'x1','x2','e1','e2','r1p','r2p','r1i','r2i'};
        for i = 1:length(flds)
            obj.(flds{i})=obj.z;
        end
        name(1)=upper(name(1));
        obj.(name)=struct();
        for t =1:length(typs)
            typ=typs{t};
            obj.(name).(typ)=obj.z;
        end

        obj.typ=Typ;
        obj.([Typ '_obs']);
    end
%-  FIT
    function LmCR(obj)
        obj.LmCR_obs();
        obj.LmCR_all();
        obj.LmCR_subj();
        obj.LmCR_each();
    end
    function LmCR_all(obj)
        obj.init('all','m');
        obj.typ='m';
        for t = 1:numel(obj.typs);
        for b = 1:obj.nBin
            obj.subjInds=[1:obj.nSubj];
            typ=obj.typs{t};

            fun=@(m) obj.objective_LmCR(m,b,typ,obj.expn);
            M=obj.fit(fun,'LmCR');
            obj.All.(typ)(:,b,:)=M;
        end
        end
    end
    function LmCR_subj(obj,Subj)
        obj.init('subj','m');
        for s = 1:obj.nSubj
        for t = 1:numel(obj.typs);
        for b = 1:obj.nBin
            obj.subjInds=s;
            typ=obj.typs{t};

            fun=@(m) obj.objective_LmCR(m,b,typ,obj.expn);
            M=obj.fit(fun,con,'LmCR');
            obj.Subj.(typ)(:,b,s)=M;
        end
        end
        end
    end
    function LmCR_each(obj)
        typs=obj.typs;
        obj.init('each','m',typs);

        for t = 1:numel(typs)
        for i = 1:obj.N
            typ=typs{t};

            con=@(m) obj.con_each_LmCR(m,i);
            fun=@(m) obj.objective_each_LmCR(m,i,typ,obj.expn);
            M=obj.fit(fun,con,'LmCR');
            obj.Each.(typ)(i)=M;
        end
        end
        obj.plotEach();

    end
    function LMCR_each(obj)
        typs={'F0'};
        obj.init('each','M',typs);

        for t = 1:numel(typs)
        for i = 1:obj.N
            typ=typs{t};

            %con=@(m) obj.con_each_LMCR(m,i);
            con=[];
            fun=@(m) obj.objective_each_LMCR(m,i,typ,obj.expn);
            M=obj.fit(fun,con,'LMCR');
            obj.Each.(typ)(i)=M;
        end
        end
        obj.plotEach('M');

    end
    function plotEach(obj,stats)
        if nargin < 2 || isempty(stats)
            stats={'M','C','R'};
        elseif ~iscell(stats)
            stats={stats};
        end
        models=fieldnames(obj.Each);
        RC=[numel(stats) numel(models)];
        Fig.new(['ExtPFit plotEach']);
        for i =1:RC(1)
        for j =1:RC(2)
            subPlot(RC,i,j);
            obj.plotEach_(models{j},stats{i},[],false);
            obj.plotEach_(models{j},stats{i},[],false);
        end
        end
    end
    function m_obs(obj)
        VarL=obj.parent.E.varLD;
        VarB=obj.parent.E.varB;

        obj.Obs=VarB./VarL;
    end
    function M_obs(obj)
        obj.Obs=obj.O.T3./obj.O.T2;;
    end
%- OBJECTIVES
    function [c,ceq]=con_each_LMCR(obj,m,ind)
        T=obj.O.T3(ind)./obj.O.T2(ind);
        obji=Ext.LMCR(obj.VarL,m,obj.covLB,T);
        [x1,r1]=obji.rho33_rho23(obj.O.Rho33(ind),true,true);
        c=[...
           x1-1 ...
          -x1-1 ...
           r1-1 ...
          -r1   ...
        ];
        %c=[...
        %  -x1-1, -x2-1,  ...
        %
        %  -r1-1, -r2-1, ...
        %];
        ceq=[];
    end
    function [c,ceq]=con_each_LmCR(obj,m,ind)
        obji=Ext.LmCR(obj.VarL,m,obj.covLB,obj.O.Rho22(ind));
        [x1,r1]=obji.rho33_rho23(obj.O.Rho33(ind),true,true);
        [x2,r2]=obji.rho23_rho33(obj.O.Rho23(ind),true,true);
        c=[...
           x1-1   x2-1  ...
          -x1-1  -x2-1  ...
           r1-1   r2-1  ...
          -r1    -r2  ...
           r1-r2 ...
        ];
        %c=[...
        %  -x1-1, -x2-1,  ...
        %
        %  -r1-1, -r2-1, ...
        %];
        ceq=[];
    end
    function out=objective_each_LMCR(obj,M,ind,typ,expn)
        x23i=0;
        x33i=0;
        bCorr=false;

        T=obj.O.T3(ind)./obj.O.T2(ind);
        r22=obj.O.Rho22(ind);
        r23=obj.O.Rho23(ind);
        r33=obj.O.Rho33(ind);

        obji=Ext.LMCR(obj.VarL,M,obj.covLB,r22);

        %[x23p,r23p]=obji.rho33_rho23(r33,bCorr,true);
        %r23i=obji.rho_interp_(x23i,'rho23',bCorr);

        switch typ
        case 'F0'
            e1=abs(r33-obji.rho33);
        case 'C0X' % X from 0
            e1=abs(x23p-x23i);
        case 'C0Y' % Y from 0
            e1=abs(r23p-r23i);
        case 'PO' % Y
            e1=abs(r23-r23p);
        case 'IO' % Y
            e1=abs(r23-r23i);
        case 'DX' % X
            out=abs(x23p-x33p);
        case 'DY' % Y
            out=abs(r33p-r23p);
        end
        out=e1.^expn;
    end
    function out=objective_each_LmCR(obj,m,ind,typ,expn)
        x23i=0;
        x33i=0;
        bCorr=false;

        r22=obj.O.Rho22(ind);
        r23=obj.O.Rho23(ind);
        r33=obj.O.Rho33(ind);
        varL=obj.E.varLD(ind);

        obji=Ext.LmCR(varL,m,obj.covLB,r22);

        if strcmp(typ,'F0')
            VarB=obj.parent.E.varB(ind);
            out=(VarB-obji.varB).^2;
            return
        end

        [x23p,r23p]=obji.rho33_rho23(r33,bCorr,true);
        [x33p,r33p]=obji.rho23_rho33(r23,bCorr,true);

        r23i=obji.rho_interp_(x23i,'rho23',bCorr);
        r33i=obji.rho_interp_(x33i,'rho33',bCorr);

        switch typ
        case 'C0X' % X from 0
            e1=abs(x23p-x23i);
            e2=abs(x33p-x33i);
        case 'C0Y' % Y from 0
            e1=abs(r23p-r23i);
            e2=abs(r33p-r33i);
        case 'PO' % Y
            e1=abs(r23-r23p);
            e2=abs(r33-r33p);
        case 'IO' % Y
            e1=abs(r23-r23i);
            e2=abs(r33-r33i);
        case 'DX' % X
            out=abs(x23p-x33p).^expn;
            return
        case 'DY' % Y
            out=abs(r33p-r23p).^expn;
            return
        end
        out=e1.^expn + e2.^expn;
    end
    function out=objective_LmCR(obj,m,binNum,typ,expn)
        % C0Y - distance between predicted and rho0
        % PO - distance between observed and predicted
        % rhop =
        if m >= 0
            out=10^10;
        end
        x=zeros(obj.N,1);
        x1i=0;
        x2i=0;
        for i = 1:obj.N
            obji=Ext.LmCR(obj.VarL,m,obj.covLB,obj.O.Rho22(i));
            [~,obj.x1(i), obj.r1p(i), obj.r1i(i)]=obji.rho33_rho23_error(obj.O.Rho33(i),x1i,true);
            [~,obj.x2(i), obj.r2p(i), obj.r2i(i)]=obji.rho23_rho33_error(obj.O.Rho23(i),x2i,true);
        end
        r1o=obj.O.Rho23;
        r2o=obj.O.Rho33;
        switch typ
        case 'C0X'  % X from 0
            e1=abs(obj.x1-x1i);
            e2=abs(obj.x2-x2i);
        case 'C0Y' % Y from 0
            e1=abs(obj.r1p-obj.r1i);
            e2=abs(obj.r2p-obj.r2i);
        case 'PO' % Y DATA
            e1=abs(r1o-obj.r1p);
            e2=abs(r2o-obj.r2p);
        case 'IO' % Y DATA
            e1=abs(r1o-obj.r1i);
            e2=abs(r2o-obj.r2i);
        case 'DX' % X
            out=sum(abs(obj.x2(:,binNum,obj.subjInds)-obj.x1(:,binNum,obj.subjInds)).^expn,[1 3]);
            return
        case 'DY' % Y
            out=sum(abs(obj.r2p(:,binNum,obj.subjInds)-obj.r1p(:,binNum,obj.subjInds)).^expn,[1 3]);
            return
        case 'DX*'
            %a=mean(cat(4,obj.x1,obj.x2),4);
            %b=mean(cat(4,r1o,2o),4);
            %out=sum(abs(a-b).^expn,[1 3]);
            return
        case 'DY*'
            %a=mean(cat(4,obj.r21p,obj.r2p),4);
            %b=mean(cat(4,r1o,2o),4);
            %out=sum(abs(a-b).^expn,[1 3]);
            return
        end
        e1=sum(e1(:,binNum,obj.subjInds).^expn,[1 3]);
        e2=sum(e2(:,binNum,obj.subjInds).^expn,[1 3]);
        out=e1+e2;
    end
%- PLOT UTIL
    function lims=getLimsDiag(obj,X,Y,bLog)
        if nargin < 3 || isempty(bLog)
            bLog=false;
        end
        X=obj.index(X);
        Y=obj.index(Y);
        if bLog
            X(X<=0)=[];
            Y(Y<=0)=[];
        end
        lims=Num.minMax([X; Y]);
    end
%- PLOT EACH
    %- C0Y
    function plotEachM_C0X(obj,bLog)
        if nargin < 2
            bLog=[];
        end
        obj.plotEach_('C0X','M',bLog);
    end
    function plotEachM_C0Y(obj,bLog)
        if nargin < 2
            bLog=[];
        end
        obj.plotEach_('C0Y','M',bLog);
    end
    %- C
    function plotEachC_C0X(obj,bLog)
        if nargin < 2
            bLog=[];
        end
        obj.plotEach_('C0X','C',bLog);
    end
    function plotEachC_C0Y(obj,bLog)
        if nargin < 2
            bLog=[];
        end
        obj.plotEach_('C0Y','C',bLog);
    end
    %- R
    function plotEachR_C0X(obj,bLog)
        if nargin < 2
            bLog=[];
        end
        obj.plotEach_('C0X','R',bLog);
    end
    function plotEachR_C0Y(obj,bLog)
        if nargin < 2
            bLog=[];
        end
        obj.plotEach_('C0Y','R',bLog);
    end

    %- Main
    function plotEach_(obj,fld,stat,bLog,bFig)
        if stat=='R'
            vstr='rho';
        else
            vstr=stat;
        end
        if nargin < 4 || isempty(bLog)
            bLog=stat~='C';
        end
        if nargin < 5 || isempty(bFig)
            bFig=true;
        end
        if bFig
            Fig.new(['ExtPFit plotEach ' vstr ' ' fld]);
        end
        hold off;

        Xr=obj.Obs;
        Yr=obj.Each.(fld);


        m=obj.typ;
        switch stat
        case 'C'
            C=cell(1,2);
            X{1}=obj.([m '2C'])(Xr,23);
            X{2}=obj.([m '2C'])(Xr,33);

            Y{1}=obj.([m '2C'])(Yr,23);
            Y{2}=obj.([m '2C'])(Yr,33);
        case 'R'
            X{1}=obj.([m '2Rho'])(Xr,23);
            X{2}=obj.([m '2Rho'])(Xr,33);

            Y{1}=obj.([m '2Rho'])(Yr,23);
            Y{2}=obj.([m '2Rho'])(Yr,33);
        case 'M'
            X{1}=Xr;
            Y{1}=Yr;
        end

        lims=obj.getLimsDiag(vertcat(X{:}),vertcat(Y{:}),bLog);
        for i = 1:numel(X)
        for b = 1:obj.nBin
        for s = 1:obj.nSubj
            x=obj.index(X{i},b,s);
            y=obj.index(Y{i},b,s);

            marker=obj.syms{s};
            color=obj.binColors{b};
            obj.scatter_(x,y,color,marker);
        end
        end
        end
        plot(lims,lims,':');
        axis square;
        Axis.format('Observed','Predicted',[ vstr ' fit: ' fld]);
        Axis.xlim(lims);
        Axis.ylim(lims);
        if bLog
            ticks=[1/7 1/3 1 3 7];
            xticks(ticks);
            yticks(ticks);
        end
        if bLog
            set(gca,'YScale','log');
            set(gca,'XScale','log');
        end
    end
    function plotFitCov(obj)
        if nargin < 2 || isempty(fld)
            %fld='IO'; % fit to data
            fld='PO'; % fit to data
            fld='C0Y'; % fit to zero
        end
        if nargin < 3 || isempty(bLog)
            bLog=false;
        end

        Fig.new(['ExtP plotFitCov']);

        X_PO=getfun('PO');
        X_C0Y=getfun('C0Y');
        X=[X_PO; X_C0Y];

        edges=obj.hist_(X,[],'k');
        hold off;
        obj.hist_(X_PO,edges,'k');
        hold on;
        obj.hist_(X_C0Y,edges,'r');

        %plot(X_PO,X_C0Y,'ko')
        if bLog
            set(gca,'YScale','log');
            set(gca,'XScale','log');
        end
        function X=getfun(fld)
            x1=zeros(obj.sz0);
            x2=zeros(obj.sz0);
            n=numel(obj.Fit.C0Y);
            subjInds=1:3;


            for i = 1:obj.N
                [~,c,s]=ind2sub(obj.sz0,i);
                M=obj.Fit.(fld)(c,s);
                obji=Ext.LmCR(obj.VarL,M,obj.covLB,obj.O.Rho22(i));
                [~,x1(i)]=obji.rho33_rho23_error(obj.O.Rho33(i),0,true);
                [~,x2(i)]=obji.rho23_rho33_error(obj.O.Rho23(i),0,true);

            end
            x1=x1(:,subjInds);
            x2=x2(:,subjInds);
            X=[x1(:); x2(:)];
        end
    end
    function plotFit(obj,bLog)
        if nargin < 2 || isempty(bLog)
            bLog=false;
        end

        X=obj.Fit.C0Y(:,1:3);
        Y=obj.Fit.PO(:,1:3);
        mx=max([X(:); Y(:)]);
        if bLog
            mn=0.25;
        else
            mn=0;
        end
        lim=[mn mx];

        Fig.new('ExtP plotFit');
        hold off;
        plotfun(1,'k');
        plotfun(2,'r');
        plot(lim,lim,':r');
        axis square;
        ylim(lim);
        xlim(lim);
        xlabel('M | E = (Rho0 - Rho predicted)^2'); % C0Y
        ylabel('M | E = (Rho predicted - Rho observed)^2'); % PO
        %ylabel('M | E = (Rho0 - Rho observed)^2'); % IO
        Axis.format();
        title('Fitted M');
        if bLog
            set(gca,'YScale','log');
            set(gca,'XScale','log');
        end
        function plotfun(bin,color)
            X=obj.Fit.C0Y(bin,1:3);
            Y=obj.Fit.PO(bin,1:3);

            plot(X(:),Y(:),['o' color]);
            hold on;

        end
    end
end
end
