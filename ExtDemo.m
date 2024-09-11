classdef ExtDemo < handle
properties
    model
    bAnimate

    covLB
    VarL
    M
    m
    T2
    T3
    Rho22
    Rho23
    Rho33

    bCorr
    FRCFlds

    SP
    figs

    iMargin
    oMargin
    oUnits
    iUnits
end
properties(Hidden)
    FRC
    frc
    sp
    fig

    xstr
    ystr
end
methods(Static)
    function P=getP_EP()
        P={...
           'fitType',0;
           'bAll',false;
        };
    end
    function P=getP()
        P={...
           'model','LmCR';
           'bAnimate',true;
            ...
           'VarL',1;
           'M',[0.8 1.5 2.6];
           'm',[0.8 1.5 2.6];
           'covLB',[];
           'Rho22',[];
           'Rho23',[];
           'Rho33',[];
           'T2',[];
           'T3',[];
           ...
           'bCorr',true;
           'FRCFlds',[];
           ...
           'oMargin',[2 .5 3 4];
           'iMargin',[1 0 1 1];
           'oUnits','char';
           'iUnits','char';
           ...
        };
    end
    function obj=fromEP(dims,SUBJS,passes,model,bNew,varargin)
        EPFigs.parser(nargin,4);
        [EPArgs,varargin]=Args.simpleLoose('args',ExtDemo.getP_EP(),varargin{:});
        fitType=EPArgs{2};
        bAll=EPArgs{4};

        if isempty(SUBJS)
            SUBJS=Figs.getSUBJS(dims(1));
        end
        nSubj=numel(SUBJS);

        if bAll
            sInd=nSubj;
        else
            sInd=1:nSubj;
        end

        S=struct(varargin{:});

        %% GET DATA
        E=EPs.ext_data(dims,SUBJS,passes,model,bNew,fitType);
        covLB=[];
        rho22=E.rho22(1,:,sInd);
        rho23=E.rho23(1,:,sInd);
        rho33=E.rho33(1,:,sInd);
        T2=E.T2;
        T3=E.T3;
        rho22=permute(rho22,[1 3 2]);
        rho23=permute(rho23,[1 3 2]);
        rho33=permute(rho33,[1 3 2]);
        T2=permute(T2,[1 3 2]);
        T3=permute(T3,[1 3 2]);

        %% rho22
        m=mean(E.rho22(:));
        st=std(E.rho22(:));
        Rho22=sort([m-st m m+st]);

        %% rho23
        Rho23=mean(rho23(:));

        %% rho33
        Rho33=mean(rho33(:));

        %% T2
        T2=mean(T2,[2 3]);

        %% T3
        T3=mean(T3,[2 3]);

        S=struct(varargin{:});
        flds={'Rho22','Rho23','Rho33','T2','T3'};

        for i = 1:length(flds);
            fld=flds{i};
            if ~isfield(S,fld)
                val=eval([fld ';']);
                varargin=[varargin fld val];
            end
        end
        obj=ExtDemo(varargin{:});
    end
end
methods
%- CON
    function obj=ExtDemo(varargin)
        Args.simple(obj,ExtDemo.getP(),varargin{:});
        if isempty(obj.FRCFlds)
            obj.FRCFlds=obj.get_sFRCFlds(true);
        end
        if isempty(obj.covLB)
            N=1000;
            c=4;
            obj.covLB=linspace(-c,c,N);
        end
    end
%- MAIN
    function d1(obj);
        if obj.bAnimate
            obj.d1_fun('rho');
        else
            obj.d1_fun('rho');
            obj.d1_fun('var');
        end
    end
    function d1_fun(obj,typ)
        %: LMCR & LmCR

        X=Num.minMax(obj.covLB);

        bRho23=~isempty(obj.Rho23) && numel(obj.Rho23) > 1;
        bRho33=~isempty(obj.Rho33) && numel(obj.Rho33) > 1;
        bRhoBoth=bRho23 & bRho33;

        obj.get_FRC();
        obj.init_SP(1,typ);

        if obj.bAnimate
            nA=6;
        else
            nA=1;
        end
        LineWidth=1.8;

        switch typ
        case 'var'
            if contains(obj.model,'m')
                m='M';
            else
                m='m';
            end
            flds={'varL','varB','varI2',m};
            Colors={'b','r','k','m'};
            LineStyle='-';
        case 'rho'
            flds={'rho22','rho23','rho33'};
            Colors={'b','m','r'};
            LineStyle='-';
        case 'm'
            flds={'M','m'};
            Colors={'b','r'};
            LineStyle='-';
        end

        plt=ExtPlotVar('bCorr',obj.bCorr,'LineWidth',LineWidth,'flds',flds,'Colors',Colors);
        for ia = 1:nA;
        for ir = 1:obj.FRC(1)
            for ii = 1:obj.FRC(2)
            for ij = 1:obj.FRC(3)
                obj.sel(ir,ii,ij);
                hold off;


                sflds=obj.get_sFRCFlds();
                args=cell(1,numel(sflds));
                [args{:}]=obj.unpackFRC();

                ext=Ext.(obj.model)(args{1:2},obj.covLB,args{3:end});
                ext

                % PLOT
                h=plt.plot(ext);

                % INTERP
                if bRho23
                    %try
                        h=obj.interpPlot(ext,'rho23','rho33',obj.Rho23,'m','r',h,ia);
                    %end
                end
                if bRho33
                    %try
                        h=obj.interpPlot(ext,'rho33','rho23',obj.Rho33,'r','m',h,ia);
                    %end
                end

                % FORMAT
                obj.get_titles(ext);
                obj.format();
                if ii ==1 && ij == 1
                    obj.legend(h,[],bRho23,bRho33,typ);
                end
            end
            end
            obj.fformat();
            if obj.bAnimate
                drawnow;
                waitforbuttonpress;
            end
        end
        end

    end
    function h=interpPlot(obj,ext,fld1,fld2,val,colo,colo2,h,ia)
        meth=[fld1 '_' fld2 '_error'];
        %meth=[fld1 '_' fld2 ];

        if obj.bCorr
            mm=[-1 1];
            X=ext.rhoLB;
        else
            mm=Num.minMax(obj.covLB);
            X=ext.covLB;
        end
        %[Xi,Ri]=ext.(meth)(val,obj.bCorr,true);
        [e,Xi,Ri,R0]=ext.(meth)(val,0,obj.bCorr);
        if (Xi < mm(1) || Xi > mm(2))
            if (ij ==1 && ii==1)
                x1=X;
                y1=[val val];
            else
                return
            end
        else
            x1=[nanmin(X) Xi];
            y1=[val val];
        end
        x2=[x1(2) Xi];
        y2=[y1(2) Ri];

        x0=[0 0];
        y0=[0 R0];


        if ~obj.bAnimate || ia==2 || ia == 3
            Plot.arrow(x1,y1,[.06 10],[ colo]);
        end
        if ~obj.bAnimate || ia==3
            Plot.arrow(x2,y2,[.03 35],[ colo]);
        end

         % 0 val
        if ~obj.bAnimate || ia == 4
             Plot.arrow(x0,y0,[.03 35],['--' colo2]);
        end

         % X error
         if ia == 5
             Plot.arrow([Xi 0],[Ri Ri],[.03 90],2,'k');
             plot([0 0],[R0 Ri],':k');
         end

         % Y error
         if  ia == 6
             Plot.arrow(x2,[Ri R0],[.06 90],2,'k');
             plot([0 Xi],[R0 R0],':k');
         end


         if  ia >= 3
             plot(Xi,Ri,'o','MarkerFaceColor',colo,'MarkerEdgeColor','none','MarkerSize',6);
         end
         if  ia >= 4
             plot(0,R0,'o','MarkerFaceColor','w','MarkerEdgeColor',colo,'MarkerSize',6);
         end

    end
%- FRC
    function sflds=get_sFRCFlds(obj,bDefault)
        switch obj.model
        case 'LMCR'
            sflds={'VarL','M','Rho22'};
        case 'LmCR'
            sflds={'VarL','m','Rho22'};
        case 'LmCT'
            sflds={'VarL','M','T2'};
        case 'TTCRR'
            sflds={'T2','T3','Rho22','Rho33'};
        end
        if nargin >= 2 && bDefault && ~strcmp(obj.model,'TTCRR')
            sflds=sflds([1 3 2]);
        end
    end
    function get_FRC(obj)
        obj.FRC=zeros(1,3);
        for i = 1:length(obj.FRCFlds);
            fld=obj.FRCFlds{i};
            obj.FRC(i)=numel(obj.(fld));
        end
        obj.frc=[0 0 0];
    end
    function varargout=unpackFRC(obj)
        bTest=false;
        varargout=cell(1,3);
        sflds=obj.get_sFRCFlds();
        for i = 1:numel(sflds)
            sfld=sflds{i};
            ind=find(strcmp(sfld,obj.FRCFlds));
            fld=obj.FRCFlds{ind};
            vals=obj.(fld);
            if ind > 3
                val=vals;
            else
                val=vals(obj.frc(ind));
            end

            %% TEST
            if bTest && all(obj.frc==[3 2 1]);
                disp([num2str(i) ' ' fld ' ' num2str(val)])
            end

            varargout{i}=val;
        end
    end
    function init_SP(obj,dnum,typ)
        if obj.bCorr
            corrstr='corr';
            obj.xstr='corr(L,B)';
        else
            corrstr='cov';
            obj.xstr='cov(L,B)';
        end
        if nargin < 3 || isempty(typ)
            typ='rho';
        end
        switch typ
        case 'rho'
            varstr='rho';
            obj.ystr='\rho';
        case 'var'
            varstr='var';
            obj.ystr='\sigma^2';
        case 'm'
            varstr='ratios';
            obj.ystr='ratios';
        end
        obj.SP=SubPlot.empty;
        for ir = 1:obj.FRC(1)
            titl=sprintf('ExtDemo %s %s %s %d-%d',obj.model,corrstr,varstr,dnum,ir);
            obj.figs(ir)=Fig.new(titl);
            obj.SP(ir)=SubPlot(obj.FRC(2:3),'iMargin',obj.iMargin,'oMargin',obj.oMargin,'iUnits',obj.iUnits,'oUnits',obj.oUnits);
            hold off;
            axis off;
        end
    end
    function sel(obj,ir,ii,ij)
        irLast=obj.frc(1);
        obj.frc=[ir,ii,ij];
        if ir ~= irLast
            figure( obj.figs(ir));
            obj.fig=obj.figs(ir);
            obj.sp =obj.SP(ir);
        end
        obj.sp.sel([ii,ij]);
    end
    function format(obj)
        obj.sp.xticks();
        obj.sp.yticks();
        obj.sp.xlabel(obj.xstr,true,[]);
        obj.sp.ylabel(obj.ystr,true,[],[]);
    end
    function fformat(obj,ir)
        obj.sp.rlabel([],[],[], 0);
        obj.sp.clabel([],[],[],-1);
        set(obj.fig,'Position',[-1312 107 1189 1061]);
    end
    function get_titles(obj,obji)
        [mStr,rStr,vStr,t2Str,r3Str,t3Str]=obji.get_titles();

        switch obj.FRCFlds{2}
        case 'VarL'
            rtitl=vStr;
        case {'M','m'}
            rtitl=mStr;
        case 'Rho22'
            rtitl=rStr;
        case 'T2'
            rtitl=t2Str;
        case 'T3'
            rtitl=t3Str;
        case 'Rho33'
            rtitl=r3Str;
        end

        switch obj.FRCFlds{3}
        case 'VarL'
            ctitl=vStr;
        case {'M','m'}
            ctitl=mStr;
        case 'Rho22'
            ctitl=rStr;
        case 'T2'
            ctitl=t2Str;
        case 'T3'
            ctitl=t3Str;
        case 'Rho33'
            ctitl=r3Str;
        end

        obj.sp.rlabels=rtitl;
        obj.sp.clabels=ctitl;
    end
    function legend(obj,h,loc,bRho23,bRho33,typ)
        if nargin < 2 || isempty(loc)
            loc='NW';
        end
        if nargin < 3 || isempty(bRho23)
            bRho23=false;
        end
        if nargin < 4 || isempty(bRho33)
            bRho33=false;
        end
        if nargin < 5 || isempty(typ)
            typ='rho';
        end


        switch typ
        case 'var'
            l={'$\sigma^2_L$','$\sigma^2_B$','$\sigma^2_I$'};
            if contains(obj.model,'m')
                m='M';
            else
                m='m';
            end
            l=[l m];
        case 'rho'
            l={'$\rho^+$','$\rho^{+*}$','$\rho^{**}$'};
            if bRho23
                l=[l '$\hat\rho^{+*}$' ];
            end
            if bRho33
                l=[l '$\hat\rho^{**}$'];
            end
        end

        warning('off','MATLAB:legend:IgnoringExtraEntries');
        l=legend(h,l,'location',loc,'Interpreter','latex');
    end
end
methods(Static)
%- SCRAP
    function obj=d2_1(covLB,T2,T3,rho22)
        % TTmCR
        N=1000;

        T2=1.25;
        T3=1.75;
        M=[1/3 1 3];
        covLB=linspace(-0.7,0.7,N);
        rho22=0.5;

        %-
        nR=numel(T2);
        nC=numel(M)+1;
        RC=[nR,nC];

        % suplot
        hold off;
        Fig.new('Ext Demo2 - Variable \sigma^2_I');

        for i = 1:nR
        for j = 1:nC
            subPlot(RC,i,j);
            if j == 1
                OBJ=Ext.TTmCR(T2(1),T3,M(1),covLB,rho22);
                OBJ.text_plotVarCov(i);
                continue
            end
            obj=Ext.TTmCR(T2(i),T3,M(j-1),covLB,rho22);
            obj.plotVarCov();
            obj.format_plotVarCov(i,j-1,RC);
        end
        end

    end
end
end
