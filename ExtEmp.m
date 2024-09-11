classdef ExtEmp < handle & ExtFitData
properties
    fig
    SP
    h
    ES

    moude
    fInds
    nFld
    fld
    flds
    f

    bLog
    marker
    markerSize
    bcolor
    bcolors={'b','r','k'};
    color
    colors
    lineWidth
    lineWidths
    lineStyle
    lineStyles
    lims

    oMargin
    iMargin
    FaceAlpha

    ext
    rho23fit
    bCorr
    bCov
    bSubjNums=false
end
properties(Hidden)
    CH
    bLog_
    ylims
    xlims
    % tmp

end
methods(Static)
    function P=getP()
        P={...
           'bZer', false;
           'bBino',false;
           'fitType', 0;

           'bAvg',true;
           'bins',[];
           'oMargin',[0 0 0 0];
           'iMargin',[0 0 0 0];
           'FaceAlpha',0.4;
           'flds',[];
           'colors',[];
           'lineWidths',[];
           'lineStyles',[];
           'K',1;
        };
    end
    function fromExtModel(EM)
        obj=ExtEmp();
        obj.subjs=EM.Corr.Data.unqiue('subjs');
        obj.bin=EM.Corr.Data.unqiue('bin');
        obj.stds=EM.Corr.Data.unqiue('stdX');
    end
end
methods
    % T
    function plot_T2(obj)
        bLog=true;
        moude='T2';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        yl=[.5 4];
        obj.plot(moude,bLog,bZer,true,yl);
    end
    function plot_T3(obj)
        bLog=true;
        moude='T3';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        cl4=obj.tmpSet('lineStyle','none');
        %lbl='A';
        lbl=[];
        yl=[.5 4];
        obj.plot(moude,bLog,bZer,true,yl);
    end
    % T DIFF
    function plot_diff(obj)
        obj.plot_Abs2Diff;
        axis square;
        ylim([0 10]);

        obj.plot_Abs3Diff;
        axis square;
        ylim([0 10]);

        obj.plot_T2Diff;
        axis square;
        obj.plot_T3Diff;
        axis square;
        ylim([0 10]);
    end
    function plot_T3Diff(obj)
    %% figure tC-1
        moude='T3';
        bAll=true;
        cl1=obj.tmpSet('bSqrt',true);
        edges=linspace(-0.4,1.6,12);
        edges
		xlims=[-0.5 1.3];
        obj.histDiff(moude,bAll,edges,xlims);
    end
    function plot_Abs3Diff(obj)
        %% fig 5C2
        moude='IE3';
        bAll=true;
        cl1=obj.tmpSet('bSqrt',true);
        %edges=linspace(-0.3,1.3,11);
        edges=linspace(-0.4,1.6,12);
		xlims=[-0.5 1.3];
        obj.histDiff(moude,bAll,edges,xlims);
    end
    function plot_T2Diff(obj)
    %% figure 6C-1
        moude='T2';
        bAll=true;
        cl1=obj.tmpSet('bSqrt',true);
        edges=linspace(-0.4,1.6,12);
		%xlims=[-0.2 0.4];
		xlims=[-0.5 1.3];
        obj.histDiff(moude,bAll,edges,xlims);
    end
    function plot_Abs2Diff(obj)
        %% fig 6C2
        moude='IE2';
        bAll=true;
        cl1=obj.tmpSet('bSqrt',true);
        edges=linspace(-0.4,1.6,12);
        %xlims=[-0.5 0.6];
		xlims=[-0.5 1.3];
        obj.histDiff(moude,bAll,edges,xlims);
    end
    % Rho
    function plot_Rho2(obj)
        bLog=false;
        moude='rho22';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',false);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        %yl=[.1 .85];
        yl=[0 1];
        obj.plot(moude,bLog,bZer,true,yl);
    end
    function plot_Rho3(obj)
        bLog=false;
        moude='rho33';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',false);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        yl=[0 1];
        obj.plot(moude,bLog,bZer,true,yl);
    end
    function plot_Rho23(obj)
        bLog=true;
        moude='rho23';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',false);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        obj.plot(moude,bLog,bZer,true);
    end
    % ratio
    function plot_Rat2(obj)
        bLog=true;
        moude='rat22';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',false);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        yl=[.3 3];
        obj.plot(moude,bLog,bZer,true,yl);
    end
    function plot_Rat3(obj)
        bLog=true;
        moude='rat33';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',false);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        yl=[.3 3];
        obj.plot(moude,bLog,bZer,true,yl);
    end
    % Rho Diff
    function plot_Rho2Diff(obj)
        moude='rho22';
        bAll=true;
        cl1=obj.tmpSet('bSqrt',false);

        %edges=linspace(-0.4,1.6,12);
		%xlims=[-1.3 1.3];
        %
        edges=linspace(-0.6,0.6,12);
		xlims=[-0.65 0.65];
        yl=[0 16];

        obj.histDiff(moude,bAll,edges,xlims,yl);
    end
    function plot_Rho3Diff(obj)
        moude='rho33';
        bAll=true;
        %edges=linspace(-0.1,1.5,7);
        bSqrt=false;
        cl1=obj.tmpSet('bSqrt',false);

        edges=linspace(-0.6,0.6,12);
		xlims=[-0.65 0.65];
        yl=[0 16];

        %edges=[];
        obj.histDiff(moude,bAll,edges,xlims,yl);
    end
    % Abs
    function plot_Abs3(obj)
        %% fig 5AB
        bLog=true;
        moude='A3';
        bZer=false;
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        lTxt={'Total','External','Internal'};;
        yl=[.18 4.5];
        obj.plot(moude,bLog,bZer,true,yl,[],lTxt);
    end
    function plot_Abs2(obj)
        %% fig 6AB
        bLog=true;
        moude='A2';
        bZer=false;
        lTxt={'Total','External','Internal'};;
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        yl=[.18 4.5];
        obj.plot(moude,bLog,bZer,true,yl,[],lTxt);
    end
    function histDiffVar_Corru(obj,varargin)
        %% fig 5D
        if nargin < 2
            moude=[];
        end
        if nargin < 3
            bAll=[];
        end
        cl1=obj.tmpSet('bSqrt',true);
        edges=linspace(-1.5,1.5,10);
        xlims=[-1.5 1.5];
        %obj.hist_fun_(3,'IE3',bAll,edges,xlims);
        %obj.hist_fun_(4,'IE3',bAll,edges,xlims);
        obj.hist_fun_(5,'IE3',bAll,edges,xlims);

        % 3 = bin 1
        % 4 = bin 2
        % 5 = bin 1:2
    end
    function histDiffVar_Flat(obj,varargin)
        %% fig 5D
        if nargin < 2
            moude=[];
        end
        if nargin < 3
            bAll=[];
        end
        cl1=obj.tmpSet('bSqrt',true);
        edges=linspace(-1.5,1.5,10);
        xlims=[-1.5 1.5];
        %obj.hist_fun_(3,'IE2',bAll,edges,xlims);
        %obj.hist_fun_(4,'IE2',bAll,edges,xlims);
        obj.hist_fun_(5,'IE2',bAll,edges,xlims);

        % 3 = bin 1
        % 4 = bin 2
        % 5 = bin 1:2
    end
    function plot_vardiff(obj)
        obj.plot_Abs3DiffVar();
        axis square;
        %Figs.setFPos(gca);

        obj.plot_Abs2DiffVar();
        axis square;
        %Figs.setFPos(gca);
    end
    function plot_Abs3DiffVar(obj)
        %% fig 5D
        moude='IE3';
        bAll=true;
        edges=linspace(-2,2,12);
        cl1=obj.tmpSet('bSqrt',true);
        xlims=[-1.6 1.6];

        % 3 = bin 1
        % 4 = bin 2
        % 5 = bin 1:2
        %obj.hist_fun_(3,moude,bAll,edges,xlims);
        %obj.hist_fun_(4,moude,bAll,edges,xlims);
        obj.hist_fun_(5,moude,bAll,edges,xlims);
    end
    function plot_Abs2DiffVar(obj)
        %% fig 6D
        moude='IE2';
        bAll=true;
        edges=linspace(-2,2,12);
        cl1=obj.tmpSet('bSqrt',true);
        xlims=[-1.6 1.6];

        % 3 = bin 1
        % 4 = bin 2
        % 5 = bin 1:2
        %obj.hist_fun_(3,moude,bAll,edges,xlims);
        %obj.hist_fun_(4,moude,bAll,edges,xlims);
        obj.hist_fun_(5,moude,bAll,edges,xlims);
    end
    % EXT
    function plot_Ext(obj)
        %% figure 7A
        bLog=false;
        moude='0';
        bZer=false;
        ylm=[-1.25 2.25];
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        lTxt={'Low','High'};
        lTxtB='Disparity Contrast';
        obj.plot(moude,bLog,bZer,false,ylm,[],lTxt,lTxtB);
    end
    function plot_Ext0(obj)
        %% figure 7C
        bLog=false;
        moude='0';
        bZer=true;
        ylm=[-1.25 2.25];
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',0);
        lTxt={'Low','High'};
        lTxtB='Disparity Contrast';
        obj.plot(moude,bLog,bZer,false,ylm,[],lTxt,lTxtB);
    end
    function plot_Ext2(obj)
        %% figure 7E
        bLog=false;
        moude='0';
        bZer=true;
        ylm=[-1.25 2.5];
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',2);
        lTxt={'Low','High'};
        lTxtB='Disparity Contrast';
        obj.plot(moude,bLog,bZer,false,ylm,[],lTxt,lTxtB);
    end
    function plot_Ext1(obj)
        bLog=false;
        moude='0';
        bZer=true;
        ylm=[-1.25 2.52];
        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',true);
        cl3=obj.tmpSet('fitType',1);
        lTxt={'Low','High'};
        lTxtB='Disparity Contrast';
        obj.plot(moude,bLog,bZer,false,ylmlTxt,[],lTxt,lTxtB);
    end
    function plot_varBDiff(obj)
        % how does varB change with bZer
        % XXX
        edges=linspace(-0.5,3,11);
        moude='varB';
        bAll=true;


        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',false);
        cl3=obj.tmpSet('fitType',0);
        %diffmoude=1;
        %obj.histDiffZero(moude,bAll,edges,xlims);
        %obj.histDiffVar(moude,bAll,edges,xlims);

        xlims=[-2.5 1.5];
        edges=linspace(-2.5,1.5,8);
        obj.hist_fun_(2.1,moude,bAll,edges,xlims);

        xlims=[-.9 .5];
        edges=linspace(-.9,.5,8);
        obj.hist_fun_(2.2,moude,bAll,edges,xlims);

    end
    function plot_ExtDiff(obj)
        % how does varB change with bZer
        % XXX
        moude='varB';
        bAll=true;


        cl1=obj.tmpSet('bSqrt',true);
        cl2=obj.tmpSet('bAvg',false);
        cl3=obj.tmpSet('fitType',0);
        %diffmoude=1;
        %obj.histDiffZero(moude,bAll,edges,xlims);
        %obj.histDiffVar(moude,bAll,edges,xlims);

        %xlims=[-2.5 1.5];
        %edges=linspace(-2.5,1.5,8);
        xlims=[-1 2.5];
        edges=linspace(-1,2,9);
        obj.hist_fun_(1,moude,bAll,edges,xlims);

        moude='varL';
        obj.hist_fun_(1,moude,bAll,edges,xlims);


    end
    function plot_Cov(obj)
        %obj.histCovAll(true,'median');
        %obj.histCovCmp(true,'median');
        bCorr=true;
        bAll=true;
        obj.histCovAll(bAll,'median',bCorr);
    end
    function plot_diag(obj)
        obj.diag(1);
    end

%- Con
    function obj=ExtEmp(dims,SUBJS,passes,model,bNew,varargin)
        if nargin < 1
            return
        end
        if nargin ==1
            if nargin >0 && numel(dims)==1
                bNew=dims;
            else
                bNew=false;
            end
            dims=[2.1 3];
            SUBJS={'JDB','DNW','AJL'};
            passes=[1 2];
            model='fix';
        else
            EPFigs.parser(nargin,4);
        end
        if isa(dims,'ExtFitData')
            obj.Data=dims;
        else
            args={bNew,'dims',dims,'subjs_',SUBJS,'passes',passes,'model',model};
            obj.parse(args{:});
        end
        obj.empInit(varargin{:});
    end

    function empInit(obj,varargin)
        Args.simple(obj,obj.getP(),varargin{:});

        obj.ES=ESubjs();
        if isempty(obj.colors)
            obj.get_def_colors;
        end
        if isempty(obj.lineWidths)
            obj.get_def_lineWidths;
        end
        if isempty(obj.lineStyles)
            obj.get_def_lineStyles;
        end
    end
    function getExt(obj)
        T2=obj.E.T2;
        T3=obj.E.T3;
        cov=0;
        rho22=obj.E.rho22;
        rho33=obj.E.rho33;
        obj.ext=Ext.TTCRR(T2,T3,cov,rho22,rho33);
    end
%- MAIN
    function diag(obj,bLog)
        if nargin < 2 || isempty(bLog)
            bLog=obj.bLog;
        end

        obj.newFig('diag');
        cl1=obj.tmpSet('bZer',false);
        cl2=obj.tmpSet('bAvg',false);
        cl3=obj.tmpSet('bKAll',false);

        for s =1:obj.nSubj
        for b =1:obj.nBin
            obj.s=s;
            obj.b=b;

            r22=obj.get_val(false,true,'rho22');
            r23=obj.get_val(false,true,'rho23');
            T2=obj.get_val(false,true,'T2');
            T3=obj.get_val(false,true,'T3');

            X=r22 .* T2;
            Y=r23 .* T3;

            plot(X(:),Y(:),...
                'Marker',obj.marker,...
                'MarkerFaceColor',obj.bcolor,...
                'MarkerEdgeColor','none',...
                'LineStyle','none'...
            );
            hold on;
        end
        end

        xl=[.1 2];
        yl=xl;
        ylim(xl);
        xlim(yl);
        plot(xl,yl,'k:');
        Axis.format();
        ylabel('\rho^{*+} T^{*}');
        xlabel('\rho^{++} T^{+}');
        if bLog
            set(gca,'YScale','log');
            set(gca,'XScale','log');
            Axis.xticksLog();
            Axis.yticksLog();
        end
        axis square;
    end
    function plotRho23(obj,bLogFit,bRmNeg)
        if nargin < 2
            bLogFit=[];
        end
        if nargin < 3
            bRmNeg=[];
        end
        if isempty(obj.ext)
            obj.getExt();
        end
        obj.fitRho(bLogFit,bRmNeg);

        obj.newFig('plotRho23');
        for b = obj.bInds
        for s = obj.sInds
            obj.b=b;
            obj.s=s;

            val=obj.get_rho23();
            valP=obj.get_rho23Ext;

            subj=obj.subjs{s};
            if strcmp(subj,'AVG')
                subj='ALL';
            end

            plot(val,valP,'Marker',obj.marker,'Color',obj.color,'MarkerFaceColor',obj.color,'LineStyle','none');
            hold on;
        end
        end

        % line
        lim=obj.get_lim();
        plot(lim,lim,'--','Color',colors{end});

        % fit
        n=length(obj.rho23fit);
        rang=1:n;
        %rang=n;
        for i = rang;
            color=colors{i};
            fit=obj.rho23fit(i);
            x=Range.log(lim(1),lim(2),100);
            y=fit.m*x + fit.b;
            plot(x,y,'Color',color);
        end


        % format
        axis square;
        xlim(lim);
        ylim(lim);

        if obj.bLog_
            set(gca,'YScale','log');
            set(gca,'XScale','log');
        end
        Axis.format('Observed','Predicted 0','\rho^{*+}');

    end
    function out=get.bCorr(obj)
        out=any(strcmp('corrLB',obj.flds));
    end
    function out=get.bCov(obj)
        out=any(ismember(obj.flds,{'covLB','corrLB'}));
    end
    function newFig(obj,caller)
        % call after setting flds
        if obj.bZer
            zerStr=' Cov0';
        else
            zerStr='';
        end
        switch obj.fitType
        case 'A'
            fitStr='';
        case 'AFitT'
            fitStr=' T';
        case 'AFitRT'
            fitStr=' RT';
        end
        titl=['ExtEmp ' caller fitStr zerStr ' ' obj.moude ];
        obj.fig=Fig.new(titl);
        obj.h=matlab.graphics.chart.primitive.Line.empty();

        obj.f=obj.fInds;
        obj.b=obj.bInds;
        obj.s=obj.sInds;
        obj.k=obj.kInds;
        obj.get_ylims();
    end
    function getSP(obj,bSwap,suplbl,llTxt,llTxtB)
        if nargin < 2 || isempty(bSwap)
            bSwap=false;
        end
        if nargin < 3
            suplbl=[];
        end
        if bSwap
            RC=[obj.nFld obj.nSubj];
        else
            RC=[obj.nBin obj.nSubj];
        end
        obj.iMargin=[0 1 0 1];
        oMargin=obj.oMargin;
        iMargin=obj.iMargin;
        ytMargin=.5;

        % XXX
        ylMargin=1;
        switch obj.moude
        case {'rat22','rat33'}
            ylTxt='\sigma_E/\sigma_I';
            yticks=[.5 1 2];
        case {'A2','A3'}
            ylTxt='Threshold Contribution';
            if obj.bLog_
                yticks=[.25 .5 1 2 4];
                ytMargin=1.5;
            end
        case {'rho22','rho33'}
            %ylTxt=['Decision Variable' newline 'Correlation'];
            ylTxt=['Correlation'];
            yticks=[0.0 .25 .5 .75 1];
            ylMargin=2;
        otherwise
            if obj.bLog_
                yticks=[];
            else
                if strcmp(obj.moude,'0')
                    yticks=[-1 0 1 2];
                end
            end
        end
        if obj.bLog_
            yscale='log';
        else
            yscale='linear';
        end



        rlTxtB='';
        rlTxt='';
        bRl=false;
        bLl=false;
        llMargin=1;
        if RC(1) > 1
            if ~isempty(llTxt) || bSwap
                bLl=true;
                if ~isempty(llTxtB)
                    llMargin=2;
                end
            else
                bRl=true;

                iMargin(3)=1;
                rlTxt=Num.toCell(obj.bins);
                if isempty(suplbl)
                    oMargin(3)=1;
                end
                rlTxtB='bin';
            end
        else
            bins=[];
            rlTxtB='';
        end

		%if ~bSwap
		%	obj.SP.rlabel();
		%end
        %
        %if obj.bLog_
        %    Axis.yticksLog();
        %end
        %if nargin >= 7 && ~isempty(ylim)
        %    ylim(ylm);
        %end
        %
        X=obj.M.X;
        if isempty(obj.xlims)
            obj.xlims=Axis.xlim(X,0.1);
        end
        xfrmt='%.2f';

        xlinds=1:2:length(X);
        ylinds=[1 3 5];
        bCov=obj.bCov;
        if bCov && ~obj.bCorr && obj.bLog_
            yyScale='log';
            yylims=[];
            yyColor='m';
        else
            yyScale='';
            yylims=[];
            yyColor='';
        end

        if obj.bSubjNums
            titls=obj.subjs;
            strrep('ALL','Average');
        else
            titls=obj.subjs;
            repinds=~ismember(obj.subjs,{'ALL','AVG'});
            titls=strrep(obj.subjs,'ALL','Average')
            titls(repinds)=cellfun(@(x) ['Observer ' num2str(x)],num2cell(1:sum(repinds)),'UniformOutput',false);
        end

        args={...
            'slOn',~isempty(suplbl),...
            'slTxt',suplbl,...
            'slFontSize',40,...
            ...
            'rlOn',bRl,...
            'rlTxtB',rlTxtB,...
            'rlCtr', false,...
            'rlTxt',rlTxt,...
            ...
            'llOn',bLl,...
            'llCtr', false,...
            'llTxt',llTxt,...
            'llTxtB',llTxtB,...
            'llMargin',llMargin,...
            ...
            'ylCtr',false,...
            'ylBy',true,...
            'ylCtr',true,...
            'ylTxt',ylTxt,...
            ...
            'xlBy',true, ...
            'xlCtr',true,...
            'xlTxt','Disparity (arcmin)',...
            ...
            'xtCtr',false,...
            ...
            'xxlOn',true,...
            'xxlTxt',titls,...
            'xxlCtr',false,...
            ...
            'yyOn',bCov,...
            'yylOn',bCov,...
            'yylBy',bCov,...
            'yylCtr',bCov,...
            'yytCtr',false,...
            'yylTxt','Corr'...
            'yyColor',yyColor,...
            'yylTxtColor',yyColor,...
            'yylims',yylims,...
            'yyScale',yyScale,...
            ...
            'xlims',obj.xlims,...
            'xt',X,...
            'xtInds',xlinds,...
            'xtFormat',xfrmt,...
            ...
            'ylims',obj.ylims,...
            'yt',yticks,...
            'yScale',yscale,...
            'ytCtr',false,...
            ...
            'ytMargin',ytMargin,...
            'ylMargin',ylMargin,...
            'xlMargin',1.5,...
            'iMargin',iMargin, ...
            'oMargin',oMargin ...
        };
        obj.SP=SubPlot(RC,args{:});
    end
    function sel(obj,b,s)
        obj.b=b;
        obj.s=s;
        obj.SP.sel([b,s]);
        hold off;
    end
    function get_ylims(obj,bAvgF)
        if nargin < 2 || isempty(bAvgF)
            bAvgF=false;
        end
        V=[];
        for f = 1:obj.nFld
            obj.f=f;
            [valf,stdf]=obj.get_val(bAvgF);
            L=stdf(1,:);
            U=stdf(2,:);
            V=[V; L(:); U(:)];
        end

        p=0.05;
        if obj.bLog_
            V(V<=0)=[];
            yl=Num.minMax(log(V));
            r=abs(diff(yl));
            m=[-1 1]*r*p;
            obj.ylims=exp(yl+m);
        else
            yl=Num.minMax(V);
            r=abs(diff(yl));
            m=[-1 1]*r*p;
            obj.ylims=yl+m;
        end
    end
    function plot(obj,moude,bLog,bZer,bSwap,ylm,suplbl,lTxt,lTxtB)
        if nargin >= 2 && ~isempty(moude)
            obj.moude=moude;
        end
        if nargin < 3 || isempty(bLog)
            bLog=obj.bLog;
        end
        if nargin >=4 && ~isempty(bZer)
            obj.bZer=bZer;
        end
        if nargin < 5 || isempty(bSwap)
            bSwap=false;
        end
        if nargin < 7 || isempty(suplbl)
            suplbl=[];
        end
        if nargin < 8 || isempty(lTxt)
            lTxt=[];
        end
        if nargin < 9 || isempty(lTxt)
            lTxtB=[];
        end
        obj.xlims=[];
        obj.bLog_=bLog;

        obj.newFig('plot');
        if nargin < 6 || isempty(ylm)
            obj.get_ylims(true);
        else
            obj.ylims=ylm;
        end
        obj.getSP(bSwap,suplbl,lTxt,lTxtB);

        Valf=cell(obj.nBin,obj.nStd,obj.nFld);
        Stdf=cell(obj.nBin,obj.nStd,obj.nFld);

        fldsI=1:obj.nFld;
        binsI=obj.bInds;
        if bSwap
            fldsI=fliplr(fldsI);
        end
        binsI=fliplr(binsI);

        for b = binsI
        for s = obj.sInds
        for f= fldsI
            obj.s=s;
            obj.b=b;
            obj.f=f;
            [Valf{b,s,f}, Stdf{b,s,f}]=obj.get_val(true);
            L=Stdf{b,s,f}(1,:,:,:,:);
            if obj.bLog_ && any(L(:)<=0)
                Stdf{b,s,f}(L<=0)=Valf{b,s,f}(L<=0)-0.0001;
            end
        end
        end
        end

        for c = 1:2
        for b = binsI
        for s = obj.sInds
            for f= fldsI
                if bSwap
                    obj.sel(f,s);
                else
                    obj.sel(b,s);
                end
                obj.f=f;


                if c==2
                    hold on;
                end
                alpha=obj.FaceAlpha+(b-1)*0.3;

                st=Stdf{b,s,f};
                v=Valf{b,s,f};
                h=obj.plot_fun(v,st,obj.color,obj.lineStyle,obj.lineWidth,alpha,c);
                hold on;
                %if c==2
                %    if bSwap && b==1 && s==1
                %        obj.h(obj.f)=h;
                %    elseif bSwap && s==1
                %        obj.h(obj.f)={obj.h(obj.f); h};
                %    elseif ~bSwap
                %        obj.h(obj.f)=h;
                %    end
                %end

                if startsWith(obj.moude,'rat') && c==1 && f==1 && ((bSwap && b==1) || (~bSwap))
                    plot(obj.xlims,[1 1],'k--','LineWidth',2);
                end
            end
            if c==2 && b==binsI(end) && s== 1
                obj.legend(bSwap);
            end
        end
        end
        end
        %obj.legend();
        Figs.setFPos(obj.SP);
        %EPFigs.setWH21([],obj.nSubj);
    end
        function h=plot_fun(obj,val,std,colr,LineStyle,LineWidth,alpha,c)
            h=[];
            if nargin < 5 || isempty(LineStyle)
                LinStyle='-';
            end
            if nargin < 6 || isempty(LineWidth)
                LineWidth=1;
            end
            if ~iscell(val)
                val={val};
                std={std};
            end
            n=numel(val);

            if ~iscell(LineStyle)
                LineStyle={LineStyle};
            end

            if numel(LineWidth) > n
                LineWidth=repmat(LineWidth,n,1);
            end
            if numel(LineStyle) > n
                LineStyle=repmat(LineStyle,n,1);
            end

            iO={'EdgeColor','none','FaceAlpha',alpha};

            X=obj.M.X;
            bCorr=any(strcmp(obj.fld,{'covLB','corrLB'}));
            for i = 1:n
                if bCorr
                    yyaxis right;
                else
                    %yyaxis left;
                end

                if c==1
                    L=std{i}(1,:);
                    U=std{i}(2,:);
                    h=Plot.interv(X,L,U,colr,iO{:}); hold on
                else

                    Y=val{i};
                    h=plot(X, Y,'Color',colr,'MarkerSize',obj.markerSize,'Marker',obj.marker,'MarkerFaceColor','white','MarkerEdgeColor','black','LineWidth',2,'LineStyle',LineStyle{i});
                    if bCorr
                        g=gca;
                        g.YAxis(2).Color=colr;
                    end
                end
                hold on;
                %yyaxis left;
            end
        end
    function hist(obj,varargin)
        obj.hist_fun_(0,varargin{:});
    end
    function histDiff(obj,varargin)
        obj.hist_fun_(1,varargin{:});
    end
    function histDiffZero(obj,varargin)
        obj.hist_fun_(2,varargin{:});
    end
    function histDiffVar(obj,varargin)
        if nargin < 2
            moude=[];
        end
        if nargin < 3
            bAll=[];
        end
        if nargin < 4
            edges=[];
        end
        obj.hist_fun_(3,varargin{:});
    end
        function hist_fun_(obj,diffmode,moude,bAll,edges,xlims,ylims)
            if nargin < 2 || isempty(diffmode)
                diffmode=1;
            end
            if nargin >= 3 && ~isempty(moude)
                obj.moude=moude;
            end
            if nargin < 4 || isempty(bAll)
                bAll=true;
            end
            if nargin < 5
                edges=[];
            end
            if nargin < 6
                xlims=[];
            end
            if nargin < 7
                ylims=[];
            end

            bEdges=~isempty(edges);

            if bAll
                N=1;
            else
                N=obj.nStd;
                X=obj.M.X;
            end
            obj.bLog_=false;
            colors=obj.colors;

            obj.newFig(['histDiff ' num2str(diffmode)]);
            if ~obj.bKAll
                clk=obj.tmpSet('k',1);
            end

            obj.b=obj.bInds;
            obj.s=obj.sInds;
            clz=obj.tmpSet('bZer',obj.bZer);
            if strcmp(obj.moude,'1') && diffmode==2
                F=2:obj.nFld;
            else
                F=1:obj.nFld;
            end

            if diffmode==1
                cl2=obj.tmpSet('bZer',obj.bZer);
            elseif diffmode==3
                cl1=obj.tmpSet('bAvg',false);
            end

            D=[];
            bNew=true;
            Ymx=~inf;
            for c= (1+bEdges):2
            for j = 1:N
            for f = 1:obj.nFld
                obj.f=f;
                if bAll
                    inds=obj.nInds;
                else
                    inds=j;
                end


                switch diffmode
                case 0
                    % just hist
                    [valf,stdf]=obj.get_val(false);
                    d=valf(inds,:,:,:);
                case 1
                    % bins diff
                    [valf,stdf]=obj.get_val(false);
                    d=valf(inds,2,:,:)-valf(inds,1,:,:);
                case {2,2.1,2.2}
                    % zer diff
                    obj.bZer=false;
                    [valf1,stdf]=obj.get_val(false);
                    obj.bZer=true;
                    [valf2,stdf]=obj.get_val(false);
                    switch diffmode
                    case 2
                        b=1:2;
                    case 2.1
                        b=1;
                    case 2.2
                        b=2;
                    end

                    d=valf2(inds,b,:,:)-valf1(inds,b,:,:);
                case {3,4,5}
                    % subtract off mean
                    s=obj.s;
                    obj.bAvg=true;
                    obj.s=obj.nSubj;
                    [vM]=obj.get_val(false);

                    obj.bAvg=false;
                    obj.s=s;
                    [valf,stdf]=obj.get_val(false);
                    if diffmode==3
                        % bin 1
                        b=1;
                    elseif diffmode==4
                        % bin 2
                        b=2;
                    elseif diffmode==5
                        % bins 3
                        b=1:2;
                    end
                    d=(valf(inds,b,:,:)-vM(inds,b,:,:));
                end


                if c==1
                    D=[D; d(:)];
                    continue
                elseif bNew
                    if ~bEdges
                        [counts,edges]=histcounts(D(:));
                    end
                    %counts=[0 counts 0];
                    %edges=Hist.padEdges(edges);
                    ctrs=Hist.edges2ctrs(edges);
                    X=linspace(min(edges),max(edges),1000);

                    bNew=false;
                end
                [counts,edges]=histcounts(d(:),edges);
                D=dstNorm(d(:));
                D.get_fit();
                D.get_pdf(X);
                Y=D.pdf;
                Y=Y/sum(Y)*numel(Y) *sum(counts)./numel(counts);

                m=D.hat.mu;
                y=max(Y);
                Ymx=max([Ymx y]);

                %m=median(d(:));
                %y=interp1(ctrs(:),counts(:),m);
                if bAll
                    [~,~,hh]=Plot.bar([],edges,counts,'Color',obj.color,'LineWidth',2); hold on
                    plot(X,Y,'k','LineWidth',2);
                    plot([m m],[0 y],'--','Color',obj.color);
                else
                    x=repmat(X(j),numel(ctrs),1);
                    hh=plot3(x,ctrs,counts,'Color',obj.color,'LineWidth',2); hold on
                    plot3([X(j) X(j)],[m m],[0 y],'--','Color',obj.color);
                end
                if j==1
                    obj.h(f)=hh(1);
                end

            end
            end
            end
            if ~isempty(xlims)
                xlim(xlims);
            end
            if ~isempty(ylims)
                ylim(ylims);
            else
                ylim([0 Ymx*1.05]);
            end
            ylabel('Count');
            if startsWith(obj.flds{1},'var')
                fld='\sigma^2';
            elseif startsWith(obj.flds{1},'std')
                if any(contains(obj.flds{1},{'T2','T3'}))
                    fld='Threshold';
                else
                    fld='Threshold Contribution';
                end

            elseif startsWith(obj.flds{1},'rho')
                fld='\rho';
            end

            xlabel(['\Delta ' fld ]);
            %obj.legend();
            obj.legendH();
            Axis.format();
            axis square;
        end
    function histCovAll(obj,bAll,loc,bCorr)
        if nargin < 2
            bAll=[];
        end
        if nargin < 3
            loc=[];
        end
        if nargin < 4
            bCorr=[];
        end

        obj.hist_cov('all',bAll,loc,bCorr);
    end
    function histCovComb(obj,bAll,loc)
        if nargin < 2
            bAll=[];
        end
        if nargin < 3
            loc=[];
        end

        obj.hist_cov('comb',bAll,loc);
    end
    function histCovCmp(obj,bAll,loc)
        if nargin < 2
            bAll=[];
        end
        if nargin < 3
            loc=[];
        end

        obj.hist_cov('cmp',bAll,loc);
    end
    function histCovCmpDiff(obj,bAll,loc,bCorr)
        if nargin < 2
            bAll=[];
        end
        if nargin < 3
            loc=[];
        end
        if nargin < 4
            bCorr=[];
        end

        obj.hist_cov('cmpDiff',bAll,loc,bCorr);
    end
        function hist_cov(obj,moude,bAll,loc,bCorr,bBar)
            if nargin < 3
                bAll=[];
            end
            if nargin < 4 || isempty(loc)
                loc=[];
            end
            if nargin < 5
                bCorr=[];
            end
            if nargin < 6 || isempty(bBar)
                bBar=[];
            end
            obj.CH=ExtEmpCov(obj,'moude',moude,'bAll',bAll,'loc',loc,'bCorr',bCorr,'bBar',bBar);
            obj.CH.plot();
        end
%- Get/set
    %% flds
    function out=get.fInds(obj);
        out=1:obj.nFld;
    end
    function out=get.nFld(obj)
        out=numel(obj.flds);
    end
    function out=get.fld(obj)
        out=obj.flds{obj.f};
    end
    function out=get.flds(obj)
        out=obj.moude2flds(obj.moude);
        if obj.bSqrt
            out=regexprep(out,'^var','std');
        end
    end
    %% width
    function out=get.lineWidth(obj)
        out=obj.lineWidths{obj.f};
    end
    function out=get.lineWidths(obj)
        flds=obj.flds;
        n=length(flds);
        widths=cell(n,1);
        for i = 1:n
            fld=flds{i};
            switch fld
            case {'varT','varE','varI'}
                widths{i}=2;
            otherwise
                widths{i}=1;
            end
        end
        out=widths;
    end
    %% style
    function out=get.lineStyle(obj)
        out=obj.lineStyles{obj.f};
    end
    function out=get.lineStyles(obj)
        flds=obj.flds;
        n=length(flds);
        styles=cell(n,1);
        for i = 1:n
            fld=flds{i};
            switch fld
            case 'varE'
                styles{i}='--';
            case 'varI'
                styles{i}=':';
            otherwise
                styles{i}='none';
            end
        end
        out=styles;
    end
    %% colors
    function out=get.bcolor(obj)
        out=obj.bcolors{obj.b};
    end
    function out=get.color(obj)
        out=obj.colors{obj.f};
    end
    function out=get.colors(obj)
        flds=obj.flds;
        n=obj.nFld;
        colors=cell(n,1);
        orange=[1 0.651 0];
        yellow=[0.5 0.5 0];
        [~,colors]=obj.moude2flds(obj.moude);
        if ~isempty(colors)
            out=colors;
            return
        end
        for i = 1:n
            fld=flds{i};
            switch fld
            case {'varL','stdL'}
                if obj.bBino
                    colors{i}='k';
                else
                    colors{i}='b';
                end
            case {'varD','stdD'}
                colors{i}=yellow;
            case {'varB','stdB'}
                colors{i}='r';
            case {'covLB','corLB'}
                colors{i}='m';

            case {'varE2','stdE2'}
                colors{i}=[.2 .2 .7];
            case {'varE3','stdE3'}
                colors{i}=[.7 .2 .2];

            case {'varI2','stdI2'}
                colors{i}=[0 0 .3];
            case {'varI3','stdI3'}
                colors{i}=[.3 0 0];

            case {'varT2','stdT2'}
                colors{i}='b';
            case {'varT3','stdT3'}
                colors{i}='r';
            case {'rho22','rat22'}
                colors{i}='b';
            case {'rho33','rat33'}
                colors{i}='r';
            case {'rho23','rat23'}
                colors{i}='m';
            otherwise
                colors{i}='k';
            end
        end
        out=colors;
    end
    %% marker
    function out=get.marker(obj)
        out=obj.ES.(obj.subj).Marker;
    end
    function out=get.markerSize(obj)
        out=obj.ES.(obj.subj).MarkerSize;
    end
    %% -
    function lims=get.lims(obj,flds)
        if nargin < 2 || isempty(flds)
            flds=obj.moude2flds(obj.moude);
        end
        if any(contains(flds,'var'))
            lims=[-1.5 4];
        elseif any(contains(flds,'rho'))
            lims=[-0.5 1];
        end
        if obj.bLog_
            lims(1)=max([lims(1) .01]);
        end
    end
%- Ext
    function fitRho(obj,bLog,bRmNeg,bFix0)
        if nargin < 2 || isempty(bLog)
            bLog=false;
        end
        if nargin < 3 || isempty(bRmNeg)
            bRmNeg=false;
        end
        if nargin < 4 || isempty(bFix0)
            bFix0=false;
        end
        if isempty(obj.ext)
            obj.getExt();
        end
        obj.rho23fit=struct();
        %sInds=1:(obj.nSubj-obj.bAvg);
        sInds=obj.sInds;
        obj.k=1;
        cl1=obj.tmpSet('k',1);
        cl2=obj.tmpSet('bKAll',false);
        for ib = 1:(obj.nBin+1)
            if ib==obj.nBin+1;
                obj.b=1:obj.nBin;
            else
                obj.b=ib;
            end
            rho=obj.get_rho23();
            rhop=obj.get_rho23Ext();

            X=rho(:);
            Y=rhop(:);
            nanind=X<=0 | Y<=0 | isnan(X) | isnan(Y);
            if bRmNeg
                X(nanind)=[];
                Y(nanind)=[];
            end
            if bLog
                X=log(X);
                Y=log(Y);
            end
            mFix=0;
            bFix0=1;
            [m,B,sigma,x,y,stats]=mydeming(X,Y,bFix0,mFix,1,0.3173);
            obj.rho23fit(ib).m=B(1);
            obj.rho23fit(ib).sigma=sigma;
            if bLog
                obj.rho23fit(ib).m=m;
                obj.rho23fit(ib).b=B;
                obj.rho23fit(ib).x=exp(x);
                obj.rho23fit(ib).y=exp(y);
            else
                obj.rho23fit(ib).m=m;
                obj.rho23fit(ib).b=B;
                obj.rho23fit(ib).x=x;
                obj.rho23fit(ib).y=y;
            end
            obj.rho23fit(ib).stats=stats;
            %SSR=sum(orthProjPoint2Line([X Y]).^2,m);
            %SST=sum( [X Y]-mean([X Y],1) ).^2;
            %R2=1-SSR./
            %obj.rho23fit(ib).SSR=SSR;
            %obj.rho23fit(ib).SST=SST;
            %obj.rho23fit(ib).R2=R2;

        end
    end
    function [valL,valp]=get_rho23Ext(obj)
        nI=obj.ext.nanInds;
        valL=nan(size(nI));
        valp=zeros(size(nI));
        valL(~nI)=obj.ext.rho23();

        b=obj.b;
        s=obj.s;
        valL=valL(:,b,s);
        valp=valp(:,b,s);
        %valp=obj.ext.rho23_s(:,b,s);
    end
    function [flds,colors]=moude2flds(obj,moude)
        if nargin < 2 || isempty(moude)
            moude=obj.moude;
        end
        if isempty(moude)
            moude=0;
        end
        if isnumeric(moude)
            moude=num2str(moude);
        end

        pink=[237/255  119/255 119/255];
        green=[43/255,111/255,70/255];
        colors=[];
        switch moude
        case '0'
            % ext
            flds={'varL','varB','covLB'};
            if obj.bBino
                flds=['varD' flds];
            end
        case '1'
            flds={'varL','varB','corLB'};
            if obj.bBino
                flds=['varD' flds];
            end
        case 'T3'
            flds={'varT3'};
            colors={[0 0 0]};
        case 'T2'
            flds={'varT2'};
            colors={[0 0 0]};
        case 'rat22'
            flds={'rat22'};
            colors={[0 0 0]};
        case 'rat33'
            flds={'rat33'};
            colors={[0 0 0]};
        case 'rho22'
            flds={'rho22'};
            colors={[0 0 0]};
        case 'rho33'
            flds={'rho33'};
            colors={[0 0 0]};
        case 'rho23'
            flds={'rho23'};
            colors={[0 0 0]};
        case 'varB'
            flds={'varB'};
            colors={[1 0 0]};
        case 'varL'
            flds={'varL'};
            colors={[0 0 1]};
        case '2'
            flds={'varT2','varE2','varT2'};
        case '3'
            flds={'varT3','varE3','varI3'};
        case 'A2'
            flds={'varT2','varE2','varI2'};
            colors={[0 0 0], green,pink};
        case 'A3'
            flds={'varT3','varE3','varI3'};
            colors={[0 0 0], green,pink};
        case 'I'
            flds={'varI2','varI3'};
        case 'IE2'
            flds={'varE2','varI2'};
            colors={green,pink};
        case 'IE3'
            flds={'varE3','varI3'};
            colors={green,pink};
        case 'E'
            flds={'varE2','varE3'};
        case 'T'
            flds={'varT2','varT3'};
        case 'rho'
            flds={'rho22','rho33'};
        otherwise
            error('invalid moude')
        end
        if obj.bSqrt
            flds=regexprep(flds,'^var','std');
        end
    end
%- FORMAT
    function leg=flds2txt(obj,flds)
        if nargin < 2 || isempty(flds)
            flds=obj.flds;
        end
        n=length(flds);
        leg=cell(n,1);
        for i = 1:n
            fld=flds{i};
            switch fld
            case 'covLB'
                leg{i}='cov(L,B)';
            case 'corLB'
                leg{i}='corr(L,B)';
            otherwise
                if startsWith(fld,'rho')
                    fld=strrep(fld,'rho','\rho');
                    fld=strrep(fld,'22','^{++}');
                    fld=strrep(fld,'33','^{**}');
                    leg{i}=strrep(fld,'23','^{+*}');
                elseif startsWith(fld,'rat')
                    fld   =strrep(fld,'rat22','^{++}\sigma_E/\sigma_I');
                    leg{i}=strrep(fld,'rat33','^{**}\sigma_E/\sigma_I');
                elseif startsWith(fld,'var')
                    fld   =regexprep(fld,'var(.)2','^+\\sigma^2_{$1}');
                    fld   =regexprep(fld,'var(.)3','^*\\sigma^2_{$1}');
                    leg{i}   =regexprep(fld,'var(.)','\\sigma^2_{$1}');
                elseif startsWith(fld,'std')
                    fld   =regexprep(fld,'std(.)3','^*\\sigma_{$1}');
                    fld   =regexprep(fld,'std(.)2','^+\\sigma_{$1}');
                    leg{i}=regexprep(fld,'std(.)','\\sigma_{$1}');
                    %leg{i}=strrep(fld,'var','var_');
                end
            end
        end
    end
    function legendH(obj)
        l=obj.flds2txt(obj.flds);
        %l=repelem(obj.flds2txt(obj.flds),2,1);
        legend(obj.h,l,'location','northwest');
    end
    function legend(obj,bSwap)
        if numel(obj.bInds) > 1 && bSwap
            l={'bin 1','bin 2','bin 1','bin 2'};
        else
            l=obj.flds2txt(obj.flds),2,1;
            %l=repelem(obj.flds2txt(obj.flds),2,1);
        end
        if obj.bCov
            h1=flipud(obj.SP.ax{obj.SP.i}.Children);

            yyaxis left;
            h2=flipud(obj.SP.ax{obj.SP.i}.Children);

            h=[h2; h1];
            h=[h2(1:(length(h2)/2)); h1(1:(length(h1)/2))];

            hold off;
        else
            h=flipud(obj.SP.ax{obj.SP.i}.Children);
        end
        legend(h,l,'location','southeast');

    end
    function format(obj,bSwap)
        % LEGEND
    end
end
end
