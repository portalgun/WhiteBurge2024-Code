classdef ImgStats < handle
properties
    P
    Blk

    PStats

    RInds
    RTable

    RCorr
end
properties(Hidden)
    p
    img
    z
    winInds
    winSzs
    statnames
end
methods
%- INIT
    function getImgStats(obj,bForce)
        if isempty(obj.P)
            obj.load_patches();
        end
        if nargin < 2 || isempty(bForce)
            bForce=[];
        end

        if bForce || isempty(obj.PStats)
            obj.get_MN();
        end
        obj.concat_blk();
        obj.corr_mn();
    end
    function load_patches(obj)
        if isempty(obj.P)
            obj.P=ptchs.load('DSP2');
        end
        %obj.Bf=Blk.load('DSP2f2');
        %[obj.Bf.blk('P').ret obj.B.blk('P').ret]

        lvlInds=[1 4 6 9 11 14 16 19 21 24]';
        obj.P.exp_init('DSP2',1,lvlInds)
        obj.P.apply_display('jburge-wheatstone.psych.upenn.edu',[])
    end
%- GET Y
    function concat_blk(obj)
        obj.init();

        stats=obj.statnames;
        %stats={'mean','max','lrg','sml'};
        BOrder={'lvlInd','blk','trl','intrvl'};
        DOrder={'lvlInd','blk','trial','aliasInd','subj','pass'};

        Blk=copy(obj.Blk);

        for i = 1:length(stats)
            Blk.add_col(stats{i},obj.PStats.(stats{i}));
        end

        % SORT
        Blk.sort(BOrder{:});
        obj.data.sort(DOrder{:});
        obj.data.order_cols(DOrder{:},':');

        aliases=obj.data.unique('aliasInd');
        subjs=obj.data.unique('subj');
        passes=obj.data.unique('pass');
        N=numel(aliases)*numel(subjs)*numel(passes);

        % Rows to columns
        Blk.rm_col({'cndInd','cmpInd','mode','intrvl'});
        Blk2=Blk.copy();
        Blk2.rm_col('P');

        Blk.rm_col(stats);
        Blk.rows2col('cmpNum','P',{'P1','P2'});

        for i = 1:length(stats)
            blk=Blk2.copy();
            stat=stats{i};
            stat1=[stat '1'];
            stat2=[stat '2'];
            rmstats=stats(~ismember(stats,stat));
            blk.rm_col(rmstats);

            % note sorting by cmpNum here in order to subtract
            blk.rows2col('cmpNum',stat,{stat1,stat2});

            [s1,s2]=blk(stat1,stat2).ret();
            Blk.add_col(stat,s1-s2);
        end

        Blk.rename_col('trl','trial');
        Blk.repelem(N);

        % reorder like 'data'
        Blk.horzcat(obj.data('aliasInd','subj','pass'));
        Blk.order_cols(DOrder{:},':');

        % JOIN
        horzcat(obj.data,Blk('P1','P2',stats{:}));

        % XXX
        %d=obj.data('d');
        for i = 1:length(stats)
            obj.data('aliasInd',2,stats{i})=0;
        end

        obj.data=obj.data('lvlInd',[1 4 6 9 11 14 16 19 21 24])

        obj.table2table();

    end
%% GX
    function get_MN(obj)
    % MN(i,k,w)
        obj.get_windows();

        N=length(obj.P.Blk.blk);
        nW=size(obj.winInds,1);

        MN=zeros(N,nW,2);
        MX=zeros(N,nW,2);
        SD=zeros(N,nW,2);

        pr=Pr(N,10);
        for i = 1:N
            pr.u();
            try
                obj.get_patch(i);
            catch
                disp(i);
                continue
            end

            winD=obj.p.win.win.dist;
            zMap=obj.z;
            for j = 1:nW
                for k = 1:2
                    % subtract win distance
                    Z=(zMap{k}(obj.winInds{j})-winD);
                    MN(i,j,k)=nanmean(Z);
                    MX(i,j,k)=nanmax(Z);
                    SD(i,j,k)=std(Z);
                end
            end
        end
        pr.u();

        obj.PStats.MN=MN;
        obj.PStats.MX=MX;
        obj.PStats.SD=SD;

        obj.get_stats();
    end
    function get_stats(obj)
        rmflds=fieldnames(obj.PStats);
        rmflds(ismember(rmflds,{'MN','MX','SD'}))=[];
        if ~isempty(rmflds)
            obj.PStats=rmfield(obj.PStats,rmflds);
        end
        MN=obj.PStats.MN;
        MN=mean(MN,3); % over k

        % SZ
        for i = 1:size(MN,2)
            j=obj.winSzs(i);
            fld=['sz' num2str(j)];
            obj.PStats.(fld)=MN(:,i);
        end

        % MEAN
        obj.PStats.mean=mean(MN,2); % over windows

        % MAX
        MX=obj.PStats.MX;
        MX=Num.maxMag(MX,[],3);
        MX=Num.maxMag(MX,[],2);
        obj.PStats.max=MX;

        %SD=obj.PStats.SD;
        %obj.PStats.sd=SD;
        %obj.PStats.std=S;
    end
    function corr_mn(obj,typ)
        if nargin < 2 || isempty(typ)
            typ='dice';
        end
        %typ='PHI';
        %typ='RT';
        %typ='SM';

        vars={'lvlInd','aliasInd','subj'};
        lvlInds=[1 4 6 9 11 14 16 19 21 24];

        stats=obj.statnames;

        DD=obj.data('lvlInd',lvlInds);
        U=obj.data.unique_rows(vars{:});
        DD.add_col('corr',nan(length(obj.data),1))
        DD.add_col('VAL',nan(length(obj.data),1))
        obj.RCorr=struct();

        N=length(U);
        n=N*numel(stats);
        pr=Pr(n,5);
        for iStat=1:length(stats)
            stat=stats{iStat};

            D=DD.copy();

            for i = 1:N
                pr.u();
                [d,a,S]=U(i).ret();
                dat=D('lvlInd',d,'aliasInd',a,'subj',S);
                UC=dat.unique('cmpd');
                for j = 1:length(UC)
                    [R,val]=dat('cmpd',UC(j),'R',stat).ret();
                    M=val>0;
                    out=simbin(R,M,typ);
                    D('lvlInd',d,'aliasInd',a,'subj',S,'cmpd',UC(j),'corr')=out;
                    D('lvlInd',d,'aliasInd',a,'subj',S,'cmpd',UC(j),'VAL')=val;
                end
                %dat
            end
            obj.RCorr.(stat)=D('lvlInd','B','aliasInd','subj','stdd','cmpd','d','corr','VAL');
        end

    end
    function imgGlm(obj)
        mdl='new2';
        obj.Glm.fit(mdl,'aliasInd',1,'d',0,'B',13);
    end
    function plot_corr_mn(obj,stat,varargin)
        %subj='DNW';
        %B=3;
        %stat='mean';
        if nargin < 3
            varargin={'aliasInd',1};
        end

        [x,y]=obj.RCorr.(stat)(varargin{:},'VAL','corr').ret();
        Fig.new(['plot_corr_mn ' stat]);
        %Fig.new(['plot_corr_mn ' stat ' ' subj ' ' num2str(B)]);
        plot(x,y,'.')
        xlabel(stat)
        ylabel('corr')
        Axis.format();
        axis square
    end
    function plot_corr_sz_all(obj,b0)
        if nargin < 2
            b0=true;
        end

        dat=obj.data('aliasInd',1,'B',1);
        L=dat.unique('lvlInd');
        subjs=dat.unique('subj');


        nS=length(subjs);
        nL=length(L);
        %tlTxt=cellstr(num2str(L));
        tlTxt={'-112.5','-9.38','-7.5','-5.62','-3.75'};
        xt=4:4:34;

        Fig.new(['plot_corr_sz_all ' num2str(b0)] );
        SP=SubPlot([nS,nL],'ylims',[0 1],'xlims',[0 34],'rlTxt',subjs,'tlTxt',tlTxt,'ylTxt','correlation','xlTxt','winSz','xlBy',true,'ylBy',true,'ytBy',true,'xtBy',true,'xt',xt,'xtInds',1:2:length(xt));
        colors={'r','y','k','b','m'};
        for i=1:nS
        for j=1:nL
            SP.sel([i,j]);

            d=dat('lvlInd',L(j));
            C=d.unique('cmpd');
            C2=d.unique('cmpD')*60;

            if b0
                obj.plot_corr_sz('-','k','d',0,'aliasInd',1,'lvlInd',L(j),'subj',subjs{i})
            else
                for c = 1:length(C)
                    obj.plot_corr_sz('-',colors{c},'cmpd',C(c),'aliasInd',1,'lvlInd',L(j),'subj',subjs{i})
                    %obj.plot_corr_sz(colors{c},'cmpd',C(c),'aliasInd',1,'lvlInd',L(j),'subj',subjs{i})
                end
            end
        end
        end
    end
    function plot_corr_sz_one(obj,varargin)
        Fig.new('plot_corr_sz_one');
        Axis.format();
        dat=obj.data('aliasInd',1,'B',13);
        L=dat.unique('lvlInd');
        nL=length(L);

        tlTxt={'-11.25','-9.38','-7.5','-5.62','-3.75'};
        xt=4:4:34;

        SP=SubPlot([1,nL],'ylims',[0 1],'xlims',[0 34],'tlTxt',tlTxt,'ylTxt','correlation','xlTxt','winSz','xlBy',true,'ylBy',true,'ytBy',true,'xtBy',true,'xt',xt,'xtInds',1:2:length(xt));

        for i=1:nL
            SP.sel(i);
            obj.plot_corr_sz('.','k','aliasInd',1,'d',0,'B',13,'lvlInd',L(i),varargin{:})
        end
        xlabel('winSize');
        ylabel('correlation');
        %zlabel('\Delta');
        %xlim([0 max(mX(:))+2]);
        %ylim([0 1]);
        axis square

    end
    function plot_corr_sz(obj,marker,color,varargin)
        if nargin < 2
            varargin={'aliasInd',1};
        end
        stats=flipud(obj.statnames);
        stats(~startsWith(stats,'sz'))=[];
        Y=[];
        X=[];
        Z=[];
        mX=[];
        mY=[];
        mZ=[];

        for i = 1:length(stats)
            stat=stats{i};
            sz=str2double(stat(3:end));
            [y,z]=obj.RCorr.(stat)(varargin{:},'corr','d').ret();

            my=median(y,1);
            %my=repmat(my,size(y,1),1);
            mY=[mY; my];

            x=repmat(sz,size(y,1),1);
            X=[X;x];
            mX=[mX; sz];

            Y=[Y;y];
            Z=[Z;z];
        end


        plot(X,Y,marker,'MarkerEdgeColor',color,'Color',color); hold on
        %plot(mX,mY,'or','MarkerFaceColor','w','MarkerEdgeColor',color);

        %plot3(X,Y,Z,'k.'); hold on
        %plot3(X,mY,zeros(size(mY)),'or','MarkerFaceColor','w');

        %Axis.format();
        %xlabel('winSize')
        %ylabel('correlation')
        %zlabel('\Delta')
        %xticks(mX);
        %xlim([0 max(mX(:))+2])
        %ylim([0 1])
        %axis square
    end
%- PLOT
    function hist_corr_mn(obj,stat,varargin)
        if nargin < 3
            varargin={'aliasInd',1};
        end
        %MN=obj.PStats.MN;
        %MN=mean(MN,3); % over k
        Fig.new('hsit_corr')
        [x,y]=obj.RCorr.(stat)(varargin{:},'VAL','corr').ret();
        hist(y)
        Axis.format();
        xlabel('correlation')
        ylabel('count')
        axis square
    end
    function validate_stats(obj)
        dat=obj.data.copy();
        dat=dat('B',13);

        x=obj.data('d').ret();
        stats=flipud(obj.statnames);
        N=length(stats);

        Fig.new('ImgStats validate_stats')
        hold off;
        SP=SubPlot([1 N],'xlTxt','d','ylTxt','c','tlTxt',stats,'ylims',[-.1 .1])
        for i = 1:N
            SP.sel(i);

            stat=stats{i};
            y=obj.data(stat).ret();
            plot(x,y,'.');
        end
        %plot('d',{'lrg','sml'},[],[],'Marker','.','LineStyle','none')
        %xl=xlim();
        %dat=dat.sort('d');
        %d=dat('d').ret()
        %plot(d)
    end
%% Helpers
    function get_windows(obj)
        obj.get_patch(1);

        SZ=obj.p.PszRC;
        obj.winSzs=(SZ:-4:4)';
        winSzs=[winSzs winSzs];

        n=size(winSzs,1);
        inds=cell(n,1);
        for i = 1:n
            sz=winSzs(i,:);
            inds{i}=obj.get_window(SZ,sz);
        end
        obj.winInds=inds;
    end
    function inds=get_window(obj,sz,PszRC,im)
        %sz=size(im)
        %sz=[32 32];
        ctr=ceil(sz./2)+0.5.*(mod(sz,2)==0);

        [x0,y0]=Rec.rectPix([0,0],PszRC(1),PszRC(2));
        AitpRC0=Set.distribute(y0(3):y0(1), x0(2):x0(1));
        subs=AitpRC0+ctr;
        inds=sub2ind(sz,subs(:,1),subs(:,2));

        %figure(99)
        %im=obj.img{1}(inds);
        %im=reshape(im,PszRC)
        %imagesc(im)


        %im(ind)
    end
%- GET
    function get_patch(obj,num)
        obj.P.get_patch(num)
    end
    function out=get.p(obj)
        out=obj.P.ptch();
    end
    function out=get.img(obj)
        out=obj.P.ptch.im.img();
    end
    function out=get.z(obj)
        % HERE
        xyz=obj.P.ptch.maps.xyzD;
        %xyz=obj.P.ptch.maps.xyz;
        out{1}=xyz{1}(:,:,3);
        out{2}=xyz{2}(:,:,3);
    end
    function out=get.Blk(obj)
        if isempty(obj.P)
            out=[];
        else
            out=obj.P.Blk.blk;
        end
    end
    function stats=get.statnames(obj)
        stats=fieldnames(obj.PStats);
        stats=stats(~ismember(stats,{'MN','MX','SD'}));
    end
%- PLOT
    function plot_im(obj,num)
        if nargin >= 2 && ~isempty(num)
            obj.get_patch(num);
        end

        m=obj.p.maps;

        flds={'xyz','xyzD','pht'};

        N=length(flds);
        for i = 1:N
            fld=flds{i};

            subplot(N,1,i);
            if startsWith(fld,'xyz')
                imagesc([m.(fld){1}(:,:,3) m.(fld){2}(:,:,3)])
            else
                imagesc([m.(fld){1} m.(fld){2}])
            end
            Axis.formatIm();
            title(fld)
            colorbar;
        end
    end
end
end
