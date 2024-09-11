classdef ExtModel < handle 
properties
    aliases
    subjs
    moude
    status
    lvlInds
    blks
    passes

    raw
    data
    splData
    splData2
    mTable
    mTableSpl
    mTableSpl2

    typ
    mdl
    name

    Rsqr
    LL

    Glm=ExtGlm()
    Cust=ExtCust()
    RNova=ExtRNova()

    Thresh
    Corr
    CorrS
    CorrSF
    Ext

    Emp
end
properties(Hidden)
    dv
    cv
    mdls

    lastOut
    lastArgs
end
methods(Static)
    function tbl=lvlsTable()
        %    stdInd; disparityInd bbin; lvl     bin
        tbl=[...
             1       1             1   -11.25  1 ;
             4       1             13  -11.25  13;
             6       2             1   -9.375  1 ;
             9       2             13  -9.375  13;
             11      3             1   -7.5    1 ;
             14      3             13  -7.5    13;
             16      4             1   -5.625  1 ;
             19      4             13  -5.625  13;
             21      5             1   -3.75   1 ;
             24      5             13  -3.75   13;
        ];
    end
end
methods
    function obj=raw_params(obj)
        tbl=ExtModel.lvlsTable();  % XXX
        lvls=tbl(:,1);
        obj.aliases={'DSP2','DSP2f2'};
        %obj.aliases={'DSP2'};
        obj.subjs={'DNW','AJL','JDB'};
        obj.moude=1;
        obj.lvlInds=lvls;
        obj.passes=[1 2];
        obj.blks=1:5;
    end
%- CON
    function init_mdls(obj,bForce)
        if nargin < 2 || isempty(bForce)
            bForce=false;
        end
        obj.mdls={'Glm','Cust','RNova'};
        for i = 1:length(obj.mdls)
            fld=obj.mdls{i};
            if isempty(obj.(fld)) || bForce
                obj.(fld)=eval(['Ext' fld ';']);
            end
            obj.(fld).parent=obj;
        end
    end
    function init(obj,bRaw)

        if nargin < 2
            bRaw=isempty(obj.raw);
        end

        if bRaw
            obj.get_raw_data();
        end
        obj.get_data();
        obj.table2table();

        obj.init_mdls();

        obj.table2split2();
        obj.load();
    end
%- SET/GET
    function out=get.Rsqr(obj)
        out=obj.get_('Rsqr');
    end
    function out=get.LL(obj)
        out=obj.get_('LL');
    end
    function out=get_(obj,Fld)
        out=struct();
        %test
        for i = 1:length(obj.mdls) % test comment
            fld=obj.mdls{i};
            s=obj.(fld).(Fld);

            flds=fieldnames(s);
            for j =1:length(flds)
                fld2=flds{i};
                out.([fld '.' fld2])=s.(fld2);
            end
        end
    end
    %% MDL
    function set.typ(obj,typ)
        if ~isprop(obj,'typ')
            error(['Invalid typ:' typ])
        end
        obj.typ=typ;
        obj.mdl=obj.(obj.typ);
    end
    function set.name(obj,name)
        if isempty(obj.typ)
            error('must set ''typ'' first');
        end
        obj.mdl.name=name;
    end
    function out=get.name(obj)
        if isempty(obj.typ)
            out=[];
        end
        out=obj.(obj.typ).name;
    end
    function out=get.mdl(obj)
        if isempty(obj.typ)
            out=[];
        end
        out=obj.(obj.typ);
    end
    function sel(obj,typ,name)
        if nargin==2
            if isprop(obj,typ)
                obj.typ=typ;
            elseif ~isempty(obj.typ) && obj.mdl.hasname(typ)
                obj.name=typ;
            else
                error('Invalid typ')
            end
        else
            obj.typ=typ;
            obj.name=name;
        end

    end
    function fit(obj,varargin)
        if nargin >= 2 && ~isempty(varargin{1})
            obj.sel(varargin{:});
        end
        obj.mdl.fit();
    end
%- EP DATA
    function load_emp(obj)
        f=Figs;
        obj.Emp=f.ExtEmp;
    end
%- DATA
    function get_raw_data(obj)
        obj.raw_params();

        sels={...
              'subj', obj.subjs, ...
              'mode','%%',unique(obj.moude), ...
              'status','>' 0, ...
              'lvlInd', obj.lvlInds, ...
              'blk', obj.blks, ...
              'pass' obj.passes, ...
        };

        Data=[];
        obj.aliases
        for i = 1:length(obj.aliases)
            E=ETable.get(obj.aliases{i});
            data=E.loadRawData(sels{:});
            data.add_col('aliasInd',i,1);
            if i==1
                Data=data;
            else
                Data=[Data; data];
            end
        end
        obj.raw=Data;
    end
    function get_data(obj)
        data=obj.raw.copy();

        rmflds={'date','iter','flags','RTime','mode','cmpInt','R'};
        data=data('iter',0); % XXX check
        data.rm_col(rmflds);

        data.rename_col('bRCmp','R');
        data.rename_col('bRCorrect','C');

        %% bin and disp from cmpX and stdX
        data.split_col();
        %data.rename_split({'stdX','cmpX'},{'Bin';'bin'});
        %data.rename_split({'stdX','cmpX'},{'Bin','bin'});
        %data.rename_split({'stdX'},{'Bin'});
        data.rename_split({'stdX','cmpX'},{'stdD','stdB','cmpD','cmpB'});
        %data
        %data.rm_col('cmpB');

        %% cmpX and stdX to deptH
        %D=data.minus('cmpD','stdD');
        %data.add_col('diff_cmpD_stdD',D,3);


        %% trlInfo to id
        [~,~,id]=data.unique_rows('lvlInd','blk','trial');
        data.add_col('id',id,1);
        [~,~,id]=data.unique_rows('subj','lvlInd','blk','trial');
        data.add_col('ids',id,1);


        %% lvlInd and aliasInd to B
        bCor=data.bfind('aliasInd',1); %NOTE MAKE SURE TO CHECK THIS
        B=data('stdB').ret();
        B(~bCor)=0;
        data.add_col('B',B);

        D=data('cmpD').ret()-data('stdD').ret();
        data.add_col('D',D);

        %[~,~,id]=data.unique_rows('aliasInd','lvlInd','trial');
        % {'aliasInd','lvlInd'}

        %% lvlInd to L
        %data.rename_col('stdB','L');

        %data.rm_col({'lvlInd','blk','trial'});
        %data.order_cols('ids','id','subj','cmpD','stdD','stdB','diff_cmpD_stdD','B','R','C','pass','aliasInd');
        data.order_cols('ids','id','aliasInd','pass','subj','lvlInd','blk','trial','cmpD','stdD','D','stdB','B','R','C');

        data.unique('B');

        obj.data=data;
        obj.gen_data_mod_();
    end
    function gen_data_mod_(obj)
        data=obj.data;

        %% TO DEPTH
        xyzF=[0 0 1];
        IPDm=0.0625; %XXX
        vrgF=XYZ.to_vrg_naive(IPDm,xyzF(3));

        flds={'cmpD','stdD'};
        locs=[5 6];
        for i = length(flds):-1:1
            fld=flds{i};
            dsp=data(fld).ret();
            vrgT=dsp+vrgF;
            dT=XYZ.to_z_naive(IPDm,vrgT);
            nfld=strrep(fld,'D','d');
            data.add_col(nfld,dT,locs);

            ldsp=log(dsp);
            nfld=['l' fld];
            data.add_col(nfld,ldsp);
        end
        d=data('cmpd').ret()-data('stdd').ret();
        data
        data.add_col('d',d);
        obj.data=data;
    end
    function get_RInds(obj)
        Raw=obj.EM.raw('aliasInd',1);
        subjs=Raw.unique('subj');
        nSubj=numel(subjs);
        nPass=max(Raw.unique('pass'));

        n=length(obj.Blk);
        inds=zeros(n*nPass*nSubj,1);
        c=0;

        [lvl,blk,trl]=obj.Blk('lvlInd','blk','trl').ret();
        for i = 1:n
            %raw=Raw('trial',trl,'lvlInd',lvl,'blk',blk);
            for iP = 1:nPass
            for iS = 1:nSubj
                c=c+1;
                %[bRCmp(i,iP,iS),bRCorr(i,iP,iS)]=raw('subj',subjs{iS},'pass',iP,'bRCmp','bRCorrect').ret();
                inds(c)=Raw.find('trial',trl(i),'lvlInd',lvl(i),'blk',blk(i),'subj',subjs{iS},'pass',iP);
            end
            end
        end
        obj.RInds=inds;
        %obj.bRCmp=bRCmp;
        %obj.bRCorrect=bRCorr;
    end
    function table2split(obj,bSubj)
        if nargin < 2
            bSubj=false;
        end
        T=obj.table2split_fun(false,bSubj);
        obj.splData=T;
    end
    function table2split2(obj,bSubj)
        if nargin < 2
            bSubj=false;
        end
        T=obj.table2split_fun(true,false);
        obj.splData2=T;
    end
    function T=table2split_fun(obj,b2,bSubj)
        % R - response
        % C - correct
        if nargin < 3
            bSubj=true;
        end
        subjs=obj.data{'subj'};
        subjs=unique(subjs);

        if b2
            A=1:2;
        else
            A=1;
        end

        if bSubj
            S=1:numel(subjs);
            idfld='id';
        else
            S=1;
            idfld='ids';
        end

        astr='';
        sstr='';
        for s = S
        for a = A
        for p = 1:2
            % strs
            pstr=num2str(p);
            if b2
                astr=num2str(a);
            end
            if bSubj
                sstr=num2str(s);
            end
            rstr=['R' astr pstr sstr];
            cstr=['C' astr pstr sstr];

            % args
            args={'pass' p};
            if b2
                args=[args 'aliasInd' a];
            end
            if bSubj
                args=[args 'subj' subjs{s}];
            end

            % get data
            t=obj.data(args{:});

            % index
            id=t(idfld).ret();
            [ids,idx]=sort(id);
            t.index(idx);

            if p==1 & a==1 && s==1;
                t.rename_col('R',rstr);
                t.rename_col('C',cstr);
                T=t;
                ids0=ids;
            else
                r=t('R').ret();
                c=t('C').ret();
                bGd=isequal(ids,ids0);
                if ~bGd
                    error('this')
                end

                T.add_col(rstr,r,lastRstr);
                T.add_col(cstr,c,lastCstr);
            end

            lastRstr=rstr;
            lastCstr=cstr;
        end
        end
        end
        T.rm_col({'stdB','pass'});
        if b2
            T.rm_col('aliasInd');
        end
        if bSubj
            T.rm_col({'subj','ids'});
        end
    end
    function subjs=get_within(obj,name)

        T=obj.splData;
        K=T.KEY;
        mtchs=regexp(K,'^R[0-9]+$','match','once');
        mtchs(cellfun(@isempty,mtchs))=[];
        subjs=zeros(numel(mtchs),1);
        for i = 1:length(mtchs)
            subjs(i)=str2double(mtchs{i}(end));
        end
    end
    function table2table(obj)
        obj.mTable=obj.data.toMTable();
        if ~isempty(obj.splData)
            obj.mTableSpl=obj.splData.toMTable();
        end
        if ~isempty(obj.splData2)
            obj.mTableSpl2=obj.splData2.toMTable();
        end
    end
    function data=toCSV(obj)
        if isempty(obj.mTable)
            obj.mTable;
        end
        name='WhiteBurge2024';

        data=obj.raw.copy();

        data.rm_col({'lvlInd','flags','iter','mode','RTime','date','bRCmp','bCorrect'})
        data.rename_col('aliasInd','exp');

        % subj data
        S=data('subj').ret();
        SNew=zeros(size(S,1),1);
        subjs={'JDB','DNW','AJL'};
        for i = 1:length(subjs)
            inds=strcmp(S,subjs{i});
            SNew(inds)=i;
        end
        data.rm_col('subj');
        data.add_col('subj',SNew,1);

        % std and cmps
        data.split_col('stdX');
        data.split_col('cmpX');
        data.rm_col('stdX(2)');
        data.rename_col('stdX(1)','stdX');
        data.rename_col('cmpX(1)','cmpX');

        % disparity contrast
        data.rename_col('cmpX(2)','CX');
        data('CX',13)=2;
        data.order_cols('subj',...
                        'exp','blk','trial','pass',...
                        'CX','stdX','cmpX',...
                        'cmpInt','R');
        dire=Env.var('Data');
        fname=[dire name];
        data.toCSV(fname);
    end
    function D=table2Input(obj,inds,gInds,varargin)
        obj.table2split2();
        Dat=obj.splData2;
        if nargin >= 4
            Dat=Dat(varargin{:});
        end

        D=struct();
        D.RCmpChs=Dat.regexp('R[0-9]*').ret();
        if nargin > 2 && ~isempty(inds)
            D.RCmpChs=D.RCmpChs(:,inds);
        end
        [subj,D.stdX,D.cmpX]=Dat.ret('subj','stdD','cmpD');
        [D.subjUnq,~,D.subjId]=unique(subj);

        if nargin > 3 && ~isempty(gInds)
            R=D.RCmpChs;
            n=max(hist(gInds,1:max(gInds)));
            N=cell(n,max(gInds));
            counts=zeros(1,max(gInds));
            for i= 1:size(R,2)
                ind=gInds(i);
                counts(ind)=counts(ind)+1;
                N{counts(ind),ind}=R(:,i);
            end
            bE=~cellfun(@isempty,N);
            nRep=sum(bE,1);
            R=cell(1,size(N,2));
            S=cell(1,size(N,2));
            C=cell(1,size(N,2));
            for i=1:size(N,2)
                R{:,i}=vertcat(N{:,i});
                S{:,i}=repmat(D.stdX,nRep(i),1);
                C{:,i}=repmat(D.cmpX,nRep(i),1);
            end
            D.RCmpChs=R;
            D.stdX=S;
            D.cmpX=C;
        end

    end
    function D=table2DVInput(obj,inds,varargin)
        D=obj.table2Input(inds,[],varargin{:});
    end
    function D=table2CurveInput(obj,inds,gInds,varargin)
        D=obj.table2Input(inds,gInds,varargin{:});
    end
%- PLOT
    function plotProb(obj)
        obj.mdl.mdl.plotResiduals('probability');
    end
    function plotFitted(obj)
        obj.mdl.mdl.plotResiduals('fitted');
    end
    function plotSym(obj)
        obj.mdl.mdl.plotResiduals('symmetry');
    end
    function plotLagged(obj)
        obj.mdl.mdl.plotResiduals('lagged');
    end
    function plotCase(obj)
        obj.mdl.mdl.plotResiduals('caseorder');
    end
    function plotHist(obj)
        obj.mdl.mdl.plotResiduals('hist');
    end
%- FIT MODELS
    function fitThresh(obj)
        nBest=5;
        nBoot=0;
        if nargin < 2
            bReset=false;
        end
        obj.thresh_fun(nBest,nBoot,bReset);
    end
    function fitCorrRR(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='RR';
        nBest=1;
        nBoot=1;
        obj.fit_corr_fun(model,nBest,nBoot,flag);
    end
    function fitCorrE(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='E';
        nBest=1;
        nBoot=1000;
        obj.fit_corr_fun(model,nBest,nBoot,flag);
    end
    function fitCorrE0(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='E0';
        nBest=50;
        nBoot=1000;
        obj.fit_corr_fun(model,nBest,nBoot,flag);
    end
    function fitCorrS(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='RR';
        nBest=1;
        nBoot=2;
        obj.fit_corr_fun(model,nBest,nBoot,flag,1);
    end
    function fitCorrSF(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='RR';
        nBest=1;
        nBoot=1000;
        obj.fit_corr_fun(model,nBest,nBoot,flag,2);
    end
    function fitCorrF(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='F';
        nBest=1;
        nBoot=1;
        obj.fit_corr_fun(model,nBest,nBoot,flag);
    end
    function fitCorrF0(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='F0';
        nBest=1;
        nBoot=1;
        obj.fit_corr_fun(model,nBest,nBoot,flag);
    end
    function fitCorrSpilt(obj,flag)
        if nargin < 2
            flag=[];
        end
        model='FB';
        nBest=1;
        nBoot=1;
        obj.fit_corr_fun(model,nBest,nBoot,flag,true);
    end
%- Thresh
    function init_thresh(obj)
        Thresh=struct();

        Key={'subj','aliasInd','bin','stdX','thresh'};
        Typs={'char','double','double','double','double'};
        Thresh.Data=Table([],Key,Typs);

        Opts=PsyCurveFit.getDefaults();
        Opts.bLinear=true;
        Opts.bLogLinear=true;
        Opts.bMuStdFix=true;
        Opts.DPCrt=1;
        Thresh.Opts=Opts;

        obj.Thresh=Thresh;
    end
    function thresh_fun(obj,nBest,nBoot,bReset)
        if isempty(obj.Thresh) || bReset
            obj.init_thresh();
        end

        T=obj.Thresh.Data;
        Opts=obj.Thresh.Opts;
        Opts.nBest=nBest;
        Opts.nBoot=nBoot;

        %% --
        inds=[1:4];
        gInds=[1 1 2 2];
        bins=obj.splData2.unique('B')';
        subjs=obj.splData2.unique('subj');
        nS=numel(subjs);

        for b = bins
        for s = 1:nS
            D=obj.table2CurveInput(inds,gInds,'subj',subjs{s},'B',b);
            nT=numel(D.RCmpChs);
            for i = 1:nT
                dat=PsyCurveData(D.stdX{i},D.cmpX{i},D.RCmpChs{i});

                obj.cv=PsyCurveFit(dat,Opts);
                obj.cv.run();
                obj.cv.plot();
                obj.cv.plotT();
                drawnow;

                T.add_row_sorted(subjs{s},i,b,dat.stdXUnq,obj.cv.tFit);
            end
        end
        end
    end
%- Corr
    function [Opts,Key,Typs]=get_corr_opts(obj,bSubj)

        Opts=DVFitter.getDefaults;
        Opts.bCtrMu=true;
        Opts.bCombineCmp=true;
        Opts.bootSeed=666;
        Opts.bootCIPrcnt=68.27;
        Opts.bObs=bSubj>0;

        if bSubj==1
            subj={'subj','subj2','E'};
            subjt={'char','char','double'};
            Opts.rhoFitInd=[1 2 2 2 2 3];
        elseif bSubj==2
            subj={'subj','subj2','E'};
            subjt={'char','char','double'};
            %Opts.rhoFitInd=[1 1 1 2 2 3 2 2 3 3];
            Opts.rhoFitInd=[0 1 1 1 1 0];
        else
            subj={'subj'};
            subjt={'char'};
            Opts.rhoFitInd=[1 2 2 2 2 3];
        end

        Key=['model','fitType', subj ,'bin',   'stdX',  'n',     'negLL'  'R1',    'R12',   'R2',    'R1B',   'R12B',   'R2B','A1','A12','A2'];
        Typs=['char','char'   , subjt,'double','double','double','double','double','double','double','double','double','double','double','double','double'];

    end
    function init_corr(obj,bSubj)
        [Opts,Key,Typs]=obj.get_corr_opts(bSubj);
        C=struct();
        C.Opts=Opts;
        C.Data=Table([],Key,Typs);


        if bSubj==1
            fld='CorrS';
        elseif bSubj==2
            fld='CorrSF';
        else
            fld='Corr';
        end
        obj.(fld)=C;
    end
    function fit_corr_fun(obj,model,nBest,nBoot,flag,bSubj)
        if isempty(flag)
            flag=0;
        end
        %
        % bSubj
        if nargin < 6 || isempty(bSubj)
            bSubj=0;
        end
        switch bSubj
        case 0
            fld='Corr';
        case 1
            fld='CorrS';
        case 2
            fld='CorrSF';
        end
        bDisp=true;

        bSplit=contains(model,'B');
        bReset=flag==-1;
        bNew=isempty(obj.(fld));

        % best and boot
        bBest=nBest > 1;
        bBoot=nBoot > 1;
        if ~bBoot && ~bBest
            typs={'test'};
        else
            typs={};
        end
        if bBest
            typs=[typs 'best'];
        end
        if bBoot
            typs=[typs 'boot'];
        end
        nT=numel(typs);

        % new or reset
        if bNew
            obj.init_corr(bSubj);
        elseif bReset
            obj.(fld).Data=obj.(fld).Data('model','~=',model,'fitType','~=',typs);
            bNew=true;
        end
        bImprove=flag==1 && ~bNew;

        % Data to save
        T=obj.(fld).Data;
        Opts=obj.(fld).Opts;

        Opts.nBest=nBest;
        Opts.nBoot=nBoot;
        Opts.modelType=model;
        keys=T.KEY;
        dKeys=keys(1:find(strcmp('stdX',keys),1,'first'));

        % inds
        inds=[1:4];

        % bins
        if bSplit
            bins=0;
            Opts.nSplit=2;
        else
            bins=obj.splData2.unique('B')';
        end
        nB=numel(bins);

        [~,rInds]=unique(Opts.rhoFitInd);
        RFlds=regexp(keys,'^R[0-9]*$','match','once');
        RFlds(cellfun(@isempty,RFlds))=[];

        % subjs
        subjs=obj.splData2.unique('subj');
        nS=numel(subjs);
        rInds=[1,2,6];
        if bSubj>0
            nS2=nS;
            nE=2;
        else
            nS2=1;
            nE=1;
        end

        %p=[];
        count=0;

        for b = bins
        for s = 1:nS
        for s2= 1:nS2
        for e = 1:nE
            subj=subjs{s};
            subj2=subjs{s2};
            if bSubj && s==s2
                continue;
            elseif ~bSubj
                subj2=subj;
            end

            % experiment number
            if bSubj>0
                aliasInds=e;
            else
                aliasInds=[1 2];
            end

            % get data
            if bSplit
                D=obj.table2DVInput(inds,'subj',subj,'B',1);
                D2=obj.table2DVInput(inds,'subj',subj,'B',13);
                D.RCmpChs=[D.RCmpChs D2.RCmpChs];
            elseif bSubj>0
                eInds=(1:2)+2.*(e-1);
                D=obj.table2DVInput(inds,'subj',subj,'B',b);
                D2=obj.table2DVInput(inds,'subj',subj2,'B',b);
                D.RCmpChs=[D.RCmpChs(:,eInds) D2.RCmpChs(:,eInds)];
            %elseif bSubj==2
                %D1=obj.table2DVInput(inds,'subj',subj,'B',b);
                %D2=obj.table2DVInput(inds,'subj',subj2,'B',b);
                %%D.RCmpChs=cat(3,D.RCmpChs,D2.RCmpChs);
                %D=struct();
                %fld='RCmpChs';
                %D.(fld)=[cat(1,D1.(fld)(:,[1,3]),D1.(fld)(:,[2,4])),...
                %         cat(1,D2.(fld)(:,[1,3]),D2.(fld)(:,[2,4]))];
                %D.stdX=[D1.stdX; D1.stdX];
                %D.cmpX=[D1.cmpX; D1.cmpX];
            else
                D=obj.table2DVInput(inds,'subj',subj,'B',b);
            end

            % conditions
            [stdXUnq,~,sC]=unique(D.stdX);
            nC=numel(stdXUnq);

            % Prog
            %if isempty(p)
                %p=Pr(nB*nS*nC,1,'Fitting Corr');
            %end

            for ic = 1:nC
                %p.u();
                stdX=stdXUnq(ic);
                ind=sC==ic;


                A1 =sum( D.RCmpChs(ind,1)==D.RCmpChs(ind,2), 1);
                A2 =sum( D.RCmpChs(ind,3)==D.RCmpChs(ind,4), 1);
                A12=sum( sum(D.RCmpChs(ind,1:2),2)==sum(D.RCmpChs(ind,3:4),2) ,1);

                for it = 1:nT
                    typ=typs{it};

                    if bSubj > 0
                        args ={Opts.modelType,typ,subj,subj2,e,b,stdX};
                    else
                        args={Opts.modelType,typ,subj,b,stdX};
                    end

                    % skip if already in
                    Args=[dKeys; args];
                    tTest=T(Args{:});
                    bnew=length(tTest) == 0;
                    if (~bImprove && ~bnew) || (bImprove && strcmp(typ,'boot'))
                        continue
                    end

                    if bSubj > 0
                        args2={Opts.modelType,typ,subj2,subj,e,b,stdX};
                        Args2=[dKeys; args2];
                        tTest2=T(Args2{:});
                        bnew2=length(tTest2)==0;
                        if (~bImprove && ~bnew2) || (bImprove && strcmp(typ,'boot'))
                            tTest2
                            s11=tTest2('subj').ret();
                            s21=tTest2('subj2').ret();
                            tTest2('subj')=s21;
                            tTest2('subj2')=s11;
                            vertcat(T,tTest2);
                            count=count+1
                            continue
                        end
                    end

                    % nlast
                    if bImprove & ~bnew
                        nlast=T(Args{:},'n').ret();
                    else
                        nlast=0;
                    end

                    % get standards
                    stds=obj.Thresh.Data{'subj',unique({subj,subj2}),'stdX',stdX,'bin',b,'aliasInd',aliasInds,'thresh'};
                    stds=stds*60;
                    Opts.stdFix=sqrt(prod(Set.distribute(stds,stds),2));
                    if bSubj==2
                        if aliasInds==1
                            Rfld='R1';
                        elseif aliasInds==2
                            Rfld='R2';
                        else
                            error('');
                        end
                        fixmdl='E';
                        df=obj.Corr.Data('fitType','best','model',fixmdl,'stdX',stdX,'bin',b);
                        df1=df('subj',subj);
                        df2=df('subj',subj2);
                        RS1=df1(1,Rfld).ret();
                        RS2=df2(1,Rfld).ret();
                        Opts.rhoFix=[RS1 RS2];
                    end

                    % RUN
                    obj.dv=DVFitter.new(D.stdX(ind,:)*60,D.cmpX(ind,:)*60,D.RCmpChs(ind,:,:), Opts);
                    obj.dv.init();


                    switch typ
                    case 'test'
                        obj.dv.run_test();
                        % orthant by comparison

                        n=nlast+1;
                        rho=Vec.row(obj.dv.DVFit.RHO(rInds,1));
                        negLL=obj.dv.DVFit.negLLAll;
                        ci=nan(numel(rInds),2);
                    case 'best'
                        obj.dv.run_best();

                        n=nlast+nBest;
                        rho=Vec.row(obj.dv.Fit.RHO.best(rInds,1));
                        negLL=obj.dv.Fit.NegLLAll.best;
                        ci=nan(numel(rInds),2);
                    case 'boot'
                        error('broken. FIX')
                        obj.dv.run_boot();
                        n=nBoot;
                        rho=obj.dv.Fit.RHO.mean(rInds,1);
                        negLL=nan;
                        cperms=[1,3,2];  %TODO GENERALIZE
                        ci=permute(obj.dv.Fit.RHO.ci(rInds,1,:),cperms);
                    end

                    % last
                    obj.lastOut=num2cell(rho);
                    obj.lastArgs=args;

                    % disp
                    if bDisp && ~bBoot
                        dispV(rho,'rhoFit');
                        if contains(model,'E')
                            [varB,covLB,varL,corrLB]=Ext.rhos2vars(rho(1),rho(2),rho(3),stds(1),stds(2));
                            dispV([stds(1) stds(2)],'stds')

                            dispV(sqrt(varB),'stdB');
                            dispV(sqrt(varL),'stdL');
                            dispV(covLB,'covLB');
                            dispV(corrLB,'corrLB');
                        end
                        dispV(negLL,'negLL');
                    end

                    % pack
                    rIn=num2cell(rho);
                    ciIn=num2cell(ci,2);
                    if bnew
                        T.add_row_sorted(args{:},n,negLL,rIn{:},ciIn{:},A1,A12,A2);
                    elseif bImprove
                        lastnLL=T(Args{:},'negLL').ret();
                        dispV(lastnLL,'lastnLL')
                        if negLL < lastnLL
                            T(Args{:},'n')=n;
                            T(Args{:},'negLL')=negLL;

                            for iRF = 1:length(RFlds)
                                T(Args{:},RFlds{iRF})=rho(iRF);
                            end
                            disp('Improved!')
                        end
                    else
                        XXX
                    end
                end

            end
        end
        end
        end
        end
        %obj.corr2ext();
    end
    function get_corr_avg(obj,model)
        ARGS={'model',model,'fiType','best'};
        T=obj.Corr.data(ARGS{:},'subj','~=','AVG');
        stds=T.unique('stdX').ret();
        for i = 1:length(stds)
            t=T('stdX',stds(i));

            n=numel(T.unique(subj));

            m=t.mean({'R1','R12','R2'});
            s=t.std( {'R1','R12','R2'});
            ci=[m-s m+s];

            args=[ARGS 'stdX',stds(i),'subj','AVG'];
            obj.Corr.data.add_row_sorted(args{:},n,nan,m(1),m(2),m(3),ci(1,:),ci(2,:),ci(3,:));
        end
    end
    function init_ext(obj)
        Ext =struct();

        inds={'model','fitType','subj','bin','stdX'};
        vars={'varB','varL','covLB','corrLB','varIB','varIL','varEB','varTB','varTL'};
        Key=[inds vars];

        Typs={'char','char','char','double','double'};
        Typs=[Typs repmat({'double'},1,numel(vars))];
        Ext.Data=Table([],Key,Typs);
        obj.Ext=Ext;
    end
    function corr2extOld(obj)
        obj.init_ext();
        ext=obj.Ext.Data;

        D=obj.Corr.Data;

        typs=D.unique('fitType');
        models=D.unique('model');
        subjs=D.unique('subj','~AVG','subj');
        bins=Vec.row(D.unique('bin'));
        stdX=Vec.row(D.unique('stdX'));

        nM=numel(models);
        nS=numel(subjs);
        nB=numel(bins);
        nL=numel(stdX);
        nT=numel(typs);
        for is=1:nS
        for ib=1:nB
        for il=1:nL
            s=subjs{is};
            b= bins(ib);
            l= stdX(il);
            tB=obj.Thresh.Data{'subj',s,'stdX',l,'bin',b,'aliasInd',1,'thresh'};
            tL=obj.Thresh.Data{'subj',s,'stdX',l,'bin',b,'aliasInd',2,'thresh'};
            for it=1:nT
            for im=1:nM
                m=models{im};
                t=typs{it};

                args={'model',m,'fitType',t,'subj',s,'bin',b,'stdX',l};

                if strcmp('typ','boot')
                    [r1,r12,r2]=D(args{:},'R1B','R12B','R2B').ret();
                else
                    [r1,r12,r2]=D(args{:},'R1','R12','R2').ret();
                end
                [varB,covLB,varL,corrLB,varIB,varIL,varEB]=Ext.rhos2vars(r1,r12,r2,tB,tL);
                if isempty(varB)
                    continue
                end
                args={m t s b l varB varL covLB corrLB varIB varIL varEB tB.^2 tL.^2};
                ext.add_row_sorted(args{:});

            end
            end
        end
        end
        end
    end

    function merge(obj)
        nBoot=1000;
        nBest=50;
        R={'R1','R2','R12'};
        RB={'R1B','R2B','R12B'};
        E={'varB','covLB','varL','corrLB','varIB','varIL','varEB'};
        RE=[R E];
        B=struct();

        T=obj.Thresh.Data.copy();
        C=obj.Corr.Data.copy();
        d=C('fitType','boot','n','>=',nBoot,'model','E');

        %% ADD STDS
        d=T.rows2col('aliasInd','thresh',{'std1','std2'});
        e=C.merge(d,{'subj','bin','stdX'});

        %% BEST
        B.B=e.notOut('fitType','best','n','>=',nBest,RB{:},'fitType');
        B.B=B.B.unique_rows();
        B.B=B.B.rename_col({'n'},{'nBest'});

        %% BOOT
        bUpper=[true false true];
        EBoot=e('fitType','boot','n','>=',nBoot);
        EBoot.le('R1','R1B',{'+','bUpper'});
        c=0;
        for fld = 'MLU'
            c=c+1;
            B.(fld)=EBoot.notOut('bUpper',bUpper(c),'bUpper',R{:},'fitType','negLL');
            B.(fld).rename_col(['n' RB],['nBoot' R]);
        end
        B.M=B.M.unique_rows();
        %% GET EXT
        for fld = 'MLUB'
            t=B.(fld);
            t.fun(@Ext.rhos2vars,5,'R1','R12','R2','std1','std2',['+',E]);
            t.rename_col({'varIB'},{'varI1'});
            t.copy_col(  {'varIL'},{'varI2'});
            t.minus('varI1','varIL',{'+','varIB'});

            t.plus('varL', 'varEB','covLB','covLB',{'+','varE1'});
            t.plus('varI1','varE1',{'+','varT1'});

            t.copy_col('varL','varE2');
            t.plus('varI2', 'varE2', {'+','varT2'});
        end
        for fld='LU'
            B.(fld)=B.(fld).rm_col('nBoot');
        end

        RE=[RE 'varI1' 'varE1' 'varT1' 'varI2','varE2','varT2'];

        %% MERGE
        Exp=B.B;
        for fld = 'MLU'
            newFld=regexprep(RE,'(.*)',['$1_' fld]);
            B.(fld)=B.(fld).rm_col({'std1','std2'});
            B.(fld)=B.(fld).rename_col(RE,newFld);
            Exp.merge(B.(fld),{'model','subj','bin','stdX'});
        end

        %% Convert
        %fun=@(x) sqrt(x.*60.^2);
        Exp=Exp.times('stdX',60);
        Exp=Exp.times({'/(var.*|covLB.*)'},3600);
        Exp=Exp.times({'std1','std2'},60);
        Exp=Exp.sqrt({'/var.*'},{'!','/var(.*)/std$1'});

        RE=regexprep(RE,'^var','std');
        obj.getAvgOld(Exp,[RE]);
        %obj.getAvg(Exp,[R E]);

        obj.Ext.Data=Exp;
    end
    function E=getAvg(obj,E,RE)
    % TODO TEST
        M=regexprep(RE,'(.*)','$1_L');
        S=regexprep(RE,'(.*)','$1_U');
        E.mean(RE,{'|','stdX','bin','model'},{'+','/$1_m'});
        E.std( RE,{'|','stdX','bin','model'},{'+','/$1_s'});
        E.plus( M,S,{'+','/([^_]*)_s/$1_u'});
        E.minus(M,S,{'+','/([^_]*)_s/$1_l'});
    end
    function getAvgOld(obj,Exp,RE)
        %% GET AVG
        models=Exp.unique('model');
        subjs=Exp.unique('subj',{'~AVG','~STD'},'subj');
        bins=Vec.row(Exp.unique('bin'));
        stdX=Vec.row(Exp.unique('stdX'));

        nM=numel(models);
        nS=numel(subjs);
        nB=numel(bins);
        nL=numel(stdX);

        REL=regexprep(RE,'(.*)','$1_L');
        REU=regexprep(RE,'(.*)','$1_U');
        REM=regexprep(RE,'(.*)','$1_M');
        ff={'std1' 'std2' 'negLL', 'nBest', 'nBoot'};

        E=Exp.copy();
        for il=1:nL
        for ib=1:nB
        for im=1:nM
            b= bins(ib);
            l= stdX(il);
            m=models{im};

            D=E('subj','~=','AVG','stdX',l,'bin',b,'model',m);
            if length(D) < 1
                continue
            end
            dm=D.mean(RE);
            ds=D.std(RE);
            dm=[dm{:}];
            ds=[ds{:}];

            lb=num2cell(minus(dm,ds));
            ub=num2cell(plus(dm,ds));
            dm=num2cell(dm);
            C=[RE; dm];
            L=[REL; lb];
            U=[REU; ub];
            M=[REM; dm];
            A=[C(:); L(:); U(:); M(:)]';

            fm=D.mean(ff);
            fs=D.std(ff);
            fm=[fm{:}];
            fs=[fs{:}];
            lb=num2cell(minus(fm,fs));
            ub=num2cell(plus(fm,fs));
            fm=num2cell(fm);
            C=[ff; fm];
            %L=[ffL; lb];
            %U=[ffRU; ub];
            %M=[ffRM; dm];
            A=[A C(:)'];

            flds=['subj','AVG','stdX',l,'bin',b,'model',m,A];
            Exp.add_row(flds{:});
        end
        end
        end
    end
    function d=getExtData(obj,model,bBoot)
        d=obj.Ext.Data('model',model);
        d=d.copy();
        d.sort('subj',{'JDB','DNW','AJL','AVG'});

        %fun=@(x) sqrt(x.*60.^2);
        %d.times({'/(var.*|covLB)'},3600);
        %d.sqrt({'/var.*'},{'!','/var(.*)/std$1'});
        %d.times('stdX',60);
    end
    function d=getCorrData(obj,model)
        D=obj.Corr.Data('model',model);
        d=D.copy();
        d.times('stdX',60);
    end
%- Plot
    function plotCorrSMagr(obj,model)
    end
    function plotCorrSBoth(obj,model,bBoot)
        if nargin < 2; model=[]; end
        if nargin < 3; bBoot=[]; end
        obj.plot_corrS_iter_('both',model,bBoot);
    end
    function plotCorrS(obj,model,bBoot)
        if nargin < 2; model=[]; end
        if nargin < 3; bBoot=[]; end
        obj.plot_corrS_iter_('',model,bBoot);
    end
    function plotCorrSVeri(obj,model,bBoot)
        if nargin < 2; model=[]; end
        if nargin < 3; bBoot=[]; end
        obj.plot_corrS_iter_('veri',model,bBoot);
    end
    function plotCorrSRatio(obj,model,bBoot)
        if nargin < 2; model=[]; end
        if nargin < 3; bBoot=[]; end
        obj.plot_corrS_iter_('ratio',model,bBoot);
    end
    function plotCorrSNegLL(obj,model,bBoot)
        if nargin < 2; model=[]; end
        if nargin < 3; bBoot=[]; end
        obj.plot_corrS_iter_('negLL',model,bBoot);
    end
    function plotCorrSAgree(obj,model,bBoot)
        if nargin < 2; model=[]; end
        if nargin < 3; bBoot=[]; end
        obj.plot_corrS_iter_('agree',model,bBoot);
    end
    function plot_corrS_iter_(obj,typ,model,bBoot)
        if nargin < 2; model=[]; end
        if nargin < 3; bBoot=[]; end
        subjs=obj.CorrS.Data.unique('subj');
        bins=[1 13];
        %for iB = 1:length(bins)
        for iS = 1:length(subjs)
            subj=subjs{iS};
            %bin=bins(iB)
            bin=[];
            obj.plotCorrS_fun_(subj,typ,model,bBoot,bin);
        end
        %end
    end
    function testcoskew(obj)
        A=Set.distribute({'S','P1','W1'},{'S','P2','W2'},{'S','P3','W3'});
        A(~any(ismember(A,{'W1','W2','W3'}),2),:)
    end
    function plotCorrS_fun_(obj,subj,ytyp,model,plotTyp,bin)
        %% decision variable correlation
        if nargin < 3 || isempty(ytyp)
            ytyp='R12';
        end
        if nargin < 4 || isempty(model)
            model='E';
        end
        if nargin < 5 || isempty(plotTyp)
            plotTyp='plot';
        end
        if nargin < 5 || isempty(bin)
            bin=[];
        end
        %plotTyp='hist';

        %plotTyp='hist'; % XXX
        %bSwap=true;


        bSF=false;
        d=obj.getCorrSData(bSF,'subj',subj,'subj2','~=',subj);

        % ytyp
        dat2=[];

        eedges=[];
        switch ytyp
        case 'both'
            dat1={'R12','R1'};
            titl='';
            eedges=linspace(0,.7,9);
        case 'rratio'
            dat1='R12';
            d.rdivide('R12','R1');
            titl='Ratio';
        case 'agree'
            dat1='A12';
            d.rdivide('A12','A1');
            titl='Agree';
        case 'ratio'
            dat1='varSE';
            titl='Ratio';
            eedges=linspace(.25,1.25,10);
        otherwise
            dat1=ytyp;
            titl=ytyp;
            eedges=linspace(0,.7,9);
        end


        % plotTyp
        if ~isempty(bin)
            crits={['(bin==' num2str(bin) ')']};
        else
            crits={'(bin==1)','(bin==13)'};
        end
        dats=Set.distribute(dat1,crits);
        dat=strcat(dats(:,1),dats(:,2));

        X={'stdX(bin==1)'};
        edges={};
        xl=[];
        xlTxt='Disparity (deg)';
        yl=[0 1];
        ylTxt='Correlation';
        switch plotTyp
        case 'boot'
            dat={dat strcat(dat,'_L') strcat(dat,'_U')};
            fun='plotInterv';
        case 'hist'
            X={};
            dat={dat};
            fun='hist';
            edges={eedges};
            xl=yl;
            xlTxt=ylTxt;
            ylTxt='count';
            yl=[0 10];
        otherwise
            dat={dat};
            fun=plotTyp;
        end

        bSwap=false;

        name=['CorrS ' plotTyp ' ' titl ' ' model ' ' subj ' ' 'bin' num2str(bin)];
        args={...
               'aliasInd','subj2',...
               'bSwapYR',bSwap,...
               'Color',{'b','r','b','r'},...
               'xlim',xl,...
               'ylim',yl,...
               'yScale','linear',...
               'xlTxt',xlTxt,...
               'ylTxt',ylTxt,...
               'slTxt',name,...
               'rlTxt',{'Exp. 1','Exp. 2'}};

        Fig.new(name);
        sp=d.(fun)(X{:},dat{:},edges{:},args{:});

    end
    function corrS(obj,bBoot)
        if nargin < 2
            bBoot=[];
        end
        obj.plotCorrS_fun2_('rhoSE','RR',true,false,bBoot)
    end
    function plotCorrS_fun2_(obj,fld,model,bSF,bHist,bBoot)
        if nargin < 2 || isempty(fld)
            fld='rhoSE';
        end
        if nargin < 3|| isempty(model)
            model='RR';
        end
        if nargin < 4|| isempty(bSF)
            bSF=false;
        end
        if nargin < 5 || isempty(bHist)
            bHist=false;
        end
        if nargin < 6 || isempty(bBoot)
            bBoot=false;
        end

        bins=[1 13];
        colors={[.6 .6 .6],'k'};
        xtInd=[1 3 5];

        fld2=[];
        fld3=[];
        xfld='Disparity (arcmin)';
        bRho=false;
        bTxt=true;

        yfld='Correlation';
        switch fld
        case {'rhoSE','rhoSEE'}
            loc='southeast';
            edges=linspace(0,1,9);
            fld2='rhoSE-M';
            fld3=fld
        case {'R12'}
            loc='northwest';
            fld2='rhoR12-SE';
            edges=linspace(0,1,18);
        case {'R11','R22'}
            loc='northwest';
            edges=linspace(0,1,18);
            bRho=true;
        case {'rhoSE1','rhoSE2'}
            yfld='Proportion shared externally driven variance';
            loc='northwest';
            %fld2='rhoMax';
            edges=linspace(0,1,9);
        case {'rhoSEM','rhoSEm','rhoSEMm'}
            yfld='normalized Correlation';
            loc='northwest';
            %fld2='rhoMax';
            edges=linspace(0,1,9);
        end

        if isa(model,'Table')
            d=model;
        else
            d=obj.getCorrSData(bSF,bBoot,'model',model);
        end
        subjs=d.unique('subj');
        subjs2=Set.choosek(2,subjs);
        str=[strcat(subjs2(:,1),{' '},subjs2(:,2))' 'Average','Bin Diff' fld3];

        obins=1;
        ibins=1;
        if bHist
            nm='Hist';
            yl=[0 6];
            xl=[0 1];
            xt=[];
            xtInd=[];
            if bHist==2
                ibins=bins;
            else
                obins=bins;
            end
            loc='northwest';
            t=yfld;
            yfld='count';
            xfld=t;
        else
            nm='Plot';
            %yl=Axis.ylim([0 1],0.1);
            yl=Axis.ylim([0 1],0.1);
            %yl=Axis.ylim([0 1],0.1);
            xt=d.unique('stdX');
            xl=Axis.xlim(xt,0.1);
            yt=[0 .25 .5 .75 1];
            ibins=bins;
        end

        SP=cell(length(obins),1);
        if bBoot
            LS='none'
        else
            LS='-';
        end
        nS=numel(subjs);
        %LS='none';

        nSP=[2,nS+2+~isempty(fld3)];
        edgesD=linspace(-0.6,0.6,12);
        edgesD2=linspace(0,1,12);
        xlD{1}=[-0.65 0.65];
        xlD{2}=yl
        ylD=[0 16];
        xl=repmat({xl},nSP);
        yl=repmat({yl},nSP);
        yt=repmat({yt},nSP);
        xt=repmat({xt},nSP);
        xtInd=repmat({xtInd},nSP);
        ytD=[0:2:16];
        xtD{1}=[-.5 0 .5];
        xtD{2}=[0 .5 1];
        xtIndD=1:3;
        K=1:-1:0
        for k1 = 1:2
        for ik2 = 1:length(K)
            k2=K(ik2);
            xl{k1,end-k2}=xlD{ik2};
            yl{k1,end-k2}=ylD;
            yt{k1,end-k2}=ytD;
            xt{k1,end-k2}=xtD{ik2};
            xtInd{k1,end-k2}=xtIndD;
        end
        end
        %xl{
        %yl

        for oB=1:length(obins)
            if bHist==1
                nm2=[' bin ' num2str(obins(oB))];
            else
                nm2='';
            end
            name=[nm ' ' fld ' ' fld2 nm2 ' ' num2str(bSF)];
            Fig.new(name);



            SP{oB}=SubPlot(nSP,...
                'llMargin',1,...
                'iMargin',[0 1 0 1],...
                'ytMargin',.5,...
                'ylMargin',1,...
                'xlim',xl,...
                'xt',xt,...
                'yt',yt,...
                'xtInd',xtInd,...
                'xxlTxt',str,...
                'xlTxt',xfld,...
                'ylTxt',yfld,...
                'ylim',yl,...
                'slTxt',name,...
                'llTxt',{'Experiment 1','Experiment 2'});

            for iS = 1:length(subjs2)
            for iA = 1:2
                SP{oB}.sel([iA,iS]);
                for iB = 1:length(ibins)
                    if bHist==1
                        bin=obins(oB);
                        col=colors{oB};
                    else
                        bin=ibins(iB);
                        col=colors{iB};
                    end
                    if bRho && iS==2
                        subj2=subjs2(iS,1);
                        subj1=subjs2(iS,2);
                    else
                        subj1=subjs2(iS,1);
                        subj2=subjs2(iS,2);
                    end
                    % XY
                    d1=d('subj',subj1,'subj2',subj2,'aliasInd',iA,'bin',bin);
                    [X,Y]=d1('stdX',fld).ret();
                    % boot
                    if bBoot
                        [L,U]=d1([fld 'L'],[fld 'U']).ret();
                    end

                    % fld2
                    if ~isempty(fld2)
                        [Y2]=d1(fld2).ret();
                    end

                    % fld 3
                    if ~isempty(fld3)
                        [Y3]=d1(fld3).ret();
                    end

                    % hist
                    if bHist==1
                        YD(:,iA,oB,iS)=Y;
                        if ~isempty(fld3)
                            YD3(:,iA,oB,iS)=Y3;
                        end
                    else
                        YD(:,iA,iB,iS)=Y;
                        if ~isempty(fld3)
                            YD3(:,iA,iB,iS)=Y3;
                        end
                    end

                    % main
                    if bHist
                        Plot.histogram(Y,edges,'Color',col,'LineWidth',2); hold on;

                        if ~isempty(fld2)
                            c='k';
                            Plot.histogram(Y2,edges,'Color',c,'LineWidth',2)
                        end
                    else
                        if bBoot
                            Plot.interv(X,U,L,'Color',col,'FaceAlpha',0.4); hold on;
                        end
                        plot(X,Y,'Color',col,'Marker','diamond','MarkerEdgeColor',col,'MarkerFaceColor','white','LineWidth',2,'LineStyle',LS,'MarkerSize',10); hold on;

                        if ~isempty(fld2)
                            plot(X,Y2,':','Color',col,'LineWidth',2)
                        end
                    end
                    if bTxt
                        if iB==1
                            tt='low';
                            ty=.1;
                        else
                            tt='high';
                            ty=.2;
                        end
                        text(4,ty,sprintf('%s mu=%f',tt,mean(Y)));
                    end

                end
            end
            end
            if bHist==1
                if oB==1
                    btxt='Low';
                else
                    btxt='High';
                end
                txt={
                    [btxt ' C_\delta'],...
                    % '',...
                    % '',...
                    % ['Max ' btxt ' C_\delta'],...
                };
            else
                txt={...
                    'Low C_\delta','Max Low C_\delta',...
                    'High C_\delta','Max High C_\delta',...
                };
            end
            SP{oB}.sel(1);
            Figs.setFPos(SP);



            % Average
            for iA = 1:size(YD,2);
                SP{oB}.sel([iA,nS+1])
                for iP = 1:2
                    for iB = size(YD,3):-1:1;
                        col=colors{iB};

                        Ym=mean(YD(:,iA,iB,:),4);
                        Ys=std(YD(:,iA,iB,:),[],4);
                        if iP==1
                            Plot.interv(X,Ym+Ys,Ym-Ys,'Color',col,'FaceAlpha',0.4); hold on;
                        elseif iP==2
                            plot(X,Ym,'Color',col,'Marker','hexagram','MarkerEdgeColor',col,'MarkerFaceColor','white','LineWidth',2,'LineStyle',LS,'MarkerSize',10);
                        end
                    end
                end
            end

            % diff
            %SP{oB}.xlims
            for iA = 1:size(YD,2);
                SP{oB}.sel([iA,nS+2])
                col=colors{1};
                Y=YD(:,iA,2,:)-YD(:,iA,1,:);
                Plot.histogram(Y(:),edgesD,'Color',col,'LineWidth',2)
                %SP{oB}.xlim=xl;
                %SP{oB}.ylim=yl;
            end
            if ~isempty(fld3)
            for iA = 1:size(YD3,2);
                SP{oB}.sel([iA,nS+3])
                for iB = 1:size(YD3,3)
                    col=colors{iB};
                    Y=YD3(:,iA,iB,:);
                    Plot.histogram(Y(:),edgesD2,'Color',col,'LineWidth',2); hold on
                end
                %SP{oB}.xlim=xl;
                %SP{oB}.ylim=yl;
            end
            end
        end

    end
    function d=getCorrSData(obj,bSF,bBoot,varargin)
        if nargin < 2
            bSF=false;
        end
        if nargin < 3
            bBoot=true;
        end
        if bSF
            cfld='CorrSF';
        else
            cfld='CorrS';
        end

        argU={'fitType','model','subj','subj2','E','bin','stdX'};
        argInds={'model','subj','subj2','E','bin','stdX'};

        % best
        d=obj.(cfld).Data.copy();
        dt=d('fitType','best');
        if isempty(dt)
            d=d('fitType','test');
        else
            d=dt;
        end
        [~,inds]=d.unique_rows(argU{:});
        d=d(inds);

        if bBoot
            db=obj.(cfld).Data.copy();
            db=db('fitType','boot');

            %% FIX FOR DUPLICATES
            bFix=true; % NOTE
            if bFix
                bflds={'R1B';'R12B';'R2B'};
                nbflds={'R1L','R1U';'R12L','R12U';'R2L','R2U'};
                nbfldsT=rot90(nbflds',2);

                % add new columns
                ve=nan(length(db),1);
                for i=1:numel(nbflds)
                    db.add_col(nbfldsT{i},ve)
                end

                [B,~,BInds]=db.unique_rows(argU{:});
                nB=length(B);
                for i = 1:nB
                    bind=find(BInds==i);
                    dbi=db.index(bind);
                    if length(dbi) == 1
                        error('only 1 match')
                    end

                    for j=1:length(bflds)
                        v=dbi(bflds{j}).ret();
                        if v(1) > v(2)
                            error('first larger than second')
                        end
                        for k = 1:2
                            db(bind(k),nbflds{j,1})=v(1);
                            db(bind(k),nbflds{j,2})=v(2);
                        end
                    end
                end
            end


            [~,inds]=db.unique_rows(argU{:});
            db=db(inds);

            db.rm_col(['negLL','A1','A2','A12','R1','R2','R12',bflds']);
            db.rename_col('n','nBoot')
            %v=repmat({'both'},length(db),1);
            db('fitType')='both';

            d.rm_col(['fitType',bflds']);
            d.rename_col('n','nBest')

            d.merge(db,argInds);
        end


        t=obj.Thresh.Data.copy();
        d.rename_col('E','aliasInd')

        subjs=d.unique('subj');
        ns=numel(subjs);

        %% APPEND Thresholds
        dn=nan(length(d),1);
        d.add_col('T1',dn,'end');
        d.add_col('T2',dn,'end');

        for ii = 1:2
            if ii == 1
                sfld='subj';
                Tfld='T1';
            else
                t.rename_col('subj','subj2');
                sfld='subj2';
                Tfld='T2';
            end
            d.newIter('aliasInd','bin','stdX',sfld);
            for i=1:d.nIter()
                [~,ind,args]=d.iter();
                T=t(args{:},'thresh').ret();
                d(ind,Tfld)=T;
            end
        end

        %% COMBINE R1 R2
        fflds={'R12','R1','R2'};
        rflds={'R12','R2','R1'};
        nF=length(fflds);

        d1=d.copy();
        d2=d.copy();
        d1.newIter('aliasInd','bin','stdX','subj','subj2');
        d2.newIter('aliasInd','bin','stdX','subj2','subj');
        for i=1:d1.nIter()
            [dd1,ind1]=d1.iter();
            [dd2,ind2]=d2.iter();
            for iF=1:nF
                ffld=fflds{iF};
                rfld=rflds{iF};
                nn1=dd1('negLL').ret();
                nn2=dd2('negLL').ret();
                f1=dd1(ffld).ret();
                f2=dd2(rfld).ret();
                if nn1 < nn2
                    f=f1;
                elseif nn1 > nn2
                    f=f2;
                else
                    f=mean([f1 f2],2);
                end
                d(ind1,ffld)=f;
                d(ind2,rfld)=f;
            end
        end
        % change units
        d.times('stdX',60);
        d.times('T1',60);
        d.times('T2',60);

        %d1=d('subj','JDB','subj2','DNW','stdX',3.75,'aliasInd',2,'bin',1);
        %d2=d('subj','DNW','subj2','JDB','stdX',3.75,'aliasInd',2,'bin',1);

        if bBoot
            n=0:2;
        else
            n=0;
        end
        for nB = n
            d=append_fun(d,nB);
        end


        % get S E P SE
        function d=append_fun(d,nBoot)
            flds={'R1','R12','R2'};
            if nBoot==0
                bChar='';
            elseif nBoot==1
                bChar='L'
            elseif nBoot==2
                bChar='U'
            end
            flds=strcat(flds,{bChar});
            [R1,R12,R2,T1,T2]=d(flds{:},'T1','T2').ret();

            % E
            E1=R1.*T1.^2;
            E2=R2.*T2.^2;
            EM=min([E1 E2],[],2);
            Em=max([E1 E2],[],2);
            ET=sqrt(E1.*E2);

            % S
            S=R12.*T1.*T2;

            % P
            P1=E1-S;
            P2=E2-S;


            SE1=S./E1;
            SE2=S./E2;
            SE =S./ET;
            SEM=S./EM;
            SEm=S./Em;
            SEMm=(S-Em)./(EM-Em);
            SEE=S./(S+P1.*P2);
            SE_M=SE./SEM;
            SE_R12=R12./SE;

            R =ET./(T1.*T2);
            RM=EM./(T1.*T2);
            Rm=Em./(T1.*T2);

            %d.add_col(['varE1' bChar], E1,'end');
            %d.add_col(['varE2' bChar], E2,'end');
            % d.add_col(['varE' bChar], ET,'end');
            %d.add_col(['varEM' bChar], EM,'end');
            %d.add_col(['varEm' bChar], Em,'end');

            %d.add_col(['varS12' bChar], S,'end');
            % d.add_col(['varP1' bChar], P1,'end');
            % d.add_col(['varP2' bChar], P2,'end');

             d.add_col(['rhoSE1' bChar], SE1,'end');
             d.add_col(['rhoSE2' bChar], SE2,'end');
              d.add_col(['rhoSE' bChar], SE,'end');
             d.add_col(['rhoSEM' bChar], SEM,'end');
            % d.add_col(['rhoSEm' bChar], SEm,'end');
            %d.add_col(['rhoSEMm' bChar], SEMm,'end');
            % d.add_col(['rhoSEE' bChar], SEE,'end');

              d.add_col(['rhoSE-M' bChar], SE_M,'end');
            d.add_col(['rhoR12-SE' bChar], SE_R12,'end');

             d.add_col(['R11' bChar], R1,'end');
             d.add_col(['R22' bChar], R2,'end');
             d.add_col(['rho' bChar], R,'end');
            d.add_col(['rhoM' bChar], RM,'end');
            d.add_col(['rhom' bChar], Rm,'end');

            %subjs=d.unique_rows('subj').ret();
            %KEY='S1','S2'
            %T=Table([],KEY,types]
            %S12=zeros(0,3,2);
            %jP12=zeros(0,3,2);
            %%1= zeros(0,3,2);
        end

    end

    function venn(obj,d)

        subjs=d.unique('subj');
        nS=length(subjs);
        %S12=zeros(:,nS,nS);
        %P12=zeros(:,nS,nS);
        %E=zeros(:,nS,nS);
        for i = 1:nS
        for j = 1:nS
            if i == j
                continue
            end
            subj1=subjs{i};
            subj2=subjs{j};
            args1={'subj1',subj1,'subj2',subj1};
            args2={'subj1',subj2,'subj2',subj2};

            S=d(args1{:},'S');
            if isempty(S)
                args=args2;
                PArg='P2'
                EArg='E2';
            else
                args=args1;
                PArg='P1'
                EArg='E1';
            end

            [s,p,e]=d(args{:},'S',PArg,EArg);
            S12(:,i,j)=s;
            P12(:,i,j)=p;
            E(:,i,j)=e;
        end
        end

        inds=1:3;
        PS=zeros(size(E));
        for i = 1:nS
        for j = 1:nS
            k=inds(~ismember(inds,[i,j]));
            PS1(:,i,j)=P12(:,i,j)-S12(:,i,k);
            PS2=P(:,i,j)-S12(:,i,k)
        end
        end

        for i = 1:nS
            PS(:,i,1)
        end

        %if nargin > 1
        %    d=d(varargin{:});
        %end
    end
    function SP=plotAbs(obj,model,num,bSwap,bBoot)
        if nargin < 2 || isempty(model)
            model='E0';
        end
        if nargin < 3 || isempty(model)
            num=1;
        end
        if nargin < 4 || isempty(bSwap)
            bSwap=true;
        end
        if nargin < 5 || isempty(bBoot)
            bBoot=true;
        end
        if bSwap
            Color={[.5 .5 .5],'k','y','r','b'};
        else
            Color={'b','r','m'};
        end
        bLog=true;
        name=['Abs ' model ' T' num2str(num) ' ' num2str(bSwap)];

        d=obj.getExtData(model,bBoot);
        if bSwap
            llTxt={'Total','External','Internal'};
        else
            llTxt={'low C\delta','high C\delta'};
        end

        if bLog
            if num==1
                yl=[.4 4.5];
            elseif num==2
                yl=[.17 2.5];
            end
        else
            yl=[0 3.7];
        end

        dat={'stdT','stdE','stdI'};
        dat=strcat(dat,num2str(num));
        if bBoot
            dat={dat strcat(dat,'_L') strcat(dat,'_U')};
        end
        Fig.new(name);
        SP=d.plotInterv('stdX',...
               dat{:},...
               'bin','subj',...
               'bSwapYR',bSwap,...
               'Color',Color,...
               'LineWidth',1,...
               'llTxt',llTxt,...
               'rlOn',false,...
               'xtInds',[1 3 5],...
               'yt',[1 2 3 4],...
                'bLog',[false true],...
               'ylim',yl,'xlTxt','Disparity (arcmin)','ylTxt','Threshold (arcmin)','slTxt',name);
        SP.sel(1,1);
        if ~bSwap
            legend({'\sigma L','\sigma B','cov LB'},'location','northwst');
        else
            legend({'low C\delta','high C\delta'},'location','northwest');
        end
    end
    function plotAbsT(obj,model,bBoot)
        if nargin < 2 || isempty(model)
            model=[];
        end
        if nargin < 3
            bBoot=false;
        end
        SP=obj.plotAbsAlt_fun(model,'T',bBoot);
    end
    function plotAbsE(obj,model,bBoot)
        if nargin < 2 || isempty(model)
            model=[];
        end
        if nargin < 3
            bBoot=true;
        end
        SP=obj.plotAbsAlt_fun(model,'E',bBoot);
    end
    function plotAbsI(obj,model,bBoot)
        if nargin < 2 || isempty(model)
            model=[];
        end
        if nargin < 3
            bBoot=true;
        end
        SP=obj.plotAbsAlt_fun(model,'I',bBoot);
    end
    function SP=plotAbsAlt_fun(obj,model,flds,bBoot)
        if nargin < 2 || isempty(model)
            model='E0';
        end
        if nargin < 3 || isempty(flds)
            flds='T';
        end
        if nargin < 4 || isempty(bBoot)
            bBoot=true;
        end
        bSwap=false;
        Color={'k',[.5 .5 .5],'r',[1 .5 .5],[0 0 .5],'b'};
        bLog=true;

        d=obj.getExtData(model,bBoot);
        if bSwap
            llTxt={'Total','External','Internal'};
        else
            llTxt={'low C\delta','high C\delta'};
        end
        llTxt='';

        if bLog
            yl=[.17 4.5];
        else
            yl=[0 3.7];
        end

        marker='o';
        if ~iscell(flds)
            switch flds
            case 'T'
                yl=[.5 4];
                marker='none';
            case 'I'
                yl=[.5 3];
            case 'E'
                yl=[.2 3];
            end

            flds={flds};
        end
        name=['AbsAlt ' model ' ' strjoin(flds)];

        dat=flds;
        ind=strcmp(dat,'T');
        if any(ind)
            dat{ind}='std';
        end
        dat(~ind)=strcat('std',dat(~ind));

        dat1=strcat(dat,num2str(1));
        dat2=strcat(dat,num2str(2));
        dat=[dat1 dat2]';
        dat1=strcat(dat,'(bin==13)');
        dat2=strcat(dat,'(bin==1)');
        dat=[dat1 dat2]';
        dat=dat(:)';
        if bBoot
            dat={dat strcat(dat,'_L') strcat(dat,'_U')};
            fun='plotInterv';
        else
            dat={dat(:)'};
            fun='plot';
        end
        Fig.new(name);
        SP=d.(fun)('stdX(bin==1)',...
               dat{:},...
               '',...
               'subj',...
               'bSwapYR',bSwap,...
               'Color',Color,...
               'Marker',marker,...
               'LineWidth',1,...
               'rlOn',false,...
               'xtInds',[1 3 5],...
               'yt',[1 2 3 4],...
               'bLog',[false true],...
               'ylim',yl,'xlTxt','Disparity (arcmin)','ylTxt','Threshold (arcmin)','slTxt',name);
        SP.sel(1,1);
        %if ~bSwap
        %    legend({'\sigma L','\sigma B','cov LB'},'location','northwst');
        %else
            legend({'NH','NL','FH','FL'},'location','northwest');
        %end
    end
    function plotAllExt(obj)
        SP(1)=obj.plotExt('E',false,true);
        SP(2)=obj.plotExt('E',true,true);
        SP(3)=obj.plotExt('E0',false,true);
        SP(4)=obj.plotExt('E0',true,true);

        drawnow;
        for i = 1:length(SP)
            SP(i).setSelWidthPix(200);
        end
    end
    function SP=plotExt(obj,model,bSwap,bBoot,bInt)
        if nargin < 2 || isempty(model)
            model='E0';
        end
        if nargin < 3 || isempty(bSwap)
            bSwap=false;
        end
        if nargin < 4 || isempty(bBoot)
            bBoot=true;
        end
        if nargin < 5 || isempty(bInt)
            bInt=false;
        end
        if bSwap
            Color={[.5 .5 .5],'k','y','r','b'};
            %Color={'b','b','r','
        else
            Color={'b','r','m'};
        end
        name=['Ext ' model ' ' num2str(bSwap)];
        b0=contains(model,'0');

        d=obj.getExtData(model);
        Fig.new(name);

        yl=[-.2 2];
        dat={'stdL','stdB'};
        if ~b0
            dat=[dat 'covLB'];
            yl(1)=-1;
        end
        if bInt
            dat=[dat 'stdIL' 'stdI1'];
            yl(2)=3.0;

            %dat=[dat 'stdIL' 'stdIB'];
            %yl=[-1.0 2.4];
        end
        tlTxt={'Observer 1','Observer 2','Observer 3','Average'};
        if bSwap
            rlTxt={'\sigma L','\sigma B'};
            if ~b0
                rlTxt=[rlTxt 'cov[L,B]'];
            end
            if bInt
                rlTxt=[rlTxt '\sigma I+','\sigma I*'];
            end
        else
            rlTxt={'low C\delta','high C\delta'};
        end
        if bBoot
            dat={dat strcat(dat,'_L') strcat(dat,'_U')};
        end
        SP=d.plotInterv('stdX',...
               dat{:},...
               'bin','subj',...
               'bSwapYR',bSwap,...
               'Color',Color,...
               'LineStyle','none',...
               'LineWidth',1,...
               'xtInds',[1 3 5],...
               'tlTxt',tlTxt,...
               'xlTxt','Disparity (arcmin)',...
               'ylim',yl,...
               'ylTxt','Threshold contribution (arcmin)',...
               'rlTxt',rlTxt,...
               'slTxt',name);
        SP.sel(1,1);
        if ~bSwap && b0
            legend({'\sigma_L','\sigma_B'});
        elseif ~bSwap
            legend({'\sigma L','\sigma B','cov LB'});
        else
            legend({'low C\delta','high C\delta'});
        end
        Figs.setFPos(SP);
    end
    function plotCmpAll(obj,model,bLog)
        if nargin < 2  || isempty(model)
            model='E0';
        end
        if nargin < 3 || isempty(bLog)
            bLog=true;
        end
        bSwap=true;
        name=['CmpAll ' model];

        if bLog
            YScale='log';
            ylim=[.3 3];
            yt=[.3 1 3];
        else
            YScale='linear';
            ylim=[0 3];
            yt=[0 1 2 3];
        end

        Color='k';

        D=obj.getExtData(model);
        D=D('subj','~=','AVG');
        s1={'std1','stdE1','stdI1'};
        s2={'std2','stdE2','stdI2'};

        dCF=cell(2,2);

        M=zeros(2,2,3);
        L=zeros(2,2,3);
        H=zeros(2,2,3);
        for i = 1:2
            if i==1
                d=D('bin',1);
            else
                d=D('bin',13);
            end
            for j = 1:2
                if j==1
                    s=s2;
                else
                    s=s1;
                end
                m=d.mean(s);
                s=d.std(s);
                m=[m{:}];
                s=[s{:}];

                M(i,j,:)=m;
                L(i,j,:)=m-s;
                H(i,j,:)=m+s;
            end
        end

        Fig.new(name);
        Opts={...
            'tlTxt',{'Total','External','Internal'};
            'xt',1:4;
            'ylTxt','Threshold (arcmin)';
            'xtTxt',{'FL','FH','NL','NH'};
            'ylim',ylim;
            'yt',yt;
            'iMargin',[0 1.5 0 0];
            'yScale',YScale;
        };
        Opts=Opts';
        SP=SubPlot([1 3],Opts{:});

        X=(1:4)';
        for k = 1:3
            l=L(:,:,k);
            h=H(:,:,k);
            SP.sel(k);
            Plot.interv(X,l(:),h(:),'EdgeColor','none');
        end
        for k = 1:3
            m=M(:,:,k);
            SP.sel(k);
            hold on;
            plot(X,m(:),'k','marker','o','MarkerFace','k');
        end

        %d.mean(

        %SP=d.(fun)('stdX',...
        %       dat{:},...
        %       'bin','subj',...
        %       'bSwapYR',bSwap,...
        %       'Color',Color,...
        %       'LineWidth',1,...
        %       'xtInds',[1 3 5],...
        %       'bLogXY',[false bLog],...
        %       'ylim',yl,'xlTxt','Disparity (arcmin)','ylTxt','Threshold (arcmin)','slTxt',name);
        %SP.sel(1,1);
    end
    function plotCVC(obj,model)
        %% decision variable correlation
        if nargin < 2 || isempty(model)
            model='E';
        end
        name=['Ext CvC ' model];
        Emp=obj.Emp;
        d=obj.getExtData(model);
        Subjs={'JDB','DNW','AJL'};
        B=[1 13];
        Fig.new(name); hold off
        SP=SubPlot([1 3],'xlim',[0 1],'ylim',[0 1],'xlTxt','Marginal fit','ylTxt','Joint fit','slTxt',name,'tlTxt',{['Correlation'],'External','Internal'});

        M='so^';
        pink=[0.9490    0.6275    0.6314];
        green=[ 0.4157    0.6039    0.4902];
        gray=[.8 .8 .8];
        vORs='std';
        bStart=true;
        for ib = 2:-1:1
            Emp.b=ib;
            b=B(ib);
            for is = 1:3
                Emp.s=is;
                s=Subjs{is};
                m=M(is);

                % RHO
                r22=Emp.get_val(false,false,'rho22');
                r33=Emp.get_val(false,false,'rho33');
                r2=r22(:,:,1,1);
                r1=r33(:,:,1,1);

                % E
                e22=Emp.get_val(false,false,[vORs 'E2']);
                e33=Emp.get_val(false,false,[vORs 'E3']);
                e2=e22(:,:,1,1);
                e1=e33(:,:,1,1);

                % I
                i22=Emp.get_val(false,false,[vORs 'I2']);
                i33=Emp.get_val(false,false,[vORs 'I3']);
                i2=i22(:,:,1,1);
                i1=i33(:,:,1,1);


                R1=d('subj',s,'R1','bin',b).ret();
                R2=d('subj',s,'R2','bin',b).ret();

                E1=d('subj',s,'stdE1','bin',b).ret();
                E2=d('subj',s,'stdE2','bin',b).ret();

                I1=d('subj',s,'stdI1','bin',b).ret();
                I2=d('subj',s,'stdI2','bin',b).ret();
                if strcmp(vORs,'var')
                    E1=E1.^2;
                    E2=E2.^2;
                    I1=I2.^2;
                    I2=I2.^2;
                end

                SP.sel(1);
                if bStart
                    lims=[.18,4.5];
                    plot(lims,lims,'k:');
                end
                plot(r2,R2,['k' m],'MarkerFaceColor',gray-(ib-1).*.4,'MarkerSize',8); hold on;
                plot(r1,R1,['k' m],'MarkerFaceColor',gray-(ib-1).*.4,'MarkerSize',8); hold on;

                SP.sel(2);
                if bStart
                    lims=[.18,4.5];
                    plot(lims,lims,'k:');
                end
                plot(e2,E2,['k' m],'MarkerFaceColor',green-(ib-1)*.3,'MarkerSize',8); hold on;
                plot(e1,E1,['k' m],'MarkerFaceColor',green-(ib-1)*.3,'MarkerSize',8); hold on;

                SP.sel(3);
                if bStart
                    lims=[0 1];
                    plot(lims,lims,'k:');
                end
                plot(i2,I2,['k' m],'MarkerFaceColor',pink-(ib-1)*.3,'MarkerSize',8); hold on;
                plot(i1,I1,['k' m],'MarkerFaceColor',pink-(ib-1)*.3,'MarkerSize',8); hold on;

                bStart=false;

            end
        end
        for i = 1:3
            SP.sel(i);
            if i  > 1
                lims=[0.18 4.5];
                xt=[.25 .5 1 2 4];

                set(gca,'YScale','Log');
                set(gca,'XScale','Log');
                xlim(lims);
                ylim(lims);
                xticks(xt);
                yticks(xt);
                if i ==2
                    yticklabels(xt);
                end
                xticklabels(xt);
            else
                legend('\rho^{**}','\rho^{++}','location','northwest');
                lims=[0 1];
                xt=[0 .5 1];
            end
        end
        Figs.setFPos(SP);
    end
    function SP=plotTvT(obj,model,bSwap,bBoot)
        if nargin < 2 || isempty(model)
            model='E0';
        end
        if nargin < 3 || isempty(bSwap)
            bSwap=true;
        end
        if nargin < 4 || isempty(bBoot)
            bBoot=true;
        end
        bDat=false;
        bLog=false;
        name=['TvT ' model num2str(bSwap) num2str(bDat)];
        d=obj.getExtData(model);
        yl=[];



        if bDat
            dat={'std1','std2'};
        else
            dat={'stdT1','stdT2'};
        end
        if bBoot
            dat={dat strcat(dat,'_L') strcat(dat,'_U')};
        else
            dat={dat};
        end

        if bSwap
            Color={[.5 .5 .5],'k','y','r','b'};
        else
            Color={'r','b','m'};
        end


        dat{:}
        Fig.new(name);

        if ~bBoot
            fun='plot';
        else
            fun='plotInterv';
        end

        SP=d.(fun)('stdX',...
               dat{:},...
               'bin','subj',...
               'bSwapYR',bSwap,...
               'Color',Color,...
               'LineWidth',1,...
               'xtInds',[1 3 5],...
               'bLogXY',[false bLog],...
               'ylim',yl,'xlTxt','Disparity (arcmin)','ylTxt','Threshold (arcmin)','slTxt',name);
        SP.sel(1,1);
        if ~bSwap
            SP.sel(1,1);
            legend('\sigma T*','\sigma T+');
        else
            legend('C \delta Low','C \delta High');
        end

    end
    function scatterExt(obj,model,bLog)
        if nargin < 2 || isempty(model)
            model='E';
        end
        if nargin < 3 || isempty(bLog)
            bLog=true;
        end
        name=['Hist IB' model ];

        %% between L and B
        name=['Hist LB' model ];
        d=obj.getExtData(model);
        X={'stdL(bin==1)','stdL(bin==13)'};
        Y={'stdB(bin==1)','stdB(bin==13)'};
        xyl=[0.3 2];
        Fig.new(name);
        d.scatter(X,Y,[],[],'Color',{'k','r'},'xylim',xyl,'slTxt',name,'bLogXY',bLog);

        %% between I and E
        name=['Hist IE' model ];
        X={'stdI1(bin==1)','stdIL(bin==13)'};
        Y={'stdEB(bin==1)','stdL(bin==13)'};
        xyl=[0.3 2.2];

        Fig.new(name);
        sp=d.scatter(X,Y,[],[],'Color',{'k','r'},'xylim',xyl,'slTxt',name,'bLogXY',bLog,'xlTxt','Internal Noise','ylTxt','External noise');
        sp.sel(1);
        legend({'bin 1','bin 13'});
    end
    function histExtModels(obj)
        models={'E','E0'};
        strs={...
            'stdL(bin==13)-stdL(bin==1)';
            'stdB(bin==13)-stdB(bin==1)';
        };
    end
    function histExt(obj,model)
        if nargin < 2
            model=[];
        end
        edges=linspace(-0.4,1.6,12);
		xlims=[-1.3 1.3];

        obj.hist_fun_(model,'ext',edges,xlims);
    end
    function histInt(obj,model)
        if nargin < 2
            model=[];
        end
        obj.hist_fun_(model,'int');
    end
    function histVar(obj,model,bBoth)
        if nargin < 2 || isempty(model)
            model='E';
        end
        if nargin < 3  || isempty(bBoth)
            bBoth=true;
        end
        v=getVar();
        switch model
        case 'E0'
            m=1;
            nbin=5;
        case 'E'
            m=2;
            nbin=[];
        end
        name=['Hist Var ' model ];

        stdL=v(:,:,m,1);
        stdB=v(:,:,m,2);

        [N,edges]=histcounts(v(:),nbin);

        Fig.new(name);

        if bBoth
            n=2;
            tlTxt={'C\delta Low','C\delta High'};
            ylim=[0,10];

        else
            n=1;
            tlTxt=[];
            ylim=[0,max(N)*1.1/4];
        end
        SP=SubPlot([1 n],'slTxt',name,'xlim',[-1 1],'xt',[-1 -.5 0 .5 1],'ylim',ylim,'ylTxt','count','xlTxt','\Delta Threshold Contribution','xlCtr',true,'tlTxt',tlTxt);

        if bBoth
            SP.sel(1);

            L=histcounts(stdL(:,1),edges);
            B=histcounts(stdB(:,1),edges);
            Plot.bar([],edges,L,'Color','b'); hold on
            Plot.bar([],edges,B,'Color','r');

            SP.sel(2);
            L=histcounts(stdL(:,2),edges);
            B=histcounts(stdB(:,2),edges);
            Plot.bar([],edges,L,'Color','b'); hold on
            Plot.bar([],edges,B,'Color','r');

        else
            SP.sel(1);
            L=histcounts(stdL(:),edges);
            B=histcounts(stdB(:),edges);
            Plot.bar([],edges,L,'Color','b'); hold on
            Plot.bar([],edges,B,'Color','r');
        end

        SP.sel(1);
        legend('','\sigma_L','','\sigma_B','location','northeast');
        Figs.setFPos(SP);

        function out=getVar()
            D=obj.Ext.Data.copy();
            models={'E','E0'};
            bins=[1 13];
            flds={'stdL','stdB'};
            out=zeros(20,length(bins),length(models),length(flds));

            for iF = 1:length(flds)
            for im = 1:length(models)
            for ib = 1:length(bins)
                d=D('model',models(im),'bin',bins(ib),flds{iF}).ret();
                out(:,ib,im,iF)=d-mean(d);

            end
            end
            end
        end
    end
    function histCov(obj)
        bSwap=false;
        d=obj.getExtData('E');
        %C=d('covLB').ret();
        %s=sign(C);
        %d('covLB')=s.*sqrt(abs(2*C));
        Y='covLB';
        bSwap=false;
        colors='k';

        name='HistCov';

        xlTxt='';
        xlims=[-1.3 1.3];

        edges=linspace(-1,1.6,18);
        nBins=[];
        Fig.new(name);

        mean(d('covLB').ret())
        std(d('covLB').ret())
        -.1125
        .2284
        SP=d.hist(Y,edges,[],[],'bSwapYR',bSwap,'slTxt',name,'FaceColor',colors,'nBins',nBins,'xlTxt',xlTxt,'xlim',xlims,'FaceColor','k','EdgeColor','k');
        Figs.setFPos(SP);
        %SP.sel(1);
        %legend(leg);
    end
    function hist_fun_(obj,model,typ,nORedges,xlims)
        if nargin < 2 || isempty(model)
            model='E0';
        end
        nBins=[];
        edges=[];
        if nargin < 4
            nORedges=[];
        end
        if numel(nORedges)==1
            nBins=nORedges;
        elseif numel(nORedges) > 1
            edges=nORedges;
        end
        if nargin < 5
            xlims=[];
        end

        bSwap=false;

        Name=['HistExt ' model];
        d=obj.getExtData(model);
        xlTxt='\Delta Threshold Contribution';
        colors={'b','r','m'};
        nBins=[];


        %% Change in factor with contrast
        switch typ
        case 'ext'
            name=[Name ' LB'];
            strs={...
                'stdL(bin==13)-stdL(bin==1)';
                'stdB(bin==13)-stdB(bin==1)';
            };
            strs2={...
                'covLB(bin==13)-covLB(bin==1)';
                %'corrLB(bin==13)-corrLB(bin==1)';
            };
            leg={'','\sigma_L','','\sigma_B'};
            if ~contains(model,'0')
                Y=[strs; strs2]';
            else
                Y=strs';
                leg=[leg 'cov[LB]'];
            end
        case 'int'
            name=[Name ' I'];
            if isempty(nORedges)
                nBins=7;
            end
            leg={'','\sigma^+_I','','\sigma^*_I'};

            Y={...
                'stdIL(bin==13) - stdIL(bin==1)';
                'stdI1(bin==13) - stdI1(bin==1)';
            }';
        end

        %% Change in noise with contrast
        colors

        Fig.new(name);
        SP=d.hist(Y,edges,[],[],'bSwapYR',bSwap,'slTxt',name,'FaceColor',colors,'nBins',nBins,'xlTxt',xlTxt,'xlim',xlims,'FaceColor','k','EdgeColor',{'b','r','m'});
        SP.sel(1);
        legend(leg);



    end
    function histCmp(obj,model)
        if nargin < 2 || isempty(model)
            model='E0';
        end
        FLDS=[];
        AB={...
            %'FH','NL';...
            %'FL','NL';...
            'FH','NH'...
        };
        for i = 1:size(AB,1)
            obj.histCmp_fun(model,AB(i,1),AB(i,2),FLDS);
        end
    end
    function histCmp_fun(obj,model,A,B,FLDS)
        bSwap=false;
        Name=['HistCmp ' model];

        d=obj.getExtData(model);
        if nargin < 3 || isempty(A)
            A={'FH'};
        end
        if nargin < 4 || isempty(B)
            A={'NL'};
        end
        if ~iscell(A)
            A={A};
        end
        if ~iscell(B)
            B={B};
        end
        if nargin < 5 || isempty(FLDS)
            FLDS={'I','E','T'};
        end
        if ~iscell(FLDS)
            FLDS={FLDS};
        end
        AB=[Vec.col(A) Vec.col(B)];
        for i = 1:length(FLDS);
            histfun(AB(:,1),AB(:,2),FLDS{i});
        end
        %% Change in noise with contrast

        function SP=histfun(A,B,fld);
            name=[Name ' ' fld];
            if numel(A)==1
                name=[name ' ' B{1} ' - ' A{1}];
            end
            Fig.new(name);
            nBins=6;
            xlTxt='\Delta Threshold Contribution';

            switch fld
            case 'I'
                bRev=false;
                colors=[0 0 1];
            case 'E'
                bRev=false;
                colors=[1 0 0];
            case 'T'
                bRev=true;
                colors=[1 1 1];
            end
            grd=linspace(.1,1,numel(A));
            colors=arrayfun(@(x) x.*colors,grd,'UniformOutput',false);
            if bRev
                colors=cellfun(@(x) 1-x,colors,'UniformOutput',false);
            end

            Y={};
            lbl={};
            for j = 1:length(A)
                [strA,lstr]=strfun(A{j},fld);
                [strB,~]=strfun(B{j},fld);
                str=[strB '-' strA];
                Y=[Y; str];
                lstr=[lstr ' ' B{j} ' - ' A{j}];
                lbl=[lbl,{''},lstr];
            end


            SP=d.hist(Y,[],[],[],'bSwapYR',bSwap,'slTxt',name,'Color',colors,'nBins',nBins,'xlTxt',xlTxt);
            SP.sel(1);
            legend(lbl);
        end
        function [YSTR,LSTR]=strfun(in,fld)
            LSTR=['\sigma_' fld];
            YSTR='std';
            if ~strcmp(fld,'T')
                YSTR=[YSTR fld];
            end
            for i = 1:2
                switch in(i)
                case 'N'
                    str='1';
                    LSTR=[LSTR '^*'];
                case 'F'
                    str='2';
                    LSTR=[LSTR '^+'];
                case 'L'
                    str='(bin==1)';
                    lstr='';
                case 'H'
                    str='(bin==13)';
                end
                YSTR=[YSTR str];
            end
        end
    end
    function histModelDiff(obj,bBoth,nBins)
        if nargin < 2
            bBoth=false;
        end
        bSwap=false;
        name=['Hist model diff'];
        if bBoth
            C='bin';
            name=[name ' both'];
        else
            C=[];
        end
        colors={'b','r'};
        Y={...
            'stdL(model==E0)-stdL(model==E)';
            'stdB(model==E0)-stdB(model==E)';
        };
        %xlTxt={'stdL','stdI1'};
        %xlTxt={'',''};
        xlTxt='\Delta Threshold Contribution';
        tlTxt={'C\delta Low','C\delta High'};
        d=obj.Ext.Data;
        xt=[-.4 -.2 0 .2 .4];
        xlim=[-.5 .5];

        %% change in noise with model
        Fig.new(name);
        SP=d.hist(Y,[],[],C,'bSwapYR',bSwap,'slTxt',name,'Color',colors,'nBins',nBins,'xlTxt',xlTxt,'xlCtr',false,'xt',xt,'xlim',xlim,'tlTxt',tlTxt,'xlCtr',true);
        SP.sel(1);
        legend('','\sigma_L','','\sigma_B','location','northeast');

        Figs.setFPos(SP);
    end
    function plotThresh(obj,model,bBoot)
        %% decision variable correlation
        if nargin < 2 || isempty(model)
            model='E0';
        end
        if nargin < 3 || isempty(bBoot)
            bBoot=false;
        end
        bSwap=true;
        d=obj.getExtData(model);
        %d=d('fitType','best');

        name=['Thresh ' model];
        Fig.new(name);

        dat={'std1','std2'};
        if bBoot
            fun='plotInterv';
            dat={dat strcat(dat,'_L') strcat(dat,'_U')};
        else
            fun='plot';
            dat={dat};
        end
        RTxt={'Natural','Flat'};
        %RTxtB='Depth Profile';
        RTxtB='';
        YTxt={'Low','Flat'};
        YTxtB='C_\delta';
        if bSwap
            color={[0.5 0.5 0.5],'k'};
            rlTxt =RTxt;
            rlTxtB=RTxtB;
            lTxt=YTxt;
        else
            color={'r','b'};
            rlTxt =YTxt;
            rlTxtB=YTxtB;
            lTxt=RTxt;
        end
        SP=d.(fun)('stdX',dat{:},'bin','subj',...
               'bSwapYR',bSwap,...
               'Marker','none',...
               'Color',color,...
               'yScale','log',...
               'ylim',[.18 4.5],...
               'yt',[.25 .5 1 2 4],...
               'xtInd',[1 3 5],...
               'rlTxt',rlTxt,...
               'rlTxtB',rlTxtB,...
               'xlTxt','Disparity (arcmin)','ylTxt','Threshold (arcmin)','slTxt',name);
        SP.sel(1);
        legend(lTxt,'location','northwest');
    end
    function plotCorr(obj,model,bBoot)
        %% decision variable correlation
        if nargin < 2 || isempty(model)
            model='E0';
        end
        if nargin < 3 || isempty(bBoot)
            bBoot=true;
        end
        d=obj.getExtData(model);
        %d=d('fitType','best');

        name=['Corr ' model];
        Fig.new(name);

        dat={'R1','R12','R2'};
        if bBoot
            dat={dat strcat(dat,'_L') strcat(dat,'_U')};
        else
            dat={dat};
        end
        if ~bBoot
            fun='plot';
        else
            fun='plotInterv';
        end
        d.(fun)('stdX',dat{:},'bin','subj',...
               'bSwapYR',true,...
               'Color',{'r','m','b'},...
               'yScale','linear',...
               'ylim',[0 1],'xlTxt','Disparity (arcmin)','ylTxt','Correlation','slTxt',name);

    end
    function plotNegLL(obj)
        d=obj.Ext.Data;
        %d=obj.getExtData(model);

        name=['NegLL models'];
        Y={...
            'negLL(model==E)';
            'negLL(model==E0)';
        };
        %Y='negLL';
        color={'m','k'};

        Fig.new(name);
        SP=d.plot('stdX(model==E)',Y,'bin','subj',...
            'xlTxt','Disparity (arcmin)',...
            'ylTxt','Negative log-likelihood',...
            'tlTxt',{'Observer 1','Observer 2','Observer 3','Average'},...
            'rlTxt',{'low C\delta','high C\delta'},...
            'slTxt',name,...
            'Color',color,...
            'xtInds',[1 3 5]);
        SP.sel(1);
        legend({'fit','fix=0'},'location','southeast');
    end
%- GEN
    function S=gen_sim2(obj)
        % INPUTS
        nTrlPerCmp=1000;

        stdX=[7.5]/60;
       % cmpX=[-3.75 -5.625 -7.5 -9.375 -11.25]/60;
        cmpX=[0.0859 0.1055 0.1250 0.1445 0.1641];
        sigma=[0.8 0.8]';
        rho=[.3];
        cr=[0]';

        S=DVCorr.genData(stdX,cmpX,[],sigma,rho,cr,nTrlPerCmp);

        Opts=DVFitter.getDefaults;
        Opts.rhoFitInd=[1];
        Opts.bCtrMu=true;
        Opts.bCombineCmp=true;

        %Opts.modelType='RRS';
        %Opts.stdFix=[];
        %   -0.2500   -0.1250         0    0.1250    0.2500
        %   -0.2772   -0.1386         0    0.1386    0.2772
        %   0.2500
        %   0.2255

        Opts.modelType='RR';
        Opts.stdFix=sigma;
        %  -0.0625   -0.0312         0    0.0312    0.0625
        %  -0.0625   -0.0312         0    0.0312    0.0625

        % -0.0781   -0.0391         0    0.0391    0.0781

        obj.DVFitter=DVFitter.new(S.stdX,S.cmpX,S.RCmpChs, Opts);
        obj.DVFitter.init();
        obj.DVFitter.run_test();

        dv=obj.DVFitter.DVFit;
        %dv.RHO
        dv.RHO
        dv.MU
        dv.STD
        obj.DVFitter.DVFit.negLLAll
        %bsxfun(@times,dv.STD,dv.MU)*60
        % 0.68 1.7: 978
        % 0.72 10   968
    end
    function gen_sim4(obj)
        Opts=DVFitter.getDefaults;

        % INPUTS
        nTrlPerCmp=100000;

        %stdX=0;
        %cmpX=0;

        %stdX=[-7.5]/60;
        %cmpX=[0.0859 0.1055 0.1250 0.1445 0.1641];
        stdX=0;
        cmpX=0;


        d=4;

        if d==4
            %% 4D
            %sigma=[0.8 0.8 0.8 0.8]';
            sigma=[1 1 1 1]';
            rho=[.15 .15 .15 .15 .15  .3]';
            cr=[0 0 0 0]';
            Opts.rhoFitInd=[1 2 2 2 2 3];
        elseif d==0
            sigma=[1 1 1 1]';
            %rho=[.15 .15 .15 .15 .15  .3]';
            rho=[0 0 0 0 0 0]';
            rho=[.3 .3 .3 .3 .3 .3]';
            cr=[0 0 0 0]';
            Opts.rhoFitInd=[1 2 2 2 2 3];
        elseif d==3
            % 3D
            sigma=[1.0 1.0 1.0]';
            %sigma=[0.8 0.8 0.8]';
            rho=[.3 .15 .15]';
            cr=[0 0 0]';
            Opts.rhoFitInd=[1 2 3];
        end

        S=DVCorr.genData(stdX,cmpX,[],sigma,rho,cr,nTrlPerCmp);

        Opts.bCtrMu=true;
        Opts.bCombineCmp=true;

        Opts.modelType='RR';
        Opts.stdFix=sigma;

        obj.DVFitter=DVFitter.new(S.stdX,S.cmpX,S.RCmpChs, Opts);
        obj.DVFitter.init();
        obj.DVFitter.run_test();

        dv=obj.DVFitter.DVFit;
        %dispV(dv.PP0.N,'NegLL')
        dispV(dv.RHO,'Rho_fit')
        dispV(dv.MU,'Mu_fit')
        dispV(dv.STD,'STD_fit')

        negLLAllE=dv.fun(rho);
        negLLAllO=obj.DVFitter.DVFit.negLLAll;
        dispV(negLLAllE,'NLL_Real');
        dispV(negLLAllO,'NLL_Fitted');

    end
    function plot_fun_gen(obj)
        dv=obj.DVFitter.DVFit;
        %rho=dv.rho;
        rho=[.3 .3 .3 .3 .3 .3]';
        [~,~,PP]=dv.fun(rho);
        P=PP.P;
        N=PP.N/sum(PP.N(:))
        d=PP.P-N;
        ind=find(abs(d) > .005)
        d(ind)
        plot(PP.N,PP.P,'.')

        [~,ind]=sort(sum(dv.F0.CRL0==0,2))
        dv.F0.CRL0=dv.F0.CRL0(ind,:)
        dv.F0.CRL0==0
    end
    function obj=gen_test0(obj)
        Opts.modelType='RRS';
        Opts.rhoFitInd=[1 2 2 2 2 3];
        Opts.bCtrMu=true;
        Opts.bCombineCmp=true;

        %inds=[1:4];
        %S=obj.table2DVInput(inds,'subj','DNW','stdD',0.0625);
        obj.DVFitter=DVFitter.new(S.stdX,S.cmpX,S.RCmpChs, Opts);
        obj.DVFitter.init();

        obj.DVFitter.run_test();

    end
    function gen_test2(obj)
        Opts=DVFitter.getDefaults;
        Opts.modelType='R';
        Opts.rhoFitInd=[1];
        Opts.bCtrMu=true;
        Opts.bCombineCmp=true;

        S=obj.table2DVInput(1:2,'subj','DNW','stdD',0.0625);
        obj.DVFitter=DVFitter.new(S.stdX,S.cmpX,S.RCmpChs, Opts);
        obj.DVFitter.init();

        obj.DVFitter.run_test();


        % 5 std * 2 std * 5 blk * 100 trl * 3 subj = 15000
    end
    function gen_test4(obj)
        Opts=DVFitter.getDefaults;
        Opts.modelType='RRS';
        Opts.rhoFitInd=[1 2 2 2 2 3];
        Opts.bCtrMu=true;
        Opts.bCombineCmp=true;

        inds=[1:2 1:2];
        %inds=[1:4];
        S=obj.table2DVInput(inds,'subj','DNW','stdD',0.0625);
        obj.DVFitter=DVFitter.new(S.stdX,S.cmpX,S.RCmpChs, Opts);
        obj.DVFitter.init();
        obj.DVFitter.run_test();


        % 5 std * 2 std * 5 blk * 100 trl * 3 subj = 15000
    end
    function gen_1(obj)
        Opts=DVFitter.getDefaults;
        Opts.modelType='RRS';
        Opts.rhoFitInd=[1 2 2 2 2 3];
        Opts.bCtrMu=true;
        Opts.bCombineCmp=true;

        Opts.stdFix=[T1 T12 T12 T];

        inds=[1:2 1:2];
        %inds=[1:4];
        S=obj.table2DVInput(inds,'subj','DNW');
        obj.DVFitter=DVFitter.new(S.stdX,S.cmpX,S.RCmpChs, Opts);
        obj.DVFitter.init();
        obj.DVFitter.run_test();


        % 5 std * 2 std * 5 blk * 100 trl * 3 subj = 15000
    end
%- SAVE
    %-PREDMODLE;
    %function fname=get_fname(obj)
    %    dire=[Env.var('DATA') 'fits' filesep];
    %    fname=[dire 'PredModel']
    %end
    %function fname=load(obj)
    %    fname=[obj.get_fname() '.mat'];
    %    load(fname);
    %    obj.RMS=RMS;
    %    obj.bRCmp=bRCmp;
    %    obj.bRCorrect=bRCorrect;
    %end
    %function save(obj)
    %    fname=obj.get_fname();
    %    RMS=obj.RMS;
    %    MN=obj.MN;
    %    bRCmp=obj.bRCmp;
    %    bRCorrect=obj.bRCorrect;
    %    save(fname,'MN','RMS','bRCorrect','bRCmp');
    %end
    function save(obj)
        dire=[Env.var('DATA') 'fits' filesep];

        %obj.save_fun('Thresh',dire);
        %obj.save_fun('Corr',  dire);
        %obj.save_fun('CorrS',  dire);
        obj.save_fun('CorrSF',  dire);
        %obj.save_fun('Ext',  dire);

    end
    function save_fun(obj,subname,dire)
        fname=[dire 'ExtModel_' subname];
        if isempty(obj.(subname)) || isempty(obj.(subname).Data)
            return
        end
        Table=obj.(subname).Data.TABLE;
        Key=  obj.(subname).Data.KEY;
        Types=obj.(subname).Data.types;
        if isfield(obj.(subname),'Opts')
            Opts= obj.(subname).Opts;
            save(fname,'Table','Key','Types','Opts');
        else
            save(fname,'Table','Key','Types');
        end
    end
%- LOAD
    function load(obj)
        dire=[Env.var('DATA') 'fits' filesep];

        obj.load_fun('Thresh',dire);
        obj.load_fun('Corr',  dire);
        %obj.load_fun('CorrS',  dire);
        obj.load_fun('CorrSF',  dire);
        obj.load_fun('Ext',   dire);
    end
    function load_fun(obj,subname,dire)
        fname=[dire 'ExtModel_' subname '.mat'];
        if Fil.exist(fname)
            D=load(fname);
            obj.(subname)
            obj.(subname).Data=Table(D.Table,D.Key,D.Types);
            if isfield(D,'Opts')
                obj.(subname).Opts=D.Opts;
            end
        else
            disp([subname ' file does not exist'])
        end
    end
    function test_fun(obj)
        %d.times('/varB/','/varL/',3600,{33},'!','stdB','-',3,{'!','stdL'},{'*','test'})
        %d.times('/varB/','/varL/',3600,33,'!','stdB','-',3,{'!','stdL'},{'*','test'});
        d2=d.times({'/varB/','/varL/'},3600)
        d3=d.times('/var[BL]',3600)
        d3
        %d2=d.times({'/var[BL]'},3600);
        %d2=d.times('/var[BL]',3600);

        %d2=d.times({'/var([BL])'},3600,{'!','//std$1'});
        dk

        d.times('/varB/','/varL/',3600,33,'!','stdB',{'!','stdL'},{'*','test'});
        d.times('varB','varL',3600)
        out=d.std('varB')
        out=d.mean('varB')
    end
end
end
