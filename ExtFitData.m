classdef ExtFitData < handle
properties
    dims
    passes
    model

    E
    M
    S
    A22
    A33
    A23

    %- DATA OPTS
    bZer=0
    bBino=0
    bAvg=0
    bKAll=0
    fitType=0
    bSqrt=0

    %- INDS
    bin
    subj

    stds
    subjs
    bins
    ks

    nStd
    nBin
    nSubj
    nK

    nInds
    bInds
    sInds
    kInds
    fname

    boot
    bootFits
    bootT
    bootR
    ciPercent=68.27
    CI
    sMean
    sStd
end
properties(Hidden)
    X22
    X23
    X33
    subjs_

    i
    b
    s
    k

    obj22
    obj33
    obj23
    A22Full
    A23Full
    A33Full
    nFits
    fitTypes={'A','AFitT','AFitRT'}
    EAll
    MAll
    SAll
    CIAll
end
methods(Static)
    function P=getP_fit()
        P={...
            'dims',[];
            'subjs_',[];
            'passes',[];
            'model',[];
        };
    end
    function obj=Get(dims,subjs,passes,model)
        bGen=0;
        if nargin ==0
            dims=[2.1 3];
            subjs={'DNW','AJL','JDB'};
            passes=[1 2];
            model='fix';
        end
        obj=ExtFitData(bGen,'dims',dims,'subjs_',subjs,'passes',passes,'model',model);
    end
    function obj=Gen(dims,subjs,passes,model)
        bGen=1;
        obj=ExtFitData(bGen,'dims',dims,'subjs_',subjs,'passes',passes,'model');
    end
    function obj=New(dims,subjs,passes,model)
        bGen=2;
        if nargin ==0
            dims=[2.1 3];
            SUBJS={'DNW','AJL','JDB'};
            passes=[1 2];
            model='fix';
        end
        obj=ExtFitData(bGen,'dims',dims,'subjs_',subjs,'passes',passes,'model');
    end
end
methods
    function obj=ExtFitData(varargin)
        if nargin < 1
            return
        end
        obj.parse(varargin{:});
    end

    function parse(obj,bGen,varargin)
        Args.simple(obj,obj.getP_fit(),varargin{:});
        if isempty(obj.subjs)
            obj.subjs=Figs.getsubjs(obj.dims(1));
        end
        obj.get(bGen);
    end
%- GET/SET
    function out=get.fname(obj)
        out=['ExtFitData.mat'];
    end
    function out=get.E(obj)
        out=obj.EAll.(obj.fitType);
    end
    function out=get.M(obj)
        out=obj.MAll.(obj.fitType);
    end
    function out=get.CI(obj)
        out=obj.CIAll.(obj.fitType);
    end
    function out=get.S(obj)
        out=obj.SAll.(obj.fitType);
    end
    function out=get.A22(obj)
        out=obj.A22Full.(obj.fitType);
    end
    function out=get.A23(obj)
        out=obj.A23Full.(obj.fitType);
    end
    function out=get.A33(obj)
        out=obj.A33Full.(obj.fitType);
    end
    function set.fitType(obj,in)
        if isnumeric(in)
            in=obj.fitTypes{in+1};
        end
        obj.fitType=in;
    end
    function out=get.nFits(obj)
        out=numel(obj.fitTypes);
    end
    %% - INDS
    % stds
    function out=get.stds(obj)
        out=obj.stds;
        if isempty(out)
            out=1:obj.nStd;
        end
    end
    function out=get.nStd(obj)
        out=size(obj.M.T2,1);
    end
    function out=get.nInds(obj)
        out=1:obj.nStd;
    end
    % bins
    function out=get.bin(obj)
        out=obj.bins(obj.b);
    end
    function out=get.bins(obj)
        out=obj.bins;
        if isempty(out)
            out=1:obj.nBin;
        end
    end
    function out=get.nBin(obj)
        out=size(obj.M.T2,2);
    end
    function out=get.bInds(obj)
        out=1:obj.nBin;
    end
    % subj
    function out=get.subj(obj)
        out=obj.subjs{obj.s};
    end
    function out=get.subjs(obj);
        out=obj.subjs_;
        if obj.bAvg
            out=[out 'ALL'];
        end
    end
    function out=get.nSubj(obj)
        out=numel(obj.subjs);
    end
    function out=get.sInds(obj)
        out=1:obj.nSubj;
    end
    % k
    function out=get.ks(obj)
        ks=num2cell(Set.distribute(obj.passes,obj.passes));
        ks=[{obj.passes} {obj.passes}; k]
        out=ks(obj.kInds,:);
    end
    function out=get.nK(obj)
        out=size(obj.M.covLB,4);
    end
    function out=get.kInds(obj)
        if obj.bKAll
            out=1;
        else
            out=2:obj.nK;
        end
    end
    function cl=tmpSet(obj,fld,val);
        oval=obj.(fld);
        obj.(fld)=val;
        cl=onCleanup(@() set_fun(obj,fld,oval));
        function set_fun(obj,fld,oval)
            obj.(fld)=oval;
        end
    end
%- GET VARS
    function [M,S]=get_val(obj,bAvgF,bKs,fld)
        if nargin < 2 || isempty(bAvgF)
            bAvgF=false;
        end
        if nargin < 3 || isempty(bKs)
            bKs=false;
        end
        if nargin < 4 || isempty(fld)
            fld=obj.fld;
        end
        b  =obj.b;
        s=obj.s;
        bKFld=any(strcmp(fld,{'covLB','corLB','varB','stdB'})) || contains(fld,'23');
        %bK=obj.bKAll && bKFld;
        if obj.bKAll || bKs
            K=2:obj.nK;
        else
            K=obj.k;
        end
        if isempty(K)
            K=1;
        end

        suf='';
        if   any(strcmp(fld,{'covLB','corLB'})) && obj.bZer
            fld='covLB0';
        elseif any(strcmp(fld,{'varB','stdB'})) && obj.bZer
            suf='B0';
        elseif any(strcmp(fld,{'varL','stdL'})) && ~obj.bBino
            suf='LD';
        elseif any(strcmp(fld,{'varD','stdD'})) && ~obj.bBino
            fld='covLB0';
        end
        if ~isempty(suf)
            fld=[fld(1:3) suf];
        end

        M=obj.M.(fld)(:,b,s,K);
        S=obj.S.(fld)(:,:,b,s,K);
        if bKs
            return
        elseif obj.bKAll
            M=median(M,4);
            S=mean(S,5);
        elseif bAvgF && size(M,4) > 1
            M=mean(M,4);
            S=mean(S,5);
        end
    end



%- LOAD/SAVE/INIT
    function get(obj,bGen)
        if nargin < 2 || isempty(bGen)
            bGen=false;
        end

        if bGen==2
            bGen=true;
            bNew=true;
        else
            bNew=false;
        end
        if bGen || ~Fil.exist(obj.fname);
            obj.gen(bNew);
            obj.unpackBoot();
            obj.fitBoot();
            obj.save();
        else
            obj.load();
        end
        mult=60;% XXX

        obj.bootFits.T=obj.bootFits.T;

        obj.unpackBoot();

        % none
        obj.CIAll.A = obj.parseBoot(obj.boot,mult);

        % T
        T=struct();
        T.RHO=obj.boot.RHO;
        T.T  =obj.bootFits.T;
        obj.CIAll.AFitT = obj.parseBoot(T,mult);

        % RT
        obj.CIAll.AFitRT = obj.parseBoot(obj.bootFits,mult);

        MAll=struct();
        SAll=struct();
        for i = 1:obj.nFits
            obj.fitType=obj.fitTypes{i};
            [MAll.(obj.fitType), SAll.(obj.fitType)]=obj.parseA(obj.A22,obj.A23,obj.A33);
        end
        obj.MAll=MAll;
        obj.SAll=SAll;
    end
    function gen(obj,bNew)
        if nargin < 2 || isempty(bNew)
            bNew=false;
        end
        obj.gen_from_EPs(bNew);
    end
    function gen_from_EPs(obj,bNew)
        if isnumeric(obj.dims) && numel(obj.dims)==2
            dims1=obj.dims(1);
            dims2=obj.dims(2);
        end
        if isnumeric(obj.passes) && numel(obj.passes)==2
            passes1=obj.passes(1);
            passes2=obj.passes(2);
        end

        [~,obj22]  =EPs.getData('ext',bNew,obj.model, obj.subjs_,dims1,obj.passes);
        [~,obj33]  =EPs.getData('ext',bNew,obj.model, obj.subjs_,dims2,obj.passes);
        [~,obj23]  =EPs.getData('ext',bNew,obj.model, obj.subjs_,obj.dims ,obj.passes);
        A22=struct();
        A33=struct();
        A23=struct();

        for i = 1:obj.nFits
            fld=obj.fitTypes{i};
            A22.(fld)=obj22.(fld);
            A33.(fld)=obj33.(fld);
            A23.(fld)=obj23.(fld);

            obj.X22=obj22.EP{1}.X;
            obj.X23=obj23.EP{1}.X;
            obj.X33=obj33.EP{1}.X;
        end

        obj.A22Full=A22;
        obj.A23Full=A23;
        obj.A33Full=A33;
        obj.obj22=obj22;
        obj.obj33=obj33;
        obj.obj23=obj23;
    end
    function save(obj)
        A22=obj.A22Full;
        A33=obj.A33Full;
        A23=obj.A23Full;
        obj22=obj.obj22;
        obj33=obj.obj33;
        obj23=obj.obj23;
        bootFits=obj.bootFits;
        save(obj.fname,'A22','A33','A23','obj22','obj23','obj33','bootFits');
    end
    function load(obj)
        load(obj.fname);
        obj.A22Full=A22;
        obj.A23Full=A23;
        obj.A33Full=A33;
        obj.obj22=obj22;
        obj.obj33=obj33;
        obj.obj23=obj23;
        obj.bootFits=bootFits;
    end
    %function unpackFits(obj)
    %    FLDS={'obj23','obj22','obj33'};
    %    flds={'RHO','T'};
    %    fldsb={'RFitBoot','TFitBoot'};
    %    %flds={'T','RHO'};
    %    %fldsb={'TFitBoot','RFitBoot'};

    %    obj1=obj.(FLDS{1}).EP{1};
    %    nBoot=obj1.Opts.dv.opts.nBoot;
    %    nStd=size(obj1.DVCORR,1);
    %    nBin=size(obj1.DVCORR,2);
    %    nK  =size(obj1.DVCORR,4);
    %    nS=length(obj.(FLDS{1}).EP);
    %    nF=length(FLDS);
    %    sz=[nBoot nStd nBin nK nS nF];

    %    boot=struct();
    %    nf=length(flds);
    %    for ff=1:nf
    %        boot.(flds{ff})=zeros(sz);
    %    end
    %    for f  = 1:nF
    %    for ff = 1:nf
    %        FLD=FLDS{f};
    %        fld=flds{ff};
    %        if strcmp(fld,'T') && strcmp(FLD,('obj23'))
    %            continue
    %        end
    %        for s  = 1:nS
    %            val=obj.(FLD).EP{s}.(fldsb{ff});

    %            if size(val,5)==1
    %                val=repmat(val,[1 1 1 1 nK]);
    %            end
    %            boot.(fld)(:,:,:,:,s,f)=permute(val(:,:,:,1,:),[1 2 3 5 4]);
    %        end
    %    end
    %    end
    %    obj.bootFits=boot;
    %end
    %function get_fits(obj)
    %    obj.get_TFit();
    %    obj.get_RFit();
    %end
    %function get_TFit(obj);
    %    flds={'obj22','obj33'};
    %    for f = 1:length(flds)
    %        OB=obj.(flds{f});
    %        for i=1:length(OB.EP)
    %            OB.EP{i}.get_TFit;
    %        end
    %    end
    %end
    %function get_RFit(obj)
    %    flds={'obj22','obj33','obj23'};
    %    for f = 1:length(flds)
    %        OB=obj.(flds{f});
    %        for i=1:numel(OB.EP)
    %            OB.EP{i}.get_RFit;
    %        end
    %    end
    %end
%- A
    function unpackBoot(obj)
        FLDS={'obj23','obj22','obj33'};
        flds={'RHO','T','varT','varE','varI'};
        boot=struct();

        KK=0;
        for f=1:length(FLDS)
            FLD=FLDS{f};
            for s=1:length(obj.(FLD).EP)
                [N,M,B,K]=size(obj.(FLD).EP{s}.DVCORR);
                KK=max(KK,K);
                for k = 1:K
                for i = 1:N
                for j = 1:M
                    B=obj.(FLD).EP{s}.DVCORR{i,j}.DVFITTER.boot;
                    for ff =1:length(flds)
                        fld=flds{ff};
                        boot.(fld)(:,i,j,k,s,f)=B.(fld)(1,:);
                    end
                end
                end
                end
                if K==1
                    for ff=1:length(flds)
                        fld=flds{ff};
                        boot.(fld)(:,:,:,:,s,f)=repmat(boot.(fld)(:,:,:,1,s,f),[1,1,1,KK,1,1]);

                    end
                end
            end
        end
        obj.boot=boot;
    end
    function fitBoot(obj)
        obj.fitBootT();
    end
    function fitBootT(obj)
        [nBoot,nStd,nBin,nK,nS,nF]=size(obj.boot.T);
        sz=[nBoot nStd nBin nK nS nF];

        flds={'T','RHO'};
        boot=struct();
        nf=length(flds);
        for ff=1:nf
            boot.(flds{ff})=zeros(sz);
        end

        f0=[0; 0];
        T=permute(obj.boot.T,[2 1 3 4 5 6]);
        R=permute(obj.boot.RHO,[2 1 3 4 5 6]);
        for m = 1:nBin
        for k = 1:nK
        for s = 1:nS
        for f = 1:nF
            YT=log(T(:,:,m,k,s,f));
            YR=    R(:,:,m,k,s,f);
            X=obj.X22(:,m);
            ffun=@(mb) (mb(1)*X+mb(2));
            for i = 1:nBoot
                yt=YT(:,i);
                yr=YR(:,i);

                fun=@(mb) mean(abs(ffun(mb) - yt).^2);
                [~,fout]=evalc('fmincon(fun,f0)');
                boot.T(i,:,m,k,s,f)=exp(ffun(fout));

                fun=@(mb) mean(abs(ffun(mb) - yr).^2);
                [~,fout]=evalc('fmincon(fun,f0)');
                boot.RHO(i,:,m,k,s,f)=ffun(fout);
            end
        end
        end
        end
        end
        obj.bootFits=boot;
        obj.bootFits.T

    end
    function out=sqrta(obj,in)
        out=sqrt(abs(in)).*sign(in);
    end
    function C=parseBoot(obj,B,mult)
        if nargin < 3  || isempty(mult)
            mult=1;
        end
        M=obj.ciPercent;
        dp=1;% XXX

        T2    =   B.T(:,:,:,:,:,2)*mult;
        T3    =   B.T(:,:,:,:,:,3)*mult;

        rho23 = B.RHO(:,:,:,:,:,1);
        rho22 = B.RHO(:,:,:,:,:,2);
        rho33 = B.RHO(:,:,:,:,:,3);

        rat22 = obj.sqrta(rho22./(1-rho22));
        rat33 = obj.sqrta(rho33./(1-rho33));

        varT2 = T2.^2.*dp;
        varT3 = T3.^2.*dp;
        varT23=obj.sqrta(varT2.*varT3);

        varE2=varT2.*rho22;
        varE3=varT3.*rho33;
        varE23=varT23.*rho23;

        varI2=varT2-varE2;
        varI3=varT3-varE3;

        % std
        stdI2=obj.sqrta(varI2);
        stdI3=obj.sqrta(varI3);

        stdE2=obj.sqrta(varE2);
        stdE3=obj.sqrta(varE3);
        stdE23=obj.sqrta(varE23);

        stdT2=obj.sqrta(varT2);
        stdT3=obj.sqrta(varT3);
        stdT23=obj.sqrta(varT23);


        % comp
        varLD =varE2;
        covLB =varE23-varLD;
        varB  =varE3-varLD-2*covLB;
        varB0 =varE3-varLD;
        varD= zeros(size(varLD)); %XXX
        varL= zeros(size(varLD)); %XXX

        corr  =covLB./sqrt(abs(varB.*varLD)); % XXX


        % std
        stdLD=obj.sqrta(varLD);
        stdB0=obj.sqrta(varB0);
        stdB =obj.sqrta(varB);
        stdD= zeros(size(varLD)); %XXX
        stdL= zeros(size(varLD)); %XXX


        %[nBoot,n,m,nK,nSubj]=size(rho23);

        C=struct();
        C.covLB=confInt(covLB, M,1);
        C.covLB0=zeros(size(C.covLB));
        C.corLB=confInt(corr,  M,1);

        C.varLD=confInt(varLD, M,1);
        C.varB0=confInt(varB0, M,1);
        C.varB =confInt(varB,  M,1);
        C.varL =confInt(varL, M,1);
        C.varD =confInt(varD, M,1);

        C.stdLD=confInt(stdLD, M,1);
        C.stdB0=confInt(stdB0, M,1);
        C.stdB =confInt(stdB,  M,1);
        C.stdD =confInt(stdD,  M,1);
        C.stdL =confInt(stdL,  M,1);

        C.varT2=confInt(varT2, M,1);
        C.varT3=confInt(varT3, M,1);
        C.varT23=confInt(varT23, M,1);
        C.varE2=confInt(varE2, M,1);
        C.varE3=confInt(varE3, M,1);
        C.varE23=confInt(varE23, M,1);
        C.varI2=confInt(varI2, M,1);
        C.varI3=confInt(varI3, M,1);

        C.T2   =confInt(T2,    M,1);
        C.T3   =confInt(T3,    M,1);
        C.rho22=confInt(rho22, M,1);
        C.rho23=confInt(rho23, M,1);
        C.rho33=confInt(rho33, M,1);

        C.rat22=confInt(rat22, M,1);
        C.rat33=confInt(rat33, M,1);

        C.stdI2=confInt(stdI2, M,1);
        C.stdI3=confInt(stdI3, M,1);

        C.stdE2=confInt(stdE2, M,1);
        C.stdE3=confInt(stdE3, M,1);
        C.stdE23=confInt(stdE23, M,1);

        C.stdT2=confInt(stdT2, M,1);
        C.stdT3=confInt(stdT3, M,1);
    end
    function [M,S]=parseA(obj,A22,A23,A33)
        E=struct();

        nK=size(A23.Rho,5);


        % int
        E.varI2=A22.vi;
        E.varI3=A33.vi;

        % ext
        E.varE2 =A22.ve;
        E.varE3 =A33.ve;

        E.varE23=A23.Rho .* sqrt(A22.vt .* A33.vt);
        %E.varE23=A23.ve;
        %E.varE23=A232.ve;
        %E.varE23=mean(cat(4,A231.ve,A232.ve),4);

        % rho
        E.rho22=A22.Rho;
        E.rho33=A33.Rho;
        E.rho23=A23.Rho;

        % ratio
        E.rat22 = obj.sqrta(E.rho22./(1-E.rho22));
        E.rat33 = obj.sqrta(E.rho33./(1-E.rho33));

        % T
        dpCrit=1; % XXX

        E.T2=obj.sqrta(A22.vt ./ dpCrit);
        E.T3=obj.sqrta(A33.vt ./ dpCrit);

        % tot
        E.varT2 =A22.vt;
        E.varT3 =A33.vt;
        E.varT23=A23.vt;

        % varLD
        E.varLD=E.varE2;
        vLD=repmat(E.varLD,[1 1 1 1 nK]);

        % cov
        E.covLB =E.varE23 - vLD;
        E.covLB0=zeros(size(E.covLB));

        % varB
        E.varB0 =E.varE3  - vLD;
        E.varB  =E.varE3  - vLD  - 2*E.covLB;

        %ind=E.varB < 0;
        %CC=E.covLB(ind) + -E.varB(ind)/2;
        %E.covLB(ind)=CC;
        %E.varB(ind)=0;

        % corr
        E.corLB =E.covLB./sqrt(abs(E.varB.*E.varLD));


        % varD
        %: y-intercept
        %: [ subj std bin ]
        E.varL=zeros(size(E.varLD));
        X=A22.X(:,1)';
        f0=[0; 0];
        for s = 1:size(E.varLD,1)
        for b = 1:size(E.varLD,3)
            Y=log(E.varLD(s,:,b));
            fun=@(mb) mean(abs((mb(1)*X+mb(2)) - Y).^2);
            [~,fout]=evalc('fmincon(fun,f0)');
            E.varL(s,:,b)=exp(fout(2));
        end
        end
        E.varD=E.varLD-E.varL;

        % std
        E.stdLD=obj.sqrta(E.varLD);
        E.stdD =obj.sqrta(E.varD);
        E.stdB =obj.sqrta(E.varB);
        E.stdB0=obj.sqrta(E.varB0);

        E.stdI2=obj.sqrta(E.varI2);
        E.stdI3=obj.sqrta(E.varI3);

        E.stdE2 =obj.sqrta(E.varE2);
        E.stdE3 =obj.sqrta(E.varE3);
        E.stdE23=obj.sqrta(E.varE23);

        E.stdT2=obj.sqrta(E.varT2);
        E.stdT3=obj.sqrta(E.varT3);


        % PERMUTE, MEAN and STDs
        %: og indeces are [subj, lvl, ind ]
        %:             -> [lvl,  ind, subj]

        flds=fieldnames(E);
        flds2={};
        M=struct();
        S=struct();
        skipflds={'X'};


        p=[2 3 4 1];
        %p=[1 2 3];
        r=[1 1 1 nK];
        for i = 1:length(flds)
            fld=flds{i};
            val=E.(fld);
            if any(strcmp(fld,skipflds))
                continue
            end
            if numel(size(val))== 5
                flds2=[flds2 fld];
                continue
            end
            val=repmat(val,r);
            % s lvl bin K
            cval=obj.CI.(fld);
            % CI lvl bin K S

            m=mean(val,1);
            s=confInt(val,obj.ciPercent,1,'approximate');
            %s=std(val,0,3);

            S.(fld)=cat(5,cval,s);
            M.(fld)=cat(1, val,m);
            M.(fld)=permute(M.(fld),p);
        end

        p23=[1 2 3 5 4];
        for i = 1:length(flds2)
            fld=flds2{i};
            val=E.(fld);
            val=permute(val,p23);
            cval=obj.CI.(fld);

            m=mean(val,1);
            s=confInt(val,obj.ciPercent,1,'approximate');
            %s=std(val,0,3);

            S.(fld)=cat(5,cval,s);
            M.(fld)=cat(1, val,m);
            M.(fld)=permute(M.(fld),p);
        end

        % lvl, bin, k, s ->
        % lvl, bin, s, k
        ps=[1 2 3 5 4];
        pa=[1 2 4 3];
        flds=fieldnames(S);
        for i = 1:length(flds)
            fld=flds{i};
            S.(fld)=permute(S.(fld),ps);
            M.(fld)=permute(M.(fld),pa);
        end
        M.X=A22.X(:,1);

    end
    function scrap(obj,A22,A33);
        % CI
        p=[2 3 1];
        varI2_SE=permute(A22.si_se, p);
        varI3_SE=permute(A33.si_se, p);

        rho22_SE=permute(A22.Rho_se, p);
        rho33_SE=permute(A33.Rho_se, p);

        varLD_SE =permute(A22.se_se, p);
        varE3_SE =permute(A33.se_se, p);
        varE23_SE=permute(A23.se_se, p23);

        varT2_SE =permute(A22.st_se, p);
        varT3_SE =permute(A33.st_se, p);


        pn=[1 1 1 5];
        vLD_SE=repmat(varLD_SE,pn);
        vE3_SE =repmat(varE3_SE,pn);
        covLB=  E.covLB(:,:,1:(end-1),:);
        varE23=E.varE23(:,:,1:(end-1),:);
        varE3=  E.varE3(:,:,1:(end-1),:);
        varLD=  E.varLD(:,:,1:(end-1),:);
        varB=    E.varB(:,:,1:(end-1),:);
        varB0=  E.varB0(:,:,1:(end-1),:);
        %covLB_SE= -covLB - ((varE23 + varE23_SE) - (varLD + vLD_SE));
        %varB_SE = -varB  - ((varE3  + vE3_SE)    - (varLD + vLD_SE) - 2*(covLB + covLB_SE));
        %varB0_SE= -varB0 - ((varE3  + vE3_SE)    - (varLD + vLD_SE));

        covLB_SE=sqrt(varE23_SE.^2 + vLD_SE.^2);
        varB_SE =sqrt(vE3_SE.^2  + vLD_SE.^2 + 2*covLB_SE.^2);
        varB0_SE=sqrt(vE3_SE.^2  + vLD_SE.^2);

        E.varLD_s=cat(3,varLD_SE, E.varLD_s);
        E.varB_s =cat(3,varB_SE,  E.varB_s);
        E.varB0_s=cat(3,varB0_SE, E.varB0_s);
        E.covLB_s=cat(3,covLB_SE, E.covLB_s);

        E.varI2_s=cat(3,varI2_SE, E.varI2_s);
        E.varI3_s=cat(3,varI3_SE, E.varI3_s);
        E.varE3_s=cat(3,varE3_SE, E.varE3_s);

        E.varT2_s=cat(3,varT2_SE, E.varT2_s);
        E.varT3_s=cat(3,varT3_SE, E.varT3_s);


        E.rho22_s=cat(3,rho22_SE, E.rho22_s);
        E.rho33_s=cat(3,rho33_SE, E.rho33_s);

        E.covLB0_s=zeros(size(E.covLB_s));

        % varB
        % covLB
        % varE23
        % varT23
        % rho23

        % XXX
        % confidence intervals
        % varL = varE
        %
        % varB0 = varE3  - varLD
        % varB  = varE3  - varLD - 2 * E.covLB
        % covLB = varE23 - varLD
    end
end
end
