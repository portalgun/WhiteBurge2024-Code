classdef Ext < handle
properties
    model
    m
    M
    varL
    varB
    covLB
    rhoLB
    nanInds

    T2
    T3

    varI2
    varI3

    varE2
    varE3
    varE23

    varT2
    varT3
    varT23

    rho22
    rho23
    rho33

    % predcitions
    rho22p
    covRho22p
    corrRho22p

    rho23p
    covRho23p
    corrRho23p

    rho33p
    covRho33p
    corrRho33p

    dpcrit=1
end
properties(Transient,Hidden)
    RCP
    rP
    cP
end
methods(Static)
    function obj=TTCRR(T2,T3,covLB,rho22,rho33)
        obj=Ext();
        obj.TTCRR_(T2,T3,covLB,rho22,rho33);
    end
    function obj=TTmCR(T2,T3,m,covLB,rho22)
        obj=Ext();
        obj.TTmCR_(T2,T3,m,covLB,rho22);
    end
    function obj=LmCR(varL,m,covLB,rho22)
        obj=Ext();
        obj.LmCR_(varL,m,covLB,rho22);
    end
    function obj=LmCT(varL,m,covLB,T2)
        obj=Ext();
        obj.LmCT_(varL,m,covLB,T2);
    end
    function obj=LMCR(varL,m,covLB,rho22)
        obj=Ext();
        obj.LMCR_(varL,m,covLB,rho22);
    end
end
methods
    function [E,covP,rho33P,rho330]=rho23_rho33_error(obj,rho23M,covP,bCorr)
        rho330=obj.cov_rho33(covP,bCorr);
        [covP,rho33P]=obj.rho23_rho33(rho23M,bCorr,true);
        E=rho33P-rho330;
    end
    function [E,covM,rho23M,rho23P]=rho33_rho23_error(obj,rho33M,covP,bCorr)
        rho23P=obj.cov_rho23(covP,bCorr);
        [covM,rho23M]=obj.rho33_rho23(rho33M,bCorr,true);
        E=rho23P-rho23M;
    end
%- covLB to rho
    function [R]=cov_rho23(obj,covq,bCorr)
        if bCorr
            X=obj.rhoLB;
        else
            X=obj.covLB;
        end
        if numel(obj.rho33)==1
            R=obj.rho23;
        else
            R=interp1(X,obj.rho23,covq);
        end
    end
    function [R]=cov_rho33(obj,covq,bCorr)
        if bCorr
            X=obj.rhoLB;
        else
            X=obj.covLB;
        end
        if numel(obj.rho33)==1
            R=obj.rho33;
        else
            R=interp1(X,obj.rho33,covq);
        end
    end
%- rho to x (to rho)
    function [X,R]=rho22_rho33(obj,Yq,bCorr,bBack)
        % Yq is rho22, R is rho33
        if nargin < 4 || isempty(bBack)
            bBack=nargout >= 2;
        end

        x=obj.covInterp_(Yq,'rho22',bCorr);;
        if bCorr
            obj.corrRho23p=x;
        else
            obj.covRho22p=x;
        end

        if bBack
            if ~isreal(x)
                obj.rho33p=nan;
            else
                obj.rho33p=obj.rho_interp_(x,'rho33',bCorr);
            end
        end

        if nargout > 1
            X=x;
            R=obj.rho22p;
        end
    end
    function [X,R]=rho23_rho33(obj,Yq,bCorr,bBack)
        % Yq is rho23, R is rho33
        if nargin < 4 || isempty(bBack)
            bBack=nargout >= 2;
        end

        x=obj.covInterp_(Yq,'rho23',bCorr);
        if bCorr
            obj.corrRho23p=x;
        else
            obj.covRho23p=x;
        end

        if bBack
            if ~isreal(x)
                obj.rho33p=nan;
            else
                obj.rho33p=obj.rho_interp_(x,'rho33',bCorr);
            end
        end

        if nargout > 0
            X=x;
            if nargout > 1
                R=obj.rho33p;
            end
        end
    end
    function [X,R]=rho33_rho23(obj,Yq,bCorr,bBack,YqObs)
        if nargin < 4 || isempty(bBack)
            bBack=nargout >= 2;
        end

        x=obj.covInterp_(Yq,'rho33',bCorr);
        if bCorr
            obj.corrRho33p=x;
        else
            obj.covRho33p=x;
        end

        if bBack
            if ~isreal(x)
                obj.rho23p=nan;
            else
                obj.rho23p=obj.rho_interp_(x,'rho23',bCorr);
            end
        end
        if ~isreal(obj.rho23p)
            obj.rho23p=nan;
        end

        if nargout > 0
            X=x;
            if nargout > 1
                R=obj.rho23p;
            end
        end
    end
    %- private
    function [out]=covInterp_(obj,Rq,fld,bCorr)
        if bCorr
            X=obj.rhoLB;
        else
            X=obj.covLB;
        end
        R=obj.(fld);

        %if any(R==0)
        %    ind=R==0;
        %    R(ind)=[];
        %    X(ind)=[];
        %end
        if numel(X)==1 || numel(R)==1
            out=X;
        else
            out=interp1(R,X,Rq,'linear','extrap');
        end
    end
    function out=rho_interp_(obj,Covq,fld,bCorr)
        if bCorr
            X=obj.rhoLB;
        else
            X=obj.covLB;
        end
        R=obj.(fld);
        if numel(R) == 1
            out=R;
        else
            out=interp1(X,R,Covq,'linear','extrap');
        end
    end
%- MODELS
    function LMCR_(obj,varL,M,covLB,rho22)
        % varT3 = M varT2
        % M = varT3/varT2

        % Tot
        varT2=varL./rho22;
        varT3=M.*varT2;
        varT23=sqrt(varT2.*varT3);

        % Int
        varI2=varT2-varL;
        varI3=varI2;

        % Ext
        varE2=varL;
        varE23=varL+covLB;
        varE3=varT3-varI3;

        % B
        varB=varE3-(varL+2.*covLB);

        % rho
        rho22=rho22;
        rho33=varE3./varT3;
        rho23=varE23./varT23;

        % Thresh
        T2=sqrt(varT2./obj.dpcrit);
        T3=sqrt(varT3./obj.dpcrit);

        % m
        m=varB./varL;

        %-
        obj.T2=T2;
        obj.T3=T3;

        %-
        obj.varL=varL;
        obj.varB=varB;
        obj.covLB=covLB;
        obj.rhoLB=obj.covLB./sqrt(obj.varB.*obj.varL);

        obj.varI2=varI2;
        obj.varI3=varI3;

        obj.varE2 =varE2;
        obj.varE3 =varE3;
        obj.varE23=varE23;

        obj.varT2 =varT2;
        obj.varT3 =varT3;
        obj.varT23=varT23;

        obj.rho22=rho22;
        obj.rho23=rho23;
        obj.rho33=rho33;
        obj.m=m;
        obj.M=M;
        obj.model='LMCR';

        obj.handle_neg_();
    end
    function BTTRR_(obj,varB,T2,T3,rho22,rho33)
        % varB = m VarL

        % tot
        varT2=T2.^2 .* obj.dpcrit;
        varT3=T3.^2 .* obj.dpcrit;
        varT23=sqrt(varT2.*varT3);


        varL=rho22.*varT2;
        varB=varB;

        % E
        varE2=rho22.*varT2;
        varE3=rho33.*varT3;
        %% varE23=varL+covLB;

        % Int
        varI2=varT2-varE2;
        varI3=varI2;

        % rho
        rho22=rho22;
        rho33=rho33;
        %% rho23=varE23./varT23;

        % vars
        varL=varL;
        %% varB=varE3-varL-2*covLB;

        % M
        M=varT3./varT2;
        m=varB./varL;

        %-
        obj.T2=T2;
        obj.T3=T3;

        %-
        obj.varL=varL;
        obj.varB=varB;
        obj.covLB=covLB;

        obj.varI2=varI2;
        obj.varI3=varI3;

        obj.varE2 =varE2;
        obj.varE3 =varE3;
        obj.varE23=varE23;

        obj.varT2 =varT2;
        obj.varT3 =varT3;
        obj.varT23=varT23;

        obj.rho22=rho22;
        obj.rho23=rho23;
        obj.rho33=rho33;
        obj.m=m;
        obj.M=M;
        obj.model='TTCRR';

        obj.handle_neg_();
    end
    function TTCRR_(obj,T2,T3,covLB,rho22,rho33)
        % varB = m VarL

        % tot
        varT2=T2.^2 .* obj.dpcrit;
        varT3=T3.^2 .* obj.dpcrit;
        varT23=sqrt(varT2.*varT3);


        varL=rho22.*varT2;

        % E
        varE2=rho22.*varT2;
        varE3=rho33.*varT3;
        varE23=varL+covLB;

        % Int
        varI2=varT2-varE2;
        varI3=varI2;

        % rho
        rho22=rho22;
        rho33=rho33;
        rho23=varE23./varT23;

        % vars
        varL=varL;
        varB=varE3-varL-2*covLB;

        % M
        M=varT3./varT2;
        m=varB./varL;

        %-
        obj.T2=T2;
        obj.T3=T3;

        %-
        obj.varL=varL;
        obj.varB=varB;
        obj.covLB=covLB;

        obj.varI2=varI2;
        obj.varI3=varI3;

        obj.varE2 =varE2;
        obj.varE3 =varE3;
        obj.varE23=varE23;

        obj.varT2 =varT2;
        obj.varT3 =varT3;
        obj.varT23=varT23;

        obj.rho22=rho22;
        obj.rho23=rho23;
        obj.rho33=rho33;
        obj.m=m;
        obj.M=M;
        obj.model='TTCRR';

        obj.handle_neg_();
    end
    function TTmCR_(obj,T2,T3,m,covLB,rho22)
        % varB = m VarL

        % totT
        varT2=T2.^2*obj.dpcrit;
        varT3=T3.^2*obj.dpcrit;

        varB=varL*m;
        varL=rho22*varT2;

        % Int
        varI2=varT2-varL;
        varI3=varI2;

        % ext
        varE2=rho22*varT2;
        varE23=varL+covLB;
        varE3=varB+varL+2.*covLB;

        % Ext
        varE2=varL;
        varE23=varL+covLB;
        varE3=varB+varL+2.*covLB;


        % Tot
        varT2=varE2+varI2;
        varT3=varE3+varI3;
        varT23=sqrt(varT2.*varT3);

        % rho
        rho22=rho22;
        rho33=varE3./varT3;
        rho23=varE23./varT23;

        % M
        M=varT3./varT2;

        %-
        obj.T2=T2;
        obj.T3=T3;

        %-
        obj.varL=varL;
        obj.varB=varB;
        obj.covLB=covLB;

        obj.varI2=varI2;
        obj.varI3=varI3;

        obj.varE2 =varE2;
        obj.varE3 =varE3;
        obj.varE23=varE23;

        obj.varT2 =varT2;
        obj.varT3 =varT3;
        obj.varT23=varT23;

        obj.rho22=rho22;
        obj.rho23=rho23;
        obj.rho33=rho33;
        obj.m=m;
        obj.M=M;
        obj.model='TTmCR';

        obj.handle_neg_();
    end
    function LmCR_(obj,varL,m,covLB,rho22)
        % varB = m VarL

        % same noise between condition0
        varB=varL*m;

        % Int
        %varI2=(varL-rho22.*varL)./rho22;
        varI2=(1./rho22-1).*varL;
        varI3=varI2;

        % Ext
        varE2=varL;
        varE23=varL+covLB;
        varE3=varB+varL+2.*covLB;


        % Tot
        varT2=varE2+varI2;
        varT3=varE3+varI3;
        varT23=sqrt(varT2.*varT3);

        % rho
        rho22=rho22;
        rho33=varE3./varT3;
        rho23=varE23./varT23;

        % Thresh
        T2=sqrt(varT2./obj.dpcrit);
        T3=sqrt(varT3./obj.dpcrit);

        % M
        M=varT3./varT2;

        %-
        obj.T2=T2;
        obj.T3=T3;

        %-
        obj.varL=varL;
        obj.varB=varB;
        obj.covLB=covLB;

        obj.varI2=varI2;
        obj.varI3=varI3;

        obj.varE2 =varE2;
        obj.varE3 =varE3;
        obj.varE23=varE23;

        obj.varT2 =varT2;
        obj.varT3 =varT3;
        obj.varT23=varT23;

        obj.rho22=rho22;
        obj.rho23=rho23;
        obj.rho33=rho33;
        obj.m=m;
        obj.M=M;
        obj.model='LmCR';

        obj.handle_neg_();
    end
    function LmCT_(obj,varL,m,covLB,T2)
        varT2=T2^2;
        varB=varL*m;

        % E
        varE2=varL;
        varE23=varL+covLB;
        varE3=varB+varL+2.*covLB;

        % I
        varI2=varT2-varE2;
        varI3=varI2;

        % T
        varT2=varT2;
        varT3=varE3+varI3;
        varT23=sqrt(varT2.*varT3);


        % Thresh
        T2=sqrt(varT2./obj.dpcrit);
        T3=sqrt(varT3./obj.dpcrit);

        % Rho
        rho22=varE2./varT2;
        rho23=varE23./varT23;
        rho33=varE3./varT3;

        % M
        m=m;
        M=T3./T2;

        %-
        obj.T2=T2;
        obj.T3=T3;

        %-
        obj.varL=varL;
        obj.varB=varB;
        obj.covLB=covLB;

        obj.varI2=varI2;
        obj.varI3=varI2;

        obj.varE2=varE2;
        obj.varE3=varE3;

        obj.varT2 =varT2;
        obj.varT3 =varT3;
        obj.varT23=varT23;

        obj.rho22=rho22;
        obj.rho23=rho23;
        obj.rho33=rho33;
        obj.model='LmCT';
        obj.m=m;
        obj.M=M;

        obj.handle_neg_();
    end
    function handle_neg_(obj)
        vflds={'varL', 'varB', ...
               'varI2','varI3',...
               'varE2','varE3','varE23'...
               'varT2','varT3','varT23' ...
               'm',...
        };
        rflds={'rho22','rho23','rho33','rhoLB'};
        tflds={'T2','T3','M'};
        oflds={'covLB'};
        inds=[];
        flds=[vflds,rflds,tflds,oflds];
        for i = 1:length(flds)

            % handle size
            if numel(obj.(flds{i}))==1;
                if isnan(obj.(flds{i}))
                    error(['nan fld ' flds{i} ]);
                else
                    continue
                end
            elseif isempty(obj.(flds{i}));
                if strcmp(flds{i},'rhoLB')
                    continue
                end
                error(['missing fld ' flds{i} ]);
            end

            % init
            if isempty(inds)
                inds=false( size( obj.(flds{i}) ));
            end

            % get inds
            inds= inds | Num.isimag(obj.(flds{i})) | isinf(obj.(flds{i}));
            if any(strcmp(flds{i},vflds)) || any(strcmp(flds{i},tflds))
                inds= inds | obj.(flds{i}) <= 0;
            elseif any(strcmp(flds{i},rflds))
                inds= inds | obj.(flds{i}) < -1 | obj.(flds{i}) > 1;
            end
        end

        % assign nan
        for i = 1:length(flds)
            if numel(obj.(flds{i}))<=1;
                continue
            end
            obj.(flds{i})(inds)=[];
            if ~isreal(obj.(flds{i})) && all(Num.isreal(obj.(flds{i})))
                obj.(flds{i})=real(obj.(flds{i}));
            end
        end
        if isempty(obj.rhoLB)
            obj.rhoLB=obj.covLB./sqrt(obj.varB .* obj.varL);
        end
        obj.nanInds=inds;

    end
%- PLOTTING
    function h=plotVarVar(obj,bCorr,LineWidth)
        P=ExtPlotVar(obj,'bCorr',bCorr,'LineWidth',LineWidth);
    end
    function h=plotVarCov(obj,bCorr,LineWidth)
        P=ExtPlotVar(obj,'bCorr',bCorr,'LineWidth',LineWidth);
        h=P.plot_var(stat);
    end
    function format_plotVarCov(obj,i,j,RC,loc,bCorr,FRCflds)
        if nargin < 5
            loc=[];
        end
        if nargin < 6
            bCorr=false;
        end
        if i == 1 && j==1
            ylabel('\rho');
            %titlef(obj);
        end
        if i == 1
            titlef(obj,j==1);
        end
        %if i == RC(1) && j==1
        %    if bCorr
        %        xlabel('corr(L,B)');
        %    else
        %        xlabel('cov(L,B)');
        %    end
        %end

        %% funcs
    end
    function text_plotVarCov(obj,i)
        % title
        txt=[sprintf('var_L= %.2f',obj.varL)];
        if numel(obj.T2)==1 && numel(obj.T3)==1
            txtT=sprintf('T    = %.2f\nT_{**} = %.2f',obj.T2,obj.T3);
            txt=[txtT newline txt ];
        end
        if i == 1
            txtR=['\rho' sprintf('    = %.2f',obj.rho22)];
            txt=[txt newline txtR ];
        end
        text(.25,.5,txt,'FontSize',12,'Interpreter','tex');
        axis off;
    end
    function [mStr,r2Str,vStr,t2Str,r3Str,t3Str]=get_titles(obj)
        % m
        mStr=obj.getMTitle();

        % rho
        if numel(obj.rho22)==1
            r2Str=['\rho^+= ' sprintf('%.3f\n',obj.rho22) ];
        else
            r2Str=[];
        end

        % rho
        if numel(obj.rho33) == 1
            r3Str=['\rho^+= ' sprintf('%.3f\n',obj.rho33) ];
        else
            r3Str='';
        end

        % VarL
        if numel(obj.varL)==1
            vStr=['var_L =' sprintf('%.3f',obj.varL)];
        else
            vStr='';
        end

        % T2
        if numel(obj.T2)==1
            t2Str=['T^+ =' sprintf('%.3f',obj.T2)];
        else
            t2Str='';
        end

        % T2
        if numel(obj.T3)==1
            t3Str=['T^* =' sprintf('%.3f',obj.T3)];
        else
            t3Str='';
        end

    end
    function mStr=getMTitle(obj)
        switch obj.model;
        case {'LmCR','LmCT'}
            m=obj.m;
        case 'LMCR'
            m=obj.M;
        otherwise
            m='';
            mStr='';
            return
        end

        frac=Num.toFrac(m);

        switch obj.model;
        case {'LmCR','LmCT'}
            mStr=['\sigma^2_B=' frac '\sigma^2_L'];
        case 'LMCR'
            mStr=['\sigma^2_{T^*}=' frac '\sigma^2_{T^+}'];
        end
    end
    function setPos(obj,f)
        if nargin < 2
            f=gcf;
        end
        set(f,'Position',[-1312 107 1189 1061]);
    end
end
methods(Static)
    function test()
        R=[ 0.3735    0.3267    0.2668];
        stdB=2;
        stdL=1;

        [varB,covLB,varL]=Ext.rhos2vars(R(1),R(2),R(3),stdB,stdL);
        [R2(1),R2(2),R2(3)]=Ext.vars2rhos(varB,covLB,varL,stdB,stdL);

    end
    function [rB,rLB,rL, corrLB]=vars2rhos(varB,covLB,varL,stdB,stdL)
        rB =(varB + varL + 2*covLB)./(stdB.^2);
        rLB=(varL + covLB)./(stdB .* stdL);
        rL =(varL)./(stdL.^2);
        corrLB=covLB./(sqrt(abs(varB)).*sqrt(abs(varL)));
    end
    function [varB,covLB,varL,corrLB,varIB,varIL,varEB]=rhos2vars(rB,rLB,rL,stdB,stdL)

        varL =rL  .* stdL.^2;
        covLB=rLB .* (stdB.*stdL)- varL;
        varB =rB  .* stdB.^2     - varL -2.*covLB;

        corrLB=covLB./(sqrt(abs(varB)).*sqrt(abs(varL)));

        varEB=(varL+varB+2*covLB);

        varIB=stdB.^2-varEB;
        varIL=stdL.^2-varL;

    end
    function ExtDemoFull()
        T2=1;
        T3=1.5;
        varL=0;
        varB=0.5;
        N=15;
        RC=[2 3];

        % cov by rho
        T23=sqrt(T2.*T3);
        %linL =linspace(0,T2,N);
        %linLC=linspace(0,T23-T2,N);
        %linC =linspace(0,T3-(2*(T23-T2))-T2,N);
        %
        linL =linspace(0,1,N);
        linLC=linspace(0,1,N);
        linC =linspace(0,1,N);

        [~,    VarL]=ndgrid(linLC,linL);
        [covLB,VarC]=ndgrid(linLC,linC);


        Fig.new('ExtDemoFull');
        %%- 1
        %[rho22,rho33,rho23]=Figs.extfuni(T2,T3,covLB,varL,VarC);
        %for i = 1:3
        %    plotfun(1,i);
        %end

        %%- 2
        %[rho22,rho33,rho23]=Figs.extfuni(T2,T3,covLB,VarL,varB);
        %for i = 1:3
        %    plotfun(2,i);
        %end


        %- 3
        RC=[1 3];
        [VarL,covLB,VarC]=meshgrid(linL,linLC,linLC);
        [rho22,rho33,rho23]=Figs.extfuni(T2,T3,covLB,VarL,VarC);
        for i = 1:3
            plotfun(1,i);
        end

        function plotfun(r,c)
            subPlot(RC,r,c);
            if c==1
                titl='\rho';
                rho=rho22;
                l=VarL(1,:,1);
                plot(l,rho(1,:,1));
                xlim([0 1]);
                ylim([0 1]);
                xlabel('varL');
                ylabel('\rho');
            elseif c==2
                titl='\rho^{*}';
                rho=rho23;
                x=VarL(:,:,1);
                y=covLB(:,:,1);
                z=rho(:,:,1);
                imagesc(x(:),y(:),z);
                xlabel('varL');
                ylabel('covLB');
                set(gca,'YDir','normal');
            elseif c==3
                titl='\rho^{**}';
                rho=rho33;
                extsurfcov(rho);
                caxis([0 1]);
            end
            %extimsagesc(rho);
            title(titl);
            axis square;
            if c==3
                cb=colorbar;
                cb.Position = cb.Position + [0.1 0 0 0];
            end
        end
        function extplotcov()
            if numel(rho22) < numel(covLB)
                rho22=repmat(rho22,1,numel(covLB));
            end
            if numel(rho33) < numel(covLB)
                rho33=repmat(rho33,1,numel(covLB));
            end
            if numel(covLB) < numel(rho33)
                covLB=repmat(covLB,1,numel(rho33));
            end
            hold off;
            plot(covLB,rho22,'b',...
                 covLB,rho33,'r',...
                 covLB,rho23,'m');
        end
        function extimagesc(rho)
            if numel(rho)==1
                rho=repmat(rho,size(covLB));
            end
            if numel(varB)==1
                varB=repmat(varB,size(covLB));
            end
            s=imagesc(covLB(:),varB(:),rho);
            xlabel('covLB');
            ylabel('varB');
        end
        function extsurfcov(rho)
            %size(VarL)
            %size(VarC)
            %size(covLB)
            %size(rho)
            %s=surf(VarL,VarC,covLB,rho);

            scatter3(VarL(:), VarC(:), covLB(:), 50, rho(:),'filled');

            %s=slice(VarL, covLB, VarC, rho, [.5],[.22 ],[.11]);
            %s.EdgeColor='none';
            xlabel('varL');
            ylabel('varB');
            zlabel('covLB');
        end
    end
end
end
