classdef ExtData < handle
properties
    aliases
    subjs
    moude
    status
    lvlInds
    blks
    passes
    hostname

    rawBlk
    Blk

    Ptchs
    patch
    Stats

    raw
    data
    splData
    splData2
    mTable
    mTableSpl
    mTableSpl2

    idsD
    idsB
end
properties(Access=protected)
    prps
    methds
    flds
end
methods(Static)
    function get_ext_data(obj)
        E=EPs.ext_data(obj.dims,obj.subjs,obj.passes,obj.model,obj.bNew);
    end
end
methods
    function out=get.patch(obj)
        out=obj.Ptchs.ptch;
    end
    function obj=ExtData()
        obj.init();
    end
    function init(obj,bRaw);
        if nargin < 2 || isempty(bRaw)
            bRaw=isempty(obj.raw);
        end
        obj.prps=fieldnames(Obj.struct(obj));
        obj.methds=methods(obj);
        obj.flds=[obj.prps; obj.methds];

        if bRaw
            obj.get_raw_data();
        end
        obj.get_data();
        obj.table2table();


        if bRaw
            obj.get_raw_blk();
        end
        obj.get_blk();
        obj.get_ptchs();
    end
    function get_stat(obj,typ,stat)
        typ='fm';
        stat='cohere';

        N=length(obj.Ptchs.fnamesB);
        out=zeros(N,1);
        for i = 1:N
            obj.Ptchs.get_patch(i);
            out(i)=obj.fft2(typ,stat);
        end
        obj.Stats.([typ '_' stat])=out;
    end
    function out=fft2(obj,typ,stat,bPlot)
        if nargin < 2 || isempty(typ)
            typ='h';
        end
        if nargin < 3 || isempty(typ)
            stat='dot';
        end
        if nargin < 4 || isempty(bPlot)
            bPlot=false;
        end

        p=obj.patch();
        sz=p.PszRC;
        ctr=sz/2;
        ctr=ctr+(mod(ctr,2)==0);
        d=max(sz./(sqrt(2)*sz/2));

        if bPlot
            C=length(p.im.img);
            RC=[6 C];
            Fig.new('test fft2');
        end

        img=p.im.img;
        for i = 1:length(img)
            im=img{i};
            fm0=fftshift(fft2(ifftshift(im)));

            % rm dc
            fm{i}=fm0;
            DC(i)=fm{i}(ctr(1),ctr(2));
            fm{i}(ctr(1),ctr(2))=0;

            % phase
            switch typ
            case 'abs'
                gm{i}=abs(fm{i});
            case 'fm'
                gm{i}=fm{i};
            case 'hc'
                gmr=real(fm{i}(ctr(1),:));
                gmi=imag(fm{i}(ctr(1),:));
                gm{i}=complex(gmr,gmi);
            case 'h'
                % no phase (hoizontal only
                gm{i}=abs(fm(ctr(1),:));
            case 'r'
                [T,A]=obj.radial_avg(abs(fm{i}),d);
                A=T';
                B=fliplr(A(2:end));
                gm{i}=[B A];
            case 'r2'
                [T,A]=obj.radial_avg(cos(angle(fm{i})),d);
                A=A';
                B=fliplr(A);
                B(1)=0;
                gi=[B A];
                gr=real(fm{i});
                gm{i}=complex(gr,gi);
            end

            if bPlot
                gmm=fftshift(ifft2(ifftshift(fm{i})));
                if size(gm,1)==1
                    am=fftshift(ifft(ifftshift(gm{i})));
                    am=repmat(am,p.PszRC(1),1);
                else
                    am=fftshift(ifft2(ifftshift(gm{i})));
                end

                plotfun(1,im);
                plotfun(2,gmm);
                plotfun(3,real(fm{i}));
                plotfun(4,imag(fm{i}));
                plotfun(5,abs(fm{i}));
                plotfun(6,am);
            end
        end

        switch stat
        case 'dot'
            out=real(dot(gm{1}(:),gm{2}(:))/real(norm(gm{1}(:)) .* norm(gm{2}(:))));
            %figure(33);
            %imagesc([gm{1} gm{2}])
            %drawnow
            %pause(.1);
        case 'cohere'
            A=gm{1}(1,:);
            B=gm{2}(1,:);
            cx=abs(A.*B).^2;
            cn=(abs(A.*conj(A)) .* abs(B.*conj(B)));
            cc=real(cx./ cn);
            cc(isnan(cc))=1;
            out=mean(cc(:));


            gmm{1}=fftshift(ifft(ifftshift(fm{1}(1,:))))';
            gmm{2}=fftshift(ifft(ifftshift(fm{2}(1,:))))';
            A=fm{1}(1,:);
            B=fm{2}(1,:);
            cx=abs(A.*B).^2;
            cn=(abs(A.*conj(A)) .* abs(B.*conj(B)));
            cc=cx./ cn;
            cc(isnan(cc))=1;
            out=mean(cc(:));


            PAi=pwelch(gmm{1});
            PA=abs(A.*conj(A));
            PB=abs(B.*conj(B));
            PBi=pwelch(gmm{2});
            [Cxy,f] = mscohere(gmm{1},gmm{2});
            figure(33)
            hold off;
            %imagesc([real(A) real(B)]);
            %plot(PA);
            %plot(PA);
            plot(PAi); hold on
            plot(PA); hold on

            %plot(Cxy); hold on
            %plot(cc);
            hold off;
            drawnow;

            %plot(cc(:),'.')

        end



        function plotfun(r,thing)
            subPlot(RC,r,i);
            imagesc(thing);
            axis image;
        end
    end
    function [Tics,Average]=radial_avg(obj,data,radial_step)
        %main axii cpecified:
        x=(1:size(data,2))-size(data,2)/2;
        y=(1:size(data,1))-size(data,1)/2;

        % coordinate grid:
        [X,Y]=meshgrid(x,y);
        % creating circular layers
        Z_integer=round(abs(X+1i*Y)/radial_step)+1;
        % % illustrating the principle:
        % % figure;imagesc(Z_integer.*data)
        % very fast MatLab calculations:
        Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
        Average=accumarray(Z_integer(:),data(:),[],@mean);
    end
%- Params
    function get_raw_data_params(obj)
        tbl=spec.lvlsTable();
        lvls=tbl(:,1);
        obj.aliases={'DSP2','DSP2f2'};
        %obj.aliases={'DSP2'};
        obj.subjs={'DNW','AJL','JDB'};
        obj.moude=1;
        obj.lvlInds=lvls;
        obj.passes=[1 2];
        obj.blks=1:5;
        obj.hostname='jburge-wheatstone.psych.upenn.edu';
    end
%- Blk
    function get_raw_blk(obj)
        obj.rawBlk=Blk.get(obj.aliases{1});
    end
    function get_blk(obj)
        B=obj.rawBlk.blk;
        B=B('mode',1,'lvlInd',obj.lvlInds);

        [obj.idsB,~,id]=B.unique_rows('lvlInd','blk','trl');
        if ~isequal(obj.idsB.ret(1),obj.idsD.ret(1))
            error('idsD and idsB are not equal')
        end
        B.add_col('id',id,1);
        %B.rm_col({'cndInd','cmpInd','lvlInd','mode','blk','trl'});
        B.sort('id','cmpNum');
        obj.Blk=obj.rawBlk.copy();
        obj.Blk.blk=B;
    end
    function get_ptchs(obj)
        p=ptchs.getRaw(obj.aliases{1},[],obj.hostname);
        p.apply_block(obj.Blk);
        obj.Ptchs=p;
    end

%- EP
%- RAW
    function get_raw_data(obj)
        obj.get_raw_data_params();

        sels={...
              'subj', obj.subjs, ...
              'mode','%%',unique(obj.moude), ...
              'status','>' 0, ...
              'lvlInd', obj.lvlInds, ...
              'blk', obj.blks, ...
              'pass' obj.passes, ...
        };

        Data=[];
        for i = 1:length(obj.aliases)
            E=ETable.get(obj.aliases{i});
            data=E.loadRawData(sels{:});
            data.add_col('aliasInd',i,1);
            if i==1
                Data=data;
            else
                Data=Data.vertcat(data);
            end
        end
        obj.raw=Data;
    end
%- Processed
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
        data.rm_col('cmpB');

        %% cmpX and stdX to depth
        D=data.minus('cmpD','stdD');
        data.add_col('D',D,3);


        %% trlInfo to id
        [obj.idsD,~,id]=data.unique_rows('lvlInd','blk','trial');
        data.add_col('id',id,1);
        [~,~,id]=data.unique_rows('subj','lvlInd','blk','trial');
        data.add_col('ids',id,1);


        %% lvlInd and aliasInd to B
        bCor=data.bfind('aliasInd',1); %NOTE MAKE SURE TO CHECK THIS
        B=data('stdB').ret();
        B(~bCor)=0;
        data.add_col('B',B);

        %[~,~,id]=data.unique_rows('aliasInd','lvlInd','trial');
        % {'aliasInd','lvlInd'}

        %% lvlInd to L
        data.rename_col('stdB','L');

        data.rm_col({'lvlInd','blk','trial'});
        data.order_cols('ids','id','subj','cmpD','stdD','D','L','B','R','C','pass','aliasInd');

        data.unique('B');

        obj.data=data;
        obj.get_data_mod_();
    end
    function get_data_mod_(obj)
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
        obj.data=data;
    end
%- Table
    function table2split(obj)
        T=obj.table2split_fun(false);
        obj.splData=T;
    end
    function table2split2(obj)
        T=obj.table2split_fun(true);
        obj.splData2=T;
    end
    function T=table2split_fun(obj,b2,bSubj)
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
        T.rm_col({'L','B','pass'});
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
    function toCSV(obj)
        if isempty(obj.mTable)
            obj.mTable;
        end
        dire=Env.var('Data');
        fname=[dire 'Ext'];
        obj.data.toCSV(fname);
    end
end
end
