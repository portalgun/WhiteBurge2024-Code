classdef ExtGen < handle & ExtMdl
properties
    Z

    DIdeal

    N

    S
    B
    L
    Inds
    nS
    nB
    nL
    sInd
    bInd
    lInd

    distMeth=0

    dist
    dcmp
    cstd
end
methods(Static)
end
methods
    function init(obj)
        dstd=obj.data('stdD');
        dcmp=obj.data('cmpD');
        obj.DIdeal=dcmp-dstd;

        switch obj.distMeth
        case 0
            obj.dist=dstd;
        case 1
            obj.dist=max([dstd dcmp]);
        case 2
            obj.dist=(dstd+dcmp)./2;
        end

        obj.N=length(obj.data);
        obj.Z=zeros(obj.N,1);

        s=obj.data{'subj'};
        b=obj.data('B').ret();
        d=obj.dist;
        [obj.S,~,obj.sInd]=unique(s);
        [obj.B,~,obj.bInd]=unique(b);
        [obj.L,~,obj.lInd]=unique(d);
        obj.nS=numel(obj.S);
        obj.nB=numel(obj.B);
        obj.nL=numel(obj.L);

        obj.Inds=Set.distribute(1:obj.nS,1:obj.nB,1:obj.nL);

    end
    function gen(obj)
        %[vL0 kl kb vs]
        vL0=params(1);
        kl=params(2);
        kb=params(3);
        vI=params(4:end);

        obj.kbfun(kb);
        vL=vL0.*exp(kl*obj.dist);
        vB=(kb*obj.dist)*vL;

        vS=obj.Z;
        for s = 1:obj.nSubj
            VS(obj.sInd==s)=vS(s);
        end
        %p=[s vL vB LB vD0 vS...];

        vT=vL+vB+vS;
        I=normrnd(0,vS,N);
        LB=0;

        R=(obj.DIdeal + vT + LB);
        R(R >= 1)=1;
        R(R  < 1)=0;
    end
    function out=kbfun(kb)
        out=kb*obj.dist;
    end
end
end
