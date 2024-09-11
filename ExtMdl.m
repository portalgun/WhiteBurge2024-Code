classdef ExtMdl < handle
properties
    mdl
    name

    Mdls

    LL
    Rsqr

end
properties(Hidden)
    data
    splData
    mTable
    mTableSpl
    mTableSpl2

    parent
end
methods(Static)
end
methods
    function set.name(obj,name)
        obj.name=name;
        if ~isfield(obj.Mdls,name)
            obj.Mdls.(name)=struct();
        end
    end
    function out=get.mdl(obj)
        out=obj.Mdls.(obj.name);
    end
    function out=set.mdl(obj,in)
        obj.Mdls.(obj.name)=in;
    end
    function out=get.data(obj)
        out=obj.parent.data;
    end
    function out=get.splData(obj)
        out=obj.parent.splData;
    end
    function out=get.mTable(obj)
        out=obj.parent.mTable;
    end
    function out=get.mTableSpl(obj)
        out=obj.parent.mTableSpl;
    end
    function out=get.mTableSpl2(obj)
        out=obj.parent.mTableSpl2;
    end
    function out=get.LL(obj)
        out=obj.get_('LL');
    end
    function out=get.Rsqr(obj)
        out=obj.get_('Rsqr');
    end
    function out=hasname(obj,name)
        out=~isempty(obj.get_formula(name));
    end
    function out=get_(obj,Fld)
        switch Fld
        case 'LL'
            F={'LogLikelihood'};
        case 'Rsqr'
            F={'Rsquared','Ordinary'};
        end

        out=struct();
        flds=fieldnames(obj.Mdls);
        for i = 1:length(flds)
            mdl=obj.Mdls.(flds{i});
            if isempty(mdl) || ~isprop(mdl,F{1})
                out.(flds{i})=nan;
            else
                if numel(F)==1
                    out.(flds{i})=mdl.(F{1});
                else
                    out.(flds{i})=mdl.(F{1}).(F{2});
                end
            end
        end
    end
end
end
