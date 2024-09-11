classdef ExtRNova < handle & ExtMdl
properties
    bWithin=true
    withinModel
end
methods(Static)
end
methods
    function obj=ExtRNova()
    end
    function args=mdl_params(obj)
        args={ ...
            'WithinModel',obj.withinModel ...
        };
    end
    function formula=get_formula(obj,name)

        num=name(end);
        if Str.Num.is(num)
            num=str2double(num);
            name=name(1:end-1);
        else
            num=1;
        end

        switch num
        case 1
            Y='R11-R22';
        case 2
            Y='R113-R223';
        otherwise
            formula=[];
        end

        switch name
        case 'base'
            X='D';
        otherwise
            formula=[];
            return
        end
        formula=[Y ' ~ ' X];

    end
    function fit(obj,name)
        if nargin >= 2 && ~isempty(name)
            obj.name=name;
        end
        formula=obj.get_formula(obj.name);
        disp(obj.name);
        disp(['  ' formula]);

        args=obj.mdl_params();
        if endsWith(obj.name,'1')
            name='splData';
            mname='mTableSpl';
        elseif endsWith(obj.name,'2')
            name='splData2';
            mname='mTableSpl2';
        end
        data=obj.(name);
        if obj.bWithin
            w=obj.parent.get_within(name);
            args=[args 'WithinDesign', w];
        end

        mdl = fitrm(obj.(mname),formula,args{:});
        obj.mdl=mdl;

        disp(obj.mdl);
        %disp(obj.Rsqr.(obj.name));

        % symmetry caseorder fitted probability
    end
end
end
