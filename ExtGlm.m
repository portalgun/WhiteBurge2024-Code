classdef ExtGlm < handle & ExtMdl
properties
    % mdl
    % name
    % Mdls

    % LL
    % Rsqr

    link='logit'
    fitMethod='Laplace'
    distribution='binomial'

    %probability
    %caseorder
    %fitted
    %histogram
    %lagged
    %symmetry
end
methods(Static)
end
methods
    function args=mdl_params(obj)
        args={ ...
            'Distribution',obj.distribution, ...
            'Link',obj.link, ...
            'FitMethod',obj.fitMethod ...
            'BinomialSize',1, ...
            %'DummyVarCoding','full' ...
            %'DispersionFlag',false ...
            %'CovariancePattern','FullCholesky'...
        };
    end
    function formula=get_formula(obj,name)

        if nargin < 2 || isempty(name)
            name='';
        end
        switch name
        case 'd'
            % 46 - 1
            formula='R ~ -1 + d ';
        case 's'
            % 53 - 3 + 1
            formula='R ~ -1 + (-1 + d | subj)';
        case 'sb'
            % 60 - 6 + 2
            formula='R ~ -1 + (-1 + d | subj) + (-1 + d | B)';
        case 'SB'
            % 67 - 9 + 1
            formula='R ~ -1 + (-1 + d | subj:B)';
        case 'SBD'
            % 72
            formula='R ~ -1 + (-1 + stdd | subj:B) + (-1 + cmpd | subj:B)'
        case 'new'
            % 68 - 24 + 2
            formula='R ~ -1 + (-1 + d | subj:B) + (-1 + d | stdd:B)'
            % 67 - 18 + 3 + 1
            formula='R ~ 1 + (d | subj:B)';
        case 'new2'
            formula='R ~ -1 + (-1 - mn | subj:B) + (-1 + d | subj:B)';
            % 60 - 9 + 1
            %formula='R ~ -1 + (-1 + mn | subj:B)';
            %% 38
            %formula='R ~ -1 + (-1 + mn | subj)';
            %% 02
            %formula='R ~ -1 + (-1 + mn | B)';
            %% 27
            %formula='R ~ -1 + (-1 + mn | B) + (-1 + mn | subj)';
            %% 16
            formula='R ~ -1 + (-1 + (d + sz24) | subj)';
            %formula='R ~ -1 + (-1 + d | subj:B:lvlInd)'
        case 'd0'
            mdl='R ~ -1 + (-1 + sz4|subj)';
        otherwise
            if ischar(name) & contains(name,'~')
                formula=name;
            else
                ormula=[];
            end
        end
        %% 60 - 24 + 2
        %formula='R ~ -1 + (-1 + d | subj:B) + (-1 + d | subj:stdd)'

    end
    function fit(obj,name,varargin)
        if nargin < 3 || isempty(varargin{1})
            tbl=obj.mTable;
        else
            dat=obj.data.copy();
            dat=dat(varargin{:});
            dat.unique('d')
            tbl=dat.toMTable();
        end

        obj.link='probit';
        obj.fitMethod='MPL';
        obj.distribution='Binomial';
        if nargin >= 2 && ~isempty(name)
            obj.name=name;
        end
        formula=obj.get_formula(obj.name);
        disp(obj.name);
        disp(['  ' formula]);

        args=obj.mdl_params();
        mdl = fitglme(tbl,formula,args{:});
        obj.mdl=mdl;

        disp(obj.mdl);
        disp(obj.Rsqr.(obj.name));

        % symmetry caseorder fitted probability
    end
end
end
