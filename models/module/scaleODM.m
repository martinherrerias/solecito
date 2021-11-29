function ODM = scaleODM(ODM,Ns,Np,unknown)
% P = SCALEODM(p,Ns,Np) - Get a One-Diode-Model parameter structure that is equivalent to
%   an Ns x Np series-parallel array of the original ODM. 
%   Fractional scales (e.g. 1/Ns and 1/Np) can be used to get sub-elements of an array:
%
%       p = SCALEODM(P,1/Ns,1/Np)  <->  P = SCALEODM(p,Ns,Np)
%
%   NOTE: ODM.name, if available, will be modified to reflect the scaling operation.
%   Fields of the structure that loose meaning during scaling (Nbpd, width, length) will be set
%   to NaN, if available. 
%
% P = SCALEODM(p,Ns,Np,KEY) - KEY = {'keep','rm','NaN'} determines what to do with fields that
%   don't belong to the standard parameter output of CHECKODM. Default is 'rm'.
%
% See also: TRANSLATEODM, CHECKODM

    narginchk(3,4);
    if nargin < 4 || isempty(unknown), unknown = 'rm'; end
    
    VOLTAGE = {'nVth','Ns','Voc_ref','Vmp_ref','Vbi','nVth_ref','di2mutau'};
    CURRENT = {'Iph','Io','Np','Isc_ref','Imp_ref','Iph_ref','Io_ref'};
    RESISTANCE = {'Rs','Rsh_ref','Rsh_0','Rsh_base'};
    POWER = {'area','Pmp_ref'};
    NDEF = {'Nbpd','width','length'}; 
        
    CONSTANT = {'AM_ref','Eg_ref','G_ref','NOCT','Rsh_exp','Tc_ref','eta_ref','isCPV',...
        'material','muEg','muIsc','muNdiode','muPmp','muVoc','nDiode','nJunct',...
        'avalanche_a','avalanche_m'};
    
    OTHER = {'name','Brev','Vbr','Rsh_model','Io_model','reverse_model'};
    
    scalefields(POWER,Ns*Np);           % Power scaling
    scalefields(CURRENT,Np);            % Current scaling
    scalefields(VOLTAGE,Ns);            % Voltage scaling
    scalefields(RESISTANCE,Ns/Np);      % Resistance scaling
    scalefields(NDEF,NaN);              % Undefined
    % scalefields(CONSTANT,1);          % Not scaled
    
    % Only modify Brev/Vbr on sets of parameters that have already been translated
    TRANSLATED = {'Iph','Io','Rs','nVth','Rsh','di2mutau','Vbi','Brev','Vbr','avalanche_a','avalanche_m'};
    if isempty(setdiff(fieldnames(ODM),TRANSLATED))
        if isfield(ODM,'Brev'), ODM.Brev = ODM.Brev*Np/Ns^2; end
        if isfield(ODM,'Vbr'), ODM.Vbr = ODM.Vbr*Ns; end
    end
    
    if isfield(ODM,'name')
    % Modify the name. Use rational expression, if possible
        [N,D] = rat([Ns,Np]);
        if all(D == 1) || ~all([Ns,Np] == N./D)
            ODM.name = sprintf('%s - scaled(%g,%g)',ODM.name,Ns,Np);
        else
            ODM.name = sprintf('%s - scaled(%d/%d,%d/%d)',ODM.name,N(1),D(1),N(2),D(2));
        end
    end
    
    ukf = setdiff(fieldnames(ODM),[VOLTAGE,CURRENT,RESISTANCE,POWER,NDEF,OTHER,CONSTANT]);
    switch lower(unknown)
        case 'keep'
        case 'rm', ODM = rmfield(ODM,ukf);
        case 'nan', scalefields(ukf,NaN);
        otherwise
            error('Expecting ''keep''/''rm''/''NaN''');
    end
       
    function scalefields(c,s)
        for j = 1:numel(c)
            if isfield(ODM,c{j})
                ODM.(c{j}) = ODM.(c{j})*s;
            end
        end
    end
end

