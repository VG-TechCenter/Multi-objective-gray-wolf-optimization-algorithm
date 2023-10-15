% The Matlab source codes to generate the boudnaries of the test instances
%   for CEC 2009 Multiobjective Optimization Competition. 
% Please refer to the report for correct one if the source codes are not
%   consist with the report.
% History:
%   v1 Sept.05 2008

function range = xboundary(name,dim)

    range = ones(dim,2);
    
    switch name
        case {'UF1','UF2','UF5','UF6','UF7','CF2'}
            range(1,1)      =  0;
            range(2:dim,1)  = -1;
        case 'UF3'
            range(:,1)      =  0;  
        case {'UF4','CF3','CF4','CF5','CF6','CF7'}
            range(1,1)      =  0;
            range(2:dim,1)  = -2;
            range(2:dim,2)  =  2; 
        case {'UF8','UF9','UF10','CF9','CF10'}
            range(1:2,1)    =  0;
            range(3:dim,1)  = -2;
            range(3:dim,2)  =  2;   
        case 'CF1'
            range(:,1)      =  0; 
        case {'CF8'}
            range(1:2,1)    =  0;
            range(3:dim,1)  = -4;
            range(3:dim,2)  =  4;             
    end
end