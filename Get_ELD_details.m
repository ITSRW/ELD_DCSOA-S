function [fobj,PD,ELD] = Get_ELD_details(F)
switch F
    case 'ELD3'
        fobj = @TD;
        PD=850;
        load('ELD3.mat');
        ELD=ELD3;
    case 'ELD13'%???
        fobj = @TD;
        PD=1800;
        load('ELD13.mat');
        ELD=ELD13;
    case 'ELD40'
        fobj = @TD;
        PD=10500;
        load('ELD40.mat');
        ELD=ELD40;
end
end

% ELD
function o = TD(mat,x,PD)
Pmat=[mat,x];

fi=Pmat(:,3).*(Pmat(:,8).^2)+Pmat(:,4).*Pmat(:,8)+Pmat(:,5)+...
    abs(Pmat(:,6).*sin(Pmat(:,7).*(Pmat(:,8)-Pmat(:,1))))+...
    (PD-sum(Pmat(:,8)))^2;

o=sum(fi);
end


