%% EIASC algorithm
function [y,yl,yr,l,r]=EIASC(Xl,Xr,...
                       Wl,Wr,needSort)

ly=length(Xl); XrEmpty=isempty(Xr);
if XrEmpty;  Xr=Xl; end
if max(Wl)==0
    yl=min(Xl); yr=max(Xr);
    y=(yl+yr)/2;  l=1; r=ly-1; return;
end
if nargin==4;  needSort=1; end

% Compute yl
if  needSort
    [Xl,index]=sort(Xl); Xr=Xr(index);
    Wl=Wl(index); Wr=Wr(index);
    Wl2=Wl; Wr2=Wr;
end
if ly==1
    yl=Xl;  l=1;
else
    yl=Xl(end); l=0;
    a=Xl*Wl'; b=sum(Wl); 
    while l<ly && yl > Xl(l+1)
        l=l+1;
        t=Wr(l)-Wl(l);
        a=a+Xl(l)*t;
        b=b+t; 
        yl=a/b;   
    end
end

% Compute yr
if ~XrEmpty && needSort==1
    [Xr,index]=sort(Xr);
    Wl=Wl2(index); Wr=Wr2(index);
end
if ly==1
    yr=Xr; r=1;
else
    r=ly; yr=Xr(1); 
    a=Xr*Wl'; b=sum(Wl);
    while r>0 && yr < Xr(r)
        t=Wr(r)-Wl(r);
        a=a+Xr(r)*t;
        b=b+t; 
        yr=a/b;  r=r-1;
    end
end
y=(yl+yr)/2;
