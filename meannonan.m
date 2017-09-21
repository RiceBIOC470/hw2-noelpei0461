function mm=meannonan(x)
ma=mean(x);
len=length(ma);
c=[]
for i=1:len
    mb=mean(x(:,i));
    if ~isnumeric(x)
        c(i)=mb
    else
        notin=isnan(x)|isinf(x);
        x(notin)=[];
        c(i)=mean(x(:,i))
    end
mm=c
end
end



