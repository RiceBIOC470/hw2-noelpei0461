function ss=protein2dna(str)
a=1;
s='';
len=length(str);
filename='codons.csv';
I=readtable(filename);
I=table2array(I(:,1:2));
I=cell2mat(I);
for ii=1:3:len-2
    ax=0;
    kx=[];
    for i=1:64
        if str(ii:ii+2)==I(i,1:3)
            ax=ax+1;
            kx(ax)=i;
        elseif str(ii:ii+2)=='End'
            break
        end
    end
    xa=randi(ax,1,1);
    s(a:a+2)=I(kx(1)-1+xa,4:6);
    a=a+3;
end
ss=s;
end