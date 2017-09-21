function ss=protein2dnaOptimized(str)
a=1
s=''
len=length(str);
k=strfind(str,'ATG');
filename='codons.csv';
I=readtable(filename);
Ix=readtable(filename);
I=table2array(I(:,1:2));
I=cell2mat(I);
Ix=table2array(Ix(:,3));
for ii=1:3:len-2
    ax=0
    axx=1
    kx=[]
    kxx=[]
    for i=1:64
        if str(ii:ii+2)==I(i,1:3)
            ax=ax+1;
            axx=axx+1;
            kx(ax)=i;
            kxx(axx)=Ix(i);
        elseif str(ii:ii+2)=='End'
            break
        end
    end
    kxx
    maxkxx=max(kxx);
    len=length(kxx)
    for iii=1:len
        if kxx(iii)==maxkxx
            axxx=iii-1;
        end
    end
    s(a:a+2)=I(kx(axxx),4:6);
    a=a+3;
end
ss=s
end