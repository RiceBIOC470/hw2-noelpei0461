function ss=dna2protein(d,x)
a=1;
s='';
len=length(d);
filename='codons.csv';
I=readtable(filename);
I=table2array(I(:,1:2));
I=cell2mat(I);
%I is the matrix containing amino acids and corresponding codons infomation.
if x~=1
    if x~=2
        if x~=3
            s(1:5)='error';
            %if the input is not 1,2,3 returns an error.
        else
            for ii=3:3:len-2
                for i=1:64
                    if d(ii:ii+2)==I(i,4:6)
                        s(a:a+2)=I(i,1:3);
                        a=a+3;
                    end
                end
            end
        end
    else
        for ii=2:3:len-2
                for i=1:64
                    if d(ii:ii+2)==I(i,4:6)
                        s(a:a+2)=I(i,1:3);
                        a=a+3;
                    end
                end
        end
    end
else
    for ii=1:3:len-2
                for i=1:64
                    if d(ii:ii+2)==I(i,4:6)
                        s(a:a+2)=I(i,1:3);
                        a=a+3;
                    end
                end
    end    
end
ss=s;
end
