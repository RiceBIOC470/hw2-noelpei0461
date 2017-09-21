function Css=ProbabilityORF(N,Nx)
    Rxx=[];
    %Rxx is an array that will be used to store the possibility values
    for xxx=3:N
        a=0;
        for i=1:1000
            s=randi(4,1,xxx);
            d='';
            for ii=1:xxx
                if s(ii)==1
                    d(ii)='A';
                elseif s(ii)==2
                    d(ii)='T';
                elseif s(ii)==3
                    d(ii)='G';
                else d(ii)='C';
                end
            end
            k=strfind(d,'ATG');
            len=length(k);
            dist=[];
            dist2=[];
            g=1;
            for i=1:len
                for ii=k(i)+3:3:xxx-2
                    if d(ii:ii+2)=='TAA'
                        m=ii;
                        dist(g)=m-k(i);
                        g=g+1;
                    elseif d(ii:ii+2)=='TAG'
                        m=ii;
                        dist(g)=m-k(i);
                        g=g+1;
                    elseif d(ii:ii+2)=='TGA'
                        m=ii;
                        dist(g)=m-k(i);
                        g=g+1;
                    end
                end
                A=min(dist);
                if isempty(A)
                    A=0;
                end
                dist2(i)=A;
            end
            distmax=max(dist2);
            if distmax>Nx
                a=a+1;
            else
                a=a;
            end
        end
        xx=a/1000;
        Rxx(N)=xx;
    end
Css=Rxx(N);
end