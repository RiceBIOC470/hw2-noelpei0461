function ff=findORF(d)
k=strfind(d,'ATG');
% find the location of all start condons 
len=length(k);
% find how many start condons are in it.
dist=[];
% dist is used to store every possible combination of a start codon to stop condons.
dist2=[];
% dist2 is used to store all the ORF length.
g=1;
for i=1:len
    for ii=k(i)+3:3:498
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
        %find all stop codons after the start codon.
        end
    end
    A=min(dist);
    % the smallest start/stop codon combination is the real ORF we are looking for.
    if isempty(A)
        A=0
    end
    % for the sequence that does not have any stop codon after their start codons, set the value of the length of ORF as 0. 
    dist2(i)=A;
    
end
distmax=max(dist2);
start=[];
stop=[];
for i=1:length(k)
    if dist2(i)==distmax
        start(1)=k(i);
        stop(1)=k(i)+distmax;
    end
end
ff=disp(distmax)|disp(start)|disp(stop)

%show the longest ORF in dist2.