function randomSeq = randdnaseq(N)
% returns a random dna sequence of length N
s=randi(4,1,x);
% set up an array with 500 random numbers, range 1 to 4 in it.
d='';
for ii=1:x
    if s(ii)==1
        d(ii)='A'
    elseif s(ii)==2
        d(ii)='T'
    elseif s(ii)==3
        d(ii)='G'
    else d(ii)='C'
    end
%For every 1 in s, input A in d; and substitue rest values as well.
end
randomSeq=d;
end
