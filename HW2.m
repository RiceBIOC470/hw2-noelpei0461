%% Problem 1. 
% in the repository you will find the meannonan.m function we discussed in
% class which produced the mean of a vector of numbers that ignores values
% of NaN or Inf. 

% Part 1. Run the following code:

xx = rand(5); % random 5x5 matrix
xx(3,2) = NaN; %put a NaN in
yy = mean(xx); 
zz = meannonan(xx);

%compare the size of yy and zz. Notice that yy produces a row vector (the
%average of the rows) and has a NaN in the column that contained the NaN whereas 
%zz is a single number (the average of all non-NaN entries in xx). Explain
%this behavior. 

% when NaN is input in the matrix, the regular mean cannot exclude it, and only do the mean for every column, and gives a row vector contains NaN in it. But for zz, it eliminates the NaN inside it, then mean function can take the average for all numbers inside. 
%Thus yy gives a row vector, while zz gives a single number.


% Part 2. Modify the meannonan code so that it behaves as the mean function
% and produces a row vector where each entry is the average of each column
% and in the column with a NaN, this NaN is ignored. 
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




%% Problem 2. ORFs using functions
% In this problem we will use functions to simplify and extend our code from HW1, prob 2 

% Part 1. Fill in the function randdnaseq.m in this repository so that it returns a random sequence
% of length N. 

function rr=randdnaseq(x)
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
rr=d;
end


% Part 2. Fill in the function findORF.m in this repository so that takes any dna
% sequence as an input and returns the length of the longest open
% reading frame and the positions of the start and stop codons. 
% Decide what your code should do when no ORF is found and
% implement this. Your function should also work whether the entered dna
% sequence is uppercase, lowercase, or some mixture. The builtin MATLAB functions
% lower and upper could be useful for this. 
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
end


% Part 3. Write another function called probabilityORF that utilizes the functions from 
% Parts 1 and 2. It should take two inputs - a sequence length (N) and an length  of an ORF (N_ORF) and
% returns the probability that that a sequence of length N contains an ORF
% of at least length N_ORF
function Css=PossibilityORF(N,Nx)
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
Css=Rxx(N)
end

PossibilityORF(N)
% Part4. Write  a final function called plotProbabilityORF.m which takes
% N_ORF as an argument and makes a plot of the probabily of having an
% ORF at least this long as a function of the dnasequence length. Decide how the
% code should determine the lengths of dna sequence to test and implement
% your decision. 
function xx=plotProbabilityORF(N_ORF)
Seql=N_ORF*10
k=[]
k(1)=0
for i=1:Seql
    k(i+1)=ProbabilityORF(i,N_ORF);
end
Seqlx=0:1:Seql;
xx=plot(Seqlx,k)

% Part 5. Write code that uses your function from part 4 to make a single
% plot with separate curves for ORF lengths 100,200,300,and 400. Make sure
% your plot has appropriate axis labels and legend. 

%% Problem 3. Codon translation and optimization

% DNA sequence gets translated into protein through a code known as the
% genetic code. Every sequence of 3 base pairs (a codon) is translated into
% 1 amino acid. The first two columns of the file codons.csv file in this repository 
% give the correspondence between codons and amino acids. 

% Part 1. Fill in the function dna2protein.m so that it takes any
% dnasequence and translates it to protein. The second argument -  frame -
% should take on values of 1,2, or 3 and should refer to whether the
% translation should start from the 1st, 2nd or 3rd base pair (that is,
% which reading frame to use). Make your code returns an error and
% appropriate message if frame isn't 1,2, or 3. 
function ss=dna2protein(d,x)
a=1
s=''
len=length(d);
k=strfind(d,'ATG');
filename='codons.csv';
I=readtable(filename);
I=table2array(I(:,1:2));
I=cell2mat(I);
if x~=1
    if x~=2
        if x~=3
            disp('error')
        else
            for ii=3:3:len-2
                for i=1:64
                    if d(ii:ii+2)==I(i,4:6)
                        s(a:a+2)=I(i,1:3);
                        a=a+3;
                    elseif d(ii:ii+2)=='TAA'
                        break
                    elseif d(ii:ii+2)=='TAG'
                        break
                    elseif d(ii:ii+2)=='TGA'
                        break
                    else 
                        break
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
                    elseif d(ii:ii+2)=='TAA'
                        break
                    elseif d(ii:ii+2)=='TAG'
                        break
                    elseif d(ii:ii+2)=='TGA'
                        break
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
                    elseif d(ii:ii+2)=='TAA'
                        break
                    elseif d(ii:ii+2)=='TAG'
                        break
                    elseif d(ii:ii+2)=='TGA'
                        break
                    end
                end
    end    
ss=s
end




% Part 2. Write code to turn your protein sequence back into DNA sequence.
% Call your function protein2dna.m
% Notice that there isn't a unique way to do this. For example, there are 4
% different codons that correspond to the amino acid Gly. For a first pass,
% choose one of these codons at random.
%
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
% Part 3. The third column of the codons.csv file contains the frequency of
% this codon's use in the human proteome in units of number of appearances per
% thousand codons. Some codons are used more than others. For example,
% for the amino acid Gln, the codon CAG is used nearly 3 times as often as
% CAA. When researchers add DNA to human cells, if it contains CAG rather than
% CAA for Gln, it will be translated  more efficiently. The process of
% taking a protein seqeunce and finding the DNA sequence that will be
% translated most efficiently in a particular organism is called codon
% optimization. Copy your function protein2dna.m from part 2 to a new
% function called protein2dnaOptimized.m that produces a codon-optimized DNA sequence using the 
% information in the third column of codons.csv. 
% In other words, for any amino acid, it always uses the codon that appears
% most frequently in the human proteome. 
%
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
ss=s;
end
