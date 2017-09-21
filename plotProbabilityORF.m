function xx=plotProbabilityORF(N_ORF)
Seql=N_ORF*2;
k=[];
k(1)=0;
for i=1:Seql
    if i<N_ORF
        k(i+1)=0;
    else
        k(i+1)=ProbabilityORF(i,N_ORF);
    end
end
Seqlx=0:1:Seql;
xx=plot(Seqlx,k)