function dag = cpdag_to_dag2(cpdags)
% dags = cpdag_to_dag(cpdags)
%
% CPDAG_TO_DAG produce a N*N matrix of a dag which instantiate cpdag.
% (also works with a cell array of cpdags, returning a cell array of dags)
% make sur that your entry is a completed PDAG
% this function can't be use instead of PDAG_TO_DAG
%
% francois.olivier.c.h@gmail.com
% 7 may 2003 - OLD version

if ~iscell(cpdags)
    cpdag=cell(1,1);
    cpdag{1}=cpdags;
else
    cpdag=cpdags;
end

for da=1:length(cpdag)

    N=length(cpdag{da});
    dag=cpdag{da}; dag2=dag;
    unprocessed = [];

    for i=1:(N-1)
        for j=(i+1):N
            if dag2(i,j)==1 & dag2(j,i)==1
                if ~myismember(i,unprocessed)
                    unprocessed = [unprocessed;i];
                end
                if ~myismember(j,unprocessed)
                    unprocessed = [unprocessed;j];
                end
            end
        end
    end

    for i=1:length(unprocessed)
        for j=1:N
            if dag(unprocessed(i),j)==1 & dag(j,unprocessed(i))==1
                dag(j,unprocessed(i))=0;
            end
        end
    end
    dags{da}=dag;
end

if ~iscell(cpdags)
    dag=dags{1};
else
    dag=dags;
end