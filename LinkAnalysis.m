clear

%% load data
% name=[1:6];
% for(i=1:6)
%     load(['graph_',num2str(name(i)),'.txt']);
%
% end

GP=load('graph_3.txt');

[Adj,V,v]=ADJ(GP);
ethon = 0.01;
[HITSAns]=HITS(v,Adj,ethon);
[PGRank]=PageRank(v,Adj,ethon);



%% SimPage
d=ones(1,v);
SI=diag(d);
S=SI;
AdjT=Adj';
% St=zeros(v,v);
StCon=1;

t=0;
while(StCon>=ethon )%   t<1
    
    for(i=1:v) %col
        for(j=i+1:v) %row
            [S(i,j)] = Sim(SI,Adj,i,j);
            S(j,i)=S(i,j);
        end
    end
    t=t+1;
    StCon=norm(S-SI);
    SI = S;
end

function [SimVal] = Sim(SI,Adj,i,j)
AdjT=Adj';
Ii=find(AdjT(i,:));
Ij=find(AdjT(j,:));
if(length(Ii)==0 || length(Ij)==0)
    SimVal=0;
else
    temp=0;
        for(i=Ii)
            for(j=Ij)
                       temp = SI(i,j)+temp;
            end
        end
    SimVal= (0.6/(length(Ii)*length(Ij))) * temp;
end
end


%% PageRank %gp1 faild %rank uncorrect
function [rank]=PageRank(v,Adj,ethon)
for(i=1:v)
    PR(i,1)=1/v;
end
AdjT=Adj'; %in-degree


CH=0;
t=1;
StopC=zeros(v,1);
flag=0;
while(flag~=1) %stop condition
    
    t=t+1;
    PR(:,t)=zeros(v,1);
    for(i=1:v) %items
        CH=find(AdjT(i,:)); % nodes(index) connected in the row
        if(CH~=0)
            for(j=1:length(CH)) %non zero element in adjT
                GCH=sum(AdjT(CH(j),:));
                if(GCH~=0)
                    PR(i,t)=PR(CH(j),t-1)/sum(AdjT(CH(j),:))+PR(i,t);
                end
                
            end
        end
        StopC(i,t)=norm(PR(i,t)-PR(i,t-1))<ethon;
    end
    %stop condition
    if(StopC(:,t)==ones(v,1))
        flag=1;
%         PR=ones(v,1);
    end
    
end

%rank
[RkPR,PgName]=sort(nonzeros(PR(:,t)),'descend');
%[~,~,rank]=unique(PR(:t));
rank=[PgName,RkPR];
rank(:,3)=1:length(PgName);
end

%% HITS
function [Ans]=HITS(v,Adj,ethon)
condition=1;

t=1;
Auth=ones(v,1);
Hub=Auth;
while(t>=1 && condition>ethon)
    
    [Auth(:,t+1),Hub(:,t+1)]=recursive(Auth(:,t),Hub(:,t),Adj);
    t=t+1;
    condition = norm(Auth(:,t)-Auth(:,t-1))+norm(Hub(:,t)-Hub(:,t-1));
    
end
Ans(:,1)=Auth(:,t);
Ans(:,2)=Hub(:,t);


    function [Auth,Hub]= recursive(Auth,Hub,Adj)
        
        Auth=Adj'*Hub; %Adj' = in-degree
        Hub=Adj*Auth; %Adj = out-degree
        Auth=Auth/norm(Auth);
        Hub=Hub/norm(Hub);
        
    end
end

%% Generate the adjacency matrix
function [Adj,V,v] = ADJ(GP)
V=unique(GP);
v=length(V);
[row,col]=size(GP);
Adj=zeros(v,v);
for(i=1:row)
    for(j=1:col)
        Adj(GP(i,1),GP(i,2))=1;
    end
end
end

