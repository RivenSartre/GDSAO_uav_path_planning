function [Best_pos,Best_score,Convergence_curve]=GDSAO_Method(nPop, MaxIt, nVar, VarMin, VarMax, CostFunction, model)
dim=nVar*3;
lb = [VarMin.r *ones(1,nVar), VarMin.psi *ones(1,nVar), VarMin.phi *ones(1,nVar)];
ub = [VarMax.r *ones(1,nVar), VarMax.psi *ones(1,nVar), VarMax.phi *ones(1,nVar)];

%% 
X=initialization_GDSAO(nPop,dim,ub,lb);

way1r=linspace(model.start(1), model.end(1), nVar+2);
way1psi=linspace(model.start(2), model.end(2), nVar+2);
way1phi=linspace(model.start(3), model.end(3), nVar+2);
X(1,:)=[way1r(2:end-1), way1psi(2:end-1), way1phi(2:end-1)];

Best_pos=zeros(1,dim);
Best_score=inf;
Objective_values = zeros(1,size(X,1));

Convergence_curve=[];
N1=floor(nPop*0.5);
Elite_pool=[];
WaitbarInter = MaxIt / 100;  
tic
h = waitbar(0, ['Completed:0%   Completeing...time:', num2str(toc)]);
Position = struct('r', [], 'psi', [], 'phi', []);


%% 
for i=1:size(X,1)

    Position.r = [X(i, 1:nVar)];
    Position.psi = [X(i, nVar+1:nVar*2)];
    Position.phi = [X(i, 2*nVar+1:end)];
    Objective_values(1,i)=CostFunction(SphericalToCart(Position,model));
    
    if i==1
        Best_pos=X(i,:);
        Best_score=Objective_values(1,i);
        %Best_sol=sol(:,i);
    elseif Objective_values(1,i)<Best_score
        Best_pos=X(i,:);
        Best_score=Objective_values(1,i);
        %Best_sol=sol(:,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end

[~,idx1]=sort(Objective_values);
second_best=X(idx1(2),:);
third_best=X(idx1(3),:);
sum1=0;
for i=1:N1
    sum1=sum1+X(idx1(i),:);
end
half_best_mean=sum1/N1;
Elite_pool(1,:)=Best_pos;
Elite_pool(2,:)=second_best;
Elite_pool(3,:)=third_best;
Elite_pool(4,:)=half_best_mean;

Convergence_curve(1) = Best_score;

for i=1:nPop
    index(i)=i;
end
std1=std(Objective_values(:)) ;
sigma(1) =1;
Na=nPop/2;
Nb=nPop/2;
DDF(1)=0.35;
DDF(2)=0.35;
b=0.6;
M(1)=0.35;

%% main loop
l=2; % 
while l<=MaxIt
    RB=randn(nPop,dim);          
    T=exp(-l/MaxIt);
    
    M(l)=DDF(l)*T;
    
    for j=1:dim
        sum1=0;
        for i=1:nPop
            sum1=sum1+X(i,j);
        end
        X_centroid(j)=sum1/nPop;
    end
    
    
    index1=randperm(nPop,Na);
    index2=setdiff(index,index1);
    
    for i=1:Na
        r1=rand;
        k1=randperm(4,1);
        for j=1:size(X,2) % in j-th 
            X(index1(i),j)= Elite_pool(k1,j)+RB(index1(i),j)*(r1*(Best_pos(j)-X(index1(i),j))+(1-r1)*(X_centroid(j)-X(index1(i),j)));
        end
    end
    
    if Na<nPop
    Na=Na+1;
    Nb=Nb-1;
    end

    
    if Nb>=1
    for i=1:Nb
        r2=2*rand-1;
        for j=1:size(X,2) % in j-th 
            X(index2(i),j)= M(l)*Best_pos(j)+RB(index2(i),j)*(r2*(Best_pos(j)-X(index2(i),j))+(1-r2)*(X_centroid(j)-X(index2(i),j)));
        end
    end
    end
    

    for i=1:size(X,1)
        for j=1:dim
            if X(i,j)>ub(j)
                X(i,j)=ub(j);
            end
            if X(i,j)<lb(j)
                X(i,j)=lb(j);
            end
        end
        


        Position.r = [X(i, 1:nVar)];
        Position.psi = [X(i, nVar+1:nVar*2)];
        Position.phi = [X(i, 2*nVar+1:end)];
        Objective_values(1,i)=CostFunction(SphericalToCart(Position, model));
        

        if Objective_values(1,i)<Best_score
            Best_pos=X(i,:);
            Best_score=Objective_values(1,i);
        end
        X_SAO(i,:) = X(i,:);
        Fit_SAO(i) = Objective_values(1,i);
    end
   
    %% NDS
    radius = pdist2(X, X_SAO, 'euclidean');        
    dist_X = squareform(pdist(X));
    r1 = randperm(nPop,nPop);

    for t=1:nPop
        neighbor(t,:) = (dist_X(t,:)<=radius(t,t));
        [~,Idx] = find(neighbor(t,:)==1);                      
        random_Idx_neighbor = randi(size(Idx,2),1,dim);

        for d=1:dim
            X_NDS(t,d) = X(t,d) + rand .*(X(Idx(random_Idx_neighbor(d)),d)- X(r1(t),d));                     
        end
        
        for j=1:dim
            if X_NDS(t,j)>ub(j)
                X_NDS(t,j)=ub(j);
            end
            if X_NDS(t,j)<lb(j)
                X_NDS(t,j)=lb(j);
            end
        end
        Position.r = [X(t, 1:nVar)];
        Position.psi = [X(t, nVar+1:nVar*2)];
        Position.phi = [X(t, 2*nVar+1:end)];
        Fit_NDS(t)=CostFunction(SphericalToCart(Position,model));

        
    end
    % Selection  
    tmp = Fit_SAO < Fit_NDS;                            
    tmp_rep = repmat(tmp',1,dim);

    tmpFit = tmp .* Fit_SAO + (1-tmp) .* Fit_NDS;
    tmpPositions = tmp_rep .* X_SAO + (1-tmp_rep) .* X_NDS;

    % Updating
    tmp = Fit_SAO <= tmpFit;                           
    tmp_rep = repmat(tmp',1,dim);

    Objective_values = tmp .* Fit_SAO + (1-tmp) .* tmpFit;
    X = tmp_rep .* X + (1-tmp_rep) .* tmpPositions;


    %% 
    [~,idx1]=sort(Objective_values);
    second_best=X(idx1(2),:);
    third_best=X(idx1(3),:);
    sum1=0;
    for i=1:N1
        sum1=sum1+X(idx1(i),:);
    end
    half_best_mean=sum1/N1;
    Elite_pool(1,:)=Best_pos;
    Elite_pool(2,:)=second_best;
    Elite_pool(3,:)=third_best;
    Elite_pool(4,:)=half_best_mean;

    Convergence_curve(l)=Best_score;
    std2=std(Objective_values(:));
    sigma(l) = std2/std1;
    if sigma(l) >2
        sigma(l) = 2;
    elseif sigma(l) <0.5
        sigma(l) = 0.5;
    end


    % 
    if mod(l, WaitbarInter) == 0
        waitbar(l / MaxIt, h, ['Completed:' num2str(l / MaxIt * 100) ...
        '%   Completeing...time:', num2str(toc),'/',num2str(toc/(l / MaxIt))])
    end

    DDF(l+1)=0.35+0.25/( 1 + exp(-10*b*(2*l/(MaxIt*sigma(l))- 1)) );
    std1=std2;
    l=l+1;
end
close(h)
end
