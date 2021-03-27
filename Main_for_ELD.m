clear;
clc;

%% Algrithm controller setting
ChaoticTrackLength_no=500;%Set length of Chaotic Track
Max_iteration=300;%Set max iteration
rate=10;%Set helix controller

Function_name='ELD3';%Input function name which is to be optimized
[fobj,PD,ELD]=Get_ELD_details(Function_name);%Get problem details from Function

ub=ELD(:,2);%get up bonds
lb=ELD(:,1);%get low bonds
dim=size(ELD,1);%get number of machines

%% Chaotic system parameters
a=10;
b=28;
c=6.2;

X0=1;
Y0=1;
Z0=1;

%% Main Function area
%Initialize readers & writers
best_fit=Inf;
current_bestfit=Inf;
current_best_Vector=[];
best_Vector=[];

fitness=ones(ChaoticTrackLength_no,1)*Inf;

%Initialize iteration logs
log_current_best_fit=zeros(Max_iteration,1);%
log_best_fit=zeros(Max_iteration,1);%

%Produce Chaotic swarm
[Track]=DLCS(X0,Y0,Z0,a,b,c,ChaoticTrackLength_no,dim)';
agents=zeros(size(Track,1),size(Track,2));
%Map Chaotic swarm into solution space
for D=1:1:dim
    agents(:,D)=((Track(:,D)-min(Track(:,D)))*(ub(D)-lb(D)))/(max(Track(:,D))-min(Track(:,D)))+lb(D);
end

%Caculate initialized fitness
for fit_index=1:1:ChaoticTrackLength_no
    fitness(fit_index)=fobj(ELD,agents(fit_index,:)',PD)';
end

%Mark initialized best solution
current_bestfit=min(fitness);
lead_index=find(fitness==current_bestfit);
current_best_Vector=agents(lead_index,:);
best_Vector=current_best_Vector;

%Main iteration area
for iter=1:1:Max_iteration
    %%  Updating swarm's positions
    %Self-adapting baseline ocsilliation swarm divide strategy(chapter 3.2)
    baseline=mean(fitness);
    core_swarm=agents(fitness<baseline,:);
    chaosdrift_swarm=agents(fitness>=baseline,:);
    
    %Chaos Drift Swarm evolution(chapter 3.3) 
    if size(chaosdrift_swarm,1)>0
        order=floor(size(chaosdrift_swarm,1)*((iter/Max_iteration)^2));
        core_swarm=[core_swarm;chaosdrift_swarm(1:1:order,:)];
        chaosdrift_swarm(1:1:order,:)=[];
        chaosdrift_swarm=chaosdrift(chaosdrift_swarm,ub,lb);
    end
    
    %Core swarm convergence in spiral equiation(chapter 3.4)
    for agent=1:1:size(core_swarm,1)
        for D=1:1:dim
            l=unifrnd(-rate,rate);
            V_D=abs(core_swarm(agent,D)-best_Vector(D))*exp(l)*cos(2*pi*l);
            core_swarm(agent,D)=best_Vector(D)+V_D;
            if core_swarm(agent,D)>ub(D) || core_swarm(agent,D)<lb(D)
                core_swarm(agent,D)=(ub(D)-lb(D))*rand+lb(D);
            end
        end        
    end
    %Completed swarm
    agents=[core_swarm;chaosdrift_swarm];
    
% % %    Temp Observe
%     figure(2)
% %     subplot(1,2,1)
% %     plot(current_best_Vector);
% %     axis([0,dim,lb,ub]);
% %     subplot(1,2,2)
%     scatter(core_swarm(:,1),core_swarm(:,2),'.','red');
%     hold on
%     scatter(chaosdrift_swarm(:,1),chaosdrift_swarm(:,2),'.','blue');
%     scatter(current_best_Vector(1),current_best_Vector(2),'*','black');
%     hold off
% %     axis([lb,ub,lb,ub]);
%     pause(0.01)

    %Mark new best solution
    for fit_index=1:1:ChaoticTrackLength_no
        fitness(fit_index)=fobj(ELD,agents(fit_index,:)',PD)';
    end
    current_bestfit=min(fitness);
    lead_index=find(fitness==current_bestfit);
    current_best_Vector=agents(lead_index,:);
    if current_bestfit<best_fit
        best_fit=current_bestfit;
        best_Vector=current_best_Vector(1,:);
    end
    disp([iter,current_bestfit-(PD-sum(best_Vector))^2,best_fit-(PD-sum(best_Vector))^2]);

    %log new solution
    log_current_best_fit(iter)=current_bestfit;
    log_best_fit(iter)=best_fit;
end

%% Plot convergence curves
figure(1);
semilogy(log_current_best_fit);
hold on
semilogy(log_best_fit);
hold off

disp('Best X£º');
disp(best_Vector);
disp(['Best fitness£º',num2str(log_best_fit(end)-(PD-sum(best_Vector))^2),'    ',num2str(sum(best_Vector))]);