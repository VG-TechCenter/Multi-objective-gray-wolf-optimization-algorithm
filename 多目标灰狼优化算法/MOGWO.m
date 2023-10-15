%___________________________________________________________________%
%  Multi-Objective Grey Wolf Optimizer (MOGWO)                      %
%  Source codes demo version 1.0                                    %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%    S. Mirjalili, S. Saremi, S. M. Mirjalili, L. Coelho,           %
%    Multi-objective grey wolf optimizer: A novel algorithm for     %
%    multi-criterion optimization, Expert Systems with Applications,%
%    in press, DOI: http://dx.doi.org/10.1016/j.eswa.2015.10.039    %       %
%                                                                   %
%___________________________________________________________________%

% I acknowledge that this version of MOGWO has been written using
% a large portion of the following code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for                                                  %
%                                                                   %
%  Multi-Objective Particle Swarm Optimization (MOPSO)              %
%  Version 1.0 - Feb. 2011                                          %
%                                                                   %
%  According to:                                                    %
%  Carlos A. Coello Coello et al.,                                  %
%  "Handling Multiple Objectives with Particle Swarm Optimization," %
%  IEEE Transactions on Evolutionary Computation, Vol. 8, No. 3,    %
%  pp. 256-279, June 2004.                                          %
%                                                                   %
%  Developed Using MATLAB R2009b (Version 7.9)                      %
%                                                                   %
%  Programmed By: S. Mostapha Kalami Heris                          %
%                                                                   %
%         e-Mail: sm.kalami@gmail.com                               %
%                 kalami@ee.kntu.ac.ir                              %
%                                                                   %
%       Homepage: http://www.kalami.ir                              %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

drawing_flag = 1;

TestProblem='UF1';
nVar=10;

fobj = cec09(TestProblem);

xrange = xboundary(TestProblem, nVar);

% Lower bound and upper bound
lb=xrange(:,1)';
ub=xrange(:,2)';

VarSize=[1 nVar];

GreyWolves_num=100;
MaxIt=1000;  % Maximum Number of Iterations
Archive_size=100;   % Repository Size

alpha=0.1;  % Grid Inflation Parameter
nGrid=10;   % Number of Grids per each Dimension
beta=4; %=4;    % Leader Selection Pressure Parameter
gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

% Initialization

GreyWolves=CreateEmptyParticle(GreyWolves_num);


for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=zeros(1,nVar);
    for j=1:nVar
        GreyWolves(i).Position(1,j)=unifrnd(lb(j),ub(j),1);
    end
    GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination(GreyWolves);

Archive=GetNonDominatedParticles(GreyWolves);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

% MOGWO main loop

for it=1:MaxIt
    a=2-it*((2)/MaxIt);
    for i=1:GreyWolves_num
        
        clear rep2
        clear rep3
        
        % Choose the alpha, beta, and delta grey wolves
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        % If there are less than three solutions in the least crowded
        % hypercube, the second least crowded hypercube is also found
        % to choose other leaders from.
        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Beta=SelectLeader(rep2,beta);
        end
        
        % This scenario is the same if the second least crowded hypercube
        % has one solution, so the delta leader should be chosen from the
        % third least crowded hypercube.
        if size(Archive,1)>2
            counter=0;
            for newi=1:size(rep2,1)
                if sum(Beta.Position~=rep2(newi).Position)~=0
                    counter=counter+1;
                    rep3(counter,1)=rep2(newi);
                end
            end
            Alpha=SelectLeader(rep3,beta);
        end
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Delta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nVar)-a;
        % Eq.(3.8) in the paper
        X1=Delta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Beta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.9) in the paper
        X2=Beta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Alpha.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.10) in the paper
        X3=Alpha.Position-A.*abs(D);
        
        % Eq.(3.11) in the paper
        GreyWolves(i).Position=(X1+X2+X3)./3;
        
        % Boundary checking
        GreyWolves(i).Position=min(max(GreyWolves(i).Position,lb),ub);
        
        GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
    end
    
    GreyWolves=DetermineDomination(GreyWolves);
    non_dominated_wolves=GetNonDominatedParticles(GreyWolves);
    
    Archive=[Archive
        non_dominated_wolves];
    
    Archive=DetermineDomination(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
    
    disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results
    
    % Results
    
    costs=GetCosts(GreyWolves);
    Archive_costs=GetCosts(Archive);
    
    if drawing_flag==1
        hold off
        plot(costs(1,:),costs(2,:),'k.');
        hold on
        plot(Archive_costs(1,:),Archive_costs(2,:),'rd');
        legend('Grey wolves','Non-dominated solutions');
        drawnow
    end
    
end


