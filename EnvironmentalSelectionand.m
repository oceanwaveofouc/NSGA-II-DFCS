function [Population,Fitness] = EnvironmentalSelectionDFCS(Population,N,Zmin)
%

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group

CV = sum(max(0,Population.cons),2);
if sum(CV==0) > N
    %% Selection among feasible solutions
    Population = Population(CV==0);
%     % Non-dominated sorting
%     [FrontNo,MaxFNo] = NDSort(Population.objs,N);
%     Next = false(1,length(FrontNo));
%     Next(FrontNo<=MaxFNo) = true;
    Next(1:size(Population,2)) = true;
    % Select the solutions in the last front
    Delete = LastSelection(Population(Next).objs,sum(Next)-N,Zmin);
    Temp = find(Next);
    Next(Temp(Delete)) = false;
    Population = Population(Next);
else
    %% Selection including infeasible solutions
    [~,rank]   = sort(CV);
    Population = Population(rank(1:N));
end

end


function Delete = LastSelection(PopObj,K,Zmin)
% Select part of the solutions in the last front
    [N,M]  = size(PopObj);
    PopObj = (PopObj-repmat(Zmin,N,1))./(repmat(max(PopObj),N,1)-repmat(Zmin,N,1));
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
      Cosine   = 1 - pdist2(PopObj,PopObj,'cosine');
      Cosine   = Cosine.*(1-eye(size(PopObj,1)));
      

    %% Environmental selection
    Delete  = false(1,N);
    % Select K solutions one by one
    while sum(Delete) < K
        [Jmin_row,Jmin_column] = find(Cosine==max(max(Cosine)));
        j = randi(length(Jmin_row));
        Temp_1 = Jmin_row(j);
        Temp_2 = Jmin_column(j);
         

        if  (   sum(PopObj(Temp_1,:))>   sum(PopObj(Temp_2,:))  ) ||  (  sum(PopObj(Temp_1,:))==sum(PopObj(Temp_2,:)) && rand<0.5)
            Delete(Temp_1) = true;
            Cosine(:,Temp_1)=0;
            Cosine(Temp_1,:)=0;
        else
            Delete(Temp_2) = true;
            Cosine(:,Temp_2)=0;
            Cosine(Temp_2,:)=0;
        end
    end
end