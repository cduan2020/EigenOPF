function [Vps1,success,maxReal,improve,best_id]=OptimizeRatio(ps,Vps0,dispratio,load_level)

% Optimize the system damping ratio by adjusting generation and demand.
% It calls the subroutine TX_opt with differnt initial stepsizes and choose the best solution
% ps: the power system data structure
% Vps0: the initial power flow solution
% disparatio: ratio of controllable load
% load_level: load factor
% Vps1: the optimized power flow solution
% success: indicator of the success of solving the problem
% maxReal: the damping ratios
% best_id: the id of the best initial stepsizes

num_job=12;
success_vec=cell(num_job,1);
maxReal_vec=cell(num_job,1);
improve_vec=zeros(num_job,1);
Vps1_vec=cell(num_job,1);
% use parfor to accelerate when possible
for ii=1:num_job
    [success_vec{ii},maxReal_vec{ii},improve_vec(ii),Vps1_vec{ii}] = TX_opt(ps,Vps0,dispratio,0.083333333333333/2*ii,load_level);
end
[~,best_id]=max(improve_vec);
success=success_vec{best_id};
maxReal=maxReal_vec{best_id};
improve=improve_vec(best_id);
Vps1=Vps1_vec{best_id};
maxReal=maxReal(end);
end
