clear all
clc
close all
%% 
SearchAgents_no=100;       
Max_iteration=300;     
Function_name=15; 
dim=10;                 
lb=-100;              
ub=100;               
fobj = @(x) cec17_func(x',Function_name);

%% GDSAO
Max_test=20;
for i=1:Max_test
    [Best_pos(i,:),Best_score(i),SAO_curve(i,:)]=GDSAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %开始优化
end

disp('-------------------------------------------------')
display(['best GDSAO: ', num2str(min(Best_score))]);
display(['mean GDSAO: ', num2str(mean(Best_score))]);
display(['worst GDSAO: ', num2str(max(Best_score))]);
display(['std GDSAO: ', num2str(std(Best_score))]);


