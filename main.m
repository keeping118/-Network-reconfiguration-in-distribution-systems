%% 清除内存空间
clc
clear
close all
yalmip('clear')
tic;
%电流kA，电压kV,功率MW
%% 系统参数
mpc = IEEE33;
nb=size(mpc.bus,1);               % 节点数
ns=size(mpc.gen,1);               % 电源节点数
nl=size(mpc.branch,1);            % 支路数
P_load=mpc.bus(:,3)/mpc.baseMVA;  % 有功负荷
Q_load=mpc.bus(:,4)/mpc.baseMVA;  % 无功负荷
r_ij=mpc.branch(:,3);             % 线路电阻
x_ij=mpc.branch(:,4);             % 线路电抗
% 最大电压
%Umax=mpc.bus(:,12).^2;
%Umin=mpc.bus(:,13).^2;
Umax=1.07^2;
Umin=0.93^2;
M=max(Umax - Umin);
% 电源功率最大值
[Lia1,gen_node_idx] = ismember(mpc.gen(:,1),mpc.bus(:,1));
P_g_max=zeros(nb,1);
P_g_max(gen_node_idx)=10000;
Q_g_max=zeros(nb,1);
Q_g_max(gen_node_idx)=10000;
% 电源功率最小值
P_g_min=zeros(nb,1);
%P_g_min(mpc.gen(:,1))=mpc.gen(:,10)/mpc.baseMVA;
Q_g_min=zeros(nb,1);
%Q_g_min(mpc.gen(:,1))=mpc.gen(:,5)/mpc.baseMVA;
% 支路功率最大值
Sij_max=1500/mpc.baseMVA;
Iij_max = mpc.branch(:,6)*1000 ./ (mpc.baseMVA*1e6./(mpc.bus(1,10)*1e3))+1;
% 流入节点的支路
branch_to_node=zeros(nb,nl);
% 流出节点的支路
branch_from_node=zeros(nb,nl);
%核心逻辑是之前节点和行数是对应的，现在不对应了
for k=1:nl
    branch_to_node(find(mpc.bus(:,1)==mpc.branch(k,2)),k)=1;
    branch_from_node(find(mpc.bus(:,1)==mpc.branch(k,1)),k)=1;
end

%寻找branch节点在bus中对应的行和列
[Lia2,from_node] = ismember(mpc.branch(:,1),mpc.bus(:,1));
[Lia3,to_node] = ismember(mpc.branch(:,2),mpc.bus(:,1));

%虚拟潮流设置：每个非根节点吸收1个单位，由根节点提供所有
Pload_vir = ones(nb, 1);
Pload_vir(gen_node_idx) = 0;
Pg_vir_max = zeros(nb, 1);
Pg_vir_max(gen_node_idx) = nb;
%% 优化变量
z_ij=binvar(nl,1);    % 支路开断情况
U_i=sdpvar(nb,1);     % 电压的平方
L_ij=sdpvar(nl,1);    % 电流的平方
P_ij=sdpvar(nl,1);    % 线路有功功率
Q_ij=sdpvar(nl,1);    % 线路无功功率
P_g=sdpvar(nb,1);     % 电源有功出力
Q_g=sdpvar(nb,1);     % 电源无功功率
Pij_vir=sdpvar(nl,1); % 虚拟线路有功功率
Pg_vir=sdpvar(nb,1);  % 虚拟发电机
u = sdpvar(1,1);     % 电压最小值辅助变量
i = sdpvar(1,1);      % 电流最大值辅助变量
y = binvar(nb, ns, 'full'); % 引入新的二元变量 y(i,k)，表示节点i是否由电源k供电
%% 目标函数
objective = sum(L_ij.*r_ij);%网损最小
%objective = i;%最小化最大电流
%objective = -u;%最大化最小电压
%% 约束条件
Constraints = [];
Constraints = [Constraints,u <= U_i]; 
Constraints = [Constraints,i >= L_ij./(Iij_max.^2)]
%% 1.潮流约束
%设置电源节点电压
Constraints = [Constraints,U_i(gen_node_idx) == mpc.gen(:,6).^2]; 

m_ij=(1-z_ij)*M;
%功率平衡约束
Constraints = [Constraints,  P_g-P_load+branch_to_node*P_ij-branch_to_node*(L_ij.*r_ij)-branch_from_node*P_ij == 0];
Constraints = [Constraints,  Q_g-Q_load+branch_to_node*Q_ij-branch_to_node*(L_ij.*x_ij)-branch_from_node*Q_ij == 0];
%电压降约束
Constraints = [Constraints,U_i(from_node)-U_i(to_node) <= m_ij + 2*r_ij.*P_ij + 2*x_ij.*Q_ij - ((r_ij.^2 + x_ij.^2)).*L_ij];
Constraints = [Constraints,U_i(from_node)-U_i(to_node) >= -m_ij + 2*r_ij.*P_ij + 2*x_ij.*Q_ij - ((r_ij.^2 + x_ij.^2)).*L_ij]; 
%功率电压二次等式约束
for k = 1:nl
    i_idx = find(mpc.bus(:,1)==mpc.branch(k,1));  % 取出起点母线索引
    Constraints = [Constraints, P_ij(k)^2 + Q_ij(k)^2  == U_i(i_idx)*L_ij(k) ];
end
Constraints = [Constraints, P_ij.^2 + Q_ij.^2 <= Sij_max'.^2.*z_ij];
Constraints = [Constraints,0 <= L_ij <= Iij_max.^2.*z_ij];
Constraints = [Constraints, Umin <= U_i <= Umax];
%% 2.拓扑约束
%约束1: 确保最终拓扑是辐射状的（总支路数 = 节点数 - 电源数）
Constraints = [Constraints , sum(z_ij) == nb-ns];
%约束2：虚拟潮流平衡约束，避免断路，确保负荷和电源相连
Constraints = [Constraints,  Pg_vir-Pload_vir+branch_to_node*Pij_vir-branch_from_node*Pij_vir == 0];
Constraints = [Constraints,  sum(Pg_vir) == nb-ns];
%虚拟潮流只能在闭合的支路上流动
M_flow = nb; % 一个足够大的数，大于等于总需求
Constraints = [Constraints, -M_flow * z_ij <= Pij_vir <= M_flow * z_ij];
Constraints = [Constraints, 0 <= Pg_vir <= Pg_vir_max];
%约束3：网络分区法，确保每个节点仅由一个电源供电
%每个节点必须且只能由一个电源供电
Constraints = [Constraints, sum(y, 2) == 1];
%约束：每个电源节点必须属于自己的分区
%gen_node_idx 是电源节点在bus矩阵中的行号索引
for k = 1:ns
    Constraints = [Constraints, y(gen_node_idx(k), k) == 1];
end
%约束：如果一条支路是闭合的，则其两端节点必须属于同一个电源分区
%from_node 和 to_node 分别是每条支路的起始节点和结束节点在bus矩阵中的行号索引
for k = 1:ns
    Constraints = [Constraints, -(1 - z_ij) <= y(from_node, k) - y(to_node, k) <= (1 - z_ij)];
end

%% 3.注入功率约束
Constraints = [Constraints, P_g_min <= P_g <= P_g_max];
Constraints = [Constraints, Q_g_min <= Q_g <= Q_g_max];
%% 设求解器
ops=sdpsettings('verbose', 2, 'solver', 'gurobi');
sol=optimize(Constraints,objective,ops);

%% 分析错误标志
if sol.problem == 0
    disp('求解成功');
else
    disp('运行出错');
    yalmiperror(sol.problem)
end
%% 结果
P_ij=value(P_ij)*mpc.baseMVA;
Q_ij=value(Q_ij)*mpc.baseMVA;
P_g=value(P_g)*mpc.baseMVA;
Q_g=value(Q_g)*mpc.baseMVA;
z_ij=value(z_ij);
U_i=value(U_i);
L_ij=value(L_ij);
i=value(i);
y=value(y);
P_loss=L_ij.*r_ij*1000*mpc.baseMVA;%(先从标幺值得到损耗，再*mpc.baseMVA转换到MW，最后*1000转换到kW)
Pg_vir=value(Pg_vir);
Pij_vir=value(Pij_vir);

result=runpf('IEEE33');
P_loss0=(result.branch(:,14)+result.branch(:,16))*1000;
V0=result.bus(:,8);
disp('******************重构前******************')
disp(['开断支路为：',num2str(find(~mpc.branch(:,11)'))])
disp(['系统网损为：',num2str(sum(P_loss0)),'kW'])
disp(['最低节点电压为：',num2str(min(V0))])
%%%%%%%差一个最大负载率
disp('******************重构后******************')
disp(['开断支路为：',num2str(find(~round(z_ij))')])
disp(['系统网损为：',num2str(value(objective)*1000*mpc.baseMVA),'kW'])
disp(['最低节点电压为：',num2str(min(sqrt(U_i)))])
disp(['最大负载率：',num2str(max(L_ij./(Iij_max.^2)))])

figure
plot(V0,'k','linewidth',1)
hold on
plot(sqrt(U_i),'k--','linewidth',1)
title('节点电压对比')
xlabel('节点')
ylabel('电压幅值/pu')
legend('重构前','重构后')
figure
plot(P_loss0,'k--','linewidth',1)
hold on
plot(P_loss,'k','linewidth',1)
title('支路有功损耗对比')
xlabel('节点')
ylabel('支路有功损耗/kW')
legend('重构前','重构后')