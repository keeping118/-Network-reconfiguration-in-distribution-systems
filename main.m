%% ����ڴ�ռ�
clc
clear
close all
yalmip('clear')
tic;
%����kA����ѹkV,����MW
%% ϵͳ����
mpc = IEEE33;
nb=size(mpc.bus,1);               % �ڵ���
ns=size(mpc.gen,1);               % ��Դ�ڵ���
nl=size(mpc.branch,1);            % ֧·��
P_load=mpc.bus(:,3)/mpc.baseMVA;  % �й�����
Q_load=mpc.bus(:,4)/mpc.baseMVA;  % �޹�����
r_ij=mpc.branch(:,3);             % ��·����
x_ij=mpc.branch(:,4);             % ��·�翹
% ����ѹ
%Umax=mpc.bus(:,12).^2;
%Umin=mpc.bus(:,13).^2;
Umax=1.07^2;
Umin=0.93^2;
M=max(Umax - Umin);
% ��Դ�������ֵ
[Lia1,gen_node_idx] = ismember(mpc.gen(:,1),mpc.bus(:,1));
P_g_max=zeros(nb,1);
P_g_max(gen_node_idx)=10000;
Q_g_max=zeros(nb,1);
Q_g_max(gen_node_idx)=10000;
% ��Դ������Сֵ
P_g_min=zeros(nb,1);
%P_g_min(mpc.gen(:,1))=mpc.gen(:,10)/mpc.baseMVA;
Q_g_min=zeros(nb,1);
%Q_g_min(mpc.gen(:,1))=mpc.gen(:,5)/mpc.baseMVA;
% ֧·�������ֵ
Sij_max=1500/mpc.baseMVA;
Iij_max = mpc.branch(:,6)*1000 ./ (mpc.baseMVA*1e6./(mpc.bus(1,10)*1e3))+1;
% ����ڵ��֧·
branch_to_node=zeros(nb,nl);
% �����ڵ��֧·
branch_from_node=zeros(nb,nl);
%�����߼���֮ǰ�ڵ�������Ƕ�Ӧ�ģ����ڲ���Ӧ��
for k=1:nl
    branch_to_node(find(mpc.bus(:,1)==mpc.branch(k,2)),k)=1;
    branch_from_node(find(mpc.bus(:,1)==mpc.branch(k,1)),k)=1;
end

%Ѱ��branch�ڵ���bus�ж�Ӧ���к���
[Lia2,from_node] = ismember(mpc.branch(:,1),mpc.bus(:,1));
[Lia3,to_node] = ismember(mpc.branch(:,2),mpc.bus(:,1));

%���⳱�����ã�ÿ���Ǹ��ڵ�����1����λ���ɸ��ڵ��ṩ����
Pload_vir = ones(nb, 1);
Pload_vir(gen_node_idx) = 0;
Pg_vir_max = zeros(nb, 1);
Pg_vir_max(gen_node_idx) = nb;
%% �Ż�����
z_ij=binvar(nl,1);    % ֧·�������
U_i=sdpvar(nb,1);     % ��ѹ��ƽ��
L_ij=sdpvar(nl,1);    % ������ƽ��
P_ij=sdpvar(nl,1);    % ��·�й�����
Q_ij=sdpvar(nl,1);    % ��·�޹�����
P_g=sdpvar(nb,1);     % ��Դ�й�����
Q_g=sdpvar(nb,1);     % ��Դ�޹�����
Pij_vir=sdpvar(nl,1); % ������·�й�����
Pg_vir=sdpvar(nb,1);  % ���ⷢ���
u = sdpvar(1,1);     % ��ѹ��Сֵ��������
i = sdpvar(1,1);      % �������ֵ��������
y = binvar(nb, ns, 'full'); % �����µĶ�Ԫ���� y(i,k)����ʾ�ڵ�i�Ƿ��ɵ�Դk����
%% Ŀ�꺯��
objective = sum(L_ij.*r_ij);%������С
%objective = i;%��С��������
%objective = -u;%�����С��ѹ
%% Լ������
Constraints = [];
Constraints = [Constraints,u <= U_i]; 
Constraints = [Constraints,i >= L_ij./(Iij_max.^2)]
%% 1.����Լ��
%���õ�Դ�ڵ��ѹ
Constraints = [Constraints,U_i(gen_node_idx) == mpc.gen(:,6).^2]; 

m_ij=(1-z_ij)*M;
%����ƽ��Լ��
Constraints = [Constraints,  P_g-P_load+branch_to_node*P_ij-branch_to_node*(L_ij.*r_ij)-branch_from_node*P_ij == 0];
Constraints = [Constraints,  Q_g-Q_load+branch_to_node*Q_ij-branch_to_node*(L_ij.*x_ij)-branch_from_node*Q_ij == 0];
%��ѹ��Լ��
Constraints = [Constraints,U_i(from_node)-U_i(to_node) <= m_ij + 2*r_ij.*P_ij + 2*x_ij.*Q_ij - ((r_ij.^2 + x_ij.^2)).*L_ij];
Constraints = [Constraints,U_i(from_node)-U_i(to_node) >= -m_ij + 2*r_ij.*P_ij + 2*x_ij.*Q_ij - ((r_ij.^2 + x_ij.^2)).*L_ij]; 
%���ʵ�ѹ���ε�ʽԼ��
for k = 1:nl
    i_idx = find(mpc.bus(:,1)==mpc.branch(k,1));  % ȡ�����ĸ������
    Constraints = [Constraints, P_ij(k)^2 + Q_ij(k)^2  == U_i(i_idx)*L_ij(k) ];
end
Constraints = [Constraints, P_ij.^2 + Q_ij.^2 <= Sij_max'.^2.*z_ij];
Constraints = [Constraints,0 <= L_ij <= Iij_max.^2.*z_ij];
Constraints = [Constraints, Umin <= U_i <= Umax];
%% 2.����Լ��
%Լ��1: ȷ�����������Ƿ���״�ģ���֧·�� = �ڵ��� - ��Դ����
Constraints = [Constraints , sum(z_ij) == nb-ns];
%Լ��2�����⳱��ƽ��Լ���������·��ȷ�����ɺ͵�Դ����
Constraints = [Constraints,  Pg_vir-Pload_vir+branch_to_node*Pij_vir-branch_from_node*Pij_vir == 0];
Constraints = [Constraints,  sum(Pg_vir) == nb-ns];
%���⳱��ֻ���ڱպϵ�֧·������
M_flow = nb; % һ���㹻����������ڵ���������
Constraints = [Constraints, -M_flow * z_ij <= Pij_vir <= M_flow * z_ij];
Constraints = [Constraints, 0 <= Pg_vir <= Pg_vir_max];
%Լ��3�������������ȷ��ÿ���ڵ����һ����Դ����
%ÿ���ڵ������ֻ����һ����Դ����
Constraints = [Constraints, sum(y, 2) == 1];
%Լ����ÿ����Դ�ڵ���������Լ��ķ���
%gen_node_idx �ǵ�Դ�ڵ���bus�����е��к�����
for k = 1:ns
    Constraints = [Constraints, y(gen_node_idx(k), k) == 1];
end
%Լ�������һ��֧·�Ǳպϵģ��������˽ڵ��������ͬһ����Դ����
%from_node �� to_node �ֱ���ÿ��֧·����ʼ�ڵ�ͽ����ڵ���bus�����е��к�����
for k = 1:ns
    Constraints = [Constraints, -(1 - z_ij) <= y(from_node, k) - y(to_node, k) <= (1 - z_ij)];
end

%% 3.ע�빦��Լ��
Constraints = [Constraints, P_g_min <= P_g <= P_g_max];
Constraints = [Constraints, Q_g_min <= Q_g <= Q_g_max];
%% �������
ops=sdpsettings('verbose', 2, 'solver', 'gurobi');
sol=optimize(Constraints,objective,ops);

%% ���������־
if sol.problem == 0
    disp('���ɹ�');
else
    disp('���г���');
    yalmiperror(sol.problem)
end
%% ���
P_ij=value(P_ij)*mpc.baseMVA;
Q_ij=value(Q_ij)*mpc.baseMVA;
P_g=value(P_g)*mpc.baseMVA;
Q_g=value(Q_g)*mpc.baseMVA;
z_ij=value(z_ij);
U_i=value(U_i);
L_ij=value(L_ij);
i=value(i);
y=value(y);
P_loss=L_ij.*r_ij*1000*mpc.baseMVA;%(�ȴӱ���ֵ�õ���ģ���*mpc.baseMVAת����MW�����*1000ת����kW)
Pg_vir=value(Pg_vir);
Pij_vir=value(Pij_vir);

result=runpf('IEEE33');
P_loss0=(result.branch(:,14)+result.branch(:,16))*1000;
V0=result.bus(:,8);
disp('******************�ع�ǰ******************')
disp(['����֧·Ϊ��',num2str(find(~mpc.branch(:,11)'))])
disp(['ϵͳ����Ϊ��',num2str(sum(P_loss0)),'kW'])
disp(['��ͽڵ��ѹΪ��',num2str(min(V0))])
%%%%%%%��һ���������
disp('******************�ع���******************')
disp(['����֧·Ϊ��',num2str(find(~round(z_ij))')])
disp(['ϵͳ����Ϊ��',num2str(value(objective)*1000*mpc.baseMVA),'kW'])
disp(['��ͽڵ��ѹΪ��',num2str(min(sqrt(U_i)))])
disp(['������ʣ�',num2str(max(L_ij./(Iij_max.^2)))])

figure
plot(V0,'k','linewidth',1)
hold on
plot(sqrt(U_i),'k--','linewidth',1)
title('�ڵ��ѹ�Ա�')
xlabel('�ڵ�')
ylabel('��ѹ��ֵ/pu')
legend('�ع�ǰ','�ع���')
figure
plot(P_loss0,'k--','linewidth',1)
hold on
plot(P_loss,'k','linewidth',1)
title('֧·�й���ĶԱ�')
xlabel('�ڵ�')
ylabel('֧·�й����/kW')
legend('�ع�ǰ','�ع���')