%% 220330124-舒晟超-电力系统分析实验1-3
clear;

n=13;           % 节点数;
nl=15;          % 支路数;
isb=1;          % 平衡母线节点号;
pr=0.00001;     % 误差精度;

% Bl矩阵：
%   [支路首端号；支路末端号；支路串联阻抗；支路对地导纳;变压器的变比；支路首端处于K侧为1，1侧为0]
Bl =[1              2       0.04642+0.09316i   0           1          0;
     2              3       0.0298+0.38375i    0           1.02381    1;
     1              4       0.08620+0.24289i   0           1          0;
     4              5       0.00911+0.37949i   0           1.02381    1;
     1              6       0.05924+0.13263i   0           1          0;
     4              6       0.05872+0.08967i   0           1          0;
     6              7       0.0072+0.18375i    0           1.02381    1;
     6              8       0.02448+0.06645i   0           1          0;
     8              10      0.00911+0.24549i   0           1          1;
     8              9       0.05612+0.08306i   0           1          0;
     10             11      0.03243+0.29142i   0           1          1;
     10             13      0.02371+0.06645i   0           1          0;
     12             13      0.00911+0.24549i   0           1.02381    1;
     2              11      0.01243+0.29142i   0           1          0;
     9              13      0.02323+0.32142i   0           1.02381    1];

% Bn矩阵：
%   [该节点发电机功率；该节点负荷功率；节点电压初始值;PV节点电压V的给定值；节点所接的无功补偿设备的容量;节点分类标号：0为平衡节点；1为PQ节点；2为PV节点]
Bn=[0               0                1.1                1.1                     0                       0;
    0               0                1                  0                       0                       1;
    0               0.343+0.21256i   1                  0                       0                       1;
    0               0                1                  0                       0                       1;
    0               0.204+0.12638i   1                  0                       0                       1;
    0               0                1                  0                       0                       1;
    0               0.306+0.18962i   1                  0                       0                       1;
    0               0.103+0.13267i   1                  0                       0                       1;
    0.5             0.372+0.18962i   1.1                1.1                     0                       2;
    0               0.343+0.21256i   1                  0                       0                       1;
    0               0                1                  0                       0                       1;
    0.1             0.204+0.12638i   1                  0                       0                       1;
    0               0                1                  0                       0                       1];

%% 求导纳矩阵
Y=zeros(n);
e=zeros(1,n);
f=zeros(1,n);
V=zeros(1,n);
sida=zeros(1,n);
S1=zeros(nl);

for i=1:nl                                              % 支路数
    if Bl(i,6)==0                                       % 左节点处于1侧
       p=Bl(i,1);q=Bl(i,2);
    else                                                % 左节点处于K侧
       p=Bl(i,2);q=Bl(i,1);
    end
    Y(p,q)=Y(p,q)-1./(Bl(i,3)*Bl(i,5));                 % 非对角元
    Y(q,p)=Y(p,q);                                      % 非对角元
    Y(q,q)=Y(q,q)+1./(Bl(i,3)*Bl(i,5)^2)+Bl(i,4)./2;    % 对角元K侧
    Y(p,p)=Y(p,p)+1./Bl(i,3)+Bl(i,4)./2;                % 对角元1侧 
end

disp('导纳矩阵 Y=');
disp(Y)

% 分解导纳阵的实部和虚部
G=real(Y);
B=imag(Y); 

%% 电压初值(直角坐标)
for i=1:n                                               % 给定各节点初始电压的实部和虚部    
    e(i)=real(Bn(i,3));
    f(i)=imag(Bn(i,3));
    V(i)=Bn(i,4);                                       % PV节点电压给定模值 
end



%% 功率初值
for i=1:n                                               % 给定各节点注入功率    
    S(i)=Bn(i,1)-Bn(i,2);                               % i节点注入功率SG-SL  
    B(i,i)=B(i,i)+Bn(i,5);                              % i节点无功补偿 
end
% 分解出各节点注入的有功和无功功率
P=real(S);
Q=imag(S);

%% Newton--Rapshon法迭代
ICT1=0;
IT2=1;
N0=2*n;
N=N0+1;
a=0;    %迭代次数ICT1,a;不满足收敛要求的节点数IT2;N0=2*n:雅可比矩阵的阶数;N=N0+1:扩展列

while IT2~=0                    
    % 迭代初始化
    IT2=0;
    a=a+1;
    for i=1:n
        %% 除平衡节点外其它节点的功率计算
        if i~=isb                   % 非平衡节点    
            C(i)=0;
            D(i)=0;
            for j1=1:n
                C(i)=C(i)+G(i,j1)*e(j1)-B(i,j1)*f(j1);          % Σ(Gij*ej-Bij*fj)
                D(i)=D(i)+G(i,j1)*f(j1)+B(i,j1)*e(j1);          % Σ(Gij*fj+Bij*ej)
            end
            P1=e(i)*(C(i))+f(i)*(D(i));                         % 节点功率P计算 eiΣ(Gij*ej-Bij*fj)+fiΣ(Gij*fj+Bij*ej)
            Q1=f(i)*(C(i))-e(i)*(D(i));                         % 节点功率Q计算 fiΣ(Gij*ej-Bij*fj)-eiΣ(Gij*fj+Bij*ej)    
            V2=e(i)^2 + f(i)^2;                                 % 电压模平方
            %%  针对非PV节点来求取功率差及Jacobi矩阵元素
            if Bn(i,6)~=2               % 非PV节点    
                DP=P(i)-P1;             % 节点有功功率差    
                DQ=Q(i)-Q1;             % 节点无功功率差  
                for j1=1:n
                    if j1~=isb & j1~=i                          % 非平衡节点 & 非对角元    
                        X1=-(G(i,j1)*e(i)+B(i,j1)*f(i));        % X1=dPi/dej    
                        X2=B(i,j1)*e(i)-G(i,j1)*f(i);           % X2=dPi/dfj    
                        X3=X2;                                  % X3=dQi/dej   X3=X2
                        X4=-X1;                                 % X4=dQi/dfj   X4=-X1
                        p=2*i-1;
                        q=2*j1-1;
                        J(p,q)=X3;J(p,N)=DQ;m=p+1;              % X3=dQ/de  J(p,N)=DQ节点无功功率差
                        J(m,q)=X1;J(m,N)=DP;q=q+1;              % X1=dP/de  J(m,N)=DP节点有功功率差
                        J(p,q)=X4;J(m,q)=X2;                    % X4=dQ/df  X2=dp/df
                    elseif j1==i & j1~=isb                      % 非平衡节点 & 对角元
                        X1=-C(i)-G(i,i)*e(i)-B(i,i)*f(i);       % X1=dPi/dei
                        X2=-D(i)+B(i,i)*e(i)-G(i,i)*f(i);       % X2=dPi/dfi 
                        X3=D(i)+B(i,i)*e(i)-G(i,i)*f(i);        % X3=dQi/dei
                        X4=-C(i)+G(i,i)*e(i)+B(i,i)*f(i);       % X4=dQi/dfi
                        p=2*i-1;q=2*j1-1;J(p,q)=X3;J(p,N)=DQ;   %扩展列△Q
                        m=p+1;
                        J(m,q)=X1;q=q+1;J(p,q)=X4;J(m,N)=DP;    %扩展列△P
                        J(m,q)=X2;
                    end
                end
            else
            %% 针对PV节点来求取Jacobi矩阵的元素
                DP=P(i)-P1;                 % PV节点有功误差
                DV=V(i)^2-V2;               % PV节点电压误差    
                for j1=1:n
                    if j1~=isb & j1~=i                          % 非平衡节点 & 非对角元    
                       X1=-(G(i,j1)*e(i)+B(i,j1)*f(i));         % X1=dPi/dej
                       X2=B(i,j1)*e(i)-G(i,j1)*f(i);            % X2=dPi/dfj   
                       X5=0;                                    % X5=dvi2/dej
                       X6=0;                                    % X6=dvi2/dfj
                       p=2*i-1;q=2*j1-1;J(p,q)=X5;J(p,N)=DV;    % PV节点电压误差
                       m=p+1;
                       J(m,q)=X1;J(m,N)=DP;q=q+1;J(p,q)=X6;     % PV节点有功误差
                       J(m,q)=X2;
                    elseif j1==i & j1~=isb                      % 非平衡节点 & 对角元    
                       X1=-C(i)-G(i,i)*e(i)-B(i,i)*f(i);        % X1=dPi/dei
                       X2=-D(i)+B(i,i)*e(i)-G(i,i)*f(i);        % X2=dPi/dfi  
                       X5=-2*e(i);                              % X5=dvi2/dei
                       X6=-2*f(i);                              % X6=dvi2/dfi
                       p=2*i-1;q=2*j1-1;J(p,q)=X5;J(p,N)=DV;  % PV节点电压误差
                       m=p+1;
                       J(m,q)=X1;J(m,N)=DP;q=q+1;J(p,q)=X6;   % PV节点有功误差
                       J(m,q)=X2;
                    end
                end
            end
        end
    end
    %% 将Jacobi矩阵化成单位矩阵
    for k=3:N0                      % N0=2*n （从第三行开始，第一、二行是平衡节点）
        k1=k+1;N1=N;                % N=N0+1 即 N=2*n+1扩展列△P、△Q 或 △U
        for k2=k1:N1                % 从k+1列的Jacobi元素到扩展列的△P、△Q 或 △U  
            J(k,k2)=J(k,k2)./J(k,k);% 用K行K列对角元素去除K行K列后的非对角元素进行规格化    
        end
            J(k,k)=1;               % 对角元规格化K行K列对角元素赋1
        %% 回代运算  
        if k~=3                     % 不是第三行  k > 3
            k4=k-1;
            for k3=3:k4                                     % 用k3行从第三行开始到当前行的前一行k4行消去
                for k2=k1:N1                                % k3行后各行上三角元素
                    J(k3,k2)=J(k3,k2)-J(k3,k)*J(k,k2);      % 消去运算（当前行k列元素消为0）
                end                                         % 用当前行K2列元素减去当前行k列元素乘以第k行K2列元素
                J(k3,k)=0;                                  % 当前行第k列元素已消为0
            end
            if k==N0                                        % 若已到最后一行
                break;
            end
        %% 前代运算  
            for k3=k1:N0                                    % 从k+1行到2*n最后一行
                for k2=k1:N1                                % 从k+1列到扩展列消去k+1行后各行下三角元素
                    J(k3,k2)=J(k3,k2)-J(k3,k)*J(k,k2);      % 消去运算
                end                                         % 用当前行K2列元素减去当前行k列元素乘以第k行K2列元素
                J(k3,k)=0;                                  % 当前行第k列元素已消为0
            end
        else                                                % 是第三行k=3
        %% 第三行k=3的前代运算      
             for k3=k1:N0                                   % 从第四行到2n行（最后一行）
                 for k2=k1:N1                               % 从第四列到2n+1列（即扩展列）
                     J(k3,k2)=J(k3,k2)-J(k3,k)*J(k,k2);     % 消去运算（当前行3列元素消为0）
                 end                                        % 用当前行K2列元素减去当前行3列元素乘以第三行K2列元素
                 J(k3,k)=0;                                 % 当前行第3列元素已消为0
             end
          end
    end

    for k=3:2:N0-1
        L=(k+1)./2;
        e(L)=e(L)-J(k,N);       % 修改节点电压实部
        k1=k+1;
        f(L)=f(L)-J(k1,N);      % 修改节点电压虚部
    end
    
    % 修改节点电压
    for k=3:N0
        DET=abs(J(k,N));
        if DET>=pr              % 电压偏差量是否满足要求
            IT2=IT2+1;          % 不满足要求的节点数加1
        end
    end
    ICT2(a)=IT2;                % 不满足要求的节点数
    ICT1=ICT1+1;                % 迭代次数
end

disp('迭代次数：');
disp(ICT1);
disp('没有达到精度要求的个数：');
disp(ICT2);
for k=1:n
    V(k)=sqrt(e(k)^2+f(k)^2);           % 计算各节点电压的模值
    sida(k)=atan(f(k)./e(k))*180./pi;   % 计算各节点电压的角度
    E(k)=e(k)+f(k)*j;                   % 将各节点电压用复数表示
end
%% 计算各输出量
disp('各节点的实际电压标幺值E为(节点号从小到大排列)：');
disp(E);    %显示各节点的实际电压标幺值E用复数表示
disp('各节点的电压大小V为(节点号从小到大排列)：');
disp(V);    %显示各节点的电压大小V的模值
disp('各节点的电压相角sida为(节点号从小到大排列)：');
disp(sida);     %显示各节点的电压相角
for p=1:n
    C(p)=0;
    for q=1:n
        C(p)=C(p)+conj(Y(p,q))*conj(E(q));  %计算各节点的注入电流的共轭值
    end
    S(p)=E(p)*C(p);     %计算各节点的功率 S = 电压 X 注入电流的共轭值
end
disp('各节点的功率S为(节点号从小到大排列)：');
disp(S);        %显示各节点的注入功率
disp('各条支路的首端功率Si为(顺序与B1一致)：');
for i=1:nl
       p=Bl(i,1);q=Bl(i,2);
       if Bl(i,6)==0
          Si(p,q)=E(p)*(conj(E(p))*conj(Bl(i,4)./2)+(conj(E(p)*Bl(i,5))...
          -conj(E(q)))*conj(1./(Bl(i,3)*Bl(i,5))));
            Siz(i)=Si(p,q);
        else
           Si(p,q)=E(p)*(conj(E(p))*conj(Bl(i,4)./2)+(conj(E(p)./Bl(i,5))...
           -conj(E(q)))*conj(1./(Bl(i,3)*Bl(i,5))));
            Siz(i)=Si(p,q);
        end
    disp(Si(p,q));
    SSi(p,q)=Si(p,q);
    ZF=['S(',num2str(p),',',num2str(q),')=',num2str(SSi(p,q))];
    disp(ZF);
end
disp('各条支路的末端功率Sj为(顺序与B1一致)：');
for i=1:nl
    p=Bl(i,1);q=Bl(i,2);
    if Bl(i,6)==0
       Sj(q,p)=E(q)*(conj(E(q))*conj(Bl(i,4)./2)+(conj(E(q)./Bl(i,5))...
       -conj(E(p)))*conj(1./(Bl(i,3)*Bl(i,5))));
        Sjy(i)=Sj(q,p);
    else
       Sj(q,p)=E(q)*(conj(E(q))*conj(Bl(i,4)./2)+(conj(E(q)*Bl(i,5))...
       -conj(E(p)))*conj(1./(Bl(i,3)*Bl(i,5))));
        Sjy(i)=Sj(q,p);
    end
    disp(Sj(q,p));
    SSj(q,p)=Sj(q,p);
    ZF=['S(',num2str(q),',',num2str(p),')=',num2str(SSj(q,p))];
    disp(ZF);
end
disp('各条支路的功率损耗DS为(顺序与B1一致)：');
for i=1:nl
        p=Bl(i,1);q=Bl(i,2);
    dS(i)=Si(p,q)+Sj(q,p);
    disp(dS(i));
    ddS(i)=dS(i);
    ZF=['dS(',num2str(p),',',num2str(q),')=',num2str(ddS(i))];
    disp(ZF);
end
