function mpc=case5_01
% MATPOWER Case Format:Version 2
mpc.version='2';
%% Power Flow Data %%
%% system MVA base 系统基准容量
mpc.baseMVA=100;    % 100MVA
%% bus data
% bus_i         设置母线编号(正整数)
% type          设置母线类型, 1 为 PQ 节点母线, 2 为 PV 节点母线, 3 为平衡(参考)节点母线，4 为孤立节点母线
% Pd Qd         设置母线注入负荷的有功功率和无功功率
% Gs Bs         设置与母线并联电导和电纳
% baseKV        设置该母线基准电压
% Vm Va         设置母线电压的幅值、相位初值
% Vmax Vmin     设置工作时母线最高、最低电压幅值
% area zone     设置电网断面号和分区号，一般都设置为 1
%         bus_i   type    Pd      Qd    Gs   Bs   area   Vm     Va    baseKV  zone Vmax   Vmin
mpc.bus=[    1      1     160     80     0    0     1     1      0     100     1    1.1   0.94;
             2      1     200     100    0    0     1     1      0     100     1    1.1   0.94;
             3      1     370     130    0    0     1     1      0     100     1    1.1   0.94;
             4      2      0       0     0    0     1     1.05   0     100     1    1.1   0.94;
             5      3      0       0     0    0     1     1.05   0     100     1    1.1   0.94;];

%% generator data
% bus       设置接入发电机(电源)的母线编号
% Pg Qg     设置接入发电机(电源)的有功功率和无功功率
% Pmax Pmin 设置接入发电机(电源)的有功功率最大、最小允许值
% Qmax Qmin 设置接入发电机(电源)的无功功率最大、最小允许值
% Vg        设置接入发电机(电源)的工作电压
% mBase     设置接入发电机(电源)的功率基准,如果为默认值,就是 baseMVA 变量的值
% status    设置发电机(电源)工作状态, 1 表示投入运行, 0 表示退出运行
%           bus     Pg      Qg    Qmax     Qmin     Vg     mBASE    status   Pmax   Pmin 
mpc.gen=[    4      500     0     99990    -9999   1.050    100        1      600     0;
             5       0      0     99990    -9999   1.050    100        1      600     0; ];

%%  branch   data
% fbus tbus             设置该支路由起始节点(母线)编号和终止节点(母线)编号
% r,x,b                 设置该支路的电阻、电抗和充电电纳
% rateA,rateB,rateC     设置该支路长期、短期和紧急允许功率
% ratio                 设置该支路的变比，如果支路元件是导线,那么 ratio 为 0；如果支路元件为变压器，则该变比为 fbus 侧母线的基准电压与 tbus 侧母线的基准电压之比
% angle                 设置支路的相位角度，如果支路元件为变压器(或移相器)，就是变压器(或移相器)的转角；如果支路元件是导线，相位角度则为 0°。
% status                设置支路工作状态，1 表示投入运行，0 表示退出运行
% angmin,angmax         设置支路相位角度的最小和最大差值
%            fbus   tbus    r      x      b      rateA     rateB    rateC   ratio   angle   status  agmin    agmax
mpc.branch=[  2      1    0.04   0.25     0.5      0         0         0      0       0        1     -360      360;
              3      1    0.1    0.35     0        0         0         0      0       0        1     -360      360;
              3      2    0.08    0.3     0.5      0         0         0      0       0        1     -360      360;
              3      5    0      0.03     0        0         0         0     1.05     0        1     -360      360;
              2      4    0      0.015    0        0         0         0     1.05     0        1     -360      360; ];
return;

