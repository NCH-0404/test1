clear
%% 参数设置
A11=
E1=
A211=
V2=
U2=
I=
tau=10^6;  %正常数 取值10^6-10^(-6)


%% 定义变量
setlmis([]);

P=lmivar(2,[]);%P(n-l)阶方阵 正定
phi=lmivar(2,[])%phi(n-rE1)阶方阵 可逆
G=lmivar(2,[])% 可逆
K=lmivar(2,[])
eps=lmivar(1,[1 0])%正常数


%L1=Kinv(G)
%% 描述方程
%11
lmiterm([1 1 1 P],A11,E1','s');
lmiterm([1 1 1 phi],A11*V2,U2,'s');

lmiterm([1 1 1 K],(-1)*E1,A211,'s');

lmiterm([1 1 1 eps],1,I,'s');

%21
lmiterm([1 2 1 P],A211,E1');
lmiterm([1 2 1 phi],A211*V2,U2);

lmiterm([1 2 1 G],-1,A211);

lmiterm([1 2 1 -K],-tau,E1');

%22
lmiterm([1 2 2 G],-tau,1,'s');

%% 结果

lmisys=getlmis;
[tmin xfeas]=feasp(lmisys)

if tmin<0
    Pmat=dec2mat(lmisys,xfeas,P) 
    phimat=dec2mat(lmisys,xfeas,phi) 
    Gmat=dec2mat(lmisys,xfeas,G) 
    Kmat=dec2mat(lmisys,xfeas,K) 
    epsmat=dec2mat(lmisys,xfeas,eps) 
    
    
    L1=K*inv(G)
end
