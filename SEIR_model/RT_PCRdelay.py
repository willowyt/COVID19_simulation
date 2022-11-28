
# ### 考虑核酸检测延迟结果1天

# In[1]:


#无检测
T=90    # 天数
N=10000         # 总人数
####各类人群
E=[0]*T         # 潜伏者人数
Ia=[5]*T        # 无症状感染者
Is=[5]*T        # 亚临床感染者
I=[0]*T         # 临床感染者
S=[N-I[0]-Ia[0]-Is[0]]*T    # 易感人数
H=[0]*T         # 隔离人数
R=[0]*T         # 康复人数

#####参数
beta=9.5/3.375# 未隔离感染者在单位时间内平均接触并传染的数量
eta_s=0.5        # 亚临床感染者的传染率
eta_A=0.5        # 无症状感染者的传染率
r=1/3            # 潜伏着转化为感染者概率(潜伏爆发率)
p=0.1           # 亚临床感染者的比例   （1-p)无症状感染者的概率
r1=1/2.5           # 亚临床感染转化为临床感染的比例
r2=1  #1/2       # 临床感染者进行隔离的比例
r3=1/5           # 隔离治疗恢复率
r4=1/7           # 无症状感染的自愈率

#无检测
for i in range(T-1):
    S[i+1]=S[i]-beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
    E[i+1]=E[i]+beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
    Is[i+1]=Is[i]+p*r*E[i]-r1*Is[i]
    Ia[i+1]=Ia[i]+(1-p)*r*E[i]-r4*Ia[i]
    I[i+1]=I[i]+r1*Is[i]-r2*I[i]
    H[i+1]=H[i]+r2*I[i]-r3*H[i]
    R[i+1]=R[i]+r3*H[i]+r4*Ia[i]

import numpy as np
I=np.array(I)
Is=np.array(Is)
Ia=np.array(Ia)
I_total=I+Is+Ia+E
print(I_total)

import pandas as pd
D={'S':S,'E':E,'Is':Is,'Ia':Ia,'I':I,'H':H,'R':R,'I_total':I_total}
D=pd.DataFrame(D)
D





###参数
T=90      # 天数
N=10000        # 总人数
####各类人群
E=[0]*T         # 潜伏者人数
Ia=[5]*T        # 无症状感染者
Is=[5]*T        # 亚临床感染者
I=[0]*T         # 临床感染者
S=[N-I[0]-Ia[0]-Is[0]]*T    # 易感人数
H=[0]*T         # 隔离人数
R=[0]*T         # 康复人数

#####参数
beta=9.5/3.375     # 未隔离感染者在单位时间内平均接触并传染的数量
eta_s=0.5        # 亚临床感染者的传染率
eta_A=0.5        # 无症状感染者的传染率
r=1/3            # 潜伏着转化为感染者概率(潜伏爆发率)
p=0.1           # 亚临床感染者的比例   （1-p)无症状感染者的概率
r1=1/2.5           # 亚临床感染转化为临床感染的比例
r2=1  #1/2       # 临床感染者进行隔离的比例
r3=1/5           # 隔离治疗恢复率
r4=1/7           # 无症状感染的自愈率
k=1           # 核酸检测时的风险系数
delta_pcrs=1     # 亚临床感染个体核酸检测准确率
delta_pcra=1      # 对无症状个体核酸检测准确率



#### 每日核酸延迟1天
Ks=[0]*T
Ka=[0]*T
i=0
S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
Is[i+1]=Is[i]+p*r*E[i]-r1*(1-delta_pcrs)*Is[i]
Ia[i+1]=Ia[i]+(1-p)*r*E[i]-(1-delta_pcra)*r4*Ia[i]
I[i+1]=I[i]+r1*(1-delta_pcrs)*Is[i]-I[i]*r2
H[i+1]=H[i]+r2*I[i]-r3*H[i]
R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)*r4*Ia[i]

i=1
S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
Is[i+1]=Is[i]+p*r*E[i]-delta_pcrs*Is[i-1]-(1-delta_pcrs)*(Is[i]-delta_pcrs*Is[i-1])*r1
Ia[i+1]=Ia[i]+(1-p)*r*E[i]-delta_pcra*Ia[i-1]-(1-delta_pcra)*(Is[i]-delta_pcrs*Is[i-1])*r4
I[i+1]=I[i]+(1-delta_pcrs)*(Is[i]-delta_pcrs*Is[i-1])*r1-r2*I[i]
H[i+1]=H[i]+delta_pcrs*Is[i-1]+delta_pcra*Ia[i-1]+r2*I[i]-r3*H[i]
R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)*(Is[i]-delta_pcrs*Is[i-1])*r4

i=2
S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
Ks[1]=Is[1]-delta_pcrs*Is[0]
Ka[1]=Ia[1]-delta_pcra*Ia[0]
Is[i+1]=Is[i]+p*r*E[i]-delta_pcrs*Ks[1]-(1-delta_pcrs)*(Is[i]-delta_pcrs*Ks[1])*r1
Ia[i+1]=Ia[i]+(1-p)*r*E[i]-delta_pcra*Ka[1]-(1-delta_pcra)*(Ia[i]-delta_pcra*Ka[1])*r4
I[i+1]=I[i]+(1-delta_pcrs)*(Is[i]-delta_pcrs*Ks[1])*r1-r2*I[i]
H[i+1]=H[i]+delta_pcrs*Ks[1]+delta_pcra*Ka[1]-r3*H[i]+r2*I[i]
R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)*(Ia[i]-delta_pcra*Ka[1])*r4

for i in range(3,T-1):
    S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
    E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
    Ks[i-1]=Is[i-1]-delta_pcrs*Ks[i-2]
    Ka[i-1]=Ia[i-1]-delta_pcra*Ka[i-2]
    Is[i+1]=Is[i]+p*r*E[i]-delta_pcrs*Ks[i-1] -(1-delta_pcrs)*(Is[i]-delta_pcrs*Ks[i-1])*r1
    Ia[i+1]=Ia[i]+(1-p)*r*E[i]-delta_pcra*Ka[i-1]-(1-delta_pcra)*(Ia[i]-delta_pcra*Ka[i-1])*r4
    I[i+1]=I[i]+(1-delta_pcrs)*(Is[i]-delta_pcrs*Ks[i-1])*r1-r2*I[i]
    H[i+1]=H[i]+delta_pcrs*Ks[i-1]+delta_pcra*Ka[i-1]-r3*H[i]+r2*I[i]
    R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)*(Ia[i]-delta_pcra*Ka[i-1])*r4

import numpy as np
Is=np.array(Is)
Ia=np.array(Ia)
I_total11=Is+Ia+I+E
#print(I_total1)

I_pro11=(sum(I_total)-sum(I_total11))/sum(I_total)
print(I_pro11)
import pandas as pd
D11={'S':S,'E':E,'Is':Is,'Ia':Ia,'I':I,'H':H,'R':R,'I_total':I_total11}
D11=pd.DataFrame(D11)
D11

#清零后累计感染人数
int(sum(D11[I_total>=1]['I_total']))



### 3天2检延迟1天
import numpy as np
T=90
arr1=np.arange(0,90,2) #测试日
arr2=arr1+1  #报告日
for i in range(T-1):
    if i in arr1:
        S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-(1-delta_pcrs)*r1*Is[i]
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-(1-delta_pcra)* r4* Ia[i]
        I[i+1]=I[i]+r1*(1-delta_pcrs)*Is[i]-r2*I[i]
        H[i+1]=H[i]+r2*I[i]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)* r4* Ia[i]
    elif i in arr2:
        S[i+1]=S[i]-beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-r1*(Is[i]-delta_pcrs*Is[i-1])-delta_pcrs*Is[i-1] 
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-r4*(Ia[i]-delta_pcra*Ia[i-1])-delta_pcra*Ia[i-1]
        I[i+1]=I[i]+r1*(Is[i]-delta_pcrs*Is[i-1])-r2*I[i]
        H[i+1]=H[i]+r2*I[i]+delta_pcrs*Is[i-1]+delta_pcra*Ia[i-1]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+r4*(Ia[i]-delta_pcra*Ia[i-1])
    else:
        S[i+1]=S[i]-beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-r1*Is[i]
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-r4*Ia[i]
        I[i+1]=I[i]+r1*Is[i]-r2*I[i]
        H[i+1]=H[i]+r2*I[i]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+r4*Ia[i]

import numpy as np
Is=np.array(Is)
Ia=np.array(Ia)
I_total22=I+Is+Ia+E
#print(I_total22)

I_pro22=(sum(I_total)-sum(I_total22))/sum(I_total)
print(I_pro22)     
D22={'S':S,'E':E,'Is':Is,'Ia':Ia,'I':I,'H':H,'R':R,'I_total':I_total22}
D22=pd.DataFrame(D22)
D22


#清零后累计感染人数
int(sum(D22[I_total>=1]['I_total']))



##3天1检 结果延迟1天
import numpy as np
arr1=np.arange(0,90,3) #测试日
arr2=arr1+1  #报告日
for i in range(T-1):
    if i in arr1:
        S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-(1-delta_pcrs)*r1*Is[i]
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-(1-delta_pcra)* r4* Ia[i]
        I[i+1]=I[i]+r1*(1-delta_pcrs)*Is[i]-r2*I[i]
        H[i+1]=H[i]+r2*I[i]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)* r4* Ia[i]
    elif i in arr2:
        S[i+1]=S[i]-beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-r1*(Is[i]-delta_pcrs*Is[i-1])-delta_pcrs*Is[i-1] 
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-r4*(Ia[i]-delta_pcra*Ia[i-1])-delta_pcra*Ia[i-1]
        I[i+1]=I[i]+r1*(Is[i]-delta_pcrs*Is[i-1])-r2*I[i]
        H[i+1]=H[i]+r2*I[i]+delta_pcrs*Is[i-1]+delta_pcra*Ia[i-1]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+r4*(Ia[i]-delta_pcra*Ia[i-1])
    else:
        S[i+1]=S[i]-beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-r1*Is[i]
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-r4*Ia[i]
        I[i+1]=I[i]+r1*Is[i]-r2*I[i]
        H[i+1]=H[i]+r2*I[i]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+r4*Ia[i]

import numpy as np
Is=np.array(Is)
Ia=np.array(Ia)
I_total33=I+Is+Ia+E
#print(I_total33)
I_pro33=(sum(I_total)-sum(I_total33))/sum(I_total)
print(I_pro33)   
D33={'S':S,'E':E,'Is':Is,'Ia':Ia,'I':I,'H':H,'R':R,'I_total':I_total33}
D33=pd.DataFrame(D33)
D33


#清零后累计感染人数
int(sum(D33[I_total>=1]['I_total']))



I_pro_delay=[I_pro11,I_pro22,I_pro33]
I_pro_delay






