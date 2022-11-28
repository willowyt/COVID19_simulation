

T=90     # 天数
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
p=0.10           # 亚临床感染者的比例   （1-p)无症状感染者的概率
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

#清零后累计感染人数
int(sum(D[I_total>=1]['I_total']))


# ### 核酸检测
# 



T=90     # 天数
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
p=0.10           # 亚临床感染者的比例   （1-p)无症状感染者的概率
r1=1/2.5           # 亚临床感染转化为临床感染的比例
r2=1  #1/2       # 临床感染者进行隔离的比例
r3=1/5           # 隔离治疗恢复率
r4=1/7           # 无症状感染的自愈率
k=1           # 核酸检测时的风险系数
delta_pcrs=1     # 亚临床感染个体核酸检测准确率
delta_pcra=1    # 对无症状个体核酸检测准确率


# ### 每日核酸

for i in range(T-1):
    S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
    E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
    Is[i+1]=Is[i]+p*r*E[i]-delta_pcrs*Is[i]-(1-delta_pcrs)*r1*Is[i]
    Ia[i+1]=Ia[i]+(1-p)*r*E[i]-delta_pcra*Ia[i]-(1-delta_pcra)* r4* Ia[i]
    I[i+1]=I[i]+r1*(1-delta_pcrs)*Is[i]-r2*I[i]
    H[i+1]=H[i]+delta_pcrs*Is[i]+delta_pcra*Ia[i]+r2*I[i]-r3*H[i]
    R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)* r4* Ia[i]

import numpy as np
Is=np.array(Is)
Ia=np.array(Ia)
I_total1=Is+Ia+I+E
print(I_total1)

print(Is)
print(Ia)
I_pro1=(sum(I_total)-sum(I_total1))/sum(I_total)
print(I_pro1)

import pandas as pd
D1={'S':S,'E':E,'Is':Is,'Ia':Ia,'I':I,'H':H,'R':R,'I_total':I_total1}
D1=pd.DataFrame(D1)
D1

#清零后累计感染人数
int(sum(D1[I_total>=1]['I_total']))


# ### 3天2检

# In[8]:


import numpy as np
T=90
arr=np.arange(0,90,2)
for i in range(T-1):
    if i in arr:
        S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-delta_pcrs*Is[i]-(1-delta_pcrs)*r1*Is[i]
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-delta_pcra*Ia[i]-(1-delta_pcra)* r4* Ia[i]
        I[i+1]=I[i]+r1*(1-delta_pcrs)*Is[i]-r2*I[i]
        H[i+1]=H[i]+delta_pcrs*Is[i]+delta_pcra*Ia[i]+r2*I[i]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)* r4* Ia[i]
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
I_total2=I+Is+Ia+E
print(Is)
print(Ia)
print(I)
print(I_total2)

I_pro2=(sum(I_total)-sum(I_total2))/sum(I_total)
print(I_pro2)  

import pandas as pd
D2={'S':S,'E':E,'Is':Is,'Ia':Ia,'I':I,'H':H,'R':R,'I_total':I_total2}
D2=pd.DataFrame(D2)
D2

#清零后累计感染人数
int(sum(D2[I_total>=1]['I_total']))



# #### 3天1检
import numpy as np
lst=np.arange(0,90,3)
for i in range(T-1):
    if i in lst:
        S[i+1]=S[i]-k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N
        E[i+1]=E[i]+k*beta*S[i]*(eta_s*Is[i]+eta_A*Ia[i]+I[i])/N-r*E[i]
        Is[i+1]=Is[i]+p*r*E[i]-delta_pcrs*Is[i]-(1-delta_pcrs)*r1*Is[i]
        Ia[i+1]=Ia[i]+(1-p)*r*E[i]-delta_pcra*Ia[i]-(1-delta_pcra)* r4* Ia[i]
        I[i+1]=I[i]+r1*(1-delta_pcrs)*Is[i]-r2*I[i]
        H[i+1]=H[i]+delta_pcrs*Is[i]+delta_pcra*Ia[i]+r2*I[i]-r3*H[i]
        R[i+1]=R[i]+r3*H[i]+(1-delta_pcra)* r4* Ia[i]
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
I_total3=I+Is+Ia+E
print(Is)
print(Ia)
print(I)
print(I_total3)
I_pro3=(sum(I_total)-sum(I_total3))/sum(I_total)
print(I_pro3) 

import pandas as pd
D3={'S':S,'E':E,'Is':Is,'Ia':Ia,'I':I,'H':H,'R':R,'I_total':I_total3}
D3=pd.DataFrame(D3)
D3

#清零后累计感染人数
int(sum(D3[I_total>=1]['I_total']))


###感染控制情况
I_all=[I_pro1,I_pro2,I_pro3]
I_pro=[1-I_pro1,1-I_pro2,1-I_pro3]
I_all


import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 微软雅黑
mpl.rcParams['font.serif'] = ['Microsoft YaHei']
test = ('每日核酸','3天2检', '3天1检')
plt.bar(test, I_all,width=0.3,color='lightcoral')
plt.title('情况对比')
plt.ylim(0,1)

from matplotlib.ticker import FuncFormatter
def to_percent(temp, position):
  return '%1.0f'%(100*temp) + '%'
plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
plt.show()






