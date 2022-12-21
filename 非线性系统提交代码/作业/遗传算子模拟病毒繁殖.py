import numpy as np
import random
from itertools import chain
import matplotlib.pyplot as plt
import matplotlib.animation as ani
#绘图设置
fig=plt.figure()
ax=fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_position(('data',0))
ax.spines['left'].set_position(('data',0))
ax.set(xlabel='generation',ylabel='fitness',title='Evolutionary trend of virus reproduction \norganisms under natural selection')
ax.grid('both',linestyle=':')
class organisms():#定义无性生物类，简化代码
    pop0x=20
    mul0x=0.01

    def __init__(self):
        pass
    def gene(self):
        x=list(np.random.normal(10,5,20))#基于正态分布设置初代基因组
        return(x)
    def innum(self):
        ininum=random.randint(20,70)
        print('随机设置初始数量{}'.format(ininum))
        return(ininum)
def inig(self):#初始个体数与基因组
    ini=self.innum()
    k=range(0,ini)
    datax={}
    for i in k:
        bgene=self.gene()
        i2={'{}'.format(i):bgene}
        datax.update(i2)
    print('设置初始基因数据{}'.format(datax))
    data=datax
    return(data,k)

def fitcode(contdata,contum):#此处采用简单版遗传算子，基于适应度获取下一代基因组
    nextpare={}
    middle=[]
    x=-1
    y=-1
    for i in contum:
        sumfitf=np.sum(contdata['{}'.format(i)])
        middle.append(sumfitf)
        middlesum=np.sum(middle)
        z = {}
        w = {}
        h = {}
        l = {}
    for i in contum:#计算个体适应度
        fitf=np.sum(contdata['{}'.format(i)])

        if fitf >= fit0x1:
            pfitf=fitf
            x=x+1
            z.update({'{}'.format(x):'{}'.format(pfitf)})
            w.update({'{}'.format(x):'{}'.format(i)})
        elif fit0x2 < fitf <fit0x1:
            y=y+1
            pfitf=fitf-100
            h.update({'{}'.format(y): '{}'.format(pfitf)})
            l.update({'{}'.format(y): '{}'.format(i)})
        xlim = range(0,x)
        ylim = range(0,y)
        add1=[]
        add2=[]
    for i in xlim:
        add1.append(float(z['{}'.format(i)]))
    for i in ylim:
        add2.append(float(h['{}'.format(i)]))
    add=np.sum(add1)+np.sum(add2)
    addave=add/len(contum)
    print('此代适应度均值为{}'.format(addave))
    print('_________________________________________________________')
    prointe1={}
    prointe2={}
    prointe3={}
    for i in xlim:
        prointe1.update({'{}'.format(w['{}'.format(i)]): float(z['{}'.format(i)])/add})
    for i in ylim:
        prointe2.update({'{}'.format(l['{}'.format(i)]): float(h['{}'.format(i)])/add})
    prointe3.update(prointe1)
    prointe3.update(prointe2)
    keys=list(prointe3.keys())
    values=[]
    for i in keys:
        values.append(prointe3['{}'.format(i)])
    return(prointe3,keys,values,addave)
def choicee(key,value,data):#设置选择方式--加权随机抽取
    middl=[]
    nextgenerate=[]
    kk=random.randint(20,70)
    print('下一代个体数{}'.format(kk))
    middl.append((random.choices(population=key,weights=value,k=kk)))
    middl2=list(chain.from_iterable(middl))
    for i in middl2:
        nextgenerate.append(data['{}'.format(i)])
    print('根据各组适应度抽取下一代{}'.format(nextgenerate))
    num=len(nextgenerate)

    midd22={}
    for i in range(0,num):
        midd22.update({'{}'.format(i):nextgenerate[i]})
    keys22=range(0,num)
    return(midd22,keys22)
def start():#程序开始设置
    organism1=organisms()
    data,k=inig(organism1)
    global fit0x1#基于初代基因组设置两组适应度
    fit0x1=np.sum(data['0'])
    global fit0x2
    fit0x2=fit0x1-30
    print('基于初始数据设置适应度标准一{}'.format(fit0x1))
    print('适应度标准二{}'.format(fit0x2))
    prointee,keys,values,ed=fitcode(data,k)
    nextgenerate,keys22=choicee(keys,values,data)
    ii=1
    endjob=[]
    ed=ed-200
    endjob.append(ed)
    while ii<=50:#循环执行算子
        ii=ii+1
        a,b,c,endd=fitcode(nextgenerate,keys22)
        d,e=choicee(b,c,nextgenerate)
        keys22=e
        nextgenerate=d
        print('第{}代{}'.format(ii,d))
        endd=endd-200
        endjob.append(endd)
    else:
        return(endjob)
a=start()
xlim=range(0,51)

ax.set(xlim=[0,55],ylim=[-50,50])
paper,=ax.plot(0,0,'o')
y0=[]
x0=[]
def update(x):#动画设置
    y01=a[x]
    y0.append(y01)
    x0.append(x)
    paper.set_data(x0,y0)
    return(paper)
anim=ani.FuncAnimation(fig=fig,func=update,frames=xlim,interval=100)
anim.save('C:\Users\Lenovo\Desktop\交大学习任务及资料\近期作业\综合集成\遗传算子动图显示.gif',fps=30)