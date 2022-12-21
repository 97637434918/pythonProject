import numpy as num
import matplotlib.pyplot as plt
import matplotlib.animation as anim
fig=plt.figure()
paper1=fig.add_subplot(221)
paper2=fig.add_subplot(222)
paper3=fig.add_subplot(223)
paper1.set(xlim=[-15,15],ylim=[-2,2])
paper2.set(xlim=[-15,15],ylim=[-2,2])
paper3.set(xlim=[-15,15],ylim=[-15,15])
paper1.spines['top'].set_visible(False)
paper2.spines['top'].set_visible(False)
paper3.spines['top'].set_visible(False)
paper1.spines['right'].set_visible(False)
paper2.spines['right'].set_visible(False)
paper3.spines['right'].set_visible(False)
paper1.spines['bottom'].set_position(('data',0))
paper2.spines['bottom'].set_position(('data',0))
paper3.spines['bottom'].set_position(('data',0))
paper1.spines['left'].set_position(('data',0))
paper2.spines['left'].set_position(('data',0))
paper3.spines['left'].set_position(('data',0))
Dx=num.linspace(-4*num.pi,4*num.pi,300)
fx1=num.sin(Dx)
fx2=num.cos(Dx)
fx3=num.tan(Dx)
ax,=paper1.plot(Dx,fx1,label='SinX')
bx,=paper2.plot(Dx,fx2,label='CosX')
cx,=paper3.plot(Dx,fx3,label='Tanx')
kax,=paper1.plot(0,0,label='k')
kbx,=paper2.plot(0,0,label='k')
kcx,=paper3.plot(0,0,label='k')
text1=paper1.text(9,1,'')
text2=paper2.text(9,1,'')
text3=paper3.text(9,8,'')
sc1,=paper1.plot(0,0,'o',color='red')
sc2,=paper2.plot(0,0,'o',color='red')
sc3,=paper3.plot(0,0,'o',color='red')
def k(x0):
    dx1 = x0 + 0.0001
    k1 = (num.sin(dx1) - num.sin(x0)) / 0.0001
    yt1 = num.sin(x0)
    b1 = yt1 - k1 * x0
    dx21 = x0 + 1
    dy21 = k1 * (dx21) + b1
    dx31 = x0 - 1
    dy31 = k1 * (dx31) + b1
    xff1 = [dx21, dx31]
    yff1 = [dy21, dy31]
    k2 = (num.cos(dx1) - num.cos(x0)) / 0.0001
    yt2 = num.cos(x0)
    b2 = yt2 - k2 * x0
    dx22 = x0 + 1
    dy22 = k2 * (dx22) + b2
    dx32 = x0 - 1
    dy32 = k2 * (dx32) + b2
    xff2 = [dx22, dx32]
    yff2 = [dy22, dy32]
    k3 = (num.tan(dx1) - num.tan(x0)) / 0.0001
    yt3 = num.tan(x0)
    b3 = yt3 - k3 * x0
    dx23 = x0 + 1
    dy23 = k3 * (dx23) + b3
    dx33 = x0 - 1
    dy33 = k3 * (dx33) + b3
    xff3 = [dx23, dx33]
    yff3 = [dy23, dy33]
    return (xff1, yff1,xff2,yff2,xff3,yff3,x0,k1,k2,k3,yt1,yt2,yt3)
def update(x):
    a,b,c,d,e,f,x0,k1,k2,k3,yt1,yt2,yt3=k(Dx[x])
    cr1=kax.set_data(a,b)
    cr2=kbx.set_data(c,d)
    cr3=kcx.set_data(e,f)
    k11 = text1.set_text(f'x={round(x0,3)}\nk={round(k1,3)}')
    k12 = text2.set_text(f'x={round(x0,3)}\nk={round(k2,3)}')
    k13 = text3.set_text(f'x={round(x0,3)}\nk={round(k3,3)}')
    sc1.set_data(x0, yt1)
    sc2.set_data(x0, yt2)
    sc3.set_data(x0, yt3)
    return(cr1,cr2,cr3,k11,k12,k13)
paper1.legend(loc=3)
paper2.legend(loc=3)
paper3.legend(loc=3)
paper1.grid('both',linestyle=':')
paper2.grid('both',linestyle=':')
paper3.grid('both',linestyle=':')
ani=anim.FuncAnimation(fig=fig,func=update,frames=range(0,200),interval=100,repeat=True)
plt.show()
