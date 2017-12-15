import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



tf = 500
dt = 1
t = np.arange(0,(tf+1),dt)
Tout = np.loadtxt('1Tout.csv',delimiter=',')
md = np.loadtxt('2md.csv',delimiter=',')
Pb = np.loadtxt('3Pb.csv',delimiter=',')
Pf = np.loadtxt('4Pf.csv', delimiter=',')
Tdist = np.loadtxt('5Tdist.csv', delimiter=',')
xdist = np.loadtxt('6xdist.csv', delimiter=',')
mdsolve, mdysolve = np.loadtxt('7mdandysolve.csv')


fig1 = plt.figure(figsize=([12.5,6.5]))
fig1.suptitle('Transient Natural Circulation of the GaTE Loop')

#Setting up settings and un-moving parts of the figures
ax1 = fig1.add_subplot(221)
ax1.set_xlim(0,1.44)
ax1.set_ylim(40,230)
ax1.plot([0.36, 0.36],[30.0, 300.0],'--k')
ax1.plot([0.36*2, 0.36*2],[30.0, 300.0],'--k')
ax1.plot([0.36*3, 0.36*3],[170.0, 300.0],'--k')
ax1.plot([0.36*3, 0.36*3],[30.0, 100.0],'--k')
ax1.set_xlim(0,1.44)
ax1.set_ylim(40,230)
ax1.set_xlabel('Spatial distribution (m)')
ax1.set_ylabel('Temperature (C)')
ax1.annotate('Core',xy=(0.13,210),xytext=(0.13,210))
ax1.annotate('Hot Plenum',xy=(0.41,210),xytext=(0.41,210))
ax1.annotate('Coolers',xy=(0.81,210),xytext=(0.81,210))
ax1.annotate('Cold Leg',xy=(1.16,210),xytext=(1.16,210))

ax2 = fig1.add_subplot(222)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Temperature (C)')
ax2.legend()
ax2.set_xlim(0,500)
ax2.set_ylim(40,230)  

ax3 = fig1.add_subplot(223)
ax3.set_xlim(0.006,0.12)
ax3.set_ylim(0,2000)
ax3.set_xlabel('Flow rate (kg/s)')
ax3.set_ylabel('Pressure (P)')
ax3.legend(loc='upper center')

ax4 = fig1.add_subplot(224)
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('Flow Rate (kg/s)')
ax4.legend()
ax4.set_xlim(0,500)
ax4.set_ylim(0,0.12)

# The actors!
line1, = ax1.plot([], [],'k', label = 'Temp Dist.(t)')
line1a, = ax1.plot([], [],'c*', label = 'Core Exit Temp')
line1b, = ax1.plot([], [],'C1*', label = 'Outlet Temp')
line1c, = ax1.plot([],[],':k', label = 'S.S. Temp Dist.')

line2, = ax2.plot([], [],'C1', label = 'Outlet Temp')
line2a, = ax2.plot([], [],'c', label = 'Core Exit Temp')

line3, = ax3.plot([], [],'g', label = 'Buoyancy(t)')
line3a, = ax3.plot([], [],'C3--', label = 'Mass Flow Solution')
line3b, = ax3.plot([], [],':g', label = 'S.S. Buoyancy')
line3c, = ax3.plot([], [],'C1', label = 'Friction')

line4, = ax4.plot([], [],'C3', label = 'Mass Flow Rate')

#legends
ax1.legend([line1,line1a,line1b,line1c], 
           [line1.get_label(),line1a.get_label(),
            line1b.get_label(), line1c.get_label()],
            loc='center right')
ax2.legend([line2,line2a,],[line2.get_label(),line2a.get_label()])
ax3.legend([line3,line3a,line3b,line3c], 
           [line3.get_label(),line3a.get_label(),
            line3b.get_label(), line3c.get_label()],
            loc='upper center')
ax4.legend([line4,],[line4.get_label()])

def init():
    line1.set_data([], [])
    line1a.set_data([], [])
    line1b.set_data([], [])
    line1c.set_data([], [])
    line2.set_data([], [])
    line2a.set_data([], [])
    line3.set_data([], [])
    line3a.set_data([], [])
    line3b.set_data([], [])
    line3c.set_data([], [])
    line4.set_data([], [])
    return line1, line1a, line1b, line1c, line2, line2a, line3, line3a, line3b, line3c, line4,

def animate(i):
    a = xdist[i]
    b = Tdist[i]
    line1.set_data(a,b)
    aa =  ([0.36, 0.36])
    bb = ([Tout[i,0], Tout[i,0]])
    line1a.set_data(aa,bb)
    aaa = ([0.36*2, 0.36*2])
    bbb = [Tout[i,40], Tout[i,40]]
    line1b.set_data(aaa,bbb)
    aaaa = (xdist[0])
    bbbb = (Tdist[0])
    line1c.set_data(aaaa,bbbb)
    
    c = t[0:i]
    d = Tout[0:i,40]
    line2.set_data(c,d)
    cc = t[0:i]
    dd = Tout[0:i,0]
    line2a.set_data(cc,dd)
    
    e = md[i]
    f = Pb[i+1]
    line3.set_data(e,f)
    ee = ([mdsolve[i+1],mdsolve[i+1]])
    ff = ([0,mdysolve[i+1]])
    line3a.set_data(ee,ff)
    eee = (md[0])
    fff = Pb[1]
    line3b.set_data(eee,fff)
    eeee = md[0]
    ffff = Pf
    line3c.set_data(eeee,ffff)
    
    x = t[0:i]
    y = mdsolve[0:i]
    line4.set_data(x,y)
    
    return line1, line1a, line1b, line1c, line2, line2a, line3, line3a, line3b, line3c, line4,

ani = FuncAnimation(fig1, animate, frames=499,
                    init_func=init, blit=False)

ani.save('Ward_ME701_Final_Animation.mp4', fps=15)
