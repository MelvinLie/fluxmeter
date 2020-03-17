import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def mpl_setup():
    #plt.style.use('my_style')
    #plt.rc('text', usetex=True)
    plt.rcParams['lines.linewidth'] = 1.5
    plt.rcParams['lines.markersize'] = 9
    #print(mpl.rcParams.keys())
    plt.rcParams['axes.labelsize'] = 'large'
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.rcParams['font.size'] = 20
    plt.rcParams['axes.grid'] = True
    plt.rcParams[ 'axes.grid.which'] = 'both'
    plt.rcParams['font.family'] = 'DejaVu Sans'

    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16

    plt.rc('text',usetex=True)
    #font = {'family':'serif','size':16}
    font = {'family':'DejaVu Sans','size':28, 'serif': ['computer modern roman']}
    
    plt.rc('font',**font)
    plt.rc('legend',**{'fontsize':22})
    mpl.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']

    #z = np.linspace(0,1,100)
    #b_3 = np.sin(z*2*np.pi)

    #fig, ax = plt.subplots()

    #ax.plot(z,b_3)
    #ax.set_ylabel('$b_n$')
    #ax.set_xlabel('$z\\rightarrow$')
    #plt.show()

    #print('stop')
