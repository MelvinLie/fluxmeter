import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

class Coil:
    def __init__(self,radius,layers,lens,widths,polarity):
        self.radius = radius
        self.layers = layers
        self.lens = lens
        self.widths = widths
        self.deltas = widths/radius
        self.lens_2 = lens/2
        self.deltas_2 = widths/radius/2
        self.polarity = polarity
        #print(self.deltas)
        
    def print(self):
        print('This is a coil with {} layers.'.format(self.layers))
        print('The radius is {} mm.'.format(self.radius*1e3))
        for i,ll in enumerate(self.lens):
            print('\tWinding {}: ( {} x {} ) mm'.format(i,ll*1e3,self.deltas[i]*radius*1e3))
            
    def compute_kn(self,n):
        kn = 0.
        if not(n == 0):
            for i,l in enumerate(self.lens):
                #print(l*np.sin(n*self.deltas_2[i]))
                kn += self.polarity[i]*l*np.sin(n*self.deltas_2[i])
                
            #print(kn)  
            kn *= 2*self.layers*self.radius**(n)/n
            
        return kn

    
    def assemble_D_fourier(self,n,z_m,K):
        
        Z_1 = max(z_m)
        Z_0 = min(z_m)
        
        L = Z_1 - Z_0
        
        M = len(z_m)
        
        D = 1j*np.zeros((M,K-1))
        D_0 = 1j*np.zeros((M,1))
        
   
        k = np.array([kk for kk in range(1,K)])
        
        kk,zz_m = np.meshgrid(k,z_m)
        

        for i,l_2 in enumerate(self.lens_2):
            D += self.polarity[i]*np.sin(n*self.deltas_2[i])*np.sin(2*np.pi*kk*l_2/L)
            #print(np.sin(n*self.deltas_2[i])*self.lens[i])
            D_0 += self.polarity[i]*np.sin(n*self.deltas_2[i])*self.lens[i]
            
        D *= 2*self.layers*self.radius*L*np.exp(2j*np.pi*kk*(Z_0-zz_m)/L)/n/kk/np.pi
        #print(D_0[0])
        D_0 *= 2*self.layers*self.radius/n     
        
        #print(D_0*self.radius**(n-1))
        
        D = np.append(np.append(np.conj(D[:,-1::-1]),D_0,axis=1),D,axis=1)
        
        return D
    
    def plot_layout(self):
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(1,1,1)
        
        for i,l_2 in enumerate(self.lens_2):
            w_2 = self.widths[i]/2*1e3
            ax.plot([-l_2*1e3,l_2*1e3,l_2*1e3,-l_2*1e3,-l_2*1e3],[-w_2,-w_2,w_2,w_2,-w_2],color=(0.72,0.45,0.2),linewidth=2)
        ax.set_xlabel(r'$z\rightarrow$ [mm]')
        ax.set_ylabel(r'$x\rightarrow$ [mm]')
        max_len = max(self.lens)
        max_width = max(self.widths)
        if max_len > max_width:
            ax.set_ylim(ax.get_xlim())
        else:
            ax.set_xlim(ax.get_ylim())
        
    def assemble_D_fluxmeter_fourier_univariate(self,z_m,K):
        
        Z_1 = max(z_m)
        Z_0 = min(z_m)
        
        L = Z_1 - Z_0
        
        M = len(z_m)
        
        D = 1j*np.zeros((M,K-1))
        D_0 = 1j*np.zeros((M,1))

        k = np.array([kk for kk in range(1,K)])
        
        kk,zz_m = np.meshgrid(k,z_m)
        
        for i,l_2 in enumerate(self.lens_2):
            D += self.polarity[i]*self.widths[i]*np.sin(2*np.pi*kk*l_2/L)*np.exp(2j*np.pi*kk*(zz_m-Z_0)/L)
            #print(np.sin(n*self.deltas_2[i])*self.lens[i])
            D_0 += self.polarity[i]*self.widths[i]*self.lens[i]
        
        D *= L/np.pi/kk
        
        D = np.append(np.append(np.conj(D[:,-1::-1]),D_0,axis=1),D,axis=1)
        
        return D      
    

    def integrate_z(self,fcn,z_m):

        int_array = np.zeros(z_m.shape)

        for m,zz_m in enumerate(z_m):

            for i,l_2 in enumerate(self.lens_2):
                int_array[m] += self.widths[i]*integrate.quad(fcn, zz_m-l_2, zz_m+l_2)[0]

        return self.layers*int_array


def transform_to_spatial_domain(u,z,Z_0,Z_1):
    
    K = np.int32(len(u)/2)
    
    k = np.array([-K+kk for kk in range(0,2*K+1)])
    
    kk,zz = np.meshgrid(k,z)
    
    return np.dot(np.exp(-2j*np.pi*kk*(zz-Z_0)/(Z_1-Z_0)),u)

def transform_to_frequency_domain(sig,z,Z_0,Z_1,K):
    

    k = np.array([-K+kk+1 for kk in range(0,2*K-1)])
    
    kk,zz = np.meshgrid(k,z)
    
    Tmp = np.exp(-2j*np.pi*kk*(zz-Z_0)/(Z_1-Z_0)) *  np.matlib.repmat(sig,2*K-1,1).T
    
    
    return np.trapz(Tmp,z,axis=0)