#Random Medium Generator
#Von Karman Correlation Function
import numpy as np
import math
#import scipy
from scipy.stats import norm
#This program will attempt to produce a realization of a 3-D random medium
#using a Von Karman correlation function. 


write_file = 1    # 1 to write file, 0 to turn writing off
is_1D = 0         # 1 for 1D OK velocity strcture, 0 for 3D heterogeneity

#Set length and spacing of domain     
N1x = 501
N1y = 501
N1z = 1001

N2x = 501
N2y = 501
N2z = 1001

Nx = N1x + N2x - 1
Ny = N1y                       # Ny = N1y = N2y
Nz = N1z                       # Nz = N1z = N2z

refine = 1                             #refinement factor 
Nx = (Nx-1)*refine + 1                 #Number of interval, which is equal to number of FFT points
Ny = (Ny-1)*refine + 1                 #Number of interval, which is equal to number of FFT points
Nz = (Nz-1)*refine + 1                 #Number of interval, which is equal to number of FFT points

hx = 0.05/refine
hy = 0.05/refine
hz = 0.05/refine
PML = 0                        #1 for PML

Lx = (Nx-2)*hx
Ly = (Ny-2)*hy
Lz = (Nz-2)*hz
#Coordinates at sample points
#[X,Y,Z] = meshgrid(-Lx : hx : Lx-hx, -Ly : hy : Ly-hy, -Lz : hz : Lz-hz);

#Set Parameters used in the PSD Function
min_wavelength = 0.5
kmax =2.0*np.pi/min_wavelength      #Cutoff Wavenumber (1/km) -- min wavelength = 200m

#Parameters used to control PSDF statistics
v = 0.3       #Hurst Exponent
a = 0.15       # Correlation distance in km
p0 = 1.0/10    # 
m = p0*(a**v)           # standard deviation of perturbations

#This is NOT the sigma parameter in the Von Karman PSDF that controls
#perturbation strength. This is the standard deviation of the prescaled
#perturbations

sigma = 1.  
E = 3.             #Euclidean Dimension

#Set up Wavenumber Domain  

kNx = np.pi/hx    #Nyquist frequency
kNy = np.pi/hy    #Nyquist frequency
kNz = np.pi/hz    #Nyquist frequency

kx = np.zeros((Nx))
ky = np.zeros((Ny))
kz = np.zeros((Nz))

#kx[1:Nx/2+1] = 2.0*[0:Nx/2]/Nx;
for i in range(0, int((Nx-1)/2+2)):
    kx[i] = (2.0*i)/(Nx-1)

    
#kx(Nx:-1:Nx/2+2) = -kx(2:Nx/2);
k0 = 0
for i in range(int(Nx-1), int((Nx-1)/2 + 1), -1):
    k0 = k0 + 1
    kx[i] = -kx[k0]


#ky(1:(Ny-1)/2+2) = 2*[0:(Ny-1)/2+1]/Ny;
for i in range(int((Ny-1)/2+2)):
    ky[i] = (2.0*i)/(Ny-1)
    

#ky(Ny:-1:(Ny-1)/2+2) = -ky(2:(Ny-1)/2+1);
k0 = 0
for i in range(int(Ny-1), int((Ny-1)/2+1), -1):
    k0 = k0 + 1
    ky[i] = -ky[k0]


#kz(1:(Nz-1)/2+2) = 2*[0:(Nz-1)/2+1]/Nz;
for i in range(int((Nz-1)/2+2)):
    kz[i] = (2.0*i)/(Nz-1)
    

#kz(Nz:-1:(Nz-1)/2+2) = -kz(2:(Nz-1)/2+1);
k0 = 0
for i in range(int(Nz-1), int((Nz-1)/2+1), -1):
    k0 = k0 + 1
    kz[i] = -kz[k0]

kx = kNx*kx
ky = kNy*ky
kz = kNz*kz

k = np.zeros((Nz,Ny,Nx))

for i in range(Nx):
    for j in range(Ny):
        for l in range(Nz):
            k[l,j,i] = np.sqrt(kx[i]**2 + ky[j]**2 + kz[l]**2)



#%Generate a Random Number for Every point in the Spatial Domain
#seed = 179; % 8, 27flipr4 (quite good), 1r16 (quite good)
#s = RandStream.create('mt19937ar','seed',seed);
#RandStream.setGlobalStream(s);
# r = randn(Nx, Ny, Nz)

seed = 2239 #179
random_generator = np.random.RandomState(seed)
#r = np.reshape(random_generator.randn(Nx*Ny*Nz), (Nx, Ny, Nz))
#np.random.seed(seed)
#r = norm.ppf(np.random.rand(Nz,Ny,Nx))
r = random_generator.randn(Nz, Ny, Nx)

#Scale by a factor of sqrt(N/L) per dimension so PSDF has unit amplitude
r = (((Nx/Lx)*(Ny/Ly)*(Nz/Lz))**(1.0/2.0))*r

#Apply Discrete Fourier Transform and scale by factor of h per dimension
r = np.fft.fftn(r)*hx*hy*hz

#Calculate PSD and check for unit amplitude
PSD = (np.abs(r)**2)/(Lx*Ly*Lz)
print('unit white noise PSD = ', (np.mean(np.mean(np.mean(PSD)))))


#Calculate Von Karmen PSDF 
VK = ((2.0*np.sqrt(np.pi)*a)**E)*math.gamma(v+E/2.0)/(math.gamma(v)*(1.0+(a*k)**2.0)**(v+E/2.0))

# p0 = sigma/a**v

#Set up step function for frequency cutoff
for i in range(Nx):
    for j in range(Ny):
        for l in range(Nz):
            if k[l,j,i] >= kmax:
                VK[l,j,i] = 0.0



#%Multiply VonKarman PSDF * Transformed noise * Cutoff Function to filter
r = r*np.sqrt(VK)
#clear VK
r[0,0,0] = 0.     #%remove k=0 component

#%Inverse Fourier Transform
r = np.fft.ifftn(r)/(hx*hy*hz)
r = np.real(r)

#%Sigma Check
sigma_check = np.sqrt(np.abs((hx*hy*hz)*np.sum(np.sum(np.sum(r**E)))/(Lx*Ly*Lz)))
ratio = sigma/sigma_check

print('calculated sigma = ', sigma_check, 'ratio = ', ratio)

#%Use the altered r matrix to create fields for
cp = np.zeros((Nz,Ny,Nx))
cs = np.zeros((Nz,Ny,Nx))

cp_avg = 6.
cs_avg = cp_avg/np.sqrt(3.)


cp = cp_avg + cp
cs = cs_avg + cs

if is_1D == 0:
  cp = cp*(1+m*r)     #%Add on average wave speed
  cs = cs*(1+m*r)


rho = 1.6612*cp-0.4721*cp**2+0.0671*cp**3-0.0043*cp**4+0.000106*cp**5

cp_check = np.reshape(cp,(Nx*Ny*Nz,1))
std_norm = np.std(cp_check)/np.mean(cp_check)
std_ratio = std_norm/m

print('std_norm = ', std_norm, 'std_ratio = ', std_ratio)

# generate material parameters
mu = rho*cs**2
lam = (cp**2)*rho-2*mu


if write_file == 1:

    
    cp_B1 = cp[:,:,0:N1x] 
    cp_B2 = cp[:,:,N1x-1:Nx] 
    cs_B1 = cs[:,:,0:N1x] 
    cs_B2 = cs[:,:,N1x-1:Nx] 
    rho_B1 = rho[:,:,0:N1x] 
    rho_B2 = rho[:,:,N1x-1:Nx] 
    
    
    mu_B1 = rho_B1*cs_B1**2
    mu_B2 = rho_B2*cs_B2**2
    lambda_B1 = (cp_B1**2)*rho_B1-2*mu_B1
    lambda_B2 = (cp_B2**2)*rho_B2-2*mu_B2

   
   
    mu_B1_R = np.reshape(mu_B1,(N1x*N1y*N1z,1))
    mu_B2_R = np.reshape(mu_B2,(N2x*N2y*N2z,1))
    lam_B1_R = np.reshape(lambda_B1,(N1x*N1y*N1z,1))
    lam_B2_R = np.reshape(lambda_B2,(N2x*N2y*N2z,1))
    rho_B1_R = np.reshape(rho_B1,(N1x*N1y*N1z,1))
    rho_B2_R = np.reshape(rho_B2,(N2x*N2y*N2z,1))
   
    # write material parameters to file
    fid_mu1 = open('mu_B1'+'_v'+str(v)+'_a'+str(a)+'_p'+str(p0), 'w')
    fid_mu1.write(mu_B1_R)
    fid_mu1.close()

    fid_mu2 = open('mu_B2'+'_v'+str(v)+'_a'+str(a)+'_p'+str(p0), 'w')
    fid_mu2.write(mu_B2_R)
    fid_mu2.close()


    fid_lam1 = open('lam_B1'+'_v'+str(v)+'_a'+str(a)+'_p'+str(p0), 'w')
    fid_lam1.write(lam_B1_R)
    fid_lam1.close()

    fid_lam2 = open('lam_B2'+'_v'+str(v)+'_a'+str(a)+'_p'+str(p0), 'w')
    fid_lam2.write(lam_B2_R)
    fid_lam2.close()


    fid_rho1 = open('rho_B1'+'_v'+str(v)+'_a'+str(a)+'_p'+str(p0), 'w')
    fid_rho1.write(rho_B1_R)
    fid_rho1.close()

    fid_rho2 = open('rho_B2'+'_v'+str(v)+'_a'+str(a)+'_p'+str(p0), 'w')
    fid_rho2.write(rho_B2_R)
    fid_rho2.close()
    

# End of random material generator

