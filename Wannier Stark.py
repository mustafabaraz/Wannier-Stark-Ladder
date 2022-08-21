# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 22:48:19 2022

@author: MUSTAFA
"""

import matplotlib.pyplot as plt
import matplotlib 
import numpy as np

# Define the position axis:
    
N=3000 # Number of lattice points.
x_start, x_end=-20, 20 # Boundries of the position.
d_x=(x_end-x_start)/(N-1) # Position step size, recall that -1 is for the endpoints.
x=np.linspace(x_start, x_end, N) # We construct the x axis.

# Define the shape of the potential:
    
p_d= -4 # Well depth in eV.
p_w= 1 # Well width in nm.
b_w= 3 # Barrier width in nm.
b_d= 0 # Barrier height in eV.

# Define the electric field values:
    
E_f= [0, 0.01] # Electric field applied in V/nm.

# Write down the potential for 10 wells and 10 barriers with E-field is on:
    
V_list= [] # Potential list.
 
for i in E_f:
    V_x=np.piecewise(x,[(x>=-b_w/2) & (x<=b_w/2),(x>=-b_w/2-p_w) & (x<-b_w/2),
    (x>=-3/2*b_w-p_w) & (x<-b_w/2-p_w), (x>=-3/2*b_w-2*p_w) & (x< -3/2*b_w-p_w),
    (x>=-5/2*b_w-2*p_w) & (x<-3/2*b_w-2*p_w), (x>=-5/2*b_w-3*p_w) & (x< -5/2*b_w-2*p_w),
    (x>=-7/2*b_w-3*p_w) & (x< -5/2*b_w-3*p_w), (x>=-7/2*b_w-4*p_w) & (x< -7/2*b_w-3*p_w),
    (x>=-9/2*b_w-4*p_w) & (x<-7/2*b_w-4*p_w), (x>= -9/2*b_w-5*p_w) & (x<-9/2*b_w-4*p_w),
    (x<-9/2*b_w-5*p_w), (x>b_w/2) & (x<= b_w/2+p_w), (x>b_w/2+p_w) & (x<= 3/2*b_w+p_w),
    (x> 3/2*b_w+p_w) & (x<=3/2*b_w+2*p_w), (x>3/2*b_w+2*p_w) & (x<= 5/2*b_w+2*p_w),
    (x>5/2*b_w+2*p_w) & ( x<= 5/2*b_w+3*p_w), (x>5/2*b_w+3*p_w) & (x<= 7/2*b_w+3*p_w),
    (x>7/2*b_w+3*p_w) & (x<= 7/2*b_w+4*p_w), (x> 7/2*b_w+4*p_w) & (x<= 9/2*b_w+4*p_w),
    (x>9/2*b_w+4*p_w) & (x<= 9/2*b_w+5*p_w), (x>9/2*b_w+5*p_w)], [lambda x: b_d-i*x, lambda x:p_d-i*x,lambda x: b_d-i*x,lambda x: p_d-i*x,lambda x: b_d-i*x, lambda x:p_d-i*x,
    lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x,lambda x: p_d-i*x,lambda x: b_d-i*x,lambda x: p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x])
    
    V_list.append(V_x)                                                              
                   
# We have to construct the discrete hamiltonian matrix:
    
H_kin= np.zeros((N,N)) # For the kinetic matrix (which includes h_bar**2/m term):

for i in range(0, N):
    if i>0:
        H_kin[i-1, i]=-0.5/d_x**2 # Superdiagonal.
    H_kin[i, i]=1/d_x**2 # Main diagonal.
    if i+1<N:
        H_kin[i+1, i]=-0.5/d_x**2 # Subdiagonal.

H_pot=([np.identity(N),np.identity(N)]) # For the diagonal potential matrix, 2 matrices for 2 different e-field values.

for k in range (0,2):
    for i in range (0, N):
        H_pot[k][i][i]= V_list[k][i] # Potential matrix must be diagonal.
              
H_list= [] # For two different hamiltonian matrices for two different e-field values.

for i in range (0,2):
    H= H_pot[i]+H_kin
    H_list.append(H)

H_array= np.array(H_list) # Make it an array.

# We can find the eigenvalues by using numpy:
    
E_list= [] # An energy list for two different e-field values.

for i in range (0,2):
    E= np.linalg.eigvalsh(H_array[i]) 
    E_list.append(E)

# Since we have 10 finite wells, we can make an index list to plot them versus eigenvalues:
   
well_index=[1,2,3,4,5,6,7,8,9,10] # index 1 is the leftmost, index 10 is the rightmost well.

# We have to pick the negative (acceptable) eigenvalues and add them to a list since our potential and e-field values allow only negative energies:
    
E_neg_zero_field=[] # For E=0.
E_neg_on_field=[] # For E=0.01.

for i in range (0,N):
    if E_list[0][i]<0:
        E_neg_zero_field.append(E_list[0][i])

    
for i in range (0,N):
    if E_list[1][i]<0:
        E_neg_on_field.append(E_list[1][i])
 
# Notice that the well with index 1 is on the outer left, thus it is on the negative axis.
# This means that the index with 1 must have the biggest energy (since -Ex contributes a positive number).
# However, the eigenvalues list is in increasing order, thus we must sort to the decreasing order for the electric field on case.

E_neg_on_field.sort()
E_neg_on_field.reverse() 

# Then, plot the results:

# Below uses LaTeX fonts:
    
matplotlib.rc('font', family='serif')
matplotlib.rc('font', size=16)
matplotlib.rc('legend', fontsize=16)
matplotlib.rc('legend', numpoints=1)
matplotlib.rc('legend', handlelength=1.5)
matplotlib.rc('legend', frameon=False)
matplotlib.rc('xtick.major', pad=7)
matplotlib.rc('xtick.minor', pad=7)
matplotlib.rc('lines', lw=1.5)
matplotlib.rc('text', usetex=True)
matplotlib.rc('text.latex', 
              preamble=[r'\usepackage[T1]{fontenc}',
                        r'\usepackage{amsmath}',
                        r'\usepackage{txfonts}',
                        r'\usepackage{textcomp}'])
 
# Plotting the eigenvalues for two different electric field values:
    
plt.plot(well_index,E_neg_zero_field, 'kx',label='$E=$'+' '+str(E_f[0])+' '+'$V/nm$')
plt.plot(well_index,E_neg_on_field, 'rx', label='$E=$'+' '+str(E_f[1])+' '+'$V/nm$')
plt.xlabel('Well Index (Unitless)')
plt.ylabel('Energy (eV)')
plt.title('Well Index (Left to Right) vs Energy ($V_{0}=$'+' '+str(p_d)+' '+ '$eV$)' ) 
plt.legend()
plt.show() 

# Now, solve for the eigenvectors using numpy again:
    
E_vect_list=[] # Eigenvector list for two different values of e-field.

for i in range (0,2):
    E_list[i], E_vect= np.linalg.eigh(H_array[i]) # Remember that E_list[i] is still in increasing order (we did not change it) in energy. 
    E_vect_list.append(E_vect) # Add the eigenvectors to the new list.

# Plotting the ground state probability densities for 10 wells, which are degenerate first, then become non-degenerate after applying e-field:
    
for i in range (0,2):
    for n in range (0,10): # Notice that n=0 will give a localized wavefunction around x=18 nm since corresponding eigenvalue is the smallest. It makes sense because of the translational symmetry breaking that e-field brings.
        plt.plot(x, abs(E_vect_list[i][:, n]/np.sqrt(d_x))**2,'k',label='$n=0$ \ $(Ground \ State)$') # We normalize the wf via dividing by sqrt(dx).
        # We take the absolute value for the ground state since it must be of the same sign for all x.
        plt.gca().set_ylim(-0.8, 1.2) # For focusing on the plot.
        plt.xlabel('Position (nm)')
        plt.ylabel(r'$|\Psi_0(x)|^{2}$ (Probability/nm)')
        plt.title('Position vs Wavefunction ($V_{0}=$'+' '+str(p_d)+' '+ '$eV$, $E=$'+' '+str(E_f[i])+' $V/nm$) ') 
        plt.legend()
        plt.show()
        
# We can plot the first excited states' probability densities by simply decreasing the p_d value and modifying the code just a little bit in terms of list indexes:

# Write down the new potential (excited states) for 10 wells and 10 barriers with E-field is on:
    
V_list_new= [] # Potential list.

p_d= -15 # in eV.

for i in E_f:
    V_x=np.piecewise(x,[(x>=-b_w/2) & (x<=b_w/2),(x>=-b_w/2-p_w) & (x<-b_w/2),
    (x>=-3/2*b_w-p_w) & (x<-b_w/2-p_w), (x>=-3/2*b_w-2*p_w) & (x< -3/2*b_w-p_w),
    (x>=-5/2*b_w-2*p_w) & (x<-3/2*b_w-2*p_w), (x>=-5/2*b_w-3*p_w) & (x< -5/2*b_w-2*p_w),
    (x>=-7/2*b_w-3*p_w) & (x< -5/2*b_w-3*p_w), (x>=-7/2*b_w-4*p_w) & (x< -7/2*b_w-3*p_w),
    (x>=-9/2*b_w-4*p_w) & (x<-7/2*b_w-4*p_w), (x>= -9/2*b_w-5*p_w) & (x<-9/2*b_w-4*p_w),
    (x<-9/2*b_w-5*p_w), (x>b_w/2) & (x<= b_w/2+p_w), (x>b_w/2+p_w) & (x<= 3/2*b_w+p_w),
    (x> 3/2*b_w+p_w) & (x<=3/2*b_w+2*p_w), (x>3/2*b_w+2*p_w) & (x<= 5/2*b_w+2*p_w),
    (x>5/2*b_w+2*p_w) & ( x<= 5/2*b_w+3*p_w), (x>5/2*b_w+3*p_w) & (x<= 7/2*b_w+3*p_w),
    (x>7/2*b_w+3*p_w) & (x<= 7/2*b_w+4*p_w), (x> 7/2*b_w+4*p_w) & (x<= 9/2*b_w+4*p_w),
    (x>9/2*b_w+4*p_w) & (x<= 9/2*b_w+5*p_w), (x>9/2*b_w+5*p_w)], [lambda x: b_d-i*x, lambda x:p_d-i*x,lambda x: b_d-i*x,lambda x: p_d-i*x,lambda x: b_d-i*x, lambda x:p_d-i*x,
    lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x,lambda x: p_d-i*x,lambda x: b_d-i*x,lambda x: p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x, lambda x:p_d-i*x, lambda x:b_d-i*x])
    
    V_list_new.append(V_x)                                                              

# Define a new potential matrix:

H_pot_excited=([np.identity(N),np.identity(N)]) # For the diagonal potential matrix, 2 matrices for 2 different e-field values.

for k in range (0,2):
    for i in range (0, N):
        H_pot_excited[k][i][i]= V_list_new[k][i] # Potential matrix must be diagonal.

# Define a
H_excited_list= [] # For two different hamiltonian matrices for two different e-field values.

for i in range (0,2):
    H= H_pot_excited[i]+H_kin
    H_excited_list.append(H)

H_array_new= np.array(H_excited_list) # Make it an array again.

# We can find the eigenvalues by using numpy:
    
E_new_list= [] # An energy list for two different e-field values.

for i in range (0,2):
    E= np.linalg.eigvalsh(H_array_new[i]) 
    E_new_list.append(E)

E_neg_zero_field_excited=[] # For E=0.
E_neg_on_field_excited=[] # For E=0.01.

for i in range (0,N):
    if E_new_list[0][i]<0:
        E_neg_zero_field_excited.append(E_new_list[0][i])

    
for i in range (0,N):
    if E_new_list[1][i]<0:
        E_neg_on_field_excited.append(E_new_list[1][i])
 

# Reverse the list for the excited state with e-field on:
    
E_neg_on_field_excited.sort()
E_neg_on_field_excited.reverse() 

# Plotting the results of eigenvalues of two states: Ground state and first excited state.:

plt.plot(well_index,E_neg_zero_field_excited[0:10], 'kx',label='$E=$'+' '+str(E_f[0])+' '+'$V/nm$ (Ground State)')
plt.plot(well_index,E_neg_on_field_excited[10:20], 'k.', label='$E=$'+' '+str(E_f[1])+' '+'$V/nm$ (Ground State)')
plt.plot(well_index,E_neg_zero_field_excited[10:20], 'rx',label='$E=$'+' '+str(E_f[0])+' '+'$V/nm$ (First Excited State)')
plt.plot(well_index,E_neg_on_field_excited[0:10], 'r.', label='$E=$'+' '+str(E_f[1])+' '+'$V/nm$ (First Excited State)')
plt.xlabel('Well Index (Unitless)')
plt.ylabel('Energy (eV)')
plt.title('Well Index (Left to Right) vs Energy ($V_{0}=$'+' '+str(p_d)+' '+ '$eV$)' ) 
plt.legend()
plt.show() 

# Now, solve for the eigenvectors:
    
E_vect_list_excited=[] # Eigenvector list for two different values of e-field.

for i in range (0,2):
    E_new_list[i], E_vect= np.linalg.eigh(H_array_new[i]) # Remember that E_list[i] is still in increasing order (we did not change it) in energy. 
    E_vect_list_excited.append(E_vect) # Add the eigenvectors to the new list.

# Plotting the ground state probability density for 10 wells, which are degenerate first, then become non-degenerate after applying e-field:
    
for i in range (0,2):
    for n in range (10,20): # Notice that n=0 will give a localized wavefunction around x=18 nm since corresponding eigenvalue is the smallest. It makes sense because of the translational symmetry breaking that e-field brings.
        plt.plot(x, abs(E_vect_list_excited[i][:, n]/np.sqrt(d_x))**2,'k',label='$n=1$ \ $(First \ Excited \ State)$') # We normalize the wf via dividing by sqrt(dx).
        # We take the absolute value for the ground state since it must be of the same sign for all x.
        plt.gca().set_ylim(-0.8, 1.2) # For focusing on the plot.
        plt.xlabel('Position (nm)')
        plt.ylabel(r'$|\Psi_1(x)|^{2}$ (Probability/nm)')
        plt.title('Position vs Wavefunction  ($V_{0}=$'+' '+str(p_d)+' '+ '$eV$, $E=$'+' '+str(E_f[i])+' $V/nm$) ') 
        plt.legend()
        plt.show()
        
# One can obtain all of the figures by changing the plot variables.
# For example, erase abs**2 to obtain the wavefunction.

# Title of the Program: Wannier-Stark Ladder
# Name: Ahmet Mustafa Baraz
# ID: 21702127