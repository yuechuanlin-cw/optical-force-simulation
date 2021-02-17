"""
    This function is used to calculate the BSC based on integral localized approximation methods 
    Reference to: 
        Ren, K. F., G. Gréhan, and G. Gouesbet. "Evaluation of laser-sheet beam shape coefficients in generalized Lorenz–Mie theory by use of a localized approximation." JOSA A 11.7 (1994): 2072-2079.
    

"""
import sys 
import numpy as np 
import math 

i = complex(0.,1.)  # imaginary unit 
def ILA_BSC_func(k,n,m,Orders = 20,w0 = [2e-6,2e-6],PositionMatrix = [0,0,0]):
    if len(w0) == 1:
        w0x = w0[0]
        w0y = w0[0]
    elif len(w0) == 2:
        w0x = w0[0]
        w0y = w0[1]
    else:
        sys.exit(['Error: Wrong Input of beam waist parameters !!'])
    
    if len(PositionMatrix) == 3:
        x0 = PositionMatrix[0]
        y0 = PositionMatrix[1]
        z0 = PositionMatrix[2] 
    else:
        sys.exit(['Error: The number of PositionMatrix should be Three!'])
    
    s = 1/(k*np.amin(w0)) 
    r = (n + 0.5)/k
    th = np.pi/2 
    
    lx = k*(w0x**2)
    ly = k*(w0y**2)
    Qx = 1/(i+(2/lx)*(r*np.cos(th)-z0))
    Qy = 1/(i+(2/ly)*(r*np.cos(th)-z0))
    A = - (r**2)*((np.sin(th))**2)*(i*Qx/(w0x**2)-i*Qy/(w0y**2))/4
    B = r*np.sin(th)*(i*Qx*x0/(w0x**2)+Qy*y0/(w0y**2))
    C = r*np.sin(th)*(i*Qx*x0/(w0x**2)-Qy*y0/(w0y**2))

    psh0 = np.sqrt(-Qx*Qy)*np.exp(-i*(Qx/(w0x**2)+Qy/(w0y**2))*(r**2)*(np.sin(th)**2)/2)*np.exp(-i*(Qx*(x0**2)/(w0x**2)+Qy*(y0**2)/(w0y**2)))
    gTMn0 = 0
    gTEn0 = 0
    gTMnm = 0
    gTEnm = 0
    tmp1 = 0 
    tmp2 = 0
 

    if m == 0:
        Zn0 = 2*n*(n+1)*i/(2*n+1)        
        for j in range(Orders):
            for q in range(Orders):
                pup2 = int(((j+1)/2))+1
                for p in range(q+pup2):
                    tmp1 = (A**(p+q))*(B**j)*(C**(j+2*q-2*p+1))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+1))
                    tmp2 = (A**(p+q))*(B**(j+2*q-2*p+1))*(C**j)/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+1))

                    gTMn0 = gTMn0 + (tmp1+tmp2)
                    gTEn0 = gTEn0 + (tmp1-tmp2)
                    tmp1 = 0
                    tmp2 = 0   
        gTMn0 = 0.5*Zn0*np.exp(i*k*z0)*psh0*gTMn0 
        gTEn0 = 0.5*Zn0*(-i)*np.exp(i*k*z0)*psh0*gTEn0 
    elif m > 0:
        Znm = (-2*i/(2*n+1))**(np.abs(m)-1)
        for j in range(Orders):
            for q in range(Orders):
                pup2 = int(((m+j-1)/2))+1
                for p in range(q+pup2):
                    gTMnm = gTMnm + (A**(p+q))*(B**(2*q-2*p+j+m-1))*(C**(j))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+m-1))
                    gTEnm = gTEnm + (A**(p+q))*(B**(2*q-2*p+j+m-1))*(C**(j))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+m-1))
        
        for j in range(Orders):
            for q in range(Orders):
                pup2 = int(((m+j+1)/2))+1
                for p in range(q+pup2):
                    gTMnm = gTMnm + (A**(p+q))*(B**(2*q-2*p+j+m+1))*(C**(j))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+m+1))
                    gTEnm = gTEnm - (A**(p+q))*(B**(2*q-2*p+j+m+1))*(C**(j))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+m+1))
    
        gTMnm = 0.5*Znm*np.exp(i*k*z0)*psh0*gTMnm 
        gTEnm = 0.5*Znm*np.exp(i*k*z0)*psh0*(-i)*gTEnm 

    else: # m<0
        Znm = (-2*i/(2*n+1))**(np.abs(m)-1)
        for j in range(Orders):
            for q in range(Orders):
                pup2 = int(((np.abs(m)+j-1)/2))+1
                for p in range(q+pup2):
                    gTMnm = gTMnm + (A**(p+q)*(B**j)*(C**(2*q-2*p+j+np.abs(m)-1)))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+np.abs(m)-1))
                    gTEnm = gTEnm + (A**(p+q)*(B**j)*(C**(2*q-2*p+j+np.abs(m)-1)))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+np.abs(m)-1))

        for j in range(Orders):
            for q in range(Orders):
                pup2= int(((np.abs(m)+j+1)/2))+1
                for p in range(q+pup2):
                    gTMnm = gTMnm + (A**(p+q)*(B**j)*(C**(2*q-2*p+j+np.abs(m)+1)))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+np.abs(m)+1))
                    gTEnm = gTEnm - (A**(p+q)*(B**j)*(C**(2*q-2*p+j+np.abs(m)+1)))/(np.math.factorial(p)*np.math.factorial(q)*np.math.factorial(j)*np.math.factorial(2*q-2*p+j+np.abs(m)+1))
        gTMnm = 0.5*Znm*np.exp(i*k*z0)*psh0*gTMnm
        gTEnm = 0.5*Znm*np.exp(i*k*z0)*psh0*i*gTEnm  
    if m == 0:
        return gTMn0, gTEn0
    else:
        return gTMnm, gTEnm  