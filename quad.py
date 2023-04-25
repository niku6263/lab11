import numpy as np

def trap(f,low,high,n):
    
    h = (high - low) / n
    
    
    integ = f(low) + f(high)
    
    for i in range(1,n):
        k = low + i*h
        integ = integ + 2 * f(k)
    
    integ = integ * h/2
    
    return integ

def simpsons(f,a,b,M):
    x = np.linspace(a,b,M);
    h = (b-a)/(M-1);
    # Explain why this defines the weights for Simpsons
    w = (h/3)*np.ones(M);
    w[1:M:2]=4*w[1:M:2];
    w[2:M-1:2]=2*w[2:M-1:2];
    I_hat = np.sum(f(x)*w);
    return I_hat;

def adaptive_quad(f,a,b,M,tol):
    maxit = 50
    left_p = np.zeros((maxit,))
    right_p = np.zeros((maxit,))
    s = np.zeros((maxit,1))
    left_p[0] = a;
    right_p[0] = b;
    # initial approx and grid
    s[0],x,_ = trap(f,a,b,M);
    # save grid
    X = []
    X.append(x)
    j = 1;
    I = 0;
    nsplit = 1;
    while j < maxit:
        # get midpoint to split interval into left and right
        c = 0.5*(left_p[j-1]+right_p[j-1]);
        # compute integral on left and right spilt intervals
        s1,x,_ = trap(f,left_p[j-1],c,M); X.append(x)
        s2,x,_ = trap(f,c,right_p[j-1],M); X.append(x)
        if np.max(np.abs(s1+s2-s[j-1])) > tol:
            left_p[j] = left_p[j-1]
            right_p[j] = 0.5*(left_p[j-1]+right_p[j-1])
            s[j] = s1
            left_p[j-1] = 0.5*(left_p[j-1]+right_p[j-1])
            s[j-1] = s2
            j = j+1
            nsplit = nsplit+1
        else:
            I = I+s1+s2
            j = j-1
            if j == 0:
                j = maxit
    return I,np.unique(X),nsplit


f = lambda x: np.sin(1/x)


print("Using trapezoidal: ", trap(f,0.1,2,5))

print("Using Simpsons: ", simpsons(f, 0.1, 2, 5))

print(adaptive_quad(f,0.1,2,5,10**-2))