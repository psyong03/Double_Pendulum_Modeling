import numpy as np

'''
n: # of variables
f: array of functions, f=np.array([f1, f2, ..., fn])
where each fi has independent variable t(scalar), x(vector)
t0: initial t
x0: initial value, x0=np.array([x1(t0), x2(t0), ..., xn(t0)])
t: the value of which we want to know x(t)=[x1(t), x2(t), ..., xn(t)]
step: # of steps
''' 
def RungeKutta4th(n, f, t0, x0, t, step):

    assert n == f.shape[0] == x0.shape[0]

    h = (t-t0)/step

    ti, xi = t0, x0
    a = np.zeros(n)

    for i in range(step):

        apply_vectorized_k1 = np.vectorize(lambda f, a: f(ti, xi))
        k1 = apply_vectorized_k1(f, a)

        apply_vectorized_k2 = np.vectorize(lambda f, a: f(ti + h/2, xi + (k1*h)/2))
        k2 = apply_vectorized_k2(f, a)

        apply_vectorized_k3 = np.vectorize(lambda f, a: f(ti + h/2, xi + (k2*h)/2))
        k3 = apply_vectorized_k3(f, a)

        apply_vectorized_k4 = np.vectorize(lambda f, a: f(ti + h, xi + k3*h))
        k4 = apply_vectorized_k4(f, a)

        ti += h
        xi += ((k1 + 2*k2 + 2*k3 + k4)*h)/6
    
    return ti, xi

'''
#a simple example

n = 3
f = np.array([lambda t, x: x[2],
              lambda t, x: x[0],
              lambda t, x: -3*x[2]+x[0]*2*np.sin(x[1])-4*x[1]+np.cos(t)])
t0 = 0
x0 = np.array([1., 0., 2.])
t = 1
h = 1000

print(RungeKutta4th(n, f, t0, x0, t, h)) # (1.0000000000000007, array([ 1.68378609,  1.50061361, -0.13909097]))
'''

m1, m2, l1, l2 = 1.0, 1.0, 1.0, 1.0
g = 9.8
n = 4
f1 = lambda t, x: x[2]
f2 = lambda t, x: x[3]
def f3(t, x):
    numer1 = -g*(2*m1+m2)*np.sin(x[0]) -m2*g*np.sin(x[0]-2*x[1])
    numer2 = -2*np.sin(x[0]-x[1])*m2*(x[3]**2*l2+x[2]**2*l1*np.cos(x[0]-x[1]))
    numer = numer1 + numer2
    denom = l1*(2*m1+m2-m2*np.cos(2*x[0]-2*x[1]))
    return numer/denom
def f4(t, x):
    numer1 = x[2]**2*l1*(m1+m2) + g*(m1+m2)*np.cos(x[0])
    numer2 = x[3]**2*l2*m2*np.cos(x[0]-x[1])
    numer = 2*np.sin(x[0]-x[1])*(numer1+numer2)
    denom = l2*(2*m1+m2-m2*np.cos(2*x[0]-2*x[1]))
    return numer/denom
f = np.array([f1, f2, f3, f4])
t0 = 0
x0 = np.array([1.0, 1.0, 0.0, 0.0])
t = 1
step = 1000
print(RungeKutta4th(n, f, t0, x0, t, step))