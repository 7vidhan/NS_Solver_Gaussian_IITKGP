from sympy import *
import numpy
from sympy.plotting import plot3d
import math
PI=3.1415
l=1
t1=60
t0=40
w=1
var('x y')
F = 0
print("n values:")
p=10
for i in range(0,p):
    n = 2*i + 1
    print(n)
    c = 4*l*(t1-t0)/(l*sinh(n*PI*w/l)*n*PI)
    F += ( c*sin(n*PI*x/l)*sinh(n*PI*y/l) )

print("plots values:")
p = None
for i in range(0,10):
    print(2*i)
    d= 2*i#difference between t and t0
    p2 = plot_implicit(F-d, (x,0,1), (y,0,1), show=False, line_color='b')
    if p:
        p.extend(p2)
    else:
        p = p2

p.show()