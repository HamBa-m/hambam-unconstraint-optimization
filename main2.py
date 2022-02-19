from multivar import gradient as gd
from functions import h,H

gd.gradient_descent(h, [60,50],plot=True, F=H, axlim=100)
gd.conjugate_gradient(h, [60,50],plot=True, F=H, axlim=100)
gd.newton_descent(h, [60,50],plot=True, F=H, axlim=100)

