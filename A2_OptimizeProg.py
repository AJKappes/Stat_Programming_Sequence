
import numpy as np
import scipy.optimize as op
import sympy as sp
from numpy import linalg as lin
from scipy.misc import factorial

#########################################################################
#                     Creating benchmark results                        #
#########################################################################

# Objective is to minimize fx "fun" wrt x = x[0] and y = x[1]

# Define function

x, y = sp.symbols('x y')

fxy = 100*(y - x**2)**2 + (1 - x)**2
f_l = sp.lambdify((x, y), fxy)

# Construct function that accepts parameter values as arguments
def f(x):

    return f_l(x[0], x[1])


dx = sp.diff(fxy, x)                            # derivative -400(y-x^2) + 2(1-x)
dy = sp.diff(fxy, y)                            # derivative 200(y-x^2)
jac_dx_l = sp.lambdify((x, y), dx)
jac_dy_l = sp.lambdify((x, y), dy)

# Construct jacobian that accepts parameter values as arguments
def jacobian(x):

    return np.array([jac_dx_l(x[0], x[1]), jac_dy_l(x[0], x[1])])


# Using base function to create benchmark results

initial = np.array([2, 5])
benchres = op.minimize(f, initial, jac=jacobian)
print('Benchmark minimization values are {}' .format(benchres.x))print()
print(); print()
# Provides successful result of min (x, y) = (1, 1) as a benchmark

# # #########################################################################
# # #                         Newton's method                               #
# # #########################################################################

# Specifying values to start with
e = 10E-5*(1+abs(np.maximum(initial[0], initial[1])))
crit = 1

# Construct the Hessian
dxx = sp.lambdify((x, y), sp.diff(fxy, x, x))   # derivatives taken from FOCs above
dxy = sp.lambdify((x, y), sp.diff(fxy, x, y))
dyy = sp.lambdify((x, y), sp.diff(fxy, y, y))

# Allow Hessian to accept numeric arguments, applying Young's theorem for dxy = dyx
def hessian(x):

    return np.array([[dxx(x[0], x[1]), dxy(x[0], x[1])],
                     [dxy(x[0], x[1]), dyy(x[0], x[1])]])

# No need to transpose gradient, jacobian already in 2x1 dimension, only renaming jacobian to gradient
gradient = jacobian

# Defining the Newton Method loop
def newt_meth(initial, crit, e):

    if crit < e:

        initial = initial
        return initial

    else:

        H = hessian(initial)
        grd = gradient(initial)

        veci = initial - lin.inv(H).dot(grd)
        crit = np.maximum(abs(veci[0] - initial[0]), abs(veci[1] - initial[1]))
        e = 10E-5*(1 + abs(np.maximum(veci[0], veci[1])))
        return newt_meth(initial=veci, crit=crit, e=e)


# Output provides minimum values (1, 1)
op_vec = newt_meth(initial, crit, e)
print('Newton\'s method f(x, y) minimization values are {}' .format(op_vec))
print()

# Optimality test
grd_op = np.maximum(abs(gradient(op_vec)[0]), abs(gradient(op_vec)[1]))
gam = 10E-5*(1 + abs(f(op_vec)))

if grd_op < gam:

    print('Newton\'s method optimality conditions met: {} < {}' .format(round(grd_op, 4), gam))

else:

    print('Optimality conditions not met')

print(); print()


#########################################################################
#                   Newton's method w/ line search                      #
#########################################################################

# "loosely" minimize lambda approach
# specify lambda guesses
lam1, lam2, lam3 = [5, 4], [3, 2], [1, 5]
lamlist = np.linspace(0, 4, 16)

sk = lin.inv(hessian(initial)).dot(gradient(initial))

lamlist = np.linspace(.5, 3.5, 12)

# Generate a functional value list for each specified lambda guess
eval_f = []
for i in lamlist:

    val = initial + i*sk
    eval_f.append(f(val))

min_ind = np.argmin(lamlist)
lamk = eval_f[min_ind]

crit = 1
e = 10E-5*(1 + np.maximum(abs(initial[0]), abs(initial[1])))


def newt_lin_search(initial, crit, e):

    if crit < e:

        initial = initial
        return initial

    else:

        veci = initial + lamk*sk
        crit = np.maximum(abs(veci[0] - initial[0]), abs(veci[1] - initial[1]))
        e = 10E-5*(1 + np.maximum(abs(veci[0]), abs(veci[1])))
        return newt_meth(initial=veci, crit=crit, e=e)

# Optimal values of (1, 1) for minimum
op_vec_ln = newt_lin_search(initial, crit, e)
print('Newton\'s w/ line search f(x, y) minimization values are {}' .format(op_vec_ln))
print()

# Optimality test
grd_op_ls = np.maximum(abs(gradient(op_vec_ln)[0]), abs(gradient(op_vec_ln)[1]))
gam_lin = 10E-5*(1 + abs(f(op_vec_ln)))

if grd_op_ls < gam_lin:

    print('Newton\'s w/ line search optimality conditions met: {} < {}' .format(round(grd_op_ls, 4), gam_lin))

else:

    print('Optimality conditions not met')

print(); print()


#########################################################################
#                           BFGS Method                                 #
#########################################################################

# Parameters
e = 10E-5
gam = 10E-5

# Constructing educated hessian guess
dxx = 800*x**2 + 200*y + 1
dxy = 200
dyy = 200
dxx_arg = sp.lambdify((x, y), dxx)

# Initial H0 Hessian
def H0(x):

    return np.array([[dxx_arg(x[0], x[1]), dxy], [dxy, dyy]])

H0 = H0(initial)

# Value construction for iterations
lamspace = np.linspace(.5, 5.5, 20)

fval = []
for i in lamspace:

    evali = initial + i*sk
    fval.append(f(evali))

lamarg = np.argmin(lamspace)
lamk = fval[lamarg]

# Param vals and update hessian guess for optimizing iterations
crit = 1
e = e*(1 + np.maximum(abs(initial[1]), abs(initial[1])))

def bfgs(initial, crit, e, H0):

    if crit < e:

        initial = initial
        return initial

    else:
        sk = -lin.inv(H0).dot(gradient(initial))
        theta_t1 = initial + lamk * sk
        zt = theta_t1 - initial
        yt = gradient(theta_t1) - gradient(initial)
        zt_T = zt[np.newaxis].T
        yt_T = yt[np.newaxis].T
        zt_T.shape = (1, 2)
        zt.shape = (2, 1)
        yt.shape = (2, 1)
        yt_T.shape = (1, 2)

        H_t1 = H0 - np.divide(np.dot(H0, zt).dot(zt_T).dot(H0), np.dot(zt_T, H0).dot(zt)) + np.divide(np.dot(yt, yt_T), np.dot(yt_T, zt))

        crit = np.maximum(abs(initial[0] - theta_t1[0]), abs(initial[1] - theta_t1[1]))
        e = e*(1 + np.maximum(abs(initial[1]), abs(initial[1])))

        return bfgs(initial=theta_t1, crit=crit, e=e, H0=H_t1)

# convergence to (0, 0) due to poor hessian guess accuracy
opvec = bfgs(initial, crit, e, H0)
print('BFGS f(x, y) minimization values are {}, showing convergence to the wrong point' .format(opvec))
print()

# Optimality test will confirm convergence to nonoptimal point
grd_op_bfgs = np.maximum(abs(gradient(opvec)[0]), abs(gradient(opvec)[1]))
gam_bfgs = gam*(1 + f(opvec))

if grd_op_bfgs < gam_bfgs:

    print('BFGS optimality conditions met: {} < {}' .format(round(grd_op, 4), gam_bfgs))

else:

    print('BFGS optimality conditions not met - Convergence to the wrong point')



