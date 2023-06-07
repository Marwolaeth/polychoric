from sympy import *

init_printing()
x, y, z, t = symbols('x y z t')
x = IndexedBase('x')
d = IndexedBase('d')
i, j, k, m, n = symbols('i j k m n', integer=True, cls=symbols)
f, g, h, P = symbols('f g h P', cls=Function)
tau, gamma = symbols('tau gamma', real=True)
D = symbols('D', integer=True)

L = Product(f(x[k])*P(d[k]), (k, 1, n))
L
log(L)

# Compute the logarithm of the product expression
l = log(L)

# Expand the logarithm using the log rule
expand_log(l, force=True)
