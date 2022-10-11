from sympy import *
u,v,w,omega = symbols("u v w omega")
R, q, u_h = symbols("R q u_h")
V = sqrt(u**2+v**2+w**2)
# q = \rho\sigma A\frac{Cl_{\alpha}}{2}
Vtip = omega*R
mu = V/Vtip
V1 = sqrt(-0.5*V**2+sqrt((0.5*V**2)**2+u_h**4))
lambdait = V1/Vtip
theta0, thetatw = symbols("theta0 thetatw")
init_printing(use_unicode=True)
T0 = (q*Vtip**2)/(1+1.5*mu**2)
TTheta = (2/3-2/3*mu**2+1.5*mu**4)*theta0
TThetaTw = 0.5*thetatw*(1-1.5*mu**2+1.5*mu**4)
TLambda = lambdait*(1-0.5*mu**2)

T = T0*(TTheta+TThetaTw+TLambda)
#print(latex(diff(T0,u)))
# print(latex(diff(T0,v)))
# print(latex(diff(T0,w)))
# print(latex(diff(T0,omega)))

#print(latex(diff(TTheta,u)))
#print(latex(diff(TTheta,v)))
#print(latex(diff(TTheta,w)))
#print(latex(diff(TTheta,omega)))

#print(latex(diff(TLambda,u)))
#print(latex(diff(TLambda,v)))
#print(latex(diff(TLambda,w)))
#print(latex(diff(TLambda,omega)))

#print(latex(diff(TThetaTw,u)))
#print(latex(diff(TThetaTw,v)))
#print(latex(diff(TThetaTw,w)))
#print(latex(diff(TThetaTw,omega)))
T = Function('T')(u,v,w,omega)
L = Function('L')(u,v,w,omega)
E = symbols('E')
# q = \frac{1}{rho*A}
ct = T/(Vtip**2)*q
alphatpp = atan2(w,u)
lambdaic = 0.5*ct/sqrt((mu*cos(alphatpp))**2+(mu*sin(alphatpp)+L)**2)
V1 = lambdaic*Vtip
C = T*V1/omega + 0.75*R/E

# print(latex(diff(lambdaic,u)))
#print(latex(diff(lambdaic,v)))
#print(latex(diff(lambdaic,w)))
#print(latex(diff(lambdaic,omega)))


#print("\\begin{equation*}")
#print("\t\\pder{C}{u}=")
#print(latex(diff(C,u)))
#print("\\end{equation*}")
#
#print("\\begin{equation*}")
#print("\t\\pder{C}{v}=")
#print(latex(diff(C,v)))
#print("\\end{equation*}")
#
#print("\\begin{equation*}")
#print("\t\\pder{C}{w}=")
#print(latex(diff(C,w)))
#print("\\end{equation*}")
#
#print("\\begin{equation*}")
#print("\t\\pder{C}{\omega}=")
#print(latex(diff(C,omega)))
#print("\\end{equation*}")

Lambda_lambda, Lambda_t, Lambda_omega, Lambda_mu, Lambda_alpha = symbols("Lambda_lambda, Lambda_t, Lambda_omega, Lambda_mu, Lambda_alpha")
T = Function('T')(u,v,w,omega)
mu = Function('mu')(u,v,w,omega)
alpha = Function('alpha')(u,v,w,omega)
LAMBDA = (1/(1-Lambda_lambda))*(Lambda_t*T+Lambda_omega*omega+Lambda_mu*mu+Lambda_alpha*alpha)

print("\\begin{equation*}")
print("\t\\pder{\\lambda}{u}=")
print(latex(diff(LAMBDA,u)))
print("\\end{equation*}")
print("\\begin{equation*}")
print("\t\\pder{\\lambda}{v}=")
print(latex(diff(LAMBDA,v)))
print("\\end{equation*}")
print("\\begin{equation*}")
print("\t\\pder{\\lambda}{w}=")
print(latex(diff(LAMBDA,w)))
print("\\end{equation*}")
print("\\begin{equation*}")
print("\t\\pder{\\lambda}{\\omega}=")
print(latex(diff(LAMBDA,omega)))
print("\\end{equation*}")