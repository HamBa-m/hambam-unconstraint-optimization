from unidim import elimination as e
from unidim import interpolation as i
from functions import g,gp,gpp

e.usfs(g,-10,10,plot=True)
e.usas(g,-10,10,plot=True)
e.bs(g,-10,10,plot=True)
e.bf(g,-10,10,plot=True)
e.halv(g,-10,10,plot=True)
e.fibo(g,-10,10,plot=True)
e.golden(g,-10,10,plot=True)

i.Newton_Rapson(g,gp,gpp,2,plot=True,)
i.Quasi_Newton(g,2,plot=True)
i.Secant_Method(g,gp,1,plot=True)