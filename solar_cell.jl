using Plots

h = 6.62607015e-34 # m^2kg/s
c = 2.99792458e8 # m/s
k = 1.380649e-23 # m^2kg/s^2k
T = 5780 # K
e = 1.602176634e-19

function EtoNu(x)
    x/h
end

function NutoX(nu)
    nu*h/k/T
end

function f(x)
    x^2/(exp(x)-1)
end

function integral(a,b,dx)
    ret = f(a)/2*dx;
    x = a+dx

    while x<b
        ret+=f(x)*dx
        x+=dx
    end

    ret+=f(x)/2*dx
end

E_min = 0.5 # eV
E_max = 3.0 # eV

E_min*=e # J
E_max*=e # J

graph_n = 100

hν_g = E_min:(E_max-E_min)/(graph_n-1):E_max
ν_g = EtoNu.(hν_g)
x_g = NutoX.(ν_g)

int_max = 15
dx = 0.001

int_x = zeros(graph_n)

for i in 1:graph_n
    int_x[i] = integral(x_g[i],int_max,dx)
end

y = zeros(graph_n)
C = h/k/T/pi^4*15*100
for i in 1:graph_n
    y[i] = C*ν_g[i]*int_x[i]
end

plot(ν_g,y)
