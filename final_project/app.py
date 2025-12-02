import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import math
import scipy as sp

#### Funciones y clases

def get_coefs_txt (txt) : 
    coef_s = txt.replace(' ', '')
    list_coefs =  coef_s.split(',')
    return list_coefs
def get_coefs_float(lista) : 
    return [float(x) for x in lista]
    
def escibir_polinomio (coefs) :
    polinomio_txt = '' 
    for a, i in zip(coefs, range(len(coefs))) : 
        if i == 0 : 
            if a == '0' : 
                continue
            else : 
                polinomio_txt += a
        elif i == 1 : 
            if a == '0' : 
                continue
            elif a[0] == '-' : 
                polinomio_txt += f'{a}x'
            else :
                polinomio_txt += f' + {a}x'
        else : 
            if a == '0' : 
                continue
            elif a[0] == '-' : 
                polinomio_txt += f'{a}x^{i}'
            else :
                polinomio_txt += f' + {a}x^{i}'

    return polinomio_txt


def coefs_dif(c1,c2) :
    n =  len(c1) - len(c2)
    if n == 0 :
        return np.array(c1) - np.array(c2)
    elif n > 0 : 
        zeros_add = np.zeros(n)
        c2_new = np.append(c2, zeros_add)
        return c1 - c2_new
    else :
        zeros_add = np.zeros(-n)
        c1_new = np.append(c1, zeros_add)
        return c1_new - c2


### Definir clases
class Polinomio (): 
    def __init__(self, C):
        self.coefs = C
    def evaluate(self, X): 
        coefs = self.coefs
        r = 0
        for c,i in zip(coefs, range(len(coefs))) : 
            r += c * (X ** i) 
        return r


class Transform() : 
    def __init__(self, a, b):
        self.a = a
        self.b = b
    def evaluate(self, x) : 
        T = (self.a + self.b ) / 2 + ((self.b - self.a ) / 2) * x
        return T
    

##3

st.title('Aproximación de centro de masa para regiones tipo I acotadas por por polinomios')
st.write('por: Daniel Vélez, Omar Neyra, Luisa Cabrera')

st.sidebar.subheader('Instrucciones')

st.sidebar.write(r'Ingresar sólo los coeficientes del polonomio. Ejemplo: para el polinomio $3x - 2^2 + x^4 $ '  )

st.sidebar.write(r'Escibir: $0, 3, -2, 0, 1$'  )
mostrar_centros = st.sidebar.toggle('Mostrar centros por región', value = False )

p1_txt = st.sidebar.text_input(r'Ingresar coeficientes  del primer polinomio ($P_1(x)$)')
p1_coefs_txt = get_coefs_txt(p1_txt)
st.sidebar.write( fr'$P_1(x) = {escibir_polinomio(p1_coefs_txt)}$ ')

p2_txt = st.sidebar.text_input(r'Ingresar coeficientes del  segundo polinomio ($P_2(x)$)')
p2_coefs_txt = get_coefs_txt(p2_txt)
st.sidebar.write( fr'$P_2(x) = {escibir_polinomio(p2_coefs_txt)}$ ')



### Trabjar y declarar los polonomios

a, b = st.slider(r'Intervalo $[a, b]$', value = [-1.0,1.0], min_value = -5.0, max_value= 5.0, step = 0.01)

X = np.linspace(-5, 5, 1000)
greater_a = X[X >= a]
interval_ab = greater_a[greater_a <= b]


p1_coefs = get_coefs_float(p1_coefs_txt)

P1 = Polinomio(p1_coefs)
p2_coefs = get_coefs_float(p2_coefs_txt)
P2 = Polinomio(p2_coefs)

P1x = P1.evaluate(X)
P2x = P2.evaluate(X)

P1ab = P1.evaluate(interval_ab)
P2ab = P2.evaluate(interval_ab)

polinomio_dif = coefs_dif(p1_coefs, p2_coefs )

### Encontrar las raíces
raices = np.roots(np.flip(polinomio_dif))
raices = np.real(raices[np.imag(raices) == 0])


## Crear los intervalos de intgeración

intervalo = np.insert(np.sort(raices), 0, a)
intervalo = np.insert(intervalo, len(intervalo) , b)

intervalo =  intervalo[ intervalo >= a]


intervalos = []
for i in range(len(intervalo) - 1) : 
    if intervalo[i] >= intervalo[i + 1] : 
        continue
    else : 
        intervalos.append((intervalo[i], intervalo[i + 1]))
    


def centro_masa (p1,p2,  intab) :
    T = Transform(intab[0], intab[1]) 
    n = max([len(p1.coefs), len(p2.coefs)]) - 1
    N = math.ceil((2*n + 1) / 2)
    A, W = np.polynomial.legendre.leggauss(N) 

    Area = 0
    AreaX = 0
    AreaY = 0

    for a_i, w_i in  zip(A,W) :
        Ta_i = T.evaluate(a_i)
        Area += w_i * (p1.evaluate(Ta_i) - p2.evaluate(Ta_i))
        AreaX += w_i *( Ta_i* (p1.evaluate(Ta_i) - p2.evaluate(Ta_i)) )
        AreaY += w_i *(p1.evaluate(Ta_i)**2/2 - p2.evaluate(Ta_i)**2/2 )


    Regla_A = Area * (intab[1] - intab[0])/2
    Regla_AreaX = AreaX * (intab[1] - intab[0])/2
    Regla_AreaY = AreaY * (intab[1] - intab[0])/2

    cX = Regla_AreaX /  Regla_A
    cY = Regla_AreaY /  Regla_A
    return (cX, cY), Regla_A 





def centro_masa_Scipy (p1,p2,  intab) :
    def P1_minus_P2(x) : 
        return p1.evaluate(x) - p2.evaluate(x)

    def x_P1_minus_P2(x) : 
        return x *(p1.evaluate(x) - p2.evaluate(x))

    def P1_sq_minus_P2_sq(x) : 
        return (1/2) *(p1.evaluate(x)**2 - p2.evaluate(x)**2)


    Area = sp.integrate.quad(P1_minus_P2, intab[0], intab[1])[0]
    AreaX = sp.integrate.quad(x_P1_minus_P2, intab[0], intab[1])[0]
    AreaY = sp.integrate.quad(P1_sq_minus_P2_sq, intab[0], intab[1])[0]

    cX = AreaX / Area
    cY = AreaY / Area
    return (cX, cY), Area



valores = np.concat((P1x, P2x)) 
valores_ab =  np.concat((P1ab, P2ab)) 




### Relizar todo el trabajo de integración:
fig, ax = plt.subplots( figsize = (15, 6))

ax.set_xlim(-5, 5)
ax.plot(X, P1x, alpha = 0.2, color = 'red', label = r'$P_1(x)$')
ax.plot(X, P2x, alpha = 0.2, color = 'blue', label = r'$P_2(x)$')
ax.fill_between(interval_ab, P1ab, P2ab, alpha = 0.7, label = 'Región', color = 'orange')

Centros = []
Regiones = []

Centros_Scipy = []
Regiones_Scipy = []

for intervalo in intervalos :
    C_i, region_i = centro_masa(P1, P2, intervalo)
    C_SP_i, region_SP_i = centro_masa_Scipy(P1, P2, intervalo)
    if mostrar_centros : 
        plt.scatter(C_i[0], C_i[1], color = 'blue', label = f'c{intervalos.index(intervalo) + 1}')

    Centros_Scipy.append(C_SP_i)
    Regiones_Scipy.append(region_SP_i)  
    Centros.append(C_i)
    Regiones.append(abs(region_i))



st.write(Regiones_Scipy)
st.write(Regiones)

Region_Total = sum(Regiones)

X = 0
Y = 0
for ci, Reg in zip(Centros, Regiones) : 
    X += ci[0] * Reg
    Y += ci[1] * Reg

CX = X / Region_Total
CY = Y / Region_Total



## Obtener las absiscas y pesos de integración Gaussiana

plt.scatter(CX,CY, s = 200, label = 'centro de masa', color = 'red', marker = '*')

ax.set_ylim(min(valores_ab) , max(valores_ab))

ax.legend()
ax.vlines(a, min(valores), max(valores), linestyles= '--', colors='black',  alpha = 0.5)
ax.vlines(b, min(valores), max(valores), linestyles= '--', colors= 'black', alpha = 0.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
st.pyplot(fig, use_container_width= True)
