import streamlit as st
import numpy as np
import matplotlib.pyplot as plt



st.title('Acertijo de los prisioneros')

with st.sidebar:
    prisioneros = st.slider('Número de prisioneros', 50, 500, 100)
    intentos = st.slider('Número de cajas por abrir', 1, prisioneros, round(prisioneros / 2))
    S =  st.slider('Número de simulaciones', 100, 10000, 1000)

    

### Definir  funciones

def ciclo_busqueda(n_p, N) : 
    ciclo = []
    c = n_p
    while True : 
        if n_p != N[c - 1] : 
           ciclo.append( N[c - 1])
           c =  N[c - 1]
        else : 
            ciclo.append( N[c - 1]) 
            break
    return ciclo

def simular_oportunidad(n_intentos, T_prisioneros, N) : 
    exito = True
    R = []
    for i in range(1, T_prisioneros + 1) : 
        r = len(ciclo_busqueda(i, N))
        R.append(r)
        if r > n_intentos : 
            exito = False
    
    return exito, max(R)

@st.cache_data
def simulacion(S, n_intentos, T_prisioneros) : 
    exitos = 0
    lista_ciclos = []
    for s in range(S) :       
        numeros = np.arange(1, T_prisioneros + 1)
        np.random.shuffle(numeros)
        resultado, ciclo_mas_largo = simular_oportunidad(n_intentos, T_prisioneros, numeros)
        exitos += resultado
        lista_ciclos.append(ciclo_mas_largo)
    return exitos, lista_ciclos

exitos, lista_ciclos = simulacion(S, intentos, prisioneros) 

D  = {i : lista_ciclos.count(i) for i in range(1, prisioneros + 1)}




fig = plt.figure(figsize=(12,5))
plt.bar(x = D.keys(), height= [ x / S for x in  D.values()])
plt.xlabel('Ciclo más largo en cada simulación')
plt.ylabel('Probabilidad de observar el ciclo')




st.pyplot(fig)


probabilidades_exito = []
for numero_cajas_por_abrir in range(1, prisioneros + 1) : 
    suma = 0
    for i in range(1, numero_cajas_por_abrir + 1):
        suma += D[i]
    probabilidades_exito.append(suma /S )

fig1 = plt.figure(figsize=(12,5))
plt.plot(probabilidades_exito)
plt.vlines(intentos, 0, 1, 'red')
plt.title('Función de densidad de probabilidad de éxito respecto al número de cajas permitidas por abrir')
plt.xlabel('Número de cajas permitidas por abrir')
plt.ylabel('Probabilidad de éxito')
st.pyplot(fig1)

st.write(f'La probabilidad de que {prisioneros} prisioneros se salven abriendo {intentos} cajas es igual a {probabilidades_exito[intentos - 1]}')
