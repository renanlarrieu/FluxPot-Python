from NewtonRaphson_fluxpot import *

S=100
c=100

regua_Barras = ['NUM', 'TIPO', 'VOLT', 'PGER', 'QGER', 'PLOAD', 'QLOAD', 'SHUNT']
dbar_4b = np.array(
             [[  1 ,      1 ,   1000 , 400/S , -60.2/S ,      0 ,     0 ,    0.0 ],
              [  2 ,      0 ,   1011 ,   0.0 ,     0.0 ,      0 ,     0 ,    0.0 ],
              [  3 ,      0 ,   1007 ,   0.0 ,     0.0 ,  200/S , 100/S ,    0.0 ],
              [  4 ,      2 ,   1072 ,   0.0 ,     0.0 ,      0 ,     0 ,    0.0 ]])



regua_LT = ['DE', 'PARA',     'R',     'X', 'MVAR', 'TAP']
dlin_4b = np.array(
        [[ 1 ,     2 ,           0 ,    1.5/c ,   0 ,    0 ],
         [ 2 ,     3 ,  0.09/(2*c) ,1.4/(2*c) ,    0,    0 ],   
         [ 3 ,     4 ,       0     ,    1.5/c ,   0 ,    0 ]])
    

regua_Barras = ['NUM', 'TIPO', 'VOLT', 'PGER', 'QGER', 'PLOAD', 'QLOAD', 'SHUNT']
dbar_9b = np.array(
               [[  1 ,     2 ,  1075 ,   0.0 ,   0.0 ,      0 ,      0 ,    0.0 ],
                [  2 ,     1 ,  1075 ,  90/S ,   0.0 ,      0 ,      0 ,     0.0],
                [  3 ,     1 ,  1075 ,  85/S ,   0.0 ,      0 ,      0 ,    0.0 ],
                [  4 ,     0 ,  1000 ,   0.0 ,   0.0 ,      0 ,      0 ,    0.0 ],
                [  5 ,     0 ,  1000 ,   0.0 ,   0.0 ,  125/S ,   50/S ,    0.0 ],
                [  6 ,     0 ,  1000 ,   0.0 ,   0.0 ,   90/S ,   30/S ,    0.0 ],
                [  7 ,     0 ,  1000 ,   0.0 ,   0.0 ,      0 ,      0 ,    0.0 ],
                [  8 ,     0 ,  1000 ,   0.0 ,   0.0 ,  100/S ,   35/S ,    0.0 ],
                [  9 ,     0 ,  1000 ,   0.0 ,   0.0 ,      0 ,      0 ,    0.0 ]])


regua_LT = ['DE', 'PARA',     'R',     'X', 'MVAR', 'TAP']
dlin_9b = np.array(
        [[ 1 ,     4 ,         0 , 5.76/c ,     0 ,    0 ],
         [ 4 ,     5 ,       1/c ,  8.5/c ,17.6/c ,    0 ],
         [ 4 ,     6 ,     1.7/c ,  9.2/c ,15.8/c ,    0 ],
         [ 6 ,     9 ,     3.9/c ,   17/c ,35.8/c ,    0 ],
         [ 5 ,     7 ,     3.2/c , 16.1/c ,30.6/c ,    0 ],
         [ 7 ,     8 ,    0.85/c ,  7.2/c ,14.9/c ,    0 ],
         [ 8 ,     9 ,    1.19/c ,10.08/c ,20.9/c ,    0 ],
         [ 3 ,     9 ,         0 , 5.86/c ,     0 ,    0 ],
         [ 2 ,     7 ,         0 , 6.25/c ,     0 ,    0 ]])

tol = 1.0e-12

simula_4b=roda_NR(tol,dbar_4b,dlin_4b,'Convergência Método NR - Caso 4 barras','n')

simula_9b=roda_NR(tol,dbar_9b,dlin_9b,'Convergência Método NR - Caso 9 barras','n)