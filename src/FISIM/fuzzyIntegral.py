from __future__ import absolute_import
from builtins import zip
from builtins import range
from copy import *

import newton


def _func(x,l):
	"""Funcion q calcula el productorio de la integral de sugeno utilizado para calcular lambda.
	Se llama en el metodo de newton-raphson
	En x van los valores q le da el metodo de newton-raphson
	En lista (q es una lista) van los valores de las g's calculadas anteriormente"""
	aux = 1
	for i in range(len(l)):
		aux = aux * (1+float(l[i])*x)
	print(aux)
	return aux - (1+x)

		
def computeLambda(l):
	"""Llama a Newton-Raphson para calcular lambda. Lambda pertenece al intervalo ]-1,inf[ 
	Busca entre ]-1,0[ dp entre ]0,inf[ y finalmente si no la encuentra entonces entre [-0.0001,0.0001]"""
	aux = newton.newtonRaphson(_func, -1, -0.000001, 1.0e-5, l)
	if aux == None:
		aux = newton.newtonRaphson(_func, 0.000001, 10000, 1.0e-5, l)
	if aux == None:
		aux = 0
	print(aux)
	return aux



def computeG(importances1,l):
	"""Funcion q calcula la lambda medida difusa para un conjunto de valores con sus g's especificadas con una lambda dada. (Ver "Fuzzy Measure Theory" Wang & Klir (libro))
	Se trata de una funcion recursiva ya q aplica la lambda medida de sugeno. Se van cogiendo los dos primeros elementos y se les aplica la expresion para crear un nuevo termino con la g combinada de dichos elementos. Se procede asi hasta q se hayan considerado todos los elementos del conjunto
	dic1 guarda los valores con las g's especificadas
	l es la lambda"""
	
	#print "Llamada", importances1
	importances = copy(importances1)
	if len(importances) == 1:
		aux = importances[0]
		del importances[0]
		return aux
	if len(importances) == 0:
		return 0
	else:
		#nom1 = (dic.keys()[0])
		#nom2 = (dic.keys()[1])
		aux1 = importances[0]
		aux2 = importances[1]
		n = len(importances)
		del importances[0]
		del importances[0]
		#aux = foo(0)
		# expresion de la integral de sugeno
		aux = aux1 + aux2 + (aux1*aux2*l)
		importances.append(aux)
		return computeG(importances,l)


def fuzzyIntegral(importances,distances,l=None):
	"""Hay que ordenar los elementos por distances (considerar pasar distances e importances en la misma variable)
	Despues calcular las importacias para los subconjuntos que se forman (anadiendo cada vez el elemento de mayor similitud)
	Finalmente calcular el maximo de los minimos de las importancias y las distancias en cada posicion"""
	if l == None: l = computeLambda(importances)

	all = list(zip(distances,importances))
	all.sort(reverse=True)
	imp = [all[i][1] for i in range(len(all))]
	dist = [all[i][0] for i in range(len(all))]
	# Calcular las importancias llamando calcularG

	G = [computeG(imp[0:i+1],l) for i in range(len(imp))]
	#print "G =",G

	minimos = [min(dist[i],G[i]) for i in range(len(G))]
	#print "Minimos =",minimos
	ret = max(minimos)
	#print "Maximo =",ret
	return ret




