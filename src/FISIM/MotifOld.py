from __future__ import division
from __future__ import print_function
from __future__ import absolute_import


# RECORDAR que la de medida de Fisher hace uso de rpy!!!


from future import standard_library
standard_library.install_aliases()
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import string
import types
from subprocess import *
from math import *
from random import *

import mosta
import stats
import tomtom
from Dirichlet import *
from numpy import *

from . import fuzzyIntegral
from . import utils


class Motif(object):
	"""Clase que representa un motivo como una PSSM.
Atributos:
	mat -> matriz (clase array) con la matriz del PSSM"""
	

	
	def __init__(self, rand=False,matrix=None, n=10, numSamples = 30, fileIn=None,fileType=None, consensus=None, name = '', refSeq = None, FSM = None,IC=None,ID=None):
		"""Inicializa un PSSM a partir de una array de shape (n,4) dado.
Si no se pasa la matriz se genera un PSSM aleatorio de longitud dada. Por defecto la longitud es 10
Parametros: 
	rand -> boolean, si es True se genera el PSSM aleatorio de longitud n
	matrix -> matriz (clase array) con la matriz del PSSM
	n -> longitud del PSSM aleatorio que se desea generar. Por defecto n = 10
	numSamples -> numero de ejemplos utilizados para la generacion del motivo. Por defecto numSamples = 30
	fileIn -> fichero abierto en modo lectura de donde se cargara el PSSM en cuestion. El formato del fichero sera el standar de TRANSFAC o JASPAR. Seria necesario especificar el numSamples si se va a trabajar con ese dato y FSM == False. Se tomara el valor que haya en la tercera columna de al lado del nombre si existe (formato que me mando Javi! Estudiarlo y ponerlo como en TRANSFAC!!!), si no se dejara al valor por defecto. Si FSM == True, numSamples se calcula automaticamente.
	fileType -> Indica el formato del fichero en cuestion. Por ahora (14-3-2008) implemento para el valor 'jaspar'
	consensus -> String de [ACGT] a partir del cual se genera el Motivo.
	name -> nombre del motivo. String
	refSeq -> codigo refSeq del motivo
	FSM -> flag que indica si se trabaja con PSM (frecuencias ralativas, FSM == False) o con FSM (frecuencias absolutas, i.e. numeros enteros, FSM = True)
Requerimientos:
	ID -> Codigo de representacion del motivo. Lo pongo xq en el transfac y en el jaspar ademas del nombre del binding site tienen un codigo del tipo MA0001 para identificar la matriz en cuestion
	Si se facilitan mas de uno de (rand,matrix,fileIn,consensus) se aplicara el primero de ellos ignorando el resto, es decir rand tiene preferencia sobre matrix y matrix sobre fileIn, etc.
	Si todos son False o None se producira un error y el PSM no se creara
	"""
		self.name = name
		self.ID = ID
		self.refSeq = refSeq
		self.exp = True  # Flag que indica si la matriz esta en forma exponencial (True) o logaritmica (False)
		self.FSM = FSM
		self.n = n	
		self.numSamples = numSamples
		self.matFSM = None
		self.matPSSM = None
		if rand:
			mat = zeros(4 * n,float)
			mat.shape = (n,4)
			cont = [0,1,2,3]		
			for i in range(n):
				shuffle(cont)
				if not FSM:
					mat[i,cont[0]] = random.random()
					for j in range(1,3):
						mat[i,cont[j]] = uniform(0,1-sum(mat[i]))
					mat[i,cont[3]] = 1-sum(mat[i])
				else:
					mat[i,cont[0]] = random.randint(0,self.numSamples)
					for j in range(1,3):
						mat[i,cont[j]] = random.randint(0,self.numSamples-sum(mat[i]))
					mat[i,cont[3]] = self.numSamples-sum(mat[i])

				if FSM:
					self.matFSM = array(mat,int)
					self.calculatePSSMFromFSM()
				else:
					self.matPSSM = array(mat,float)
					self.calculateFSMFromPSSM(numSamples)
		
		elif(matrix != None):
			if self.FSM == None:
				self.FSM = utils.decideFSMorNot(matrix)

			
			if self.FSM == True:
				self.numSamples = sum(matrix[0,:])
				self.matFSM = array(matrix,int)
				self.n = matrix.shape[0]
				self.calculatePSSMFromFSM()
			elif self.FSM == False:
				self.matPSSM = array(matrix,float)
				self.n = matrix.shape[0]
				self.calculateFSMFromPSSM(numSamples)

		elif(fileIn != None):
			self.readFromFile(fileIn,fileType=fileType) # 6-3-2008 Pongo lo de especial para que me coja el nombre bien del fichero de jaspar que me he bajado. Es decir ignore el codigo de la matriz y ponga una concatenacion del nombre y la familia (segundo y tecer campo de la primera fila de cada motivo)
		elif(consensus != None):
			self.calculateMotifFromConsensus(consensus)
		else:
			print("Wrong usage for PSSM.__init__() -> mat, rand, fileIn and consensus cannot be all false")

	def calculateMotifFromConsensus(self,consensus):
		self.n = len(consensus)
		mat = zeros(self.n*4, float)
		mat.shape = (self.n, 4)
		for i in range(len(consensus)):
			if consensus[i] == 'A':
				mat[i][0] = 1.0
			elif consensus[i] == 'C':
				mat[i][1] = 1.0
			elif consensus[i] == 'G':
				mat[i][2] = 1.0
			elif consensus[i] == 'T':
				mat[i][3] = 1.0
			elif consensus[i] == 'N':
				mat[i][0] = 0.25
				mat[i][1] = 0.25
				mat[i][2] = 0.25
				mat[i][3] = 0.25
			else:
				print("Wrong value for consensus at position",str(i))
				print("Found", consensus[i], "expected [ACGTN]")
				sys.exit(-1)
		self.matFSM = array(mat,int)
		self.matPSSM = array(mat,float)	

	def calculatePSSMFromFSM(self):
		self.matPSSM = array(self.matFSM,float)
		self.matPSSM /= self.numSamples

	def calculateFSMFromPSSM(self,numSamples=None):
		if numSamples == None:
			numSamples = self.numSamples

		self.numSamples = numSamples
		self.matFSM = array(self.matPSSM*self.numSamples).astype(int)
		for i in range(self.n):
			s = sum(self.matFSM[i,:])
			ind = self.matFSM[i,:].tolist().index(max(self.matFSM[i,:]))
			self.matFSM[i,ind] += self.numSamples - s
		


	def pseudoCountFSM(self):
		if not all(self.matFSM):
			self.matFSM += 1
			self.numSamples += 4
			self.calculatePSSMFromFSM()

	def pseudoCountPSSM(self):
		if not all(self.matPSSM):
			self.matPSSM += 0.000001
			self.calculateFSMFromPSSM()

						 

	def logForm(self):
		"""Transforma la matriz de representacion del PSSM en logaritmica para temas de precision
Modifica:
	self.matFSM: Forma logaritimica
	self.matPSSM: Forma logaritimica
	sel.exp: Pasa a ser False para indicar que la matriz esta en forma logaritmica"""
		
		if(self.exp):
			self.pseudoCount()
			self.matPSSM = log(4.0*self.matPSSM)
			self.calculateFSMFromPSSM()

		self.exp = False

	def expForm(self):
		"""Transforma la matriz de representacion del PSSM y FSM en exponencial (la que tenia originalmente) para deshacer el cambio logaritmico.
Modifica:
	self.matFSM: Forma exponencial
	self.matPSSM: Forma exponencial
	sel.exp: Pasa a ser True para indicar q			for i in range(self.matFSM.shape[0]):
				for j in range(self.matFSM.shape[1]-1):
					fileOut.write(str(self.matFSM[i,j])+"\t")
				fileOut.write(str(self.matFSM[i,self.matFSM.shape[1]-1])+"\n")ue la matriz NO esta en forma logaritmica"""
		if(self.exp == False):
			self.matPSSM[:] = exp(self.matPSSM[:]) / 4.0
			self.calculateFSMFromPSSM()
		self.exp = True
		
	def calculateRevMatrixPSSM(self):
		"""Devuelve un array de shape (n,4) con la representacion matricial del reverse opposite del PSSM original"""
		ret = array(zeros(4*self.matPSSM.shape[0],float))
		ret.shape = (self.matPSSM.shape[0],4)
		for i in range(self.matPSSM.shape[0]):
				for j in range(4):
					ret[i,j] =  self.matPSSM[self.matPSSM.shape[0]-1-i,3-j];

		return ret;

	def calculateRevMatrixFSM(self):
		"""Devuelve un array de shape (n,4) con la representacion matricial del reverse opposite del FSM original"""
		ret = array(zeros(4*self.matFSM.shape[0],float))
		ret.shape = (self.matFSM.shape[0],4)
		for i in range(self.matFSM.shape[0]):
				for j in range(4):
					ret[i,j] =  self.matFSM[self.matFSM.shape[0]-1-i,3-j];

		return ret;


	def ICarray(self,start=None,end=None):
		"""Devuelve en un array (Numeric module) el Information Content de cada posicion del PSSM.
Parametros:
	start: Desde que posicion del PSSM queremos calcular el IC (debe estar entre 0 y la longitud del PSSM-1). Si no se especifica se empieza por la primera posicion (0)
	end: Hasta que posicion del PSSM queremos calcular el IC (debe estar entre start+1 y la longitud del PSSM). Si no se especifica se calcula hasta la ultima posicion (longitud del PSSM)
Esta implementado en forma vectorial para mejorar la eficiencia pero la formula que implementa es la del IC tipica"""
		if end and (end > self.n or end <= start):
			print("Motif.py -> ICarray -> Wrong value for end paramater. It should be in ["+ str(start+1)+","+str(self.n)+"]. Found: "+str(end))
			sys.exit(-1)

		if not start: start = 0 # Si no se especifica el comienzo se empieza por el 0
		if (start > self.n-1 or start < 0):
			print("Motif.py -> ICarray -> Wrong value for start paramater. It should be in [0,"+str(self.n-1)+"]. Found: "+str(start))
			sys.exit(-1)

		self.pseudoCountPSSM()
		#print self.matPSSM
		
		matAux = self.matPSSM * (old_div(log(self.matPSSM),log(2)))
		if end: matAux = matAux[start:end]
		else: matAux = matAux[start:]
		ic = zeros(len(matAux))
		for i,v in enumerate(matAux):
			ic[i] = 2 + sum(v)
		
		return ic


	def IC(self,start=None,end=None):
		"""Devuelve el valor de la suma del Information Content las posiciones del PSSM.
Esta implementado en forma vectorial para mejorar la eficiencia pero la formula que implementa es la del IC tipica
Parametros:
	start: Desde que posicion del PSSM queremos calcular el IC (debe estar entre 0 y la longitud del PSSM-1). Si no se especifica se empieza por la primera posicion (0)
	end: Hasta que posicion del PSSM queremos calcular el IC (debe estar entre start+1 y la longitud del PSSM). Si no se especifica se calcula hasta la ultima posicion (longitud del PSSM)
Esta implementado en forma vectorial para mejorar la eficiencia pero la formula que implementa es la del IC tipica"""
		if end and (end > self.n or end <= start):
			print("Motif.py -> ICarray -> Wrong value for end paramater. It should be in ["+ str(start+1)+","+str(self.n)+"]. Found: "+str(end))
			sys.exit(-1)

		if not start: start = 0 # Si no se especifica el comienzo se empieza por el 0
		if (start > self.n-1 or start < 0):
			print("Motif.py -> ICarray -> Wrong value for start paramater. It should be in [0,"+str(self.n-1)+"]. Found: "+str(start))
			sys.exit(-1)

		self.pseudoCountPSSM()
		#print self.matPSSM
		
		matAux = self.matPSSM * (old_div(log(self.matPSSM),log(2)))
		if end: matAux = matAux[start:end]
		else: matAux = matAux[start:]
		ic = zeros(len(matAux))
		for i,v in enumerate(matAux):
			ic[i] = 2 + sum(v)
		
		return sum(ic)
		
		return sum(ic)
		
	def insertMotif(self,motif, numSamples=30, start = -1):
		"""Inserta un motivo a partir de una dada del PSSM original.
Si no se le facilita la posicion se genera aleatoriamente.
Si no cabe el motivo a partir de la posicion dada se devuelve -1 y acaba.
Si todo va bien, substituye en la matriz del PSSM los valores originales por los generados a partir de hacer la media de numSamples muestras de la distribucion de Dirichlet con los parametros indicados por cada fila en motif
Parametros:
	motif -> PSSM con los pesos para cada base en cada posicion del motivo
	numSamples -> numeros de muestras que se realizaran para el computo de los valores a insertar
	start -> posicion a partir de la cual se inserta el nuevo motivo. Si es -1 se elige aleatoriamente

Modifica: 
	self.matPSSM convenientemente
	self.matFSM convenientemente
Devuelve:
	-1 si el motif a insertar es mas largo que el el original
	0 si funciona correctamente"""

		if(motif.matPSSM.shape[0] > self.matPSSM.shape[0]):
			print("Error. PSSM.insertMotif: El motivo a insertar es mas largo que el PSSM")
			return -1
		elif start == -1:
			pos = arange(0,self.matPSSM.shape[0]-motif.matPSSM.shape[0]+1,1)
			start = random.choice(pos)
			#print "Start = ", start

		else:
			if (self.matPSSM.shape[0] - start < motif.matPSSM.shape[0]):
				print("Error. PSSM.insertMotif: El motivo no se puede insertar a partir de la posicion " + str(start) + ". No caberia")
				return -1
		
		#self.pseudoCount()
		
		motif.pseudoCountPSSM()
		
		insert = motif.createDirichletSample(numSamp=numSamples)

		self.matPSSM[start:start + motif.matPSSM.shape[0],:] = insert[:,:]

		self.calculateFSMFromPSSM()

		return 0
		
			

	def createDirichletSample(self,numSamp=30):
		"""Realiza numSamp muestras de la distribucion de Dirichlet subyacente a la matriz de representacion del PSSM y devuelve la media de cada posicion en un array de shape (n,4).
Hace uso del modulo Dirichlet.py
Cuanto mayor sea numSamp, mas se asemejara la matriz devuelta al PSSM"""
		tempList =[]
		for base in self.matPSSM:
			#print base
			dist = Dirichlet(base)
			currDist = []
			for i in range(numSamp):
				currDist.append(dist.sample())

			tempList.append(average(currDist,axis=0))
			#print "anadiendo", tempList
		#print tempList
		return array(tempList)
 
	def consensus(self):
		"""Devuelve en un string la secuencia de consenso (ACGT) del PSSM. Si una posicion no tiene ninguna base con una frecuencia mayor que 0.35 se le asigna la letra N."""
		letters = 'ACGT'
		ret = []
		for i in range(self.matPSSM.shape[0]):
			maximum = max(self.matPSSM[i])
			if (maximum > 0.35): 
				pos = self.matPSSM[i].tolist().index(maximum)
				ret.append(letters[pos])
			else:
				ret.append('N')
		
		return string.join(ret,'')


	def distanceFuzzy(self, curPSSM, bigger=False):
		"""Calcula la distancia difusa entre dos PSSM. La distancia esta definida basandose en la que aparece en Fuzzy Logic in Medicine and Bioinformatics, A Torres, JJ Nieto - Journal of Biomedicine and Biotechnology. Se computa esa medida para cada columna y la distancia total sera la suma de las columnas
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#print "Correlation of " , self.consensus() , " and " , curPSSM.consensus()
		
		
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceFuzzy(self,True)
				
		#print self.matPSSM
		#print curPSSM.matPSSM
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		minDist = 9999.0
		
		while (curTime < times):
			curDist = 0.0

			for i in range(minSize):
				num = 0.0
				den = 0.0
				for j in range(4):
					num += abs(float(self.matPSSM[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(self.matPSSM[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
					#print num,den
				curDist += old_div(num, den)
                       
			#curDist = num / den
			#print curDist
			#raw_input()

			
			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist	
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0


		while (curTime < times):
			curDist = 0.0
		
			for i in range(minSize):
				num = 0.0
				den = 0.0
				for j in range(4):
					num += abs(float(revMatrix[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(revMatrix[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
           			curDist += old_div(num, den)
			#curDist = num / den
			#print curDist
			#raw_input()

			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist
				bestRev = True
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return minDist, startPosition, bestRev

	def distanceFuzzyWeighted(self, curPSSM, bigger=False):
		"""Calcula la distancia difusa entre dos PSSM. Es como la distanceFuzzyNew que aplica la distancia definida en Fuzzy Logic in Medicine and Bioinformatics, A Torres, JJ Nieto - Journal of Biomedicine and Biotechnology, pero Multiplica la distancia de cada posicion de las columnas por el maximo de los dos valores para dar mas importancia a las distancias de las posiciones mas consevadas
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#print "Correlation of " , self.consensus() , " and " , curPSSM.consensus()
		
		
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceFuzzyWeighted(self,True)
				
		#print self.matPSSM
		#print curPSSM.matPSSM
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		minDist = 9999.0
		
		while (curTime < times):
			curDist = 0.0


			for i in range(minSize):
				num = 0.0
				den = 0.0
				for j in range(4):
					num += abs(float(self.matPSSM[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))*max(float(self.matPSSM[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(self.matPSSM[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
				curDist += old_div(num, den)

			
			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist	
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0


		while (curTime < times):
			curDist = 0.0


			for i in range(minSize):
				num = 0.0
				den = 0.0
				for j in range(4):
					num += abs(float(revMatrix[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))*max(float(revMatrix[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(revMatrix[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
				curDist += old_div(num, den)

			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist
				bestRev = True
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return minDist, startPosition, bestRev

	def distanceIntegralFuzzy(self, curPSSM, bigger=False):
		"""Calcula la integral difusa entre dos PSSM. Se basa en el algoritmo propuesto en "Fuzy Measure Theory pag 332"... Seguir comentando!
Devuelve: (minDist, startPosition, bestRev)
PROGRAMAR!!"""
		
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceIntegralFuzzy(self,True)
				
		#print self.matPSSM
		#print curPSSM.matPSSM
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxDist = -1
		
		while (curTime < times):
			curDist = 0.0


			for i in range(minSize):
				#print self.matPSSM
				#print curPSSM.matPSSM
				importances = [max(value) for value in zip(self.matPSSM[i],curPSSM.matPSSM[i+curTime])]
				distances = [1-abs(value[0]-value[1]) for value in zip(self.matPSSM[i],curPSSM.matPSSM[i+curTime])]
				#print "Importances["+str(i)+"] ", importances
				#print "Distances["+str(i)+"] ", distances

				curDist += fuzzyIntegral.calculaIntegralDifusa(importances, distances)
				#print curDist
				#raw_input()
			if(curDist > maxDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist	
				maxDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0


		while (curTime < times):
			curDist = 0.0
			for i in range(minSize):
				#print self.matPSSM
				#print curPSSM.matPSSM
				importances = [max(value) for value in zip(revMatrix[i],curPSSM.matPSSM[i+curTime])]
				distances = [1-abs(value[0]-value[1]) for value in zip(revMatrix[i],curPSSM.matPSSM[i+curTime])] # En realidad se trata de similitudes!
				#print importances
				#print distances
				#raw_input()
				curDist += fuzzyIntegral.calculaIntegralDifusa(importances, distances)

			if(curDist > maxDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist
				bestRev = True
				maxDist = curDist
				startPosition = curTime
		
			curTime += 1
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
		#print "minsize=",minSize
		return old_div(maxDist,minSize), startPosition, bestRev

	def distanceIntegralChoquet(self, curPSSM, bigger=False):
		"""Calcula la integral difusa de CHOQUET entre dos PSSM. ... Seguir comentando!
Devuelve: (minDist, startPosition, bestRev)
PROGRAMAR!!"""
		
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceIntegralChoquet(self,True)
				
		#print self.matPSSM
		#print curPSSM.matPSSM
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxDist = -1
		
		while (curTime < times):
			curDist = 0.0


			for i in range(minSize):
				#print self.matPSSM
				#print curPSSM.matPSSM
				importances = [max(value) for value in zip(self.matPSSM[i],curPSSM.matPSSM[i+curTime])]
				distances = [1-abs(value[0]-value[1]) for value in zip(self.matPSSM[i],curPSSM.matPSSM[i+curTime])]
				#print importances
				#print distances
				#raw_input()
				curDist += fuzzyIntegral.calculaIntegralChoquet(importances, distances)

			
			if(curDist > maxDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist	
				maxDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0


		while (curTime < times):
			curDist = 0.0
			for i in range(minSize):
				#print self.matPSSM
				#print curPSSM.matPSSM
				importances = [max(value) for value in zip(revMatrix[i],curPSSM.matPSSM[i+curTime])]
				distances = [1-abs(value[0]-value[1]) for value in zip(revMatrix[i],curPSSM.matPSSM[i+curTime])]
				#print importances
				#print distances
				#raw_input()
				curDist += fuzzyIntegral.calculaIntegralChoquet(importances, distances)

			if(curDist > maxDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist
				bestRev = True
				maxDist = curDist
				startPosition = curTime
		
			curTime += 1
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return old_div(maxDist,minSize), startPosition, bestRev

			

	def distanceFuzzyOriginal(self, curPSSM, bigger=False):
		"""Calcula la distancia difusa entre dos PSSM. La distancia esta definida como en Fuzzy Logic in Medicine and Bioinformatics, A Torres, JJ Nieto - Journal of Biomedicine and Biotechnology. La diferencia con la que yo propongo es que hace da todas las columnas un solo array en vez de considerarlas separadamente
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#print "Correlation of " , self.consensus() , " and " , curPSSM.consensus()
		
		
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceFuzzy(self,True)
				
		#print self.matPSSM
		#print curPSSM.matPSSM
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		minDist = 9999.0
		
		while (curTime < times):
			curDist = 0.0
			num = 0.0
			den = 0.0
			for i in range(minSize):

				for j in range(4):
					num += abs(float(self.matPSSM[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(self.matPSSM[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
					#print num,den

                       
			curDist = old_div(num, den)
			#print curDist
			#raw_input()

			
			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist	
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0


		while (curTime < times):
			curDist = 0.0
			num = 0.0
			den = 0.0
			for i in range(minSize):
				for j in range(4):
					num += abs(float(revMatrix[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(revMatrix[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
           	
			curDist = old_div(num, den)
			#print curDist
			#raw_input()

			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist
				bestRev = True
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
			
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return minDist, startPosition, bestRev

	

	def distanceTomtom(self, curPSSM, dirTomtom = '/home/fernan/articuloFuzzyMotifs/relatedSoftware/metameme_3.4/bin/'):
		"""Calcula la distancia (Evalue) segun el programa TOMTOM entre dos PSSM. Cuanto mas pequena sea la distancia (Evalue), mas cercanos estan los PSSM. La distancia esta definida con los valores por defecto (distancia euclidea) del software descargado del supplementary material del paper de Gupta et al: Quantifying similarity between motifs, Genome Biology 2007, 8:R24.
	Internamente se crea un fichero temporal con los PSSM a comparar y se rescatan los valores de la distancia (Evalue) de los ficheros que genera el programa. Finalmente se borran los ficheros generados.
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
	
		minDist, startPosition, bestRev = tomtom.distanceTomtom(self,curPSSM,dirTomtom=dirTomtom)
		return minDist, startPosition, bestRev


	def distanceMosta(self, curPSSM, dirMosta = '/home/fernan/articuloFuzzyMotifs/relatedSoftware/mosta_src/'):
		"""Documentar fijandose en la distanceTomtom!!!
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
	
		minDist, startPosition, bestRev = mosta.distanceMosta(self,curPSSM,dirMosta=dirMosta)
		return minDist, startPosition, bestRev

	def distanceMostaMax(self, curPSSM, dirMosta = '/home/fernan/articuloFuzzyMotifs/relatedSoftware/mosta_src/'):
		"""Documentar fijandose en la distanceTomtom!!!
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
	
		minDist, startPosition, bestRev = mosta.distanceMostaMax(self,curPSSM,dirMosta=dirMosta)
		return minDist, startPosition, bestRev

			
	def distanceChoquet(self, curPSSM, bigger=False):
		T1 = self.matPSSM[0]
		T2 = curPSSM.matPSSM[0]
		dif = [math.sqrt((i[0]-i[1])*(i[0]-i[1])) for i in zip(T1,T2)]
		print("Dist Euclidea = ", sum(dif))
		veces = 1
		for i in range(len(dif)):
			if dif.count(dif[i]) > 1:
				dif[i] += 0.000001 * veces
				veces += 1

		d = dict.fromkeys(dif)
		k = 0
		for i,j in zip(T1,T2):
			d[dif[k]] = max(i,j)
			k += 1
		dif.append(0.0)
		dif.sort(reverse=True)
		choquet = 0
		for i in range(len(dif)-1):
			choquet += d[dif[i]] * (dif[i] - dif[i+1])
		print(T1)
		print(T2)
		print(dif)
		print(d)
		return choquet
		
	def distanceChoquetSuma(self, curPSSM, bigger=False):
		choquet = 0
		for indice in range(len(self.matPSSM)):
			T1 = self.matPSSM[indice]
			T2 = curPSSM.matPSSM[indice]
			dif = [math.sqrt((i[0]-i[1])*(i[0]-i[1])) for i in zip(T1,T2)]
			print("Dist Euclidea = ", sum(dif))
			veces = 1
			for i in range(len(dif)):
				if dif.count(dif[i]) > 1:
					dif[i] += 0.000001 * veces
					veces += 1

			d = dict.fromkeys(dif)
			k = 0
			for i,j in zip(T1,T2):
				d[dif[k]] = max(i,j)
				k += 1
			dif.append(0.0)
			dif.sort(reverse=True)
			choqTemp = 0
			for i in range(len(dif)-1):
				choqTemp += d[dif[i]] * (dif[i] - dif[i+1])
				choquet += d[dif[i]] * (dif[i] - dif[i+1])
			print(T1)
			print(T2)
			#print dif
			#print d
			print("Distancia Choquet temporal = ", choqTemp, "Acumulada = ", choquet) 
			
			
		return choquet
			
	def distanceChoquetLineal(self, curPSSM, bigger=False):
		dim = self.matPSSM.shape[0]
		
		self.matPSSM.shape = (4*dim,1)
		curPSSM.matPSSM.shape = (4*dim,1)
		
		T1 = self.matPSSM[:]
		T2 = curPSSM.matPSSM[:]
		dif = [math.sqrt((i[0]-i[1])*(i[0]-i[1])) for i in zip(T1,T2)]
		print("Dist Euclidea = ", sum(dif))
		veces = 1
		for i in range(len(dif)):
			if dif.count(dif[i]) > 1:
				dif[i] += 0.000001 * veces
				veces += 1

		d = dict.fromkeys(dif)
		k = 0
		for i,j in zip(T1,T2):
			d[dif[k]] = max(i,j)
			k += 1
		dif.append(0.0)
		dif.sort(reverse=True)
		choquet = 0
		for i in range(len(dif)-1):
			choquet += d[dif[i]] * (dif[i] - dif[i+1])
		print("T1 = ", T1)
		print(T2)
		print(dif)
		print(d)
		
		self.matPSSM.shape = (dim,4)
		curPSSM.matPSSM.shape = (dim,4)
		return choquet
		
	def distanceChoquetLineal2(self, curPSSM, bigger=False):
		
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceChoquetLineal2(self,True)
				
		#print self.matPSSM
		#print curPSSM.matPSSM
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		minDist = 9999.0
		
		dimSelf = self.matPSSM.shape[0]
		dimCur = curPSSM.matPSSM.shape[0]
		
		self.matPSSM.shape = (4*dimSelf,1)
		curPSSM.matPSSM.shape = (4*dimCur,1)
		
		while (curTime < times):
			curDist = 0.0
		
		
			
		
			T1 = self.matPSSM[:]
			T2 = curPSSM.matPSSM[curTime*4:(curTime+dimSelf)*4]
			#print T1.shape, T2.shape, curTime*4, (curTime+dimSelf)*4
			#raw_input()
			dif = [math.sqrt((i[0]-i[1])*(i[0]-i[1])) for i in zip(T1,T2)]
			#print "Dist Euclidea = ", sum(dif)
			veces = 1
			for i in range(len(dif)):
				if dif.count(dif[i]) > 1:
					dif[i] += 0.000001 * veces
					veces += 1

			d = dict.fromkeys(dif)
			k = 0
			for i,j in zip(T1,T2):
				d[dif[k]] = max(i,j)
				k += 1
			dif.append(0.0)
			dif.sort(reverse=True)
			#choquet = 0
			for i in range(len(dif)-1):
				curDist += d[dif[i]] * (dif[i] - dif[i+1])
			#print "T1 = ", T1
			#print "T2 = ", T2
			#print "dif = ", dif
			#print "d = ", d
			
			if(curDist[0] < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist[0]
				minDist = curDist[0]
				startPosition = curTime
			
			curTime += 1
		
		self.matPSSM.shape = (dimSelf,4)
		curPSSM.matPSSM.shape = (dimCur,4)
		
		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0
		
		self.matPSSM.shape = (4*dimSelf,1)
		curPSSM.matPSSM.shape = (4*dimCur,1)
		
		print("Probando con el revOpp")
		
		while (curTime < times):
			curDist = 0.0
			
			T1 = self.matPSSM[:]
			T2 = curPSSM.matPSSM[curTime*4:(curTime+dimSelf)*4]
			#print T1.shape, T2.shape, curTime*4, (curTime+dimSelf)*4
			#raw_input()
			dif = [math.sqrt((i[0]-i[1])*(i[0]-i[1])) for i in zip(T1,T2)]
			#print "Dist Euclidea = ", sum(dif)
			veces = 1
			for i in range(len(dif)):
				if dif.count(dif[i]) > 1:
					dif[i] += 0.000001 * veces
					veces += 1

			d = dict.fromkeys(dif)
			k = 0
			for i,j in zip(T1,T2):
				d[dif[k]] = max(i,j)
				k += 1
			dif.append(0.0)
			dif.sort(reverse=True)
			#choquet = 0
			for i in range(len(dif)-1):
				curDist += d[dif[i]] * (dif[i] - dif[i+1])
			#print "T1 = ", T1
			#print "T2 = ", T2
			#print "dif = ", dif
			#print "d = ", d
			
			if(curDist[0] < minDist):
				
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist[0]
				bestRev = True
				minDist = curDist[0]
				startPosition = curTime
			
			curTime += 1
		
		self.matPSSM.shape = (dimSelf,4)
		curPSSM.matPSSM.shape = (dimCur,4)
		
		return minDist, startPosition, bestRev

	
		
		
	# SEGUIR!!! 
	# He cambiado la forma de calcular la medida (suma de columnas)... Discutirla!
	def distanceFuzzyNew(self, curPSSM, bigger=False):
		"""Calcula la distancia difusa entre dos PSSM. La distancia esta definida como en Fuzzy Logic in Medicine and Bioinformatics, A Torres, JJ Nieto - Journal of Biomedicine and Biotechnology.
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#print "Correlation of " , self.consensus() , " and " , curPSSM.consensus()
		
		
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceFuzzyNew(self,True)
				
		#print self.matPSSM
		#print curPSSM.matPSSM
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		minDist = 9999.0
		
		while (curTime < times):
			curDist = 0.0
			num = 0.0
			den = 0.0

			for i in range(minSize):
				num = 0.0
				den = 0.0
				for j in range(4):
					num += abs(float(self.matPSSM[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(self.matPSSM[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
				curDist += old_div(num, den)
				#print num, den, curDist
                       
			#curDist = num / den

			
			if(curDist < minDist):
				
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0


		while (curTime < times):
			curDist = 0.0
			num = 0.0
			den = 0.0

			for i in range(minSize):
				num = 0.0
				den = 0.0
				for j in range(4):
					num += abs(float(revMatrix[i][j]) - float(curPSSM.matPSSM[i+curTime,j]))
					den += max(float(revMatrix[i][j]),float(curPSSM.matPSSM[i+curTime,j]))
				curDist += old_div(num, den)
				#print num, den, curDist
			#curDist = num / den

			if(curDist < minDist):
				bestRev = True
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		

		
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return minDist, startPosition, bestRev

	def	distanceCorrelation(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceCorrelation(self,True)
				
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxCorr = -9999.0
				
		while (curTime < times):
			curCorr = 0.0
			curX = 0; # X = self (its length will be always <= curPSMM)
			curY = 0; # Y = curPSSM
	
			for i in range(minSize):
				
				sumX = float(sum(self.matPSSM[i,:]))
				sumY = float(sum(curPSSM.matPSSM[i+curTime,:]))
				avX = sumX / 4.0
				avY = sumY / 4.0

				
				numerator = 0.0
				denominator1 = 0.0
				denominator2 = 0.0
				for j in range(4):
					T1 = ((old_div(self.matPSSM[i,j],sumX)) - avX)
					T2 = ((old_div(curPSSM.matPSSM[i+curTime,j],sumY)) - avY)
					numerator += T1 * T2
					denominator1 += pow(T1,2)
					denominator2 += pow(T2,2)

				curCorr += old_div(numerator, sqrt((denominator1 * denominator2))) 

			
			if(curCorr > maxCorr):
				maxCorr = curCorr
				startPosition = curTime
			
			curTime += 1
		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0

		while (curTime < times):
			curCorr = 0.0
			curX = 0 # X = self (its length will be always <= curPSMM)
			curY = 0 # Y = curPSSM
	
			for i in range(minSize):


				sumX = float(sum(revMatrix[i,:]))
				sumY = float(sum(curPSSM.matPSSM[i+curTime,:]))
				avX = sumX / 4.0
				avY = sumY / 4.0
				
				numerator = 0.0
				denominator1 = 0.0
				denominator2 = 0.0
				for j in range(4):
					T1 = ((old_div(revMatrix[i,j],sumX)) - avX)
					T2 = ((old_div(curPSSM.matPSSM[i+curTime,j],sumY)) - avY)
					numerator += T1 * T2
					denominator1 += pow(T1,2)
					denominator2 += pow(T2,2)
				
				curCorr += old_div(numerator, sqrt((denominator1 * denominator2)))
			
			
			if(curCorr > maxCorr):
				bestRev = True
				maxCorr = curCorr
				startPosition = curTime
			
			curTime += 1
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return maxCorr, startPosition, bestRev

	def	distanceALLR(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceALLR(self,True)
				

		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxALLR = -9999.0
				
		while (curTime < times):
			curALLR = 0.0
	
			for i in range(minSize):
				
				numerator1 = 0.0
				numerator2 = 0.0
				denominator = 0.0

				for j in range(4):
					numerator1 += (float(self.matPSSM[i,j]) * math.log((float(curPSSM.matPSSM[i+curTime,j])/float(sum(curPSSM.matPSSM[i+curTime,:])))/0.25))
					numerator2 += float(curPSSM.matPSSM[i+curTime,j]) * math.log((float(self.matPSSM[i,j])/float(sum(self.matPSSM[i,:])))/0.25)
					denominator += float(self.matPSSM[i,j] + curPSSM.matPSSM[i+curTime,j])

				curALLR += old_div((numerator1 + numerator2), denominator)
			
			
			if(curALLR > maxALLR):
				maxALLR = curALLR
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0

		while (curTime < times):
			curALLR = 0.0
	
			for i in range(minSize):
				
				numerator1 = 0.0
				numerator2 = 0.0
				denominator = 0.0

				for j in range(4):
					numerator1 += (float(revMatrix[i,j]) * math.log((float(curPSSM.matPSSM[i+curTime,j])/float(sum(curPSSM.matPSSM[i+curTime,:])))/0.25))
					numerator2 += (float(curPSSM.matPSSM[i+curTime,j]) * math.log((float(revMatrix[i,j])/float(sum(revMatrix[i,:])))/0.25))
					denominator += float(revMatrix[i,j] + curPSSM.matPSSM[i+curTime,j])

				
				curALLR += old_div((numerator1 + numerator2), denominator)
			
			
			if(curALLR > maxALLR):
				bestRev = True
				maxALLR = curALLR
				startPosition = curTime
			
			curTime += 1
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return maxALLR, startPosition, bestRev

	def	distanceKLD(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceKLD(self,True)
				

		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		minKLD = 9999.0
				
		while (curTime < times):
			curKLD = 0.0
	
			for i in range(minSize):
				
				numerator1 = 0.0
				numerator2 = 0.0
				try:
					for j in range(4):
						sumX = float(sum(self.matPSSM[i,:]))
						sumY = float(sum(curPSSM.matPSSM[i+curTime,:]))
						probX = old_div(self.matPSSM[i,j], sumX)
						probY = old_div(curPSSM.matPSSM[i+curTime,j], sumY)
						numerator1 += probX * math.log(old_div(probX, probY))
						numerator2 += probY * math.log(old_div(probY, probX))

					curKLD += (numerator1 + numerator2)
				except:
					curKLD = 9999.0
			
			
			if(curKLD < minKLD):
				minKLD = curKLD
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0

		while (curTime < times):
			curKLD = 0.0
	
			for i in range(minSize):
				
				numerator1 = 0.0
				numerator2 = 0.0
				try:
					for j in range(4):
						sumX = float(sum(revMatrix[i,:]))
						sumY = float(sum(curPSSM.matPSSM[i+curTime,:]))
						probX = old_div(revMatrix[i,j], sumX)
						probY = old_div(curPSSM.matPSSM[i+curTime,j], sumY)
						numerator1 += probX * math.log(old_div(probX, probY))
						numerator2 += probY * math.log(old_div(probY, probX))

					curKLD += (numerator1 + numerator2)
				except:
					curKLD = 9999.0
			
			
			if(curKLD < minKLD):
				bestRev = True
				minKLD = curKLD
				startPosition = curTime
			
			curTime += 1
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return minKLD/2.0, startPosition, bestRev # La division por 2.0 es para implementar la formula del KLD tal y como viene en la literatura pero en realidad no afecta a nuestros propositos

	def	distanceEuclidean(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceEuclidean(self,True)
				

		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		minDist = 9999.0
				
		while (curTime < times):
			curDist = 0.0
	
			for i in range(minSize):
				T1 = 0.0
				for j in range(4):
					sumX = float(sum(self.matPSSM[i,:]))
					sumY = float(sum(curPSSM.matPSSM[i+curTime,:]))
					probX = old_div(self.matPSSM[i,j], sumX)
					probY = old_div(curPSSM.matPSSM[i+curTime,j], sumY)
					T1 += pow(probX-probY,2)

				curDist += sqrt(T1)
			
			
			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist
				minDist = curDist
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0

		while (curTime < times):
			curDist = 0.0
	
			for i in range(minSize):
				T1 = 0.0
				for j in range(4):
					sumX = float(sum(revMatrix[i,:]))
					sumY = float(sum(curPSSM.matPSSM[i+curTime,:]))
					probX = old_div(revMatrix[i,j], sumX)
					probY = old_div(curPSSM.matPSSM[i+curTime,j], sumY)
					T1 += pow(probX-probY,2)

				curDist += sqrt(T1)
			
			
			if(curDist < minDist):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curDist
				bestRev = True
				minDist = curDist
				startPosition = curTime
			
			curTime += 1
		
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return minDist, startPosition, bestRev

	def	distanceChi2(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: Float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#self.pseudoCount()
		#curPSSM.pseudoCount()
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matFSM.shape[0] <= curPSSM.matFSM.shape[0]):
			minSize = self.matFSM.shape[0]
			times = curPSSM.matFSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceChi2(self,True)
		#print self.matFSM
		#print curPSSM.matFSM
		#print self.numSamples, curPSSM.numSamples	

		#raw_input()
		#print selfFSM.matrix
		#print curFSM.matrix
		#print selfFSM.numSamples, curFSM.numSamples	
		#raw_input()
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxChi = 0.0
				
		while (curTime < times):
			curChi = 0
			window = []
			for i in range(minSize):	
				chiX = 0.0
				chiY = 0.0
				x = array(self.matFSM[i,:],float)
				y = array(curPSSM.matFSM[i+curTime,:],float)
				xy = x+y
				Nx = sum(x)
				Ny = sum(y)
				N = sum(xy)
				Nex = old_div(Nx*xy, N)
				Ney = old_div(Ny*xy, N)
				for index,v in enumerate(Nex):
					if v == 0: Nex[index] = 0.0000001
				for index,v in enumerate(Ney):
					if v == 0: Ney[index] = 0.0000001

				chiX += sum(old_div(power(x-Nex,2), Nex))
				chiY += sum(old_div(power(y-Ney,2), Ney))
				
				window.append(stats.lchisqprob(chiX+chiY,3))
				
				#DEBUG
				#print x,y,xy
				#print Nx,Ny,N
				#print Nex,Ney
				#print chiX,chiY
				#print sum(chiX + chiY)
				#raw_input()
				#print stats.lchisqprob(sum(chiX + chiY),3)
				#raw_input()
				#END DEBUG
			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				maxChi = curChi
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixFSM()
		curTime = 0

		while (curTime < times):
			curChi = 0.0
			window = []
			for i in range(minSize):
				chiX = 0.0
				chiY = 0.0
				x = revMatrix[i,:]
				y = curPSSM.matFSM[i+curTime,:]
				xy = x+y
				Nx = sum(x)
				Ny = sum(y)
				N = sum(xy)
				Nex = old_div(Nx*xy, N)
				Ney = old_div(Ny*xy, N)
				for index,v in enumerate(Nex):
					if v == 0: Nex[index] = 0.0000001
				for index,v in enumerate(Ney):
					if v == 0: Ney[index] = 0.0000001
				
				chiX += sum(old_div(power(x-Nex,2), Nex))
				chiY += sum(old_div(power(y-Ney,2), Ney))
				window.append(stats.lchisqprob(chiX+chiY,3))

			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				bestRev = True
				maxChi = curChi
				startPosition = curTime


			
			curTime += 1
		
		#print maxChi
		
		#print str(maxChi)
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return maxChi, startPosition, bestRev

	def	distanceChi2Float(self, curPSSM, bigger=False):
		"""Calcula la el valor de la Chi2 entre dos PSSM. Se diferencia de la Chi2 normal en que utiliza las PSSM en lugar de los FSM, esto es trabaja con las frecuencias (valores reales entre 0 y 1) en lugar de con las cuentas (valores enteros con el numero de ocurrencias de cada nucleotido).
Devuelve: (maxChi2, startPosition, bestRev)
	maxChi2: Float con el mejor pvalue entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia optima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#self.pseudoCount()
		#curPSSM.pseudoCount()
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matPSSM.shape[0] <= curPSSM.matPSSM.shape[0]):
			minSize = self.matPSSM.shape[0]
			times = curPSSM.matPSSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceChi2(self,True)
		#print self.matFSM
		#print curPSSM.matFSM
		#print self.numSamples, curPSSM.numSamples	

		#raw_input()
		#print selfFSM.matrix
		#print curFSM.matrix
		#print selfFSM.numSamples, curFSM.numSamples	
		#raw_input()
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxChi = 0.0
				
		while (curTime < times):
			curChi = 0
			window = []
			for i in range(minSize):	
				chiX = 0.0
				chiY = 0.0
				x = array(self.matPSSM[i,:],float)
				y = array(curPSSM.matPSSM[i+curTime,:],float)
				xy = x+y
				Nx = sum(x)
				Ny = sum(y)
				N = sum(xy)
				Nex = old_div(Nx*xy, N)
				Ney = old_div(Ny*xy, N)
				for index,v in enumerate(Nex):
					if v == 0: Nex[index] = 0.0000001
				for index,v in enumerate(Ney):
					if v == 0: Ney[index] = 0.0000001

				chiX += sum(old_div(power(x-Nex,2), Nex))
				chiY += sum(old_div(power(y-Ney,2), Ney))
				
				window.append(stats.lchisqprob(chiX+chiY,3))
				
				#DEBUG
				#print x,y,xy
				#print Nx,Ny,N
				#print Nex,Ney
				#print chiX,chiY
				#print sum(chiX + chiY)
				#raw_input()
				#print stats.lchisqprob(sum(chiX + chiY),3)
				#raw_input()
				#END DEBUG
			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				maxChi = curChi
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0

		while (curTime < times):
			curChi = 0.0
			window = []
			for i in range(minSize):
				chiX = 0.0
				chiY = 0.0
				x = revMatrix[i,:]
				y = curPSSM.matPSSM[i+curTime,:]
				xy = x+y
				Nx = sum(x)
				Ny = sum(y)
				N = sum(xy)
				Nex = old_div(Nx*xy, N)
				Ney = old_div(Ny*xy, N)
				for index,v in enumerate(Nex):
					if v == 0: Nex[index] = 0.0000001
				for index,v in enumerate(Ney):
					if v == 0: Ney[index] = 0.0000001
				
				chiX += sum(old_div(power(x-Nex,2), Nex))
				chiY += sum(old_div(power(y-Ney,2), Ney))
				window.append(stats.lchisqprob(chiX+chiY,3))

			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				bestRev = True
				maxChi = curChi
				startPosition = curTime


			
			curTime += 1
		
		#print maxChi
		
		#print str(maxChi)
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return maxChi, startPosition, bestRev



	def	distanceChiPru(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: Float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#self.pseudoCount()
		#curPSSM.pseudoCount()
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matFSM.shape[0] <= curPSSM.matFSM.shape[0]):
			minSize = self.matFSM.shape[0]
			times = curPSSM.matFSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceChi2(self,True)
		#print self.matFSM
		#print curPSSM.matFSM
		#print self.numSamples, curPSSM.numSamples	

		#raw_input()
		#print selfFSM.matrix
		#print curFSM.matrix
		#print selfFSM.numSamples, curFSM.numSamples	
		#raw_input()
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxChi = 0.0
				
		while (curTime < times):
			curChi = 0
			window = []
			for i in range(minSize):	
				chiX = 0.0
				chiY = 0.0
				x = array(self.matPSSM[i,:],float)
				y = array(curPSSM.matPSSM[i+curTime,:],float)
				xy = x+y
				Nx = sum(x)
				Ny = sum(y)
				N = sum(xy)
				Nex = old_div(Nx*xy, N)
				Ney = old_div(Ny*xy, N)
				for index,v in enumerate(Nex):
					if v == 0: Nex[index] = 0.0000001
				for index,v in enumerate(Ney):
					if v == 0: Ney[index] = 0.0000001

				chiX += sum(old_div(power(x-Nex,2), Nex))
				chiY += sum(old_div(power(y-Ney,2), Ney))
				
				window.append(stats.lchisqprob(chiX+chiY,3))
				
				#DEBUG
				#print x,y,xy
				#print Nx,Ny,N
				#print Nex,Ney
				#print chiX,chiY
				#print sum(chiX + chiY)
				#raw_input()
				#print stats.lchisqprob(sum(chiX + chiY),3)
				#raw_input()
				#END DEBUG
			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				maxChi = curChi
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixPSSM()
		curTime = 0

		while (curTime < times):
			curChi = 0.0
			window = []
			for i in range(minSize):
				chiX = 0.0
				chiY = 0.0
				x = revMatrix[i,:]
				y = curPSSM.matPSSM[i+curTime,:]
				xy = x+y
				Nx = sum(x)
				Ny = sum(y)
				N = sum(xy)
				Nex = old_div(Nx*xy, N)
				Ney = old_div(Ny*xy, N)
				for index,v in enumerate(Nex):
					if v == 0: Nex[index] = 0.0000001
				for index,v in enumerate(Ney):
					if v == 0: Ney[index] = 0.0000001
				
				chiX += sum(old_div(power(x-Nex,2), Nex))
				chiY += sum(old_div(power(y-Ney,2), Ney))
				window.append(stats.lchisqprob(chiX+chiY,3))

			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				bestRev = True
				maxChi = curChi
				startPosition = curTime


			
			curTime += 1
		
		#print maxChi
		
		#print str(maxChi)
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return maxChi, startPosition, bestRev


	def	distanceChi2New(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: Float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""
		#self.pseudoCount()
		#curPSSM.pseudoCount()
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matFSM.shape[0] <= curPSSM.matFSM.shape[0]):
			minSize = self.matFSM.shape[0]
			times = curPSSM.matFSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceChi2(self,True)
		#print self.matFSM
		#print curPSSM.matFSM
		#print self.numSamples, curPSSM.numSamples	

		#raw_input()
		#print selfFSM.matrix
		#print curFSM.matrix
		#print selfFSM.numSamples, curFSM.numSamples	
		#raw_input()
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxChi = 0.0
				
		while (curTime < times):
			curChi = 0
			window = []
			for i in range(minSize):	
				
				x = array(self.matFSM[i,:],float)
				y = array(curPSSM.matFSM[i+curTime,:],float)
				Nx = sum(x)
				Ny = sum(y)
				den = Nx*Ny*(x+y)
				
				for i,v in enumerate(den):
					if v == 0: den[i] = 1
				
				chi2 = sum(old_div((power((Ny*x -Nx*y),2)), (den)))
				
				window.append(stats.lchisqprob(chi2,3))
				
				#DEBUG
				#print x,y,xy
				#print Nx,Ny,N
				#print Nex,Ney
				#print chiX,chiY
				#print sum(chiX + chiY)
				#raw_input()
				#print stats.lchisqprob(sum(chiX + chiY),3)
				#raw_input()
				#END DEBUG
			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				maxChi = curChi
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixFSM()
		curTime = 0

		while (curTime < times):
			curChi = 0.0
			window = []
			for i in range(minSize):
				x = array(revMatrix[i,:],float)
				y = array(curPSSM.matFSM[i+curTime,:],float)

				Nx = sum(x)
				Ny = sum(y)
				den = Nx*Ny*(x+y)
				
				for i,v in enumerate(den):
					if v == 0: den[i] = 1
				
				chi2 = sum(old_div((power((Ny*x -Nx*y),2)), (den)))
				
				window.append(stats.lchisqprob(chi2,3))

			#print "Window = ", window
			curChi = utils.geometricMean(window)
			if(curChi > maxChi):
				#print curTime, "Vez. La distancia se mejora. Pasa a ser:", curChi
				bestRev = True
				maxChi = curChi
				startPosition = curTime


			
			curTime += 1
		
		#print maxChi
		maxPv = stats.lchisqprob(maxChi,3)
		print(str(maxChi),maxPv)
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return maxChi, startPosition, bestRev

	def	distanceFisher(self, curPSSM, bigger=False):
		"""Calcula la correlacion de Pearson entre dos PSSM. La correlacion se calcula sumando las correlaciones de cada una de las columnas de las matrices que definen los PSSM's.
Devuelve: (minDist, startPosition, bestRev)
	minDist: Float con la distancia minima entre los dos PSSM
	startPosition: Entero con la posicion a partir de la cual el mas pequeno se parece en grado mayor al mas grande
	bestRev: Boolean que indica si la distancia minima se obtuvo haciendo el revOpp del PSSM original
Parametros:
	curPSSM: PSSM que queremos comparar"""

		self.pseudoCount()
		curPSSM.pseudoCount()
		
		#print "Correlation of " , self.consensus() , " and " + curPSSM.consensus()
		if (self.matFSM.shape[0] <= curPSSM.matFSM.shape[0]):
			minSize = self.matFSM.shape[0]
			times = curPSSM.matFSM.shape[0] - minSize + 1
		else:
			return curPSSM.distanceFisher(self,True)
				
		
		#print "Times -> ", times
		
		bestRev = False
		curTime = 0
		startPosition = 0
		maxPv = -1.0
				
		while (curTime < times):
			curPv = 0.0
			for i in range(minSize):
				x = self.matFSM[i,:]
				y = curPSSM.matFSM[i+curTime,:]
				
				curPv += fisherTest(x,y)
			
			
			if(curPv > maxPv):
				maxPv = curPv
				startPosition = curTime
			
			curTime += 1		

		revMatrix = self.calculateRevMatrixFSM()
		curTime = 0

		while (curTime < times):
			curPv = 0.0
			for i in range(minSize):
				
				x = revMatrix[i,:]
				y = curPSSM.matFSM[i+curTime,:]
				curPv += fisherTest(x,y)
			

			if(curPv > maxPv):
				bestRev = True
				maxPv = curPv
				startPosition = curTime
			
			curTime += 1
		
		#print str(minChi)
		# Comento esto para que me devuelva siempre la posicion a partir de la cual el mas
		# pequeno esta contenido en el mas grande	
		#if(bigger):
		#	startPosition = 0
	
		return maxPv, startPosition, bestRev
	

	
	def matrixForImage(self):
		"""Devuelve la matriz que necesita el seqlogo para generar las imagenes.
Es usada por Motif.createImage"""
		self.expForm()

		ret = 100 * self.matPSSM[:]

		for i in range(self.matPSSM.shape[0]):
			s = sum(ret[i,:])
			ind = ret[i,:].tolist().index(ret[i,:].max())
			ret[i,ind] += 100 - s

		return ret


	def matrixRevForImage(self):
		"""Devuelve la matriz Reverse Opposite que necesita el seqlogo para generar las imagenes.
Es usada por PSSM.createImage"""
		self.expForm()

		ret = array(zeros(4*self.matPSSM.shape[0]))
		ret.shape = (self.matPSSM.shape[0],4)

		ret[:] = 100 * self.calculateRevMatrixPSSM()[:]

		for i in range(self.matPSSM.shape[0]):
			s = sum(ret[i,:])
			ind = ret[i,:].tolist().index(ret[i,:].max())
			ret[i,ind] += 100 - s

		return ret

	def createImage(self,name, rev = False):
		"""Crea la(s) imagen(es) PNG del PSSM utilizando el programa seqlogo.
El programa seqlogo se encuentra en /home/fernan/homer/weblogo/seqlogo.
La(s) imagen(es) se crean en el path facilitado en name.
Si se indica se crea la imagen del PSSM reverse opposite
Parametros:
	name: Path donde se crea(n) la(s) imagen(es)
	rev: Si True se crea la imagen rev opp en el path nameR, es decir se anade una R al nombre facilitado. Por defecto rev = False"""
		matAux = self.matrixForImage()
		matRev = self.matrixRevForImage()
		letters = 'ACGT'
		
		f = open("tmp",'w')
		if rev:
			fR = open("tmpR",'w')
		for k in range(100):
			for i in range(self.matPSSM.shape[0]):
				for j in range(4):
					if (matAux[i,j] > 0):
						f.write(letters[j])
						matAux[i,j] -= 1
						break
					
			f.write("\n")

			if rev:
				for i in range(self.matPSSM.shape[0]):
					for j in range(4):
						if (matRev[i,j] > 0):
							fR.write(letters[j])
							matRev[i,j] -= 1
							break
					
				fR.write("\n")
				
		f.close()
		if rev:
			fR.close()

		cmd = "/home/fernan/articuloFuzzyMotifs/weblogo/seqlogo -Sc -f tmp -F PNG -o " + name
		(status, output) = getstatusoutput(cmd)
		if(status):
			print("Error while trying to execute ", cmd)
		os.remove("tmp")

		if rev:
			cmd = "/home/fernan/homer/articuloFuzzyMotifs/seqlogo -Sc -f tmpR -F PNG -o " + name + "R";
			(status, output) = getstatusoutput(cmd)
			if(status):
				print("Error while trying to execute ", cmd)
			os.remove("tmpR")


				
	def readFromFile(self,fileIn, fileType):# 6-3-2008 Pongo lo de especial para que me coja el nombre bien del fichero de jaspar que me he bajado. Es decir ignore el codigo de la matriz y ponga una concatenacion del nombre y la familia (segundo y tecer campo de la primera fila de cada motivo)
		if type(fileIn) == bytes: fileIn = open(fileIn)
		header = fileIn.readline().split()
		if fileType in ['jaspar','JASPAR','Jaspar']:
			name = header[1] + '-' + header[2]
			self.ID = header[0].replace('>','')
		else:  name = header[0].split('/')[0][1:]  #La cabecera tendra '/' o no dependiendo si es un fichero de Jaspar o de Transfac de los q me dio Chris

		# 6-3-2008 Pongo lo de especial para que me coja el nombre bien del fichero de jaspar que me he bajado. Es decir ignore el codigo de la matriz y ponga una concatenacion del nombre y la familia (segundo y tecer campo de la primera fila de cada motivo)
		mat = []
		col = fileIn.readline().split()
		pos = fileIn.tell()
		while(col and col[0][0] != '>'):
			mat.append([float(s) for s in col])
			pos = fileIn.tell()
			col = fileIn.readline().split()
			
		
		mat = array(mat,float)
		#utils.checkFSM(mat)  # Esto lo pongo por si hago algo de comprobar la bondad del motivo. P.ej ver que todas las columnas sumen lo mismo		
	
		#mat /= sum(mat[0,:]) # Paso de matriz de frecuencias a PSSM
		if (col):
			fileIn.seek(pos)
		
		if self.FSM == None:
			self.FSM = utils.decideFSMorNot(mat)
			
		if self.FSM == True:
			self.numSamples = sum(mat[0,:])
			self.matFSM = mat.astype(int)
			self.n = mat.shape[0]
			self.calculatePSSMFromFSM()
			self.name = name
		elif self.FSM == False:
			self.matPSSM = array(mat,float)
			self.n = mat.shape[0]
			self.calculateFSMFromPSSM()
			self.name = name

	def core(self,n=5):
		"""Devuelve un Motif que representa al core del motivo original. El core se calcula como las n posiciones consecutivas mas conservadas dentro del motivo original. Esta definicion es la que dan en "Similarity of position frequency matrices for transcription factor binding sites" que es la que dicen que toma el TRANSFAC. Habria que estudiar si se pueden definir,  cores de manera mas "inteligente", por ejemplo mediante las lambda medidas difusas (para ver la importancia)
	Recibe:
	n: Longitud del core
	Devuelve:
	Un Motif con el core del motivo original"""

		if n > self.n: n = self.n # Si se busca un core de longitud mayor que el motivo

		inicioCore = 0
		bestConservation = 0
		for i in range(self.n - n):
			curConservation = 0.0
			for j in range(n):
				curConservation += max(self.matPSSM[i+j,:])

			if curConservation > bestConservation: 
				bestConservation = curConservation
				inicioCore = i
		
		core = Motif(matrix = self.matPSSM[inicioCore:inicioCore+n,:],name=self.name+'_core',ID=self.ID+'_core')
		return core
	
	
	def toFile(self,fileOut,FSM=True):
		if type(fileOut) == bytes: fileOut = open(fileOut,'w')
		if self.name == '':
			name = self.consensus()
		else:
			name = self.name
		#print name
		fileOut.write(">"+name + "\t" + str(self.IC()) + "\n")
		#print self.matFSM
		#raw_input()
		if FSM:
			for i in range(self.matFSM.shape[0]):
				for j in range(self.matFSM.shape[1]-1):
					fileOut.write(str(self.matFSM[i,j])+"\t")
				fileOut.write(str(self.matFSM[i,self.matFSM.shape[1]-1])+"\n")
		else:
			for i in range(self.matPSSM.shape[0]):
				for j in range(self.matPSSM.shape[1]-1):
					fileOut.write(str(self.matPSSM[i,j])+"\t")
				fileOut.write(str(self.matPSSM[i,self.matPSSM.shape[1]-1])+"\n")

	def distances(self,curPSSM, f=None):
		print("Euclidea = ", self.distanceEuclidean(curPSSM))
		print("Chi2 = ", self.distanceChi2(curPSSM))
		print("Fuzzy = ", self.distanceFuzzy(curPSSM))
		print("Choquet = ", self.distanceChoquetLineal2(curPSSM))
		if(f!=None):
			if type(f) == bytes: f = open(f,'w')
			self.toFile(f)
			curPSSM.toFile(f)
			f.write("Euclidea = " + str(self.distanceEuclidean(curPSSM))+ "\n")
			f.write("Chi2 = " + str(self.distanceChi2(curPSSM))+ "\n")
			f.write("Fuzzy = " + str(self.distanceFuzzy(curPSSM))+ "\n")
			f.write("Choquet = " + str(self.distanceChoquetLineal2(curPSSM))+ "\n")
			f.write("\n")
		
	
