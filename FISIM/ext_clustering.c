#include <Python.h>              /* Python as seen from C */
#include </usr/local/include/python2.7/numarray/numarray.h> /* NumPy  as seen from C */
#include <math.h>
#include <stdio.h>               /* for debug output */
#include <stdlib.h>
#include "NumPy_macros.h"

//#define DEBUG

double distance(PyArrayObject *mat, long ii, long jj)
{
	int n = mat->dimensions[0];	// Number of objects
  	int d = mat->dimensions[1];	// Number of dimensions

	#ifdef DEBUG
	if(ii >= n || jj >= n)
	{
		fprintf(stderr, "Problem in distance.\n");
		return -1.0;
	}
	#endif
	
	double acum = 0.0;
	int k;
	
	for(k = 0; k < d; k++)
	{
		double dif = IND2(mat, ii, k) - IND2(mat, jj, k);
		acum += dif * dif;
	}
	
	return acum;
}

double distance_to_prototype(PyArrayObject *mat, long ii, double *prototype)
{
	int n = mat->dimensions[0];	// Number of objects
  	int d = mat->dimensions[1];	// Number of dimensions

	#ifdef DEBUG
	if(ii >= n)
	{
		fprintf(stderr, "Problem in distance.\n");
		return -1.0;
	}
	#endif
	
	double acum = 0.0;
	int k;
	
	for(k = 0; k < d; k++)
	{
		double dif = IND2(mat, ii, k) - prototype[k];
		acum += dif * dif;
	}
	
	return acum;
}

static PyObject *cmeans(PyObject *self, PyObject *args)
{
	PyObject *_mat;
	int k, n_it, n_rep;

	// Parsing arguments: mat, k, n_it, n_rep
	if(!PyArg_ParseTuple(args, "Oiii:cmeans", &_mat, &k, &n_it, &n_rep))
		return NULL;
	
	PyArrayObject *mat = NA_InputArray(_mat, tFloat64, NUM_C_ARRAY);
	if(mat == NULL)
	{
		PyErr_Format(PyExc_ValueError, "mat is not an array");
		return NULL;
	}
			
	// Checking array: 2 dimensional and double
	if (mat->nd != 2 || mat->descr->type_num != tFloat64)
	{
    	PyErr_Format(PyExc_ValueError, "mat array is %d-dimensional or not of type Float", mat->nd);
    	return NULL;
  	}
	
	if (k < 0 || n_it < 0 || n_rep < 0)
	{
		PyErr_Format(PyExc_ValueError, "k or n_it or n_rep are negative (%d%d%d)", k, n_it, n_rep);
		return NULL;
	}

	int n = mat->dimensions[0];	// Number of objects
  	int d = mat->dimensions[1];	// Number of dimensions
  	int c, p, i, j;
  	
  	// A clustering for each rep and a next_clustering
  	int **clustering = (int **)malloc(n_rep * sizeof(int*));
  	for(i = 0; i < n_rep; i++)
  		clustering[i] = 0;
  	
  	int *next_clustering = (int *)malloc(n * sizeof(int));
  	
  	// Prototypes
  	double **prototypes = (double **)malloc(k * sizeof(double*));
  	for(i = 0; i < k; i++)
  		prototypes[i] = (double *)malloc(d * sizeof(double));
  	
  	// Costs
  	double *costs = (double *)malloc(n_rep * sizeof(double));
  	
  	int rep, it;
  	int converge = 1;
  	
  	for(rep = 0; rep < n_rep && converge; rep++)
  	{
  		#ifdef DEBUG
  			fprintf(stderr, "rep = %d\n", rep);
  		#endif
	  	clustering[rep] = (int *)malloc(n * sizeof(int));
	  	
	  	// Random clustering
  		for(i = 0; i < n; i++)
  			clustering[rep][i] = rand() % k;
  	
  		double cost = 0.0;
  		int clustersize = 0;
  		
  		for(it = 0; it < n_it; it++)
  		{
	  		#ifdef DEBUG
  				fprintf(stderr, "it = %d\n", it);
  			#endif
  			for(c = 0; c < k; c++)	// For each cluster
  			{
  				// Compute cluster prototype
  				clustersize = 0;
  				
  				for(i = 0; i < d; i++)
  					prototypes[c][i] = 0.0;
  				
  				for(p = 0; p < n; p++)
  					if(clustering[rep][p] == c)
  					{
  						for(i = 0; i < d; i++)
  							prototypes[c][i] += IND2(mat, p, i);
  						
  						clustersize++;
  					}
  				
  				for(i = 0; i < d; i++)
  					prototypes[c][i] /= clustersize;
  			}
  			
  			cost = 0.0;
  			
  			// For each point
  			for(p = 0; p < n; p++)
  			{
  				double min_dist = distance_to_prototype(mat, p, prototypes[0]);
  				int min_index = 0;
  				double aux;
  				
  				// For each cluster: computing min distance
  				for(c = 1; c < k; c++)
  					if((aux = distance_to_prototype(mat, p, prototypes[c])) < min_dist)
  					{
  						min_dist = aux;
  						min_index = c;
  					}
  				
  				// Update point in next_clustering
  				next_clustering[p] = min_index;
  				cost += min_dist;
  			}
  			
  			int *aux;
  			aux = clustering[rep];
  			clustering[rep] = next_clustering;
  			next_clustering = aux;
	  		
	  		#ifdef DEBUG
  				fprintf(stderr, "cost = %f\n", cost);
  			#endif
  			
  			// If next_clustering equals previous one I break the iterations
  			int equals = 1;
  			for(p = 0; p < n; p++)
  				if(clustering[rep][p] != next_clustering[p])
  				{
  					equals = 0;
  					break;
  				}
  			
  			if(equals) break;
  		}
  		
  		if(it == n_it) converge = 0;
  		#ifdef DEBUG
  			fprintf(stderr, "No more iterations\ncost = %f\n", cost);
  		#endif
  		
  		costs[rep] = cost;
  	}
  	
  	// Search for min cost
  	int min = 0;
  	for(i = 0; i < rep; i++)
  		if(costs[i] < costs[min]) min = i;
  	
  	// Construct python list to return
  	PyObject *l = PyTuple_New(n);
  	for(i = 0; i < n; i++)
  		PyTuple_SetItem(l, i, Py_BuildValue("i", clustering[min][i]));
  	
  	// Freeing memory
  	for(i = 0; i < rep; i++)
  		if(clustering[i] != 0) free(clustering[i]);  		
  	free(clustering);
  	free(next_clustering);
  	
  	for(i = 0; i < k; i++)
  		free(prototypes[i]);
  	free(prototypes);
  	
  	free(costs);
  	Py_DECREF(mat);
  	
	return l;
}


static PyObject *cmeans_from_distances(PyObject *self, PyObject *args)
{
	PyObject *_distances;
	int k, n_it, n_rep;

	// Parsing arguments: mat, k, n_it, n_rep
	if(!PyArg_ParseTuple(args, "Oiii:cmeans", &_distances, &k, &n_it, &n_rep))
		return NULL;
	
	PyArrayObject *distances = NA_InputArray(_distances, tFloat64, NUM_C_ARRAY);
	if(distances == NULL)
	{
		PyErr_Format(PyExc_ValueError, "distances is not an array");
		return NULL;
	}
			
	// Checking array: 2 dimensional and double
	if (distances->nd != 2 || distances->descr->type_num != tFloat64)
	{
    	PyErr_Format(PyExc_ValueError, "mat array is %d-dimensional or not of type Float", distances->nd);
    	return NULL;
  	}
	
	if (k < 0 || n_it < 0 || n_rep < 0)
	{
		PyErr_Format(PyExc_ValueError, "k or n_it or n_rep are negative");
		return NULL;
	}

	int n = distances->dimensions[0];	// Number of objects
  	
  	if ( n != distances->dimensions[1] )
	{
		PyErr_Format(PyExc_ValueError, "distances is not a square matrix");
		return NULL;
	}
  	
  	int c, p, i, j;
  	
  	// A clustering for each rep and a next_clustering
  	int **clustering = (int **)malloc(n_rep * sizeof(int*));
  	for(i = 0; i < n_rep; i++)
  		clustering[i] = 0;
  	
  	int *next_clustering = (int *)malloc(n * sizeof(int));
  	
  	// Distances to prototypes [i][j] -> distance from prototype i to point j
  	double **prototypesD = (double **)malloc(k * sizeof(double*));
  	for(i = 0; i < k; i++)
  		prototypesD[i] = (double *)malloc(n * sizeof(double));
  	
  	// Costs
  	double *costs = (double *)malloc(n_rep * sizeof(double));
  	
  	int rep, it;
  	int converge = 1;
  	
  	for(rep = 0; rep < n_rep && converge; rep++)
  	{
  		#ifdef DEBUG
  			fprintf(stderr, "rep = %d\n", rep);
  		#endif
	  	clustering[rep] = (int *)malloc(n * sizeof(int));
	  	
	  	// Random clustering
  		for(i = 0; i < n; i++)
  			clustering[rep][i] = rand() % k;
  	
  		double cost = 0.0;
  		int clustersize = 0;
  		
  		for(it = 0; it < n_it; it++)
  		{
	  		#ifdef DEBUG
  				fprintf(stderr, "it = %d\n", it);
  			#endif
  			for(c = 0; c < k; c++)	// For each cluster
  			{
  				// Compute distances to prototypes
  				clustersize = 0;
  				
  				for(p = 0; p < n; p++)
  				{
  					prototypesD[c][p] = 0.0;
  					int e;
  					
  					for(e = 0; e < n; e++)
  						if(clustering[rep][e] == c)
  						{
							prototypesD[c][p] += IND2(distances, e, p); // Acumulamos en las distancias del prototipo del cluster c al punto p la distancia del punto perteneciente al cluster c, e  al punto p. Esto viene de la idea de que la distancia de un punto al centroide del cluster se puede calcular como la media de las distancias de los puntos pertenecientes a ese cluster hasta el punto en cuestion.
  						
  							clustersize++;
  						}
  				}
  				for(p = 0; p < n; p++)
  				{
  					prototypesD[c][p] /= clustersize;
  					
//  					if(prototypesD[c][p] < 0.0)
//  						fprintf(stderr, "Negative distance to prototype.\n");
  				}
  			}
  			
  			cost = 0.0;
  			
  			// For each point
  			for(p = 0; p < n; p++)
  			{
  				double min_dist = prototypesD[0][p];
  				int min_index = 0;
  				double aux;
  				
  				// For each cluster: computing min distance
  				for(c = 1; c < k; c++)
  					if((aux = prototypesD[c][p]) < min_dist)
  					{
  						min_dist = aux;
  						min_index = c;
  					}
  				
  				// Update point in next_clustering
  				next_clustering[p] = min_index;
  				cost += min_dist;
  			}
  			
  			int *aux;
  			aux = clustering[rep];
  			clustering[rep] = next_clustering;
  			next_clustering = aux;
	  		
	  		#ifdef DEBUG
  				fprintf(stderr, "cost = %f\n", cost);
  			#endif
  			
  			// If next_clustering equals previous one I break the iterations
  			int equals = 1;
  			for(p = 0; p < n; p++)
  				if(clustering[rep][p] != next_clustering[p])
  				{
  					equals = 0;
  					break;
  				}
  			
  			if(equals) break;
  		}
  		
  		if(it == n_it){
			printf("It did not converge.\n");
		       	converge = 0;
		}
  		#ifdef DEBUG
  			fprintf(stderr, "No more iterations\ncost = %f\n", cost);
  		#endif
  		
  		costs[rep] = cost;
  	}
  	
  	// Search for min cost
  	int min = 0;
  	for(i = 0; i < rep; i++)
  		if(costs[i] > 0) min = i;
  	
  	for(; i < rep; i++)
  		if(costs[i] < costs[min] && costs[i] > 0) min = i;  	
  	
  	// Construct python list to return
  	PyObject *l = PyTuple_New(n);
  	
  	for(i = 0; i < n; i++)
  		PyTuple_SetItem(l, i, Py_BuildValue("i", clustering[min][i]));
  	
  	// Freeing memory
  	for(i = 0; i < rep; i++)
  		if(clustering[i] != 0) free(clustering[i]);  		
  	free(clustering);
  	free(next_clustering);
  	
  	for(i = 0; i < k; i++)
  		free(prototypesD[i]);
  	free(prototypesD);
  	
  	free(costs);
  	Py_DECREF(distances);
  	
	return l;
}

static PyObject *cmeans_from_similarities(PyObject *self, PyObject *args)
{
	PyObject *_sims;
	int k, n_it, n_rep;

	// Parsing arguments: mat, k, n_it, n_rep
	if(!PyArg_ParseTuple(args, "Oiii:cmeans", &_sims, &k, &n_it, &n_rep))
		return NULL;
	
	PyArrayObject *sims = NA_InputArray(_sims, tFloat64, NUM_C_ARRAY);
	if(sims == NULL)
	{
		PyErr_Format(PyExc_ValueError, "similarities is not an array");
		return NULL;
	}
			
	// Checking array: 2 dimensional and double
	if (sims->nd != 2 || sims->descr->type_num != tFloat64)
	{
    	PyErr_Format(PyExc_ValueError, "mat array is %d-dimensional or not of type Float", sims->nd);
    	return NULL;
  	}
	
	if (k < 0 || n_it < 0 || n_rep < 0)
	{
		PyErr_Format(PyExc_ValueError, "k or n_it or n_rep are negative");
		return NULL;
	}

	int n = sims->dimensions[0];	// Number of objects
  	
  	if ( n != sims->dimensions[1] )
	{
		PyErr_Format(PyExc_ValueError, "similarities is not a square matrix");
		return NULL;
	}
  	
  	int c, p, i, j;
  	
  	// A clustering for each rep and a next_clustering
  	int **clustering = (int **)malloc(n_rep * sizeof(int*));
  	for(i = 0; i < n_rep; i++)
  		clustering[i] = 0;
  	
  	int *next_clustering = (int *)malloc(n * sizeof(int));
  	
  	// Distances to prototypes [i][j] -> distance from prototype i to point j
  	double **prototypesD = (double **)malloc(k * sizeof(double*));
  	for(i = 0; i < k; i++)
  		prototypesD[i] = (double *)malloc(n * sizeof(double));
  	
  	// Costs
  	double *costs = (double *)malloc(n_rep * sizeof(double));
  	
  	int rep, it;
  	int converge = 1;
  	
  	for(rep = 0; rep < n_rep && converge; rep++)
  	{
  		#ifdef DEBUG
  			fprintf(stderr, "rep = %d\n", rep);
  		#endif
	  	clustering[rep] = (int *)malloc(n * sizeof(int));
	  	
	  	// Random clustering
  		for(i = 0; i < n; i++)
  			clustering[rep][i] = rand() % k;
  	
  		double cost = 0.0;
  		int clustersize = 0;
  		
  		for(it = 0; it < n_it; it++)
  		{
	  		#ifdef DEBUG
  				fprintf(stderr, "it = %d\n", it);
  			#endif

  			for(c = 0; c < k; c++)	// For each cluster
  			{
  				// Término común para este cluster
				// K_c = sum( [kernel[l, m] for l in c for m in c] )
				clustersize = 0;
				
				double K_c = 0.0;
				int l, m;
				for(l = 0; l < n; l++)
					if(clustering[rep][l] == c)
					{
						for(m = 0; m < n; m++)
							if(clustering[rep][m] == c)
								K_c += IND2(sims, l, m);
						
						clustersize++;
					}
				K_c /= (clustersize * clustersize);
				
  				// Compute distances from p to prototypes
  				for(p = 0; p < n; p++)
  				{
  					int e;
  					double K_p = 0.0;
  					
  					for(e = 0; e < n; e++)
  						if(clustering[rep][e] == c)
							K_p += IND2(sims, p, e);

					K_p /= clustersize;
					
  					prototypesD[c][p] = IND2(sims, p, p) + K_c - 2 * K_p;
  					//if(prototypesD[c][p] < 0)
  					//	fprintf(stderr, "Negative distance to prototype.\n");
					if(prototypesD[c][p] < 0)
						//fprintf(stderr, "Negative distance to prototype.\n");
  				}
  			}
  			
  			cost = 0.0;
  			
  			// For each point
  			for(p = 0; p < n; p++)
  			{
  				double min_dist = prototypesD[0][p];
  				int min_index = 0;
  				double aux;
  				
  				// For each cluster: computing min distance
  				for(c = 1; c < k; c++)
  					if((aux = prototypesD[c][p]) < min_dist)
  					{
  						min_dist = aux;
  						min_index = c;
  					}
  				
  				// Update point in next_clustering
  				next_clustering[p] = min_index;
  				cost += min_dist;
  			}
  			
  			int *aux;
  			aux = clustering[rep];
  			clustering[rep] = next_clustering;
  			next_clustering = aux;
	  		
	  		#ifdef DEBUG
  				fprintf(stderr, "cost = %f\n", cost);
  			#endif
  			
  			// If next_clustering equals previous one I break the iterations
  			int equals = 1;
  			for(p = 0; p < n; p++)
  				if(clustering[rep][p] != next_clustering[p])
  				{
  					equals = 0;
  					break;
  				}
  			
  			if(equals) break;
  		}
  		
  		if(it == n_it) converge = 0;
  		#ifdef DEBUG
  			fprintf(stderr, "No more iterations\ncost = %f\n", cost);
  		#endif
  		
  		costs[rep] = cost;
  	}
  	
  	// Search for min cost
  	int min = 0;
  	for(i = 0; i < rep; i++)
  		if(costs[i] < costs[min]) min = i;
  	
  	// Construct python list to return
  	PyObject *l = PyTuple_New(n);
  	
  	for(i = 0; i < n; i++)
  		PyTuple_SetItem(l, i, Py_BuildValue("i", clustering[min][i]));
  	
  	// Freeing memory
  	for(i = 0; i < rep; i++)
  		if(clustering[i] != 0) free(clustering[i]);  		
  	free(clustering);
  	free(next_clustering);
  	
  	for(i = 0; i < k; i++)
  		free(prototypesD[i]);
  	free(prototypesD);
  	
  	free(costs);
  	Py_DECREF(sims);
  	
	return l;
}

/* doc strings: */
static char cmeans_doc[] = \
  "c = cmeans(mat, k, n_it, rep)\n\
   where:\n\
   - mat: a 2-d numarray.array with objects we want to cluster.\n\
   - k\n\
   - n_it: max number of iterations in each execution.\n\
   - rep: number of executions.\n\
   returns:\n\
   - a 1-dimensional numarray.array where c[i] = j where j is the cluster where i belongs.\n";

static char cmeans_from_distances_doc[] = \
  "c = cmeans_from_distances(distances, k, n_it, rep)\n\
   where:\n\
   - disances: a squared 2-d numarray.array with the distances among objects.\n\
   - k\n\
   - n_it: max number of iterations in each execution.\n\
   - rep: number of executions.\n\
   returns:\n\
   - a 1-dimensional numarray.array where c[i] = j where j is the cluster where i belongs.\n";

static char cmeans_from_similarities_doc[] = \
  "c = cmeans_from_similarities(similarities, k, n_it, rep)\n\
   where:\n\
   - disances: a squared 2-d numarray.array with the distances among objects.\n\
   - k\n\
   - n_it: max number of iterations in each execution.\n\
   - rep: number of executions.\n\
   returns:\n\
   - a 1-dimensional numarray.array where c[i] = j where j is the cluster where i belongs.\n";

static char module_doc[] = \
  "module ext_clustering:\n\
   c = cmeans(mat, k, n_it, rep)\n\
   c = cmeans_from_distances(distances, k, n_it, rep)\n\
   c = cmeans_from_similarities(similarities, k, n_it, rep)\n\
  ";

/* 
   The method table must always be present - it lists the 
   functions that should be callable from Python: 
*/
static PyMethodDef ext_clustering_methods[] = {
  {"cmeans",    /* name of func when called from Python */
   cmeans,      /* corresponding C function */
   METH_VARARGS,   /* ordinary (not keyword) arguments */
   cmeans_doc}, /* doc string for gridloop1 function */
  {"cmeans_from_distances",    /* name of func when called from Python */
   cmeans_from_distances,      /* corresponding C function */
   METH_VARARGS,   /* ordinary (not keyword) arguments */
   cmeans_from_distances_doc}, /* doc string for gridloop1 function */
  {"cmeans_from_similarities",    /* name of func when called from Python */
   cmeans_from_similarities,      /* corresponding C function */
   METH_VARARGS,   /* ordinary (not keyword) arguments */
   cmeans_from_similarities_doc}, /* doc string for gridloop1 function */
  {NULL, NULL}     /* required ending of the method table */
};

PyMODINIT_FUNC initext_clustering()
{
  /* Assign the name of the module and the name of the
     method table and (optionally) a module doc string:
  */
  Py_InitModule3("ext_clustering", ext_clustering_methods, module_doc);
  /* without module doc string: 
  Py_InitModule ("ext_gridloop", ext_gridloop_methods); */

  import_libnumarray();   /* required NumPy initialization */
}
