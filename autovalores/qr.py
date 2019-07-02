# Escreva e teste um programa que usa o metodo QR para transformar uma matriz simetrica A, dada como entrada, em uma matriz DIAGONAL.
# O programa deve ler a matriz e uma tolerancia a partir de um arquivo e imprimir: a matriz de entrada, a matriz diagonal e a matriz de acumulada Q.
# Teste se as colunas da matriz de Q final sao automotores da matriz de entrada.

import numpy as np
import math

# matrix = np.array([[4,1,-2,2],[1,2,0,1],[-2,0,3,-2],[2,1,-2,-1]])

def reafFile():

	with open('data2.txt', 'r') as f:
		content = f.read().splitlines()
		
		# Lendo n
		n = int(content[0])

		# Lendo a matriz
		matrix = []
		for x in xrange(1,n+1):
			line = content[x].split(",")
			line = [int(i) for i in line]
			matrix.append(line)
		matrix = np.asarray(matrix)
	
		# Lendo o erro
		error = float(content[n+1])
	return n, matrix, error

def qr(matrix, error, matrix_order):
	original_Q = np.identity(matrix_order)
	modified_matrix = matrix
	norm = float("infinity")

	while norm > error:

		Q, R = obtain_QR(modified_matrix,matrix_order)
		
		new_modified_matrix = R.dot(Q)
		
		original_Q = original_Q.dot(Q)

		modified_matrix = new_modified_matrix

		#criterio de parada
		# retorna os indices dos elementos abaixo da diagonal principal
		lower_matrix_indexes = np.tril_indices(matrix_order,-1)
		lower_triangular_matrix = new_modified_matrix[lower_matrix_indexes] # matriz triangular modificada
		norm = math.sqrt(sum(np.square(lower_triangular_matrix))) # soma quadrada de todos os elementos
	return modified_matrix, original_Q


def obtain_QR(matrix, matrix_order):
	original_Q = np.identity(matrix_order)
	# da coluna 0 ate a coluna n-1
	for j in xrange(0,matrix_order-1):
		# vai ate a coluna n
		for i in xrange(j+1,matrix_order):
		
			ij = matrix[i,j]
			jj = matrix[j,j]

			modified_Q = mount_Qij(ij, jj, i, j,matrix_order)
			modified_matrix = (modified_Q.T).dot(matrix)
			original_Q = original_Q.dot(modified_Q)
			matrix = modified_matrix

	return original_Q, modified_matrix

def mount_Qij(ij, jj, i, j,matrix_order):
	original_Q = np.identity(matrix_order)
	if jj == 0:
		theta = math.pi/2
	else:
		theta = math.atan(ij/jj)
	original_Q[j,j] = math.cos(theta)
	original_Q[i,i] = math.cos(theta)
	original_Q[j,i] = -math.sin(theta)
	original_Q[i,j] = math.sin(theta)

	return original_Q

def main():
	n, matrix ,error = reafFile()

	np.set_printoptions(suppress=True, precision=4)
	
	autovalor,autovetor = qr(matrix, error, n)

	print "Matriz Original: \n"
	print matrix
	print "\nMatriz Diagonal: \n"
	print autovalor
	print "\nMatriz Acumulada: \n"
	print autovetor
	# with open('outputHH.txt',"w") as f:	
	# 	f.write("Matrix Original: \n")
	# 	f.write(" \n".join(map(str, matrix)))
	# 	f.write("\n\n")
	# 	f.write("Matrix Triangular: \n")
	# 	f.write(" \n".join(map(str, t_matrix)))
	# 	f.write("\n\n")
	# 	f.write("Matrix de Householder: \n")
	# 	f.write(" \n".join(map(str, hh_matrix)))
	# 	f.close()

if __name__ == "__main__": main()
