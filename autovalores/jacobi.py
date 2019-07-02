# Escreva e teste um programa que usa o metodo de Jacobi para transformar uma matriz simetrica A, dada como entrada, em uma matriz DIAGONAL.
# O programa deve ler a matriz e uma tolerancia a partir de um arquivo e imprimir: 
# a matriz de entrada, a matriz diagonal e a matriz de Jacobi. 
# Teste se as colunas da matriz de Jacobi sao automotores da matriz de entrada.

import numpy as np
import math

# matrix = np.array([[4,2,2,1],[2,-3,1,1],[2,1,3,1],[1,1,1,2]])

# matrix = np.array([[4,1,-2,2],[1,2,0,1],[-2,0,3,-2],[2,1,-2,-1]])
#matrix = np.array([[1,2],[3,4]])


def jacobi(matrix, error, matrix_order):
	modified_matrix = matrix
	jacobi_matrix = np.identity(matrix_order)
	e = float("infinity")

	while e > error:
		# da coluna 0 ate a coluna n-1
		for j in range(0,matrix_order-1):
			# vai ate a coluna n
			for i in range(j+1,matrix_order):
				# utiliza os elementos ajj, aij, aii da matrix modificada a
				jj = modified_matrix[j,j]
				ij = modified_matrix[i,j]
				ii = modified_matrix[i,i]
				J_matrix = construct_J(jj, ij, ii, i, j, matrix_order)

				# print("MATRIZ H")
				# print(J_matrix)
				# print("\n")

				modified_matrix = (J_matrix.T).dot(modified_matrix).dot(J_matrix)

				# print("MATRIZ A")
				# print(modified_matrix)
				# print("\n")

				jacobi_matrix = jacobi_matrix.dot(J_matrix)

				# print("MATRIZ ACUMULADA")
				# print(J_matrix)
				# print("\n===============================================")
		
		#criterio de parada
		# retorna os indices dos elementos abaixo da diagonal principal
		lower_matrix_indexes = np.tril_indices(matrix_order,-1)
		lower_triangular_matrix = modified_matrix[lower_matrix_indexes] # matriz triangular modificada
		e = np.linalg.norm(lower_triangular_matrix) # soma quadrada de todos os elementos

	return modified_matrix, jacobi_matrix


def construct_J(jj, ij, ii, i, j, matrix_order):
	J_matrix = np.identity(matrix_order)
	theta = 0
	if ii == jj:
		theta = math.pi/4
	else:
		theta = (math.atan(2*ij/(jj-ii)))/2

	J_matrix[j,j] =	math.cos(theta)
	J_matrix[j,i] = -math.sin(theta)
	J_matrix[i,j] = math.sin(theta)
	J_matrix[i,i] = math.cos(theta)

	return J_matrix


def main():
	A = 3
	B = 7
	C = 7
	D = 7
	E = 3
	F = 1
	matrix = np.matrix([[30+A+F , A , B , C , D],[A, 10 + B + E, E , F , A + B],[B , E , 50 + C + D, B + C , C + D]])

	erro = 0.01

	np.set_printoptions(suppress=True, precision=4)
	print("Matriz Original: \n")
	print( matrix)
	AAt = matrix*matrix.transpose()
	print("AAt:\n")
	print(AAt)	
	diagonal, U = jacobi(AAt, erro, 3)
	print ("\nMatriz Diagonal: \n")
	print (diagonal)
	print ("\nMatriz Jacobi: \n")
	print (U)

	S = np.sqrt(np.absolute(diagonal))
	print(S)

	AtA = matrix.transpose()*matrix
	print("AAt:\n")
	print(AAt)
	diagonal, V = jacobi(AtA, erro, 5)
	print ("\nMatriz Diagonal: \n")
	print (diagonal)
	print ("\nMatriz Jacobi: \n")
	print (V)

	Sigma = np.append(S,[[0,0],[0,0],[0,0]],axis=1)
	print(Sigma)
	print("U*S*Vt")
	print(U*Sigma*V.transpose())
	print(matrix)

	# USVt = U*S*

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

main()
