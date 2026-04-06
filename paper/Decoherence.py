#qutip imports
from qutip import *
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from qutip.tensor import flatten

#standard python imports
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.constants as sc


N=2
def Decoherence1():
    #offd_vector1 = fock(8, 2)
    #offd_vector1_t = fock(8,1).trans()
    #offd_vector2 = fock(8, 1)
    #offd_vector2_t = fock(8,2).trans()
    #lower_matrix_decoherence1 = offd_vector1 * offd_vector1_t
    #lower_matrix_decoherence2 = offd_vector2 * offd_vector2_t
    #offd_decoherence = lower_matrix_decoherence1 + lower_matrix_decoherence2
    #first_decoherence_matrix = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix = Qobj([[0,-1],[1,0]])
    #third_decoherence_matrix = Qobj([[0,1],[-1,0]])
    #first_decoherence_matrix1 = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix1 = Qobj([[0,1],[1,0]])
    #third_decoherence_matrix1 = Qobj([[0,-1],[1,0]])
    #first_decoherence_matrix2 = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix2 = Qobj([[0,1],[-1,0]])
    #third_decoherence_matrix2 = Qobj([[0,1],[1,0]])
    #offd_decoherence1 = tensor(first_decoherence_matrix, second_decoherence_matrix, third_decoherence_matrix)
    #offd_decoherence2 = tensor(first_decoherence_matrix, sigmax(), sigmax())
    #offd_decoherence1 = tensor(first_decoherence_matrix1, second_decoherence_matrix1, third_decoherence_matrix1)
    #offd_decoherence2 = tensor(first_decoherence_matrix2, second_decoherence_matrix2, third_decoherence_matrix2)
    #offd_decoherence = 0.5 * (offd_decoherence1+offd_decoherence2)
    projection = (tensor(fock(N,0),fock(N,1), fock(N,0))) * (tensor(fock(N,0),fock(N,1), fock(N,0))).dag()
    decay = tensor(fock(N,0),fock(N,1),fock(N,0)) * tensor(fock(N,0),fock(N,0),fock(N,1)).trans()

    return decay

def Decoherence2():
    #offd_vector1 = fock(8, 2)
    #offd_vector1_t = fock(8,1).trans()
    #offd_vector2 = fock(8, 1)
    #offd_vector2_t = fock(8,2).trans()
    #lower_matrix_decoherence1 = offd_vector1 * offd_vector1_t
    #lower_matrix_decoherence2 = offd_vector2 * offd_vector2_t
    #offd_decoherence = lower_matrix_decoherence1 + lower_matrix_decoherence2
    #first_decoherence_matrix = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix = Qobj([[0,-1],[1,0]])
    #third_decoherence_matrix = Qobj([[0,1],[-1,0]])
    #first_decoherence_matrix1 = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix1 = Qobj([[0,1],[1,0]])
    #third_decoherence_matrix1 = Qobj([[0,-1],[1,0]])
    #first_decoherence_matrix2 = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix2 = Qobj([[0,1],[-1,0]])
    #third_decoherence_matrix2 = Qobj([[0,1],[1,0]])
    #offd_decoherence1 = tensor(first_decoherence_matrix, second_decoherence_matrix, third_decoherence_matrix)
    #offd_decoherence2 = tensor(first_decoherence_matrix, sigmax(), sigmax())
    #offd_decoherence1 = tensor(first_decoherence_matrix1, second_decoherence_matrix1, third_decoherence_matrix1)
    #offd_decoherence2 = tensor(first_decoherence_matrix2, second_decoherence_matrix2, third_decoherence_matrix2)
    #offd_decoherence = 0.5 * (offd_decoherence1+offd_decoherence2)
    projection = (tensor(fock(N,0),fock(N,0), fock(N,1))) * (tensor(fock(N,0),fock(N,0), fock(N,1))).dag()


    return projection

def Decoherence3():
    #offd_vector1 = fock(8, 2)
    #offd_vector1_t = fock(8,1).trans()
    #offd_vector2 = fock(8, 1)
    #offd_vector2_t = fock(8,2).trans()
    #lower_matrix_decoherence1 = offd_vector1 * offd_vector1_t
    #lower_matrix_decoherence2 = offd_vector2 * offd_vector2_t
    #offd_decoherence = lower_matrix_decoherence1 + lower_matrix_decoherence2
    #first_decoherence_matrix = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix = Qobj([[0,-1],[1,0]])
    #third_decoherence_matrix = Qobj([[0,1],[-1,0]])
    #first_decoherence_matrix1 = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix1 = Qobj([[0,1],[1,0]])
    #third_decoherence_matrix1 = Qobj([[0,-1],[1,0]])
    #first_decoherence_matrix2 = Qobj([[1,0],[0,0]])
    #second_decoherence_matrix2 = Qobj([[0,1],[-1,0]])
    #third_decoherence_matrix2 = Qobj([[0,1],[1,0]])
    #offd_decoherence1 = tensor(first_decoherence_matrix, second_decoherence_matrix, third_decoherence_matrix)
    #offd_decoherence2 = tensor(first_decoherence_matrix, sigmax(), sigmax())
    #offd_decoherence1 = tensor(first_decoherence_matrix1, second_decoherence_matrix1, third_decoherence_matrix1)
    #offd_decoherence2 = tensor(first_decoherence_matrix2, second_decoherence_matrix2, third_decoherence_matrix2)
    #offd_decoherence = 0.5 * (offd_decoherence1+offd_decoherence2)
    projection = (tensor(fock(N,1),fock(N,0), fock(N,0))) * (tensor(fock(N,1),fock(N,0), fock(N,0))).dag()

    return projection
