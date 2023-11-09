import ASCsoft.bla as bla
import numpy as np
import pytest
import pickle

@pytest.fixture
def simple_matrix():
  """Gives the matrix

  0 1 2
  1 2 3
  2 3 4
  """
  m = bla.Matrix(3,3)
  for p in [(i,j,i+j) for i in range(3) for j in range(3)]:
    (r,c,v) = p
    m[r,c] = v
  return m

@pytest.fixture
def regular_matrix():
  """Gives the matrix

  0 0 0 1
  0 0 1 0
  0 1 0 0
  1 0 0 0
  """
  m = bla.Matrix(4,4)
  for p in [(i,j,1 if i == 3-j else 0) for i in range(4) for j in range(4)]:
    (r,c,v) = p
    m[r,c] = v
  return m

@pytest.fixture
def vector4():
  """Gives the vector

  1
  2
  3
  4
  """
  x = bla.Vector(4)
  for i in range(4):
    x[i] = i+1
  return x

def assert_matrix(bla_matrix, array):
  (h, w) = bla_matrix.shape
  for p in [(i,j) for i in range(h) for j in range(w)]:
    (i,j) = p
    assert bla_matrix[i,j] == array[i][j]

def test_scalar_multiplication(simple_matrix):
  print(simple_matrix)
  scaled = 3*simple_matrix
  print(scaled)

  assert_matrix(
    scaled,
    [[0,3,6],
     [3,6,9],
     [6,9,12]]
  )

def test_matrix_times_vector(regular_matrix, vector4):
  print(regular_matrix)
  multiplied = regular_matrix * vector4
  print(vector4)
  print(multiplied)

  assert multiplied[0] == 4
  assert multiplied[1] == 3
  assert multiplied[2] == 2
  assert multiplied[3] == 1

def test_matrix_multiplication(simple_matrix):
  multiplied = simple_matrix*simple_matrix

  assert_matrix(
    multiplied,
    [[5,8,11],
     [8,14,20],
     [11,20,29]]
  )

def test_vector_slice(vector4):
  print(vector4)
  vector4[1:2] = 100
  print(vector4)
  assert vector4[0] == 1
  assert vector4[1] == 100
  assert vector4[2] == 100
  assert vector4[3] == 4

def test_pickle(simple_matrix):
  pickled_matrix = pickle.loads(pickle.dumps(simple_matrix))

  assert_matrix(
    pickled_matrix,
    [[0,1,2],
     [1,2,3],
     [2,3,4]]
  )

def test_inverse(regular_matrix):
  inverse = regular_matrix.inverse()
  print(inverse)

  assert_matrix(
    inverse,
    [[0,0,0,1],
     [0,0,1,0],
     [0,1,0,0],
     [1,0,0,0]]
  )