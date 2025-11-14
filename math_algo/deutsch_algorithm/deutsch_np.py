import numpy as np

zero = np.array([1, 0])
one = np.array([0, 1])

def x_gate(x):
    return np.matmul(np.array([[0, 1], [1, 0]]), x)

def cx_gate(x):
    cx = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
    ])
    return np.matmul(cx, x)

def h_gate(x):
    h = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]])
    return np.matmul(h, x)


x = h_gate(zero)
y = h_gate(one)

def f(x):
    return x

def oracle_f(x, y):
    result = np.sqrt(1/2) * ((-1)**f(x)) * np.kron(x, (zero - one))

    return result

def deutsch_algorithm():
    xy = oracle_f(x, y)

    x_final = h_gate(xy[:2])

    return x_final
                            

result = deutsch_algorithm()
print("Result:", result)