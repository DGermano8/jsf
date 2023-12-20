cdef extern from "../src/hello.h":
    long fibonacci(int n)

def py_fibonacci(int n):
    return fibonacci(n)
