def lu_decomposition(A):
    n = len(A)
    for row in A:
        if len(row) != n:
            raise ValueError("Матрица должна быть квадратной")
 
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]
 
    for i in range(n):
        # i-я строка матрицы U
        for j in range(i, n):
            s = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = A[i][j] - s
 
        if abs(U[i][i]) < 1e-12:
            raise ZeroDivisionError(
                f"u[{i + 1}][{i + 1}] = 0: разложение без перестановки строк невозможно"
            )
 
        # i-й столбец матрицы L
        L[i][i] = 1.0
        for j in range(i + 1, n):
            s = sum(L[j][k] * U[k][i] for k in range(i))
            L[j][i] = (A[j][i] - s) / U[i][i]
 
    return L, U


def multiply(X, Y):
    n, m, p = len(X), len(Y), len(Y[0])
    return [[sum(X[i][k] * Y[k][j] for k in range(m)) for j in range(p)] for i in range(n)]

def lu_solve(L, U, b):
    """
    Решает систему L * U * x = b.
    1) прямая подстановка:   L * y = b
    2) обратная подстановка: U * x = y
    Возвращает решение x в виде списка.
    """
    n = len(L)
    if len(b) != n:
        raise ValueError("Длина b должна совпадать с размером матрицы")
 
    # Прямая подстановка (сверху вниз)
    y = [0.0] * n
    for i in range(n):
        s = sum(L[i][k] * y[k] for k in range(i))
        y[i] = (b[i] - s) / L[i][i]
 
    # Обратная подстановка (снизу вверх)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(U[i][k] * x[k] for k in range(i + 1, n))
        x[i] = (y[i] - s) / U[i][i]
 
    return x




