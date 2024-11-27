import random

def generate_spanning_tree_test_case(N, M):
    if M < N - 1 or M > N * (N - 1) // 2:
        raise ValueError("Number of edges M must be in the range [N-1, N*(N-1)/2]")

    edges = set()
    for i in range(1, N):
        A = i
        B = random.randint(0, i - 1)
        C = random.randint(1, 100)
        edges.add((A, B, C))

    while len(edges) < M:
        A = random.randint(0, N - 1)
        B = random.randint(0, N - 1)
        if A != B:
            C = random.randint(1, 100)
            if (A, B, C) not in edges and (B, A, C) not in edges:
                edges.add((A, B, C))

    edges = list(edges)
    random.shuffle(edges)

    print(N, M)
    for A, B, C in edges:
        print(1 + A, 1 + B, C)

if __name__ == "__main__":
    N = 10000
    M = 30000
    generate_spanning_tree_test_case(N, M)
