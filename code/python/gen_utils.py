import numpy as np

def generate(K: int, P: int, seed: int, write=False, path="D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\data\\utils\\"):
    np.random.seed(seed)
    base_price = np.random.uniform(50, 200, P)
    utils = np.transpose([np.random.normal(mu, 30, K) for mu in base_price])
    if write:
        with open(path + 'up_k%d_p%d_s%d.txt' % (K, P, seed), 'w') as f:
            f.write("%d %d \n" % (K, P))
            for line in utils:
                f.write(" ".join([str(int(x)) for x in line]) + "\n")