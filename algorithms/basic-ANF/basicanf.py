import numpy as np
import tqdm


class BasicANF():
    def __init__(self, edge_list, r, k):
        self.edge_list = edge_list
        self.r = r
        self.k = k
        self.rng = np.random.default_rng(seed=103)

        self.nodes = list(set(np.array(edge_list).flatten()))

        self.n = len(self.nodes)
        self.bitmask_length = int(np.log2(self.n) + self.r)

    def generate_bit_mask(self):
        mask = np.full(self.bitmask_length, False)

        odds = 0.5**(np.arange(0, self.bitmask_length) + 1)
        index = self.rng.choice(a=np.arange(0, self.bitmask_length + 1), p=[*odds, 1 - np.sum(odds)])

        if index != len(mask):
            mask[index] = True

        return mask

    def generate_concatinated_bitmasks(self):
        k_masks = np.full(self.k * self.bitmask_length, False)
        # Create k masks and concat
        for k in range(self.k):
            k_masks[k * self.bitmask_length: (k + 1) * self.bitmask_length] = self.generate_bit_mask()

        return k_masks

    def compute(self, max_distance):
        M = {}
        IN = {}

        # Generate k bitmask for each node
        M[0] = {node: self.generate_concatinated_bitmasks() for node in self.nodes}
        IN[0] = {node: 0 for node in self.nodes}

        for h in range(1, max_distance):
            # Set the value of M[h] to previous values
            M[h] = {node: mask for node, mask in M[h - 1].items()}

            # For every edge pair x y set the value of M[h][x] to the current + previous from y
            for x, y in self.edge_list:
                np.logical_or(M[h][x], M[h - 1][y], out=M[h][x])

            # For each node compute the approximate distance based on the k masks
            IN[h] = {node: self.approximate_distance_from_mask(M[h][node]) for node in self.nodes}

        # Return mean neighbourhood size for every distance
        return {h: np.mean(list(IN[h].values())) for h in range(0, max_distance)}, IN, M

    def approximate_distance_from_mask(self, combined_mask):
        return 2**self.get_mean_least_zero_bit(combined_mask) / 0.77351

    def get_mean_least_zero_bit(self, combined_mask):
        # Find the least significant bit that is set to 0

        s = 0
        for k in range(self.k):
            s += self.get_index_of_0(combined_mask[k * self.bitmask_length: (k + 1) * self.bitmask_length])
        return s / self.k

    def get_index_of_0(self, mask):
        # Find the index of the least significant bit that is set to 0
        args = np.argwhere(mask == False)
        if len(args) == 0:
            return self.bitmask_length

        return args[0][0]

    def _test_generate_bitmask(self):
        masks = {b: 0 for b in range(self.bitmask_length + 1)}
        N_samples = 1_000_00
        for _ in tqdm.tqdm(range(N_samples)):
            mask = self.generate_bit_mask()
            b = np.argwhere(mask == True).reshape(-1)

            if len(b) == 0:
                masks[self.bitmask_length] += 1
            else:
                masks[b[0]] += 1

        for k, v in masks.items():
            print(k, v / N_samples)
