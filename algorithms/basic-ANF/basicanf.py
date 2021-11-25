import numpy as np

rng = np.random.default_rng()


class BasicANF():
    def __init__(self, edge_list, r, k):
        self.edge_list = edge_list
        self.r = r
        self.k = k

        self.nodes = list(set(np.array(edge_list).flatten()))

        self.n = len(self.nodes)
        self.bitmask_length = int(np.log2(self.n) + self.r)

    def generate_bit_mask(self):
        mask = []
        for i in range(self.bitmask_length):
            # Check whether or not we assign bit to node
            if rng.random() < 0.5**(i + 1):
                mask.append('1')
                break
            else:
                mask.append('0')

        # Append the rest of the bits
        return "".join(mask) + "0" * (self.bitmask_length - len(mask))

    def generate_concatinated_bitmasks(self):
        k_masks = ''
        # Create k masks and concat
        for _ in range(self.k):
            k_masks += self.generate_bit_mask()

        return k_masks

    def compute(self, max_distance):
        M = {}
        IN = {}

        # Generate k bitmask for each node
        M[0] = {node: self.generate_concatinated_bitmasks() for node in self.nodes}

        for h in range(1, max_distance):
            # Set the value of M[h] to previous values
            M[h] = {node: mask for node, mask in M[h - 1].items()}

            # For every edge pair x y set the value of M[h][x] to the current + previous from y
            for x, y in self.edge_list:
                M[h][x] = self.bitwise_or_for_mask(M[h][x], M[h - 1][y])

            # For each node compute the approximate distance based on the k masks
            IN[h] = {node: self.approximate_distance_from_mask(M[h][node]) for node in self.nodes}

        # Return mean neighbourhood size for every distance
        return {h: np.mean(list(IN[h].values())) for h in range(1, max_distance)}, IN, M

    def get_mean_least_zero_bit(self, combined_mask):
        # Find the least significant bit that is set to 0
        splits = self.split_by_k(combined_mask, self.bitmask_length)
        mean = np.mean([self.get_index_of_0(split) for split in splits])
        return mean

    def get_index_of_0(self, mask):
        if '0' in mask:
            return mask.index('0')
        return len(mask)

    def approximate_distance_from_mask(self, combined_mask):
        return 2**self.get_mean_least_zero_bit(combined_mask) / 0.77351

    def split_by_k(self, seq, n):
        '''A generator to divide a sequence into chunks of n units.'''
        while seq:
            yield seq[:n]
            seq = seq[n:]

    def bitwise_or_for_mask(self, mask1, mask2):
        r = ''
        for i in range(len(mask1)):
            if mask1[i] == "1" or mask2[i] == "1":
                r += "1"
            else:
                r += "0"
        return r
