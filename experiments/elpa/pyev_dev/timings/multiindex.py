import numpy as np

class MultiIndex:
    def __init__(self, values:list[list[str]]):
        self.values = values
        self.ndims = len(values)
        self.dims = np.zeros(self.ndims, dtype=int)
        for idim in range(self.ndims):
            self.dims[idim] = len(values[idim])
        self.indx = np.zeros_like(self.dims)

    def value(self):
        result = []
        for idim in range(self.ndims):
            result.append(self.values[idim][self.indx[idim]])
        return result

    def __len__(self):
        return np.prod(self.dims)

    def increment(self):
        for d in range(self.ndims):
            self.indx[d] += 1
            if self.indx[d] == self.dims[d]:
                self.indx[d] = 0
                continue
            else:
                break
        return self.indx

    def __str__(self):
        return f"{str(self.indx)} -> {str(self.value())}"

if __name__ == '__main__':
    criteria = [
        [16,25,36],
        ['s','e'],
        ['A']
    ]
    indx = MultiIndex(criteria)
    print(indx)
    print(len(indx))
    for i in range(len(indx)):
        print(indx.value())
        indx.increment()
    