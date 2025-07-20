import numpy as np

class MultiIndex:
    def __init__(self, dims):
        self.dims = np.array(dims)
        self.ndims = len(dims)
        self.indx = np.zeros_like(dims)
        self.indx[0] = -1
    
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
        return str(self.indx)

if __name__ == '__main__':
    criteria =[
        ('nranks',[16,25,36]),
        ('ba',['s','e']),
    ]
    indx = MultiIndex(criteria)
    print(indx)
    print(len(indx))
    for i in range(len(indx)):
        # print(indx.increment())
        print(indx.next())
    