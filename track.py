
import math
import copy


class Track():

    def __init__(self, _dim=None):
        self.dim = _dim
        self.values = None
        self.ans = None
        self.children = []

    def __create_child(self, d_val):
        for i in range(d_val):
            self.children.append(Track())

    def create_dimension(self, dim):
        if len(dim) == 0:
            self.ans = math.inf
        else:
            d_val = dim.pop(0)
            self.__create_child(d_val)
            i = 0
            for t in self.children:
                t.values = i
                t.create_dimension(copy.copy(dim))
                i += 1

    def get_value(self, div):
        if len(self.children) == 0:
            return self.ans
        else:
            val = div.pop(0)
            for t in self.children:
                if t.values == val:
                    return t.get_value(copy.copy(div))

    def set_value(self, div, val):
        if len(self.children) == 0:
            self.ans = val
        else:
            tempVal = div.pop(0)
            for t in self.children:
                if t.values == tempVal:
                    t.set_value(copy.copy(div), val)

    # def get_missing_count(self, dim_index, f_dim):
        # d = self.dim[dim_index]
        # print(d)
        # cnt = 0
        # for i in range(d):
        #     print(self.get_value([i, f_dim]))
        #     if self.get_value([i, f_dim]) == math.inf:
        #         cnt = cnt + 1
        # return cnt
