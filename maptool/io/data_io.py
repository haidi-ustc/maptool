import numpy as np
from   monty.io import zopen
from monty.json import MontyEncoder,MontyDecoder
from json import dumps,loads
from copy import deepcopy

class DataIO():
    def __init__(self,data,row_head=None,col_head=None,
                 indent=4,fmt_all=None,fmt_single="%15.6f"):
        if row_head:
            assert isinstance(row_head,list)
            for i in row_head:
                assert isinstance(i,str)
        self._row_head=row_head 

        if isinstance(col_head,list):
           self._col_head=" ".join(col_head)
        elif isinstance(col_head,str):
           self._col_head=col_head
        else:
           self._col_head=''

        self._data=data
        self._str_data=None
        self._data_dim=None
        self.fmt_all=fmt_all
        self.fmt_single=fmt_single

        if isinstance(data,str):
            self._str_data=deepcopy(data)
            self._data_dim=0
        if isinstance(data,np.ndarray):
            if data.ndim==2:
               self._data_dim=self._data.shape
               self.array2str()
            else:
                raise RuntimeError("Unsupported data")
        if isinstance(data,dict):
            self._data_dim=0
            self.dict2str(indent=indent)

    def get_data_shape(self):
        return self._data_dim

    def get_raw_data(self):
        return self._data

    def get_str(self):
        return self._str_data

    def dict2str(self,indent=4):
        self._str_data=dumps(self._data,cls=MontyEncoder,indent=indent)      

    def array2str(self):
        ret=''
        row,col=self.get_data_shape()
        if self._row_head:
            assert row==len(self._row_head)
        ret+=self._col_head+'\n'

        if self.fmt_all:
            fmt=self.fmt_all
        else:
            fmt=self.fmt_single*col+'\n'
        i=0
        for data in self._data.tolist():
            if self._row_head:
                ret += "{} ".format(self._row_head[i])+fmt % tuple(data)
            else:
                ret += fmt % tuple(data)
            i+=1
        self._str_data=ret

    def write(self,file_name):
        with zopen(file_name, "wt") as f:
             f.write(self.get_str())

    @classmethod
    def read_data(cls,file_name,fmt=None):
        with zopen(file_name, "r") as f:
            data=f.read()
        if fmt=="json":
            return cls(loads(data))
        else:
            return cls(data)

    def plot(self):
        pass

if __name__=='__main__':
    import numpy as np
    import os
    rd=np.random.rand(3,4)
    ret=DataIO(rd)
    print(ret.get_str())

    ret=DataIO(rd*100)
    print(ret.get_str())

    ret=DataIO(rd*1000,fmt_single="%15.3e")
    print(ret.get_str())

    ret=DataIO(rd*100,fmt_all ="%10d %10d %10d % 10d\n")
    print(ret.get_str())

    ret=DataIO(rd*100,fmt_all ="%10d %10.5f %10.4f % 10.3e\n")
    print(ret.get_str())
    
    col_head="{0:^10s}{1:^10s}{2:^10s}{3:^10s}".format("A","BDD","XXX","NAME")
    ret=DataIO(rd*100, col_head=col_head, fmt_all ="%10d %10.5f %10.4f % 10.3e\n")
    print(ret.get_str())

    col_head="{0:^10s}{1:^10s}{2:^10s}{3:^10s}".format("A","BDD","XXX","NAME")
    row_head=["A","B","C"]
    ret=DataIO(rd*100, row_head=row_head,col_head=col_head, fmt_all ="%10d %10.5f %10.4f % 10.3e\n")
    print(ret.get_str())

    d={'species': [{'element': 'Ta', 'occu': 1}],
             'abc': [0.162300004, 0.329299995, 0.751900039],
              'lattice': {'@module': 'pymatgen.core.lattice',
                    '@class': 'Lattice',
                      'matrix': [[14.4218997955, 0.0, 0.0],
                             [0.0, 5.7512998581, 0.0],
                                [0.0, 0.0, 5.0816001892]]},
                       '@module': 'pymatgen.core.sites',
                        '@class': 'PeriodicSite',
                         'properties': {}}

    ret=DataIO(d)
    print(ret.get_str())
    ret.write('t.json')

    ret=DataIO("\nhahahahah\njust test")
    print(ret.get_str())
 
    if os.path.exists('t.json'):
       ret=DataIO.read_data(file_name='t.json')
       print(ret.get_str())
