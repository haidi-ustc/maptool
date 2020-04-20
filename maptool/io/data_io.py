import os
import numpy as np
from monty.io import zopen
from monty.json import MontyEncoder,MontyDecoder
from json import dumps,loads
from copy import deepcopy

class DataIO():
    """
    A basic class to operate array data, string or dict data
    and to convert  into format type, to write data into file.

    :param data: (2D numpy.ndarray or string or dict) data to be formatted
    :param row_head: (str or list of str) head line for output data
    :param col_head: (list) first column for output data, then length
                     must equals to length of np.array
    :param indent: (int) indent for json data
    :fmt_all: (str) the format string for every column, must match with
                    data in every column
    :fmt_single: (str) every column use the same format 

    Example:
    >>> rd=np.random.rand(3,4)
    >>> ret=DataIO(rd)
    >>> print(ret.get_str())

    0.691986       0.507859       0.770485       0.956558
    0.525326       0.317083       0.518152       0.486881
    0.612499       0.889347       0.587508       0.160017

    >>> ret=DataIO(rd*100)
    >>> print(ret.get_str())

    69.198561      50.785903      77.048507      95.655789
    52.532644      31.708269      51.815226      48.688147
    61.249944      88.934738      58.750788      16.001678

    >>> ret=DataIO(rd*1000,fmt_single="%15.3e")
    >>> print(ret.get_str())

    6.920e+02      5.079e+02      7.705e+02      9.566e+02
    5.253e+02      3.171e+02      5.182e+02      4.869e+02
    6.125e+02      8.893e+02      5.875e+02      1.600e+02

    >>> ret=DataIO(rd*100,fmt_all ="%10d %10d %10d % 10d\n")
    >>> print(ret.get_str())

        69         50         77         95
        52         31         51         48
        61         88         58         16

    >>> ret=DataIO(rd*100,fmt_all ="%10d %10.5f %10.4f % 10.3e\n")
    >>> print(ret.get_str())

    69   50.78590    77.0485  9.566e+01
    52   31.70827    51.8152  4.869e+01
    61   88.93474    58.7508  1.600e+01

    >>> col_head="{0:^10s}{1:^10s}{2:^10s}{3:^10s}".format("A","BDD","XXX","NAME")
    >>> ret=DataIO(rd*100, col_head=col_head, fmt_all ="%10d %10.5f %10.4f % 10.3e\n")
    >>> print(ret.get_str())

    A        BDD       XXX       NAME   
    69   50.78590    77.0485  9.566e+01
    52   31.70827    51.8152  4.869e+01
    61   88.93474    58.7508  1.600e+01

    >>> col_head="{0:^10s}{1:^10s}{2:^10s}{3:^10s}".format("A","BDD","XXX","NAME")
    >>> row_head=["A","B","C"]
    >>> ret=DataIO(rd*100, row_head=row_head,col_head=col_head, fmt_all ="%10d %10.5f %10.4f % 10.3e\n")
    >>> print(ret.get_str())

              A        BDD       XXX       NAME   
    A         69   50.78590    77.0485  9.566e+01
    B         52   31.70827    51.8152  4.869e+01
    C         61   88.93474    58.7508  1.600e+01

    >>> ret.write('test.txt')

    """

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
