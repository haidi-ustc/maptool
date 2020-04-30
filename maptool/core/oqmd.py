import os
import json
import requests
from pymatgen import Structure
from monty.serialization import dumpfn
from maptool.util.utils import wait_sep,wait,procs,sepline
#from pymatgen.entries.computed_entries import ComputedEntry
#from pymatgen.analysis.phase_diagram import *

web="http://oqmd.org"
Rester_ENDPOINT = os.environ.get('oqmd_url', web)

class QMPYRester(object):
    def __init__(self, endpoint=Rester_ENDPOINT):
        self.preamble = endpoint
        self.session = requests.Session()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.session.close()

    def _make_requests(self, sub_url, payload=None, method='GET'):
        response = None
        url = self.preamble + sub_url

        if method == 'GET':
            response = self.session.get(url, params=payload, verify=True)
            
            if response.status_code in [200, 400]:
                data = json.loads(response.text)
                return data

    def get_oqmd_phases(self, verbose=True, **kwargs):
        """
        Input:
            verbose: boolean
            **kwargs: dict
        Output:
            dict
        """

        # URL paramters
        url_args = []
        kwargs_list = ['composition', 'icsd', 'filter',
                       'sort_by', 'desc', 'sort_offset',
                       'limit', 'offset']

        # Attributes for filters
        filter_args = []
        filter_list = ['element_set', 'element', 'spacegroup'
                       'prototype', 'generic', 'volume',
                       'natoms', 'ntypes', 'stability',
                       'delta_e', 'band_gap']

        for k in kwargs.keys():
            if k in kwargs_list:
                url_args.append('%s=%s' %(k, kwargs[k]))
            elif k in filter_list:
                if '>' in kwargs[k] or '<' in kwargs[k]:
                    filter_args.append('%s%s' %(k, kwargs[k]))
                else:
                    filter_args.append('%s=%s' %(k, kwargs[k]))
            elif k == 'fields':
                if '!' in kwargs[k]:
                    url_args.append('fields!=%s' %kwargs[k].replace('!',''))
                else:
                    url_args.append('fields=%s' %kwargs[k])

        if filter_args != []:
            filters_tag = ' AND '.join(filter_args)
            url_args.append('filter='+filters_tag)
            
        if verbose:
            print("Your filters are:")
            if url_args == []:
                print("   No filters?")
            else:
                for arg in url_args:
                    print("   ", arg)

            ans = input('Proceed? [Y/n]:')

            if ans not in ['Y', 'y', 'Yes', 'yes']:
                return

        _url = '&'.join(url_args)
        self.suburl = _url

        return self._make_requests('/oqmdapi/formationenergy?%s'%_url)

    def get_oqmd_phase_space(self, space, **kwargs):
        """
        Input:
            space: str
                e.g. 'Al-O', 'Pd-Se-Te-Pb'
            verbose: boolean
            **kwargs: dict
        Output:
            dict
        """
        _url = 'composition={}'.format(space)
        _url += '&fields=name,delta_e,stability'
        _url += '&limit=200'
        output_data = self._make_requests('/oqmdapi/formationenergy?%s'%_url)

        while output_data['links']['next']:
            new_data = self._make_requests(output_data['links']['next'].replace(self.preamble, ''))
            output_data['data'].extend(new_data['data'])
            output_data['links']['next'] = new_data['links']['next']
            output_data['meta']['data_returned'] += new_data['meta']['data_returned']

        output_data['meta']['more_data_available'] = False

        return output_data

    def get_oqmd_phase_by_id(self, fe_id, fields=None):
        if fields:
            if '!' in fields:
                ex_fields = fields.replace('!', '')
                return self._make_requests('/oqmdapi/formationenergy/%d?fields!=%s'%(fe_id, ex_fields))
            else:
                return self._make_requests('/oqmdapi/formationenergy/%d?fields=%s'%(fe_id, fields))

        return self._make_requests('/oqmdapi/formationenergy/%d'%fe_id)

    def get_optimade_structures(self, verbose=True, **kwargs):
        """
        Input:
            verbose: boolean
            **kwargs: dict
        Output:
            dict
        """

        # URL paramters
        url_args = []
        kwargs_list = ['limit', 'offset', 'filter']

        # Attributes for filters
        filter_args = []
        filter_list = ['elements', 'nelements',
                       'chemical_formula', 'formula_prototype', 
                       '_oqmd_volume', '_oqmd_spacegroup',
                       '_oqmd_natoms', '_oqmd_prototype', 
                       '_oqmd_stability', '_oqmd_delta_e', 
                       '_oqmd_band_gap']

        for k in kwargs.keys():
            if k in kwargs_list:
                url_args.append('%s=%s' %(k, kwargs[k]))
            elif k in filter_list:
                if '>' in kwargs[k] or '<' in kwargs[k]:
                    filter_args.append('%s%s' %(k, kwargs[k]))
                else:
                    filter_args.append('%s=%s' %(k, kwargs[k]))
            elif k == 'fields':
                if '!' in kwargs[k]:
                    url_args.append('fields!=%s' %kwargs[k].replace('!',''))
                else:
                    url_args.append('fields=%s' %kwargs[k])


        if filter_args != []:
            filters_tag = ' AND '.join(filter_args)
            url_args.append('filter='+filters_tag)
            
        if verbose:
            print("Your filters are:")
            if url_args == []:
                print("   No filters?")
            else:
                for arg in url_args:
                    print("   ", arg)

            ans = input('Proceed? [Y/n]:')

            if ans not in ['Y', 'y', 'Yes', 'yes']:
                return

        _url = '&'.join(url_args)
        self.suburl = _url

        return self._make_requests('/optimade/structures?%s'%_url)

    def get_optimade_structure_by_id(self, id, fields=None):
        """
        #optimade format
        elements: the set of elements that the compound must have, e.g. Si,O
        nelements: number of elements types in the compound, e.g. 2, <3
        chemical_formula: compostion of the materials, e.g. Al2O3
        formula_prototype: chemical formula abstract, e.g. AB, AB2
        _oqmd_natoms: number of atoms in the supercell, e.g. 2, >5
        _oqmd_volume: volume of the supercell, e.g. >10
        _oqmd_spacegroup: the space group of the structure, e.g. Fm-3m
        _oqmd_prototype: structure prototype of that compound, e.g. Cu, CsCl
        _oqmd_stability: hull distance of the compound, e.g. 0, <-0.1,
        _oqmd_delta_e: formation energy of that compound, e.g. <-0.5,
        _oqmd_band_gap: band gap of the materials, e.g. 0, >2
        fields: return subset of fields, e.g. 'elements,chemical_formula', '!_oqmd_sites'
        filter: customized filters, e.g. 'elements=O AND ( _oqmd_stability<-0.1 OR _oqmd_delta_e<-0.5 )'
        limit: number of data return at once
        offset: the offset of data return
        """
        if fields:
            if '!' in fields:
                ex_fields = fields.replace('!', '')
                return self._make_requests('/optimade/structures/%d?fields!=%s'%(id, ex_fields))
            else:
                return self._make_requests('/optimade/structures/%d?fields=%s'%(id, fields))

        return self._make_requests('/optimade/structures/%d'%id)

def get_oqmd_structure():
    print('{} >>> {}'.format('1','get a structure by oqmd-ID'))
    print('{} >>> {}'.format('2','get structures by fomular'))
    print('{} >>> {}'.format('3','get structures by elements'))
    print('{} >>> {}'.format('4','get structures by filters'))
    wait_sep()
    in_str=wait()
    choice=in_str
    qr=QMPYRester()
    if choice=="1":
       print("input the oqmd-ID")
       wait_sep()
       in_str=wait()
       oqmd_id=int(in_str)
       step_count=1
       proc_str="Reading Data From "+ web +" ..."
       procs(proc_str,step_count,sp='-->>')
       data=qr.get_optimade_structure_by_id(oqmd_id)
       if data is None:
           print("Unknown OQMD-ID, check the OQMD-ID")
           return None
       else:
           struct=data_to_structure(data)
       filename=str(oqmd_id)+'.vasp'
       step_count+=1
       proc_str="Writing Data to "+filename +" File..."
       procs(proc_str,step_count,sp='-->>')
       struct.to(filename=filename,fmt='POSCAR')
       return True
    elif choice=="2":
       print("Not supported now!")
       os._exit(0)
    elif choice=="3":
       print("Not supported now!")
       os._exit(0)
    elif choice=="4":
       print("Not supported now!")
       os._exit(0)
       #print("elements=O AND ( _oqmd_stability<-0.1 OR _oqmd_delta_e<-0.5)")
    else:
       print('Unknow choice')
       return None

def data_to_structure(data):
    '''
    exctract structure data from oqmd request data
    '''
    lattice=data['_oqmd_unit_cell']
    species=[];coords=[]
    for d in data['_oqmd_sites']:
        _d=d.split("@")
        species.append(_d[0].strip())
        coords.append(list(map(float,_d[1].split())))
    return Structure(lattice=lattice,species=species,coords=coords)

if __name__=="__main__":
   from monty.serialization import loadfn,dumpfn 
   with QMPYRester() as q:
        data = q.get_oqmd_phase_space('Pd-O',verbose=False)
        print(data)
        dumpfn(data,'test.json',indent=4) 
