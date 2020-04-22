#!/usr/bin/env python
# coding: utf-8
import os
import numpy as np
from maptool import mlog
from pymatgen import Structure
from pymatgen.io.vasp import Poscar
from monty.serialization import dumpfn,loadfn
from maptool.util.utils import wait,wait_sep,multi_structs,warn_tip
from maptool.io.read_structure import read_structures,ase2pmg
from maptool.io.data_io import DataIO
from maptool.core.substrate import match_substrate
from maptool.core.mpdb import get_mp_film_substrate
from maptool.core.selection import atom_selection

def move_to_zcenter(struct):
    '''
    move slab into center along z direction

    Args:
        struct: Structure obj.

    Returns:
        struct: Structure obj.
    '''
    frac=struct.frac_coords
    zcenter=np.mean(frac[:,2])
    zshift=0.5-zcenter
    frac[:,2]=frac[:,2]+zshift
    return Structure(struct.lattice,struct.species,frac)

def ripple(structs,fnames,supercell,strain_range,atom_index,auto_fix,margin_dist=1.0):
    '''
    build ripple structures and strained ones 

    Args:
        structs: list of Structure obj.
        fnames: list of str for input file name.
        supercell: list, used to build supercell along x or y direction.
        strain_range: list of strain 
        atom_index: int , specify which atom will be constrained
        auto_fix: boolean, whether decide which atoms to be constrained automaticly
        margin_dist: float,
    Returns:
        None
    Output:
        0.04000_wo_xxxx.vasp  
        0.04000_w_xxxx.vasp
    '''
    for struct,fname in zip(structs,fnames):
        if not struct.lattice.is_orthogonal:
            warn_tip(0,"{} is not orthogonal, skip it !!!".format(fname))
            continue
        struct_sc=struct.copy()
        struct_sc.make_supercell(supercell)
        natom=struct_sc.num_sites
        min_z=np.min(struct_sc.cart_coords[:,2])
        tmp_coords=np.ones((struct_sc.num_sites,1))*min_z
        cart_coords=struct_sc.cart_coords;
        cart_coords[:,2]=cart_coords[:,2]-tmp_coords.T+0.01
        frac_coords_new=np.dot(cart_coords,np.linalg.inv(struct_sc.lattice.matrix.copy()))
        for i_strain in strain_range:
            new_lat_matrix=struct_sc.lattice.matrix.copy()
            new_lat_matrix[direction,direction]=struct_sc.lattice.matrix[direction,direction]*(1-i_strain)
            # structure only applied with in-plan strain
            filename="%10.5f"%(i_strain)+'_wo_'+fname+'.vasp'
            struct_wo_ripple=Structure(new_lat_matrix,struct_sc.species,frac_coords_new)
            struct_wo_ripple.to(filename=filename.strip(),fmt='poscar')
            frac_coords_new_cp=struct_wo_ripple.frac_coords.copy()
            cart_coords_new_cp=struct_wo_ripple.cart_coords.copy()
            nz=0
            selective_dynamics=[[True for col in range(3)] for row in range(natom)]  
            z_shift=np.zeros((natom,3))
            for i_atom in range(natom):
                #To do : z_shit should support more function 
                z_shift[i_atom,2]=40*(2*i_strain-10*i_strain**2)*np.sin(cart_coords_new_cp[i_atom,direction]*np.pi/new_lat_matrix[direction,direction]) 
                if cart_coords_new_cp[i_atom,direction]< nz or cart_coords_new_cp[i_atom,direction] > new_lat_matrix[direction,direction] -nz:
                   z_shift[i_atom,2]=0.0
                   
                if auto_fix:
                   if struct_wo_ripple[i_atom].coords[direction]<margin_dist or \
                      struct_wo_ripple[i_atom].coords[direction] > new_lat_matrix[direction,direction]-margin_dist:
                      selective_dynamics[i_atom]=[False,False,False]
                else:     
                   if i_atom in atom_index:
                      selective_dynamics[i_atom]=[False,False,False]
            
            struct_w_ripple=Structure(new_lat_matrix,struct_sc.species,cart_coords_new_cp+z_shift,coords_are_cartesian=True,\
                                     site_properties={'selective_dynamics':selective_dynamics})
            # structure  applied with in-plan strain  and ripple    
            filename="%10.5f"%(i_strain)+'_w_'+fname+'.vasp' 
            struct_w_ripple.to(filename=filename.strip(),fmt='poscar')        

def multi_layers(structs,fnames,layer_number,layer_distance):
    '''
    build ripple structures and strained ones 

    Args:
        structs: list of Structure obj.
        fnames: list of str for input file name.
        layer_number: int, number of layers .
        layer_distance: float, distance between two layers.
    Returns:
        None
    Output:
        layer_num_xxxx.vasp  
    '''
    for struct,fname in zip(structs,fnames):
        species=[]
        new_struct=move_to_zcenter(struct)
        struct_thickness=np.max(new_struct.cart_coords[:,2])-np.min(new_struct.cart_coords[:,2])
        natom=new_struct.num_sites
        new_cart_coords=np.zeros((natom*layer_number,3))
        for i in range(layer_number):
            new_cart_coords[i*natom:(i+1)*natom,0:2]=new_struct.cart_coords[:,0:2]
            new_cart_coords[i*natom:(i+1)*natom,2]=new_struct.cart_coords[:,2]+i*(layer_distance+struct_thickness)
            species.extend(new_struct.species)
        new_lat=new_struct.lattice.matrix.copy()
        new_lat[2,2]=new_lat[2,2]+layer_distance*layer_number
        tmp_struct=Structure(new_lat,species,new_cart_coords,coords_are_cartesian=True)
        tmp1_struct=move_to_zcenter(tmp_struct)
        tmp2_struct=tmp1_struct.get_sorted_structure()
        tmp2_struct.to(filename='layer_'+str(layer_number)+'_'+fname+'.vasp',fmt='poscar')

def split(structs,fnames,atom_index,in_str,SplitDistance,NumberSplitSite,ProperDist=3.5,DenseFrac=0.75):
    # len(structs)==1
    '''
    split one or several layers from bulk or structures 

    Args:
        structs: list of Structure obj.
        fnames: list of str for input file name.
        atom_index: list, which atom will be split .
        SplitDistance: float, the largest distance for splitting.
        NumberSplitSite: int, how many structures to be constructured during split 

    Returns:
        None
    Output:
        split_xxxx.dat  
        Num_xxx.vasp
    '''
    col_head="#%(key1)+12s  %(key2)+12s"%{'key1':'index','key2':'distance/Ang'}
    for struct,fname in zip(structs,fnames):
        DensityN=int(NumberSplitSite*DenseFrac)
        SparseN=NumberSplitSite-DensityN+1
        dist=ProperDist/(DensityN-1)
        SplitDistanceArray=np.zeros(NumberSplitSite+1)        
        for Nsite in range(DensityN):
            SplitDistanceArray[Nsite]=(Nsite)*dist
        
        dist=(SplitDistance-ProperDist)/SparseN
        for  Nsite in range(SparseN):
             SplitDistanceArray[Nsite+DensityN]=ProperDist+(Nsite+1)*dist

        coords=struct.cart_coords
        for Nsite in range(NumberSplitSite+1):
            coords=struct.cart_coords
            for atom in atom_index:
                coords[atom,2]=coords[atom,2]+SplitDistanceArray[Nsite]
            tmp_struct=Structure(struct.lattice,struct.species,coords,coords_are_cartesian=True) 
            filename=str(Nsite)+'_'+fname+'.vasp'
            tmp_struct.to(filename=filename,fmt='poscar')
        data=np.zeros((NumberSplitSite+1,2))
        for i,j in enumerate(SplitDistanceArray):
            data[i][0]=i
            data[i][1]=j
        ret=DataIO(data,col_head=col_head,fmt_all ="%12d %12.6f"+'\n')
        filename='split_'+fname+'.dat'
        ret.write(filename)

def resize_vacuum(structs,fnames,nvac_layer_thickness):
    for struct,fname in zip(structs,fnames):
        if not struct.lattice.is_orthogonal:
            warn_tip(0,"{} is not orthogonal, skip it !!!".format(fname))
            continue
        new_struct=move_to_zcenter(struct)
        struct_thickness=np.max(new_struct.cart_coords[:,2])-np.min(new_struct.cart_coords[:,2])
        vac_layer_thickness=new_struct.lattice.c-struct_thickness

        new_lat=new_struct.lattice.matrix.copy()
        new_lat[2,2]=nvac_layer_thickness
        tmp_struct=Structure(new_lat,new_struct.species,new_struct.cart_coords,coords_are_cartesian=True)
        center_struct=move_to_zcenter(tmp_struct)
        center_struct.to(filename='new_vacuum_'+fname+'.vasp',fmt='poscar')

def constrain(structs,fnames,atom_index,in_str=" "):
    #len(structs)==1
    mlog.debug("len(structs): {} ".format(len(structs)))
    mlog.debug("fnames: {} ".format(fnames[0]))
    for struct,fname in zip(structs,fnames):
        mlog.debug(fname)
        mlog.debug(struct)
        natom=struct.num_sites
        selective_dynamics=[[True for col in range(3)] for row in range(natom)]  
        for i in range(natom):
            if i in atom_index:
                selective_dynamics[i]=[False,False,False]
        tmp_struct=Structure(struct.lattice,struct.species,struct.frac_coords,site_properties={'selective_dynamics':selective_dynamics})
        poscar=Poscar(tmp_struct)
        poscar.comment=poscar.comment+' |--> '+in_str
        filename='Fixed_'+fname+'.vasp'
        poscar.write_file(filename)

def circle_strain(structs,fnames,C11,C12,C22,C66,sigma,nps=37):
    for struct,fname in zip(structs,fnames):
        if not struct.lattice.is_orthogonal:
            warn_tip(0,"{} is not orthogonal, skip it !!!".format(fname))
            continue
        new_struct=move_to_zcenter(struct)
        orig_struct=new_struct.copy()
        new_struct =new_struct.copy()
        natom=orig_struct.num_sites
        lat=orig_struct.lattice.matrix.copy()
        pos=orig_struct.frac_coords
        phi=np.linspace(0,360,nps)*np.pi/180
        vzz=C12/C22
        temp_num=(C11*C22-C12**2)/(C22*C66)
        d1=C11/C22 +1.0-temp_num;
        d2=-(2.0*C12/C22-temp_num);
        d3=C11/C22;
        F=sigma*C22/(C11*C22-C12**2.0); 
        Poisson=(vzz*(np.cos(phi))**4.0-d1*(np.cos(phi))**2.0*(np.sin(phi))**2.0+vzz*(np.sin(phi))**4.0)/\
                ((np.cos(phi))**4.0+d2*(np.cos(phi))**2.0*(np.sin(phi))**2.0+d3*(np.sin(phi))**4.0)
        
        eps_theta=F*((np.cos(phi))**4+d2*(np.cos(phi))**2.0*(np.sin(phi))**2.0+d3*(np.sin(phi))**4.0) 
        t = sympy.Symbol('t', real=True)
        e = sympy.Symbol('e', real=True)
        v = sympy.Symbol('v', real=True)
        eprim=sympy.Matrix([[e+1,0],[ 0,1-e*v]])
        R=sympy.Matrix([[sympy.cos(t),-sympy.sin(t)],[ sympy.sin(t),sympy.cos(t)]])
        eps_mat=R*eprim*R.adjugate()
        for k in range(len(phi)):
            cur__phi=phi[k]*180/np.pi 
            Rot=eps_mat.subs({e:eps_theta[k],v:Poisson[k],t:phi[k]})
            filename=str(k)+"_"+fname+'.vasp'
            final_lat=np.matrix(np.eye(3))
            final_lat[0,0]=Rot[0,0]
            final_lat[0,1]=Rot[0,1]
            final_lat[1,0]=Rot[1,0]
            final_lat[1,1]=Rot[1,1]
            lat_new=lat*final_lat
            tmp_struct=Structure(lat_new,new_struct.species,pos)
            tmp_struct.to(filename=filename,fmt='poscar') 

def twod_operation(choice):
    assert choice in ["1","2","3","4","5","6","7"] 
    if choice=="1":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        print("input the supercell supercell factor")
        print('for x direction can be: 10 1 1')
        print('for y direction can be: 1 10 1')
        wait_sep()
        scale_str=wait()
        supercell=[int(x) for x in scale_str.split()]
        if (len(supercell)!=3):
            print('Unknown format')
            os._exit(0)

        if supercell[0]>=supercell[1]:
            direction=0  # for x direction
        else:
            direction=1  # for y direction

        print("input the strain range")
        print("example: 0.02:0.1:10 ")
        wait_sep()
        strain_str=wait()
        tmp=[float(x) for x in strain_str.split(":")]
        strain=[]
        if len(tmp)==3:
             strain_range=np.linspace(tmp[0],tmp[1],int(tmp[2]))
        else:
            print("Unknown format")
            os._exit(0)
        print("input the index of atom need to be fixed")
        print("example: 1 10 11 20 ")
        print("0 means fix the atom automatically")
        wait_sep()
        atom_index_str=wait()
        if len(atom_index_str.split())>1:
            atom_index=[int(x) for x in atom_index_str.split("")]
            auto_fix=False
        else:
            atom_index=[int(atom_index_str)]
            auto_fix=True
        ripple(structs,fnames,supercell,direction,strain_range,atom_index,auto_fix)
        return True

    elif choice=="2":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        print("Input the number of layers")
        wait_sep()
        in_str=wait()
        layer_number=int(in_str)

        print("Input the layer distance")
        wait_sep()
        in_str=wait()
        layer_distance=float(in_str)
        multi_layers(structs,fnames,layer_number,layer_distance)

        return True
 
    elif choice=="3":
        """
          splitting sep
          .... .  .  .   .   .
        """
        structs,fnames=read_structures()
        if len(structs)>1:
           print("Splitting dont support multi-structures!!!")
           os._exit(0)
        atom_index,in_str=atom_selection(structs[0])
        print("Input the splitting distance, 10 Ang is enough!")
        wait_sep()
        in_str=wait()
        SplitDistance=float(in_str)
        print("Numbers of splitting site, 50 sites are enough!")
        wait_sep()
        in_str=wait()
        NumberSplitSite=int(in_str)
        split(structs,fnames,atom_index,in_str,SplitDistance,NumberSplitSite,ProperDist=3.5,DenseFrac=0.75)
        return True

    elif choice=="4":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        print("Input the new value of vacuum layer thickness")
        wait_sep()
        in_str=wait()
        wait_sep()
        nvac_layer_thickness=float(in_str)
        resize_vacuum(structs,fnames,nvac_layer_thickness)
        return True

    elif choice=="5":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        for struct,fname in zip(structs,fnames):
            new_struct=move_to_zcenter(struct)
            new_struct.to(filename='z-center_'+fname+'.vasp',fmt='poscar')
        return True

    elif choice==6:
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        try:
           import sympy
        except ImportError:
           print("You have to install sympy module")
           os._exit(0)

        print("Input the elastic of material by order : C11 C12 C22 C66")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip().split()
        elastic_constant=[float(x) for x in in_str] 
        if len(elastic_constant) !=4:
           print("You have to input C11 C12 C22 C66")
           return None
        C11=elastic_constant[0]
        C12=elastic_constant[1]
        C22=elastic_constant[2]
        C66=elastic_constant[3]
        print("Input applied force: e.x. 1.0 GPa nm")
        wait_sep()
        in_str=wait()
        sigma=float(in_str)
        circle_strain(structs,fnames,C11,C12,C22,C66,sigma)
        return True

    elif choice=="7":
        structs,fnames=read_structures()
        if len(structs)>1:
           print("Constrain doesnot support multi-structures!!!")
           os._exit(0)
        atom_index,in_str=atom_selection(structs[0])
        mlog.debug("constrained atom index")
        mlog.debug(' '.join(map(str,atom_index)))
        constrain(structs,fnames,atom_index,in_str)
        return True

    elif choice=="8":
        print('your choice ?')
        print('{} >>> {}'.format('1','input 2D structure from local disk'))
        print('{} >>> {}'.format('2','get 2D structure online'))
        wait_sep()
        in_str=wait()
        _choice=int(in_str)

        if _choice=="1":
           mpid=None
           structs,fnames=read_structures()
           if len(structs)>1:
               print("Matching dont support multi-structures!!!")
               os._exit(0)
        else:
           print("Input the mp-id for your structure")
           wait_sep()
           in_str=wait()
           mpid=in_str
           struct=None
          
        film,substrates=get_mp_film_substrate(mpid=mpid,struct=struct)
        df=match_substrate(film,substrates)
        dumpfn(df.to_dict(),'substrate_'+fnames[0]+'.json',indent=4)
        #df.to_csv('substrate.csv', sep=',', header=True, index=True)
        return True
    else:
        print("Unkonw choice")
        return None


if __name__=="__main__":
   # structs,fnames=read_structures()
   # if len(structs)>1:
   #    print("Constrain doesnot support multi-structures!!!")
   #    os._exit(0)
   # atom_index,in_str=atom_selection(structs[0])
   # mlog.debug("constrained atom index")
   # mlog.debug(' '.join(map(str,atom_index)))
   # constrain(structs,fnames,atom_index,in_str="")
    st=Structure.from_file("POSCAR")
    constrain([st],['POSCAR'],atom_index=[0,1,2],in_str="")  
