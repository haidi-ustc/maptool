#!/usr/bin/env python
import math
from pymatgen import Structure
from monty.serialization import dumpfn,loadfn
from maptool.util.utils import wait,wait_sep,multi_structs
from maptool.io.read_structure import read_structures,ase2pmg


def move_to_zcenter(struct):
    frac=struct.frac_coords
    zcenter=np.mean(frac[:,2])
    zshift=0.5-zcenter
    frac[:,2]=frac[:,2]+zshift
    return Structure(struct.lattice,struct.species,frac)
    
def twod_operation(choice):
    assert choice in ["1","2","3","4","5","6","7"] 
    if choice=="1":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        margin_dist=1.0 # 
        assert(struct.lattice.is_orthogonal)
        print("input the supercell scaling factor")
        print('for x direction can be: 10 1 1')
        print('for y direction can be: 1 10 1')
        wait_sep()
        scale_str=wait()
        scaling=[int(x) for x in scale_str.split()]
        if (len(scaling)!=3):
            print('unknow format')
        else:
            if scaling[0]>=scaling[1]:
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
        #     print(strain_range)
        else:
            print("unknow format")
            return
        print("input the index of atom need to be fixed")
        print("example: 1 10 11 20 ")
        print("0 means fix the atom automatically")
        wait_sep()
        atom_index_str=wait()
        try:
            atom_index=[int(x) for x in strain_str.split("")]
            auto_fix=False
        except:
            atom_index=[int(atom_index_str)]
            auto_fix=True
        
        struct_sc=struct.copy()
        struct_sc.make_supercell(scaling)
        natom=struct_sc.num_sites
        min_z=np.min(struct_sc.cart_coords[:,2])
        tmp_coords=np.ones((struct_sc.num_sites,1))*min_z
        cart_coords=struct_sc.cart_coords;
        cart_coords[:,2]=cart_coords[:,2]-tmp_coords.T+0.01
        frac_coords_new=np.dot(cart_coords,np.linalg.inv(struct_sc.lattice.matrix))
        for i_strain in strain_range:
            new_lat_matrix=struct_sc.lattice.matrix.copy()
            new_lat_matrix[direction,direction]=struct_sc.lattice.matrix[direction,direction]*(1-i_strain)
            fname="%10.5f"%(i_strain)+'_wo.vasp' # structure only applied with in-plan strain
            struct_wo_ripple=Structure(new_lat_matrix,struct_sc.species,frac_coords_new)
            struct_wo_ripple.to(filename=fname.strip(),fmt='poscar')
            frac_coords_new_cp=struct_wo_ripple.frac_coords.copy()
            cart_coords_new_cp=struct_wo_ripple.cart_coords.copy()
            nz=0
            selective_dynamics=[[True for col in range(3)] for row in range(natom)]  
            z_shift=np.zeros((natom,3))
            for i_atom in range(natom):
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
            fname="%10.5f"%(i_strain)+'_w.vasp' # structure  applied with in-plan strain  and ripple    
            struct_w_ripple.to(filename=fname.strip(),fmt='poscar')        
        return True

    elif choice=="2":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        print("input the number of layers")
        wait_sep()
        in_str=wait()
        layer_number=int(in_str)

        print("input the layer distance")
        wait_sep()
        in_str=wait()
        species=[]
        layer_distance=float(in_str)

        new_struct=move_to_zcenter(struct)
        struct_thickness=np.max(new_struct.cart_coords[:,2])-np.min(new_struct.cart_coords[:,2])
        natom=new_struct.num_sites
        new_cart_coords=np.zeros((natom*layer_number,3))
        for i in range(layer_number):
            new_cart_coords[i*natom:(i+1)*natom,0:2]=new_struct.cart_coords[:,0:2]
            new_cart_coords[i*natom:(i+1)*natom,2]=new_struct.cart_coords[:,2]+i*(layer_distance+struct_thickness)
            species.extend(new_struct.species)
        new_lat=new_struct.lattice.matrix
        new_lat[2,2]=new_lat[2,2]+layer_distance*layer_number
        tmp_struct=Structure(new_lat,species,new_cart_coords,coords_are_cartesian=True)
        tmp1_struct=move_to_zcenter(tmp_struct)
        tmp2_struct=tmp1_struct.get_sorted_structure()
        tmp2_struct.to(filename='layer_'+str(layer_number)+'.vasp',fmt='poscar')
 
    elif choice=="3":
        ProperDist=3.5  # Ang
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        (atom_index,in_str)=atom_selection(struct)
        print("input the splitting distance, 10 Ang is enough!")
        wait_sep()
        in_str=wait()
        SplitDistance=float(in_str)
        print("numbers of splitting site, 50 sites are enough!")
        wait_sep()
        in_str=wait()
        NumberSplitSite=int(in_str)
        DensityN=int(NumberSplitSite*0.75)
        SparseN=NumberSplitSite-DensityN+1
#        print(DensityN,SparseN)
        dist=ProperDist/(DensityN-1)
        SplitDistanceArray=np.zeros(NumberSplitSite+1)        
        for Nsite in range(DensityN):
            SplitDistanceArray[Nsite]=(Nsite)*dist
        
        dist=(SplitDistance-ProperDist)/SparseN
        for  Nsite in range(SparseN):
             SplitDistanceArray[Nsite+DensityN]=ProperDist+(Nsite+1)*dist

#        print(SplitDistanceArray)
        coords=struct.cart_coords
        for Nsite in range(NumberSplitSite+1):
            coords=struct.cart_coords
            for atom in atom_index:
                coords[atom,2]=coords[atom,2]+SplitDistanceArray[Nsite]
            tmp_struct=Structure(struct.lattice,struct.species,coords,coords_are_cartesian=True) 
            fname=str(Nsite)+'.vasp'
            tmp_struct.to(filename=fname,fmt='poscar')
        data=np.zeros((NumberSplitSite+1,2))
        for i,j in enumerate(SplitDistanceArray):
            data[i][0]=i
            data[i][1]=j
        head_line="#%(key1)+12s  %(key2)+12s"%{'key1':'index','key2':'distance/Ang'}
        fmt="%12d %12.6f"+'\n'
        write_col_data('split.dat',data,head_line,sp_fmt=fmt) 
        return

    elif choice=="4":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        new_struct=move_to_zcenter(struct)
        struct_thickness=np.max(new_struct.cart_coords[:,2])-np.min(new_struct.cart_coords[:,2])
        vac_layer_thickness=new_struct.lattice.c-struct_thickness
        print("current vacuum layer thickness is %6.3f Ang"%(vac_layer_thickness))
        print("input the new value of vacuum layer thickness")
        wait_sep()
        in_str=wait()
        nvac_layer_thickness=float(in_str)

        assert(nvac_layer_thickness>struct_thickness)
        new_lat=new_struct.lattice.matrix
        new_lat[2,2]=nvac_layer_thickness
        tmp_struct=Structure(new_lat,new_struct.species,new_struct.cart_coords,coords_are_cartesian=True)
        center_struct=move_to_zcenter(tmp_struct)
        center_struct.to(filename='new_vacuum.vasp',fmt='poscar') 
        return

    elif choice=="5":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        new_struct=move_to_zcenter(struct)
        new_struct.to(filename='z-center.vasp',fmt='poscar')
        return

    elif choice==6:
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        assert(struct.lattice.is_orthogonal)
        try:
           import sympy
        except:
           print("you must install sympy module")
           return None

        new_struct=move_to_zcenter(struct)
        print("input the elastic of material by order : C11 C12 C22 C66")
        wait_sep()
        in_str=""
        while in_str=="":
              in_str=input().strip().split()
        elastic_constant=[float(x) for x in in_str] 
        if len(elastic_constant) !=4:
           print("you must input C11 C12 C22 C66")
           return None
        C11=elastic_constant[0]
        C12=elastic_constant[1]
        C22=elastic_constant[2]
        C66=elastic_constant[3]
        print("input applied force: e.x. 1.0 GPa nm")
        wait_sep()
        in_str=wait()
        sigma=float(in_str)
        #for struct,fname in zip(structs,fnames):
        orig_struct=new_struct.copy()
        new_struct =new_struct.copy()
        natom=orig_struct.num_sites
        lat=orig_struct.lattice.matrix
        pos=orig_struct.frac_coords
        nps=36
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
            fname=str(k)+'.vasp'
            final_lat=np.matrix(np.eye(3))
            final_lat[0,0]=Rot[0,0]
            final_lat[0,1]=Rot[0,1]
            final_lat[1,0]=Rot[1,0]
            final_lat[1,1]=Rot[1,1]
            lat_new=lat*final_lat
            tmp_struct=Structure(lat_new,new_struct.species,pos)
            tmp_struct.to(filename=fname,fmt='poscar') 
        return True

    elif choice=="7":
        structs,fnames=read_structures()
        multi_structs(structs,fnames)
        for struct,fname in zip(structs,fnames):
            natom=struct.num_sites
            atom_index,in_str=atom_selection(struct)
            selective_dynamics=[[True for col in range(3)] for row in range(natom)]  
            for i in range(natom):
                if i in atom_index:
                    selective_dynamics[i]=[False,False,False]
            tmp_struct=Structure(struct.lattice,struct.species,struct.frac_coords,site_properties={'selective_dynamics':selective_dynamics})
            poscar=Poscar(tmp_struct)
            poscar.comment=poscar.comment+' |--> '+in_str
            fname='Fixed_'+fname+'.vasp'
            poscar.write_file(fname)
        return True

    elif choice=="8":
        print('your choice ?')
        print('{} >>> {}'.format('1','input 2D structure from local disk'))
        print('{} >>> {}'.format('2','get 2D structure online'))
        wait_sep()
        in_str=wait()
        _choice=int(in_str)

        if _choice==1:
           mpid=None
           structs,fnames=read_structures()
           multi_structs(structs,fnames)
        else:
           print("input the mp-id for your structure")
           wait_sep()
           in_str=wait()
           mpid=in_str
           struct=None
          
        film,substrates=make_connect(mpid=mpid,struct=struct)
        df=get_subs(film,substrates)
        dumpfn(df.to_dict(),'substrate.json',indent=4)
        #df.to_csv('substrate.csv', sep=',', header=True, index=True)
        return True
    else:
        print("unkonw choice")
        return None

