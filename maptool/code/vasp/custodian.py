from maptool.util.utils import sepline,box_center
error_msgs = {
    "tet": [
        "Tetrahedron method fails for NKPT<4",
        "Fatal error detecting k-mesh",
        "Fatal error: unable to match k-point",
        "Routine TETIRR needs special values",
        "Tetrahedron method fails (number of k-points < 4)",
    ],
    "inv_rot_mat": [
        "inverse of rotation matrix was not found (increase " "SYMPREC)"
    ],
    "brmix": ["BRMIX: very serious problems"],
    "subspacematrix": ["WARNING: Sub-Space-Matrix is not hermitian in " "DAV"],
    "tetirr": ["Routine TETIRR needs special values"],
    "incorrect_shift": ["Could not get correct shifts"],
    "real_optlay": ["REAL_OPTLAY: internal error", "REAL_OPT: internal ERROR"],
    "rspher": ["ERROR RSPHER"],
    "dentet": ["DENTET"],
    "too_few_bands": ["TOO FEW BANDS"],
    "triple_product": ["ERROR: the triple product of the basis vectors"],
    "rot_matrix": ["Found some non-integer element in rotation matrix"],
    "brions": ["BRIONS problems: POTIM should be increased"],
    "pricel": ["internal error in subroutine PRICEL"],
    "zpotrf": ["LAPACK: Routine ZPOTRF failed"],
    "amin": ["One of the lattice vectors is very long (>50 A), but AMIN"],
    "zbrent": ["ZBRENT: fatal internal in", "ZBRENT: fatal error in bracketing"],
    "pssyevx": ["ERROR in subspace rotation PSSYEVX"],
    "eddrmm": ["WARNING in EDDRMM: call to ZHEGV failed"],
    "edddav": ["Error EDDDAV: Call to ZHEGV failed"],
    "grad_not_orth": ["EDWAV: internal error, the gradient is not orthogonal"],
    "nicht_konv": ["ERROR: SBESSELITER : nicht konvergent"],
    "zheev": ["ERROR EDDIAG: Call to routine ZHEEV failed!"],
    "elf_kpar": ["ELF: KPAR>1 not implemented"],
    "elf_ncl": ["WARNING: ELF not implemented for non collinear case"],
    "rhosyg": ["RHOSYG internal error"],
    "posmap": ["POSMAP internal error: symmetry equivalent atom not found"],
    "point_group": ["Error: point group operation missing"],
}
box_center(ch='_',fill='_',sp=' ')
box_center(ch="Vasp Errors",fill=' ',sp='|')
box_center(ch="",fill=' ',sp='|')
for k in  error_msgs.keys():
    box_center(ch=k,fill='-',sp="|")
    print('\n'.join(error_msgs[k]))
box_center(ch='_',fill='_',sp='|')
