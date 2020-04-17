* Pymatgen Structure/Molecule instance are used as core structure for different routing, 
if you have other class to parse your structure, the final structure should be converted to Pymatgen one.
* All of building and analysing modulea are in maptool/core
* The basic data output should use io/DataIO class, the json or yaml format is also accept.
* Try to avoid `goto` command, unless you have no choice.
* The final functional function should only accept the structure and other parameters, no keyboard input function should be include in it.
