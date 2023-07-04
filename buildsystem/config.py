import re 
import sys
import os


def get_options(fin):
    options = {}
    for line in fin:
        s = line.split()
        if(len(s)>0):
            if(s[0][0] != "#"):
                val = s[0].split("=")
                if len(val) > 1:
                    options[val[0]] =  val[1]
                else:
                    options[val[0]] = None
    
    return options
    
def out1(options, fname):
    f = open(fname, "w")
    
    keys = list(options.keys())
    keys.sort()
    
    for key in keys:
        if options[key] is None:
            f.write("#define " + key + "\n")
        else:
            f.write("#define " + key + " " + options[key] + "\n")
            
def out2(options, fname):
    f = open(fname, "w")
    
    str = """
#include \"gadgetconfig.h\"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include \"data/dtypes.h\"
#include \"data/allvars.h\"
#include \"main/main.h\"
void output_compile_time_options(void)\n{
printf(
"""
    f.write(str)
    
    keys = list(options.keys())
    keys.sort()
    
    for key in keys:
        if options[key] is None:
            f.write("\"    " + key + "\\n\"\n")
        else:
            f.write("\"    " + key + "=" + options[key] + "\\n\"\n")
            
    str = """);
}"""
    f.write(str)
    
def out3(options, fname):
    f = open(fname, "w")
    
    str = """
#include \"gadgetconfig.h\"
#include <mpi.h>
#include <stdio.h>
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include \"data/constants.h\"
#include \"data/dtypes.h\"
#include \"data/macros.h\"
#include \"io/io.h\"
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);
herr_t my_H5Tclose(hid_t type_id);

void IO_Def::write_compile_time_options_in_hdf5(hid_t handle)
{
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);\n
"""
    f.write(str)
    
    keys = list(options.keys())
    keys.sort()
    
    for key in keys:
        f.write("hdf5_dataspace = my_H5Screate(H5S_SCALAR);")
            
        if options[key] is None:
            f.write("hdf5_attribute = my_H5Acreate(handle, \"" + key + "\" , atype, hdf5_dataspace, H5P_DEFAULT);\n")
            f.write("my_H5Awrite(hdf5_attribute, atype, \"\", \"" + key + "\");\n")
        else:
            f.write("hdf5_attribute = my_H5Acreate(handle, \"" + key + "\" , H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);\n")
            f.write("val = " + options[key] + ";\n")
            f.write("my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, \"" + key + "\");\n")
            
        f.write("my_H5Aclose(hdf5_attribute, \"" + key + "\");\n")
        f.write("my_H5Sclose(hdf5_dataspace, H5S_SCALAR);\n\n")

    f.write("my_H5Tclose(atype);\n")
    f.write("}\n")
    f.write("\n")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python check.py <Config.sh> <build_dir> <curr_dir> <src_dir>")
        exit(1)
        
       
    fin = open(sys.argv[1],'r')
    fout = sys.argv[2]
    
    options = get_options(fin)
    
    out1(options, sys.argv[2] + "/gadgetconfig.h")
    out2(options, sys.argv[2] + "/compile_time_info.cc")
    out3(options, sys.argv[2] + "/compile_time_info_hdf5.cc")
