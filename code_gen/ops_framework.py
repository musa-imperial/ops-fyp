from interface import *

def global_params():
    code = """
"""
    for p in indepedant_constants:
        code += f"""{p[0]} {p[1]} = {p[2]};\n"""

    for p in dependent_constants:
        code += f"""{p[0]} {p[1]} = {p[2]};\n"""

    return code 



def header():

    code = """
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
"""

    code += global_params()


    code += """
#define OPS_2D
#include <ops_seq_v2.h>
"""

    code += f"""#include "{code_name}_kernels.h" """

    code +="""
int main(int argc, const char** argv)
{
    ops_init(argc, argv,1);
"""

    return code


def datasets():
    code ="""
    """

    for name in data_variables:
        code += f"""ops_dat d_{name}  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp, "double", "{name}");
    """

    for name in data_variables_copy:
        code += f"""ops_dat d_{name}  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp, "double", "{name}");
    """

    return code


def build_datasets():
    code = """
    // block
    ops_block block = ops_decl_block(2, "2D_grid");
    
    // Block params
    int size[] = {Nx, Ny};
    int base[] = {0,0};
    int d_m[] =  {-1,-1};
    int d_p[] =  {1,1};
    double* temp = NULL;
"""
    code += datasets()

    return code

def decl_global_constants():

    code = """
    """

    for p in indepedant_constants:
        code += f"""ops_decl_const("{p[1]}",1,"{p[0]}",&{p[1]});
    """
    code += """\n"""
    for p in dependent_constants:
        code += f"""ops_decl_const("{p[1]}",1,"{p[0]}",&{p[1]});
    """
    
    return code


def finish():
    code = """
    //Finalising the OPS library

    ops_exit();
    free(temp);
    """
    return code

def loop_ranges():
    code = """

    //loop ranges
    int bottom_left[] = {-1, 0, -1, 0};

    int bottom_right[] = {Nx, Nx+1, -1, 0};

    int top_left[] = {-1, 0, Ny, Ny+1};

    int top_right[] = {Nx, Nx+1, Ny, Ny+1};

    int bottom[] = {0, Nx, -1, 0};

    int top[] = {0, Nx, Ny, Ny+1};

    int left[] = {-1, 0, 0, Ny};

    int right[] = {Nx, Nx+1, 0, Ny};

    int interior[] = {0, Nx, 0, Ny};

    int all[] = {-1, Nx+1, -1, Ny+1};
"""
    return code

def set_zero():

    code = """
    """

    for name in data_variables: 
        code += f"""
    ops_par_loop(set_zero, "set_zero", block, 2, all,
            ops_arg_dat(d_{name}, 1, S2D_00, "double", OPS_WRITE));\n"""

    for name in data_variables_copy: 
        code += f"""
    ops_par_loop(set_zero, "set_zero", block, 2, all,
            ops_arg_dat(d_{name}, 1, S2D_00, "double", OPS_WRITE));\n"""
    
    return code

def copy():
    code = """
    """

    for index, name in enumerate(data_variables):
        code += f"""
    ops_par_loop(copy, "copy", block, 2, all,
        ops_arg_dat(d_{name},    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_{data_variables_copy[index]}, 1, S2D_00, "double", OPS_READ));
"""
    return code
    
def FOR(i, start, finish, inc):

    code = """
"""

    code += f"""for (double {i} = {start}; {i} < {finish}; {i} += {inc}) {{"""

    return code

def END_SCOPE():
    return """
    }"""

def results():
    code = """
    """
    for p in data_variables:
        code += f"""
        ops_print_dat_to_txtfile(d_{p}, "{p}_results.txt");   
        """

    return code