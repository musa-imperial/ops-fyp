code_name = "reaction_diffusion"

indepedant_constants = [
("double", "dt", 0.001),
("double", "T", 100),
("int", "Nx", 251),
("int", "Ny", 251),
("double", "a", 0.75),
("double", "b", 0.06),
("double", "mu1", 5.0),
("double", "mu2", 0.0),
("double", "eps", 13.0),
("double", "dx", 1.0),
("double", "dy", 1.0)]

dependent_constants = [
("double", "Lx", "dx*(Nx-1)"),
("double", "Ly", "dy*(Ny-1)"),
("double", "h", "dx"),
("double", "hmu1dt", "mu1/h/h*dt"),
("double", "hmu2dt", "mu2/h/h*dt"),
("double", "div_a", "1/a")]

data_variables = ["u", "v"]
data_variables_copy = ["u_calc", "v_calc"]

kernel_names = [
    "bottom_left", 
    "bottom_right",
    "top_left",
    "top_right",
    "left",
    "right",
    "top",
    "bottom",
]


def generate_c_code(c_code, file_path):
    with open(file_path, 'w') as file:
        file.write(c_code)


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

def stencils():
    code = """
    int s2d_00[] = {0,0};
    ops_stencil S2D_00 = ops_decl_stencil(2,1,s2d_00,"0,0");

    //corner bc stencils
    //point x = 0,  y = 0 stencil (bottom_left)
    int s2d_bottom_left[] = {0,0,   1,0,  0,1};
    ops_stencil S2D_BOTTOM_LEFT = ops_decl_stencil(2,3,s2d_bottom_left,"bottom_left");

    //point x = 0,  y = Ny-1 stencil (top_left)
    int s2d_top_left[] = {0,0,  1,0, 0,-1};
    ops_stencil S2D_TOP_LEFT = ops_decl_stencil(2,3,s2d_top_left,"top_left");

    //point x = Nx-1, y = 0 stencil (bottom_right)
    int s2d_bottom_right[] = {0,0,  -1,0, 0,1};
    ops_stencil S2D_BOTTOM_RIGHT = ops_decl_stencil(2,3,s2d_bottom_right,"bottom_right");

    //point x = Nx-1, y = Ny-1 stencil (top_right)
    int s2d_top_right[] = {0,0, -1,0,  0,-1};
    ops_stencil S2D_TOP_RIGHT = ops_decl_stencil(2,3,s2d_top_right,"top_right");

    //bc stencils
    int s2d_left[] = {0,0, 1,0, 0,1, 0,-1};
    ops_stencil S2D_LEFT = ops_decl_stencil(2,4,s2d_left,"left");

    int s2d_right[] = {0,0, -1,0, 0,1, 0,-1};
    ops_stencil S2D_RIGHT = ops_decl_stencil(2,4,s2d_right,"right");

    int s2d_top[] = {0,0, 1,0, -1,0, 0,-1};
    ops_stencil S2D_TOP = ops_decl_stencil(2,4,s2d_top,"top");

    int s2d_bottom[] = {0,0, 1,0, -1,0, 0,1};
    ops_stencil S2D_BOTTOM = ops_decl_stencil(2,4,s2d_bottom,"bottom");

    int s2d_interior[] = {0,0, 1,0, -1,0, 0,1, 0,-1};
    ops_stencil S2D_INTERIOR = ops_decl_stencil(2,5,s2d_interior,"interior");
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

def initial_condition():

    code = """
"""

    for name in data_variables:
        code += f"""
    ops_par_loop({name}_initcond_stencil, "{name}_initcond_stencil", block, 2, all,
            ops_arg_dat(d_{name},    1, S2D_00, "double", OPS_WRITE),
            ops_arg_idx());
"""

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

def boundary():

    code = """
    
    """
    for k in kernel_names:
        for index, name in enumerate(data_variables):
            code += f"""ops_par_loop({k}_{name}, "{k}_{name}", block, 2, {k},
                ops_arg_dat(d_{name},    1,   S2D_{k.upper()}, "double", OPS_READ),
                ops_arg_dat(d_{data_variables_copy[index]}, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_{data_variables[index-1]},    1,   S2D_00, "double", OPS_READ));

    """

    return code

def interior():
    code = """
    """
    for index, name in enumerate(data_variables):
        code += f"""ops_par_loop(interior_stencil_{name}, "interior_stencil_{name}", block, 2, interior,
        ops_arg_dat(d_{name},    1,   S2D_INTERIOR, "double", OPS_READ),
        ops_arg_dat(d_{data_variables_copy[index]}, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_{data_variables[index-1]},    1,   S2D_00, "double", OPS_READ));

    """
    return code

def results():
    code = """
    """
    for p in data_variables:
        code += f"""
        ops_print_dat_to_txtfile(d_{p}, "{p}_results.txt");   
        """

    return code
    
def finish():
    code = """
    //Finalising the OPS library

    ops_exit();
    free(temp);
    """
    return code

def main(code):

    code += build_datasets()
    code += decl_global_constants()
    code += stencils()

    code += """
    ops_partition("");"""

    code += loop_ranges()

    code += set_zero()

    code += initial_condition()

    code +=FOR("t", "0", "T", "dt")

    code += boundary()

    code += interior()

    code += copy()

    code +=END_SCOPE()

    code += results()

    code += finish()

    code += """
    return 0;
}
"""
    return code



if __name__ == "__main__":
    output_file = "reaction_diffusion.cpp"

    c_code = header()
    c_code = main(c_code)
    generate_c_code(c_code, output_file)
    print(f"C code generated and saved to {output_file}")

