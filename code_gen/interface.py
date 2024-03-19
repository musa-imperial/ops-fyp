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

#ops stencil def, stencil name, no. stencil points
stencils_def = [
("{0,0}", "00", 1),
("{0,0,   1,0,  0,1}", "bottom_left", 3),
("{0,0,  1,0, 0,-1}", "top_left", 3),
("{0,0,  -1,0, 0,1}", "bottom_right", 3),
("{0,0, -1,0,  0,-1}", "top_right", 3),
("{0,0, 1,0, 0,1, 0,-1}", "left", 4),
("{0,0, -1,0, 0,1, 0,-1}", "right", 4),
("{0,0, 1,0, -1,0, 0,-1}", "top", 4),
("{0,0, 1,0, -1,0, 0,1}", "bottom", 4),
("{0,0, 1,0, -1,0, 0,1, 0,-1}", "interior", 5),
]

def stencils():
    code = """
"""
    for p in stencils_def:
        code +=f"""
    int s2d_{p[1]}[] = {p[0]}; 
    ops_stencil S2D_{p[1].upper()} = ops_decl_stencil(2,{p[2]},s2d_{p[1]},"{p[1]}");

    """
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