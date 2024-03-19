from interface import *
from ops_framework import *



def generate_c_code(c_code, file_path):
    with open(file_path, 'w') as file:
        file.write(c_code)

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

