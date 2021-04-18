#!/usr/bin/env python

f90_file = "input_parameter.f90"

input_parameter = [
    {
        "name": "calculation",
        "variables": [
            {"name": "theory", "type": "character(256)", "default": "'perturb_dielec'"},
        ]
    },
    {
        "name": "control",
        "variables": [
            {"name": "sysname", "type": "character(256)", "default": "'test'"},
            {"name": "base_directory", "type": "character(256)", "default": "'./'"},
            {"name": "gs_directory", "type": "character(256)", "default": "'./'"},
            {"name": "read_sbe_gs_bin", "type": "character(256)", "default": "'n'"},
            {"name": "write_sbe_gs_bin", "type": "character(256)", "default": "'y'"},
        ]
    },
    {
        "name": "system",
        "variables": [
            {"name": "al", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "al_vec1", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "al_vec2", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "al_vec3", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "nstate", "type": "integer", "default": "0"},
            {"name": "nelec", "type": "integer", "default": "0"},
            {"name": "nstate_sbe", "type": "integer", "default": "0"},
        ]
    },
    {
        "name": "kgrid",
        "variables": [
            {"name": "nkgrid", "type": "integer, dimension(3)", "default": "(/0, 0, 0/)"},
        ]
    },
    {
        "name": "tgrid",
        "variables": [
            {"name": "nt", "type": "integer", "default": "1000"},
            {"name": "dt", "type": "real(8)", "default": "1.0d-2"},
        ]
    },
    {
        "name": "emfield",
        "variables": [
            {"name": "e_impulse", "type": "real(8)", "default": "0.0d0"},
            {"name": "ae_shape1", "type": "character(256)", "default": "'none'"},
            {"name": "ae_shape2", "type": "character(256)", "default": "'none'"},
            {"name": "epdir_re1", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "epdir_re2", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "epdir_im1", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "epdir_im2", "type": "real(8), dimension(3)", "default": "(/0.0, 0.0, 0.0/)"},
            {"name": "phi_cep1", "type": "real(8)", "default": "0.0d0"},
            {"name": "phi_cep2", "type": "real(8)", "default": "0.0d0"},
            {"name": "E_amplitude1", "type": "real(8)", "default": "0.0d0"},
            {"name": "E_amplitude2", "type": "real(8)", "default": "0.0d0"},
            {"name": "I_wcm2_1", "type": "real(8)", "default": "0.0d0"},
            {"name": "I_wcm2_2", "type": "real(8)", "default": "0.0d0"},
            {"name": "tw1", "type": "real(8)", "default": "0.0d0"},
            {"name": "tw2", "type": "real(8)", "default": "0.0d0"},
            {"name": "omega1", "type": "real(8)", "default": "0.0d0"},
            {"name": "omega2", "type": "real(8)", "default": "0.0d0"},
            {"name": "t1_t2", "type": "real(8)", "default": "0.0d0"},
            {"name": "t1_start", "type": "real(8)", "default": "0.0d0"},
        ]
    },
    {
        "name": "analysis",
        "variables": [
            {"name": "nenergy", "type": "integer", "default": "1000"},
            {"name": "de", "type": "real(8)", "default": "1.0d-3"},
            {"name": "gamma", "type": "real(8)", "default": "5.0d-3"},
        ]
    },
]


template = """! This file is automatically created by input_parameter.py
module input_parameter
    use salmon_file, only: open_filehandle
    implicit none

{DEF_VARIABLE}

contains

    subroutine read_input()
        implicit none
        integer :: ret, fh
        character(256) :: tmp

{DEF_NAMELIST}

{VARIABLE_DEFAULT}

        fh = open_filehandle('.namelist.tmp')
        do while (.true.)
            read(*, '(a)', iostat=ret) tmp
            if (ret < 0) exit ! End of file
            tmp = adjustl(tmp)
            if (tmp(1:1) .ne. '!') write(fh, '(a)') trim(tmp)
        end do
{READ_NAMELIST}

        close(fh)

{VAR_DUMP}
    end subroutine read_input
end module input_parameter
"""

import os
import sys

def_variable = ""
def_namelist = ""

# DEF_VARIABLES
for group in input_parameter:
    for variable in group['variables']:
        def_variable += (" " * 4) + "%s :: %s\n" % (variable['type'], variable['name'])

# DEF_NAMELIST
def_namelist = ""
for group in input_parameter:
    def_namelist +=  (" " * 8) + "namelist/%s/" % group["name"]
    for variable in group['variables']:
        def_namelist += ", &\n" + (" " * 8) + "& " + variable["name"]
    def_namelist += "\n\n"
def_namelist = def_namelist.replace("/,", "/")

# VARIABLE_DEFAULT
variable_default = ""
for group in input_parameter:
    for variable in group['variables']:
        variable_default += (" " * 8) + "%s = %s\n" % (variable['name'], variable['default'])

# READ_NAMELIST
read_namelist = ""
for group in input_parameter:
    read_namelist += (" " * 8) + "rewind(fh); read(fh, nml=%s, iostat=ret)\n" % group['name']

# VARIABLE_DEFAULT
var_dump = ""
for group in input_parameter:
    var_dump += (" " * 8) + "write(*, '(a)') '# %s'\n" % group["name"]
    for variable in group['variables']:
        if "int" in variable["type"]:
            var_dump += (" " * 8) + "write(*, '(a, 9(i9))') '# %s=', %s\n" % (variable["name"], variable["name"])
        elif "real" in variable["type"]:
            var_dump += (" " * 8) + "write(*, '(a, 9(es23.15e3))') '# %s=', %s\n" % (variable["name"], variable["name"])
        elif "char" in variable["type"]:
            var_dump += (" " * 8) + "write(*, '(a, a)') '# %s=', trim(%s)\n" % (variable["name"], variable["name"])



with open(f90_file, "w") as fh:
    fh.write(template.format(
        DEF_VARIABLE=def_variable,
        DEF_NAMELIST=def_namelist,
        READ_NAMELIST=read_namelist,
        VARIABLE_DEFAULT=variable_default,
        VAR_DUMP=var_dump
    ))

