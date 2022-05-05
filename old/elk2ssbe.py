#!/usr/bin/env python3
import numpy as np

# Parse Elk's output file and convert to SALMON results

def read_geometory(filename="GEOMETRY.OUT"):
    avec = np.empty([3, 3])

    print("# Retrieve %s ..." % filename)
    with open(filename, "r") as fh:
        # Skip 14 lines
        for i in range(14):
            line = fh.readline()
        # Read vector
        for i in range(3):
            line = fh.readline()
            tmp = line.split()
            avec[i, 0] = float(tmp[0])
            avec[i, 1] = float(tmp[1])
            avec[i, 2] = float(tmp[2])
            print("# i=%d avec=%+.3f,%+.3f,%+.3f"
                % (i, avec[i, 0], avec[i, 1], avec[i, 2]))
    # Results
    return avec


# Read "KPOINTS.OUT"
def read_kpoints(filename="KPOINTS.OUT"):
    print("# Retrieve %s ..." % filename)
    with open(filename, "r") as fh:
        # Line 1
        line = fh.readline()
        nkpt = int(line.split(":")[0])
        print("# nkpt=%d" % nkpt)

        vkl = np.empty([nkpt, 3], dtype=float)
        wkpt = np.empty([nkpt], dtype=float)

        for ik in range(nkpt):
            line = fh.readline()
            tmp = line.split()
            jk = int(tmp[0])
            vkl[ik, 0] = float(tmp[1])
            vkl[ik, 1] = float(tmp[2])
            vkl[ik, 2] = float(tmp[3])
            wkpt[ik] = float(tmp[4])
            # print("# ik=%d jk=%d vkl=%+.3f,%+.3f,%+.3f wkpt=%.3f"
            #     % (ik, jk, vkl[ik, 0], vkl[ik, 1], vkl[ik, 2], wkpt[ik]))
    # Results
    return nkpt, vkl, wkpt


# Read "EIGVAL.OUT"
def read_eigval(filename="EIGVAL.OUT"):
    print("# Retrieve %s ..." % filename)
    with open(filename, "r") as fh:
        # Line 1
        line = fh.readline()
        nkpt = int(line.split(":")[0])
        print("# nkpt=%d" % nkpt)
        # Line 2
        line = fh.readline()
        nstsv = int(line.split(":")[0])
        print("# nstsv=%d" % nstsv)

        vkl = np.empty([nkpt, 3], dtype=float)
        eig = np.empty([nkpt, nstsv], dtype=float)
        occ = np.empty([nkpt, nstsv], dtype=float)

        for ik in range(nkpt):
            # Skip single line
            line = fh.readline()
            # Header
            line = fh.readline()
            tmp = line.split()
            jk = int(tmp[0])
            vkl[ik, 0] = float(tmp[1])
            vkl[ik, 1] = float(tmp[2])
            vkl[ik, 2] = float(tmp[3])
            # print("# ik=%d jk=%d kpt=%+.3f,%+.3f,%+.3f"
            #     % (ik, jk, vkl[ik, 0], vkl[ik, 1], vkl[ik, 2]))
            # Skip single line
            line = fh.readline()
            for ib in range(nstsv):
                # Eigenenergy
                line = fh.readline()
                tmp = line.split()
                jb = int(tmp[0])
                eig[ik, ib] = float(tmp[1])
                occ[ik, ib] = float(tmp[2])
                # print("# ib=%d jb=%d eig=%+.3f occ=%.3f"
                #     % (ib, jb, eig[ik, ib], occ[ik, ib]))
            # Skip single line
            line = fh.readline()
    # Results    
    return nkpt, nstsv, vkl, eig, occ


def read_pmat(nkpt, nstsv, filename="PMAT.OUT"):
    vkl = np.empty([nkpt, 3])
    pmat = np.empty([nkpt, nstsv, nstsv, 3], dtype=np.complex128)

    print("# Retrieve %s ..." % filename)
    with open("PMAT.OUT", "rb") as fh:
        for ik in range(nkpt):
            vkl[ik, :] = np.fromfile(fh, dtype=np.float64, count=3)
            nstsv_tmp = np.fromfile(fh, dtype=np.int32, count=1)
            pmat_tmp = np.fromfile(fh, dtype=np.complex128, count=nstsv*nstsv*3)
            pmat[ik, :, :, :] = pmat_tmp.reshape([nstsv, nstsv, 3])
    # Results
    return vkl, pmat


def calc_nelec(nkpt, nstsv, occ, wkpt):
    tmp = 0.0
    for ik in range(nkpt):
        tmp += wkpt[ik] * sum(occ[ik, :]) 
    nelec = int(np.rint(tmp / sum(wkpt)))
    return nelec


header_k = """# k-point distribution
# ik: k-point index
# kx,ky,kz: Reduced coordinate of k-points
# wk: Weight of k-point
# 1:ik[none] 2:kx[none] 3:ky[none] 4:kz[none] 5:wk[none]"""

header_eigen = """#esp: single-particle energies (eigen energies)
#occ: occupation numbers, io: orbital index
# 1:io, 2:esp[a.u.], 3:occ
k=     1,  spin=     1"""

header_tm = """# #Transition Moment between occupied and unocupied orbitals in GS
# # (Separated analysis tool is available)
 #<u_nk|p_j|u_mk>  (j=x,y,z)"""

template_inp = """
&calculation
    theory = 'sbe_pulse'
/

&control
    sysname = 'elk'
    base_directory = "./"
    gs_directory = "./"
/

&system
    al_vec1(1:3) = {AX1}, {AY1}, {AZ1}
    al_vec2(1:3) = {AX2}, {AY2}, {AZ2}
    al_vec3(1:3) = {AX3}, {AY3}, {AZ3}
    nstate = {NSTSV}
    nelec = {NELEC}
    nstate_sbe = {NSTSV}
/

&kgrid
    nkgrid(1:3) = {NKPT}, 1, 1
/
"""


def main():
    avec = read_geometory()
    nkpt, vkl, wkpt = read_kpoints()
    nkpt, nstsv, vkl, eig, occ = read_eigval()
    vkl, pmat = read_pmat(nkpt, nstsv)

    print("# Generate elk_k.data")
    with open("elk_k.data", "w") as fh:
        fh.write(header_k)
        for ik in range(nkpt):
            fh.write("%d %e %e %e %e\n" % (ik+1, vkl[ik, 0], vkl[ik, 1], vkl[ik, 2], wkpt[ik]))

    print("# Generate elk_eigen.data")
    with open("elk_eigen.data", "w") as fh:
        fh.write(header_eigen)
        for ik in range(nkpt):
            fh.write("k=     %d,  spin=     1\n" % (ik + 1))
            for ib in range(nstsv):
                fh.write("%d %e %e\n" % (ib + 1, eig[ik, ib], occ[ik, ib]))

    print("# Generate elk_tm.data")
    with open("elk_tm.data", "w") as fh:
        fh.write(header_tm)
        for ik in range(nkpt):
            for ib in range(nstsv):
                for jb in range(nstsv):
                    fh.write("%d %d %d %e %e %e %e %e %e\n" % (
                        ik + 1, ib + 1, jb + 1,
                        np.real(pmat[ik, ib, jb, 0]), np.imag(pmat[ik, ib, jb, 0]),
                        np.real(pmat[ik, ib, jb, 1]), np.imag(pmat[ik, ib, jb, 1]),
                        np.real(pmat[ik, ib, jb, 2]), np.imag(pmat[ik, ib, jb, 2]),
                    ))
        fh.write("# nonlocal component is not exist\n")
        for ik in range(nkpt):
            for ib in range(nstsv):
                for jb in range(nstsv):
                    fh.write("%d %d %d %e %e %e %e %e %e\n" % (
                        ik + 1, ib + 1, jb + 1,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                    ))


                

    nelec = calc_nelec(nkpt, nstsv, occ, wkpt)


    print(template_inp.format(
        NKPT=nkpt, NSTSV=nstsv, NELEC=nelec,
        AX1 = avec[0, 0], AY1 = avec[0, 1], AZ1 = avec[0, 2],
        AX2 = avec[1, 0], AY2 = avec[1, 1], AZ2 = avec[1, 2],
        AX3 = avec[2, 0], AY3 = avec[2, 1], AZ3 = avec[2, 2],
    ))