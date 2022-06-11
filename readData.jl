using HDF5

function loadHDF5(filename="MHD_256_Mag.h5")
    """Load magnetic field from HDF5 data"""
    f = h5open(filename, "r")
    mm = read(f, "maxmin")'
    Bx = read(f, "BX") * ((mm[1, 1] - mm[2, 1]) / 256 + mm[2, 1])
    By = read(f, "BX") * ((mm[1, 2] - mm[2, 2]) / 256 + mm[2, 2])
    Bz = read(f, "BX") * ((mm[1, 3] - mm[2, 3]) / 256 + mm[2, 3])
    close(f)

    return(Bx, By, Bz)
end