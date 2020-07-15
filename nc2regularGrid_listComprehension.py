from dfm_tools.get_nc import get_netdata, get_ncmodeldata
from dfm_tools.get_nc_helpers import get_ncvardimlist
from dfm_tools.regulargrid import scatter_to_regulargrid
import os
import numpy as np
from netCDF4 import Dataset
import time as tm


"""original model located at
    p:\11202428-hisea\03-Model\Greece_model\waq_model\
"""


def regularGrid_to_netcdf(fp_in, nx, ny, treg, lreg):
    dir_output = os.path.abspath(os.path.join(os.path.dirname(__file__),'..', 'output'))
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    file_nc = fp_in
    input_nc = Dataset(file_nc, 'r', format='NetCDF4')
    time_old = input_nc.variables['time'][:]
    if treg != 'all':
        time_old = np.take(time_old, treg)

    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    df = vars_pd
    """
    ####################################################################################################################
    #   Regularise all files with 3 dimensions (time, nFaces, layers). 
    #   This will be equal to four dimensions in the regular grid format since nFaces is the x- and y- dimension.
    ####################################################################################################################
    """
    df2 = df.loc[df['ndims'] == 3]
    data_frommap_x = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x')
    data_frommap_y = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y')
    time = get_ncmodeldata(file_nc=file_nc, varname='time', timestep=treg)
    outname = '%s_regular.nc' % os.path.split(fileNC)[1][0:-3]
    file_nc_reg = os.path.join(dir_output, outname)
    root_grp = Dataset(file_nc_reg, 'w', format='NETCDF4')
    root_grp.description = 'Example simulation data'
    first_read = True
    i = 0
    for index, row in df2.iterrows():

        if row['dimensions'][1] == 'mesh2d_nEdges':
            continue
        data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=row['nc_varkeys'], timestep=treg, layer=lreg)
        data_frommap_var = data_frommap_var.filled(np.nan)
        field_array = np.empty((data_frommap_var.shape[0], ny, nx, data_frommap_var.shape[-1]))
        tms = data_frommap_var.shape[0]
        lrs = data_frommap_var.shape[-1]
        trange = range(0, tms)
        lrange = range(0, lrs)

        A = np.array([scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny,
                                             values=data_frommap_var[t, :, l].flatten(), method='linear') for t in
                      trange for l in lrange])
        x_grid = A[0][0]
        y_grid = A[0][1]
        A = A[:, 2, :, :]
        A = np.moveaxis(A, [0], [2])
        subs = np.split(A, tms, axis=2)

        field_array[:, :, :, 0:lrs] = [subs[tn] for tn in trange]
        field_array = np.ma.masked_invalid(field_array)
        print('done with variable %s' % row['nc_varkeys'])

        if first_read:
            unout = 'seconds since 2015-01-01 00:00:00'
            lon = x_grid[0, :]
            lat = y_grid[:, 0]
            # create dimensions
            root_grp.createDimension('time', None)
            root_grp.createDimension('lon', lon.shape[0])
            root_grp.createDimension('lat', lat.shape[0])
            root_grp.createDimension('layer', lrs)
            lonvar = root_grp.createVariable('lon', 'float32', 'lon')
            lonvar[:] = lon
            latvar = root_grp.createVariable('lat', 'float32', 'lat')
            latvar[:] = lat

            layervar = root_grp.createVariable('layer', 'float32', 'layer')
            layervar[:] = range(0, lrs)

            timevar = root_grp.createVariable('time', 'float64', 'time')
            timevar.setncattr('units', unout)
            timevar[:] = time_old

        fieldName = row['nc_varkeys']
        fieldvar = root_grp.createVariable(fieldName, 'float32', ('time', 'lat', 'lon', 'layer'), fill_value=-999)
        key = fieldName
        for ncattr in input_nc.variables[key].ncattrs():
            if ncattr != "_FillValue":
                root_grp.variables[fieldName].setncattr(ncattr, input_nc.variables[key].getncattr(ncattr))

        fieldvar[:] = field_array
        first_read = False
        i += 1
        if i>2:
            break

    """
    ####################################################################################################################
    #   Regularise all files with 2 dimensions (time, nFaces, layers).
    #   This will be equal to 3 dimensions in the regular grid format since nFaces is the x- and y- dimension.
    ####################################################################################################################
    """
    print('STARTING 2D')
    df2 = df.loc[df['ndims'] == 2]
    excludeList = ['edge', 'face', 'x', 'y']
    for index, row in df2.iterrows():
        test = any(n in str(row['nc_varkeys']) for n in excludeList)
        if not test:
            if row['dimensions'][1] == 'mesh2d_nEdges':
                continue
            ntimes = row['shape'][0]
            data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=row['nc_varkeys'], timestep=treg)
            data_frommap_var = data_frommap_var.filled(np.nan)
            field_array = np.empty((data_frommap_var.shape[0], ny, nx))
            trange = range(0, data_frommap_var.shape[0])
            tms = data_frommap_var.shape[0]
            A = np.array([scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny,
                                                 values=data_frommap_var[t, :].flatten(), method='linear') for t in
                          trange])

            A = A[:, 2, :, :]
            field_array[:, :, :] = A
            field_array = np.ma.masked_invalid(field_array)

            """write data to new netcdf"""
            fieldName = row['nc_varkeys']
            fieldvar = root_grp.createVariable(fieldName, 'float32', ('time', 'lat', 'lon'), fill_value=-999)
            key = fieldName
            for ncattr in input_nc.variables[key].ncattrs():
                if ncattr != "_FillValue":
                    root_grp.variables[fieldName].setncattr(ncattr, input_nc.variables[key].getncattr(ncattr))
            fieldvar[:] = field_array
    root_grp.close()


if __name__ == '__main__':
    time_start = tm.time()
    # dir_input = os.path.abspath(os.path.join('P:\\', '11202428-hisea', '03-Model', 'Greece_model','waq_model', 'run0_20200603', 'DFM_OUTPUT_tttz_waq'))
    dir_input = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'input', 'DFM_OUTPUT_tttz_waq'))
    fileNC = os.path.join(dir_input, 'tttz_waq_0000_map.nc')
    xpts = 100
    ypts = 80
    tms = np.arange(0,4,2)
    lrs = np.arange(0,3)
    # tms = 'all'
    # lrs = 'all'
    regularGrid_to_netcdf(fileNC, xpts, ypts, tms, lrs)
    time_elapsed = tm.time() - time_start
    print('Duration: %f s' %time_elapsed)

