import numpy
import gdspy
from phidl import Device, Path, CrossSection
from phidl import Group
import phidl.geometry as pg
import phidl.path as pp

def generate_aper_list(num_cavity_holes,
    num_mirror_holes_middle,
    num_mirror_holes_end,
    aper_mir_nm,
    aper_cav_nm,
    num_PhC_per_beam):
    aper_list = []
    print(num_mirror_holes_middle)

    middle_mir_aper_list = numpy.ones(num_mirror_holes_middle)*aper_mir_nm
    ends_mir_aper_list = numpy.ones(num_mirror_holes_end)*aper_mir_nm
    nmax = (num_cavity_holes - 1) / 2.0
    aper_list_idx = numpy.linspace(0 - nmax, nmax, num_cavity_holes)
    aper_a = (aper_mir_nm - aper_cav_nm) / (nmax * nmax - 0.25)
    aper_c = (aper_mir_nm - 4 * aper_cav_nm * nmax * nmax) / (1 - 4 * nmax * nmax)
    cav_aper_list = [aper_a * x * x + aper_c for x in aper_list_idx]

    #there are PhCparams['num_PhC_per_wg'] cavities to define

    #add left mirror aper list
    aper_list.extend(ends_mir_aper_list)

    #add cavity set of holes
    aper_list.extend(cav_aper_list)

    #for however many more than 1 number of cavities per beam,
    for num_more_1_cav in range(0,num_PhC_per_beam-1):

        # add middle mirror set
        aper_list.extend(middle_mir_aper_list)
        #add cavity set (last in for loop will be rightmost cavity on this beam)
        aper_list.extend(cav_aper_list)

    #add final set of end mirrors to aper list
    aper_list.extend(ends_mir_aper_list)

    return aper_list

def generate_PhC_holes(
    PhCparams,
    # this_PhC_set_scale_min,
    # this_PhC_set_scale_max,
    num_cavity_holes,
    num_mirror_holes_middle,
    num_mirror_holes_end
):
    # scale_list_phc = numpy.linspace(this_PhC_set_scale_min, this_PhC_set_scale_max, PhCparams['num_PhC_per_GC'])  # param2['num_phc'] #SRP: this is where PhC scale list gets defined by sweeping over n total cavities - n/2 on the bottom waveguide and n/2 on the top waveguide

    #get aper_list
    aper_list = generate_aper_list(num_cavity_holes,num_mirror_holes_middle,num_mirror_holes_end,PhCparams['aper_mir'],PhCparams['aper_cav'],PhCparams['num_PhC_per_wg'])

    # for i in range(0, len(param['rad_list_mat'][:, 1])):
    #     hole_center_x = cell_edge_x + param['aper_list'][i] / 2.0  # SRP: aper_list is used to set center_x of ellipse
    #     hole_center_y_new += param['aper_list'][i]  # SRP: aper_list is used to set center_y of ellipse
    single_beam_of_holes = Device()
    return [single_beam_of_holes,aper_list]