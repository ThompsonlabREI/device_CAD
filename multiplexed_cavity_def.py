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
    num_PhC_per_beam,
    param_sweep_scaling,
    baseline_hx,
    baseline_hy):
    aper_list = []
    ellipse_dims_x=[]
    ellipse_dims_y=[]

    middle_mir_aper_list = numpy.ones(num_mirror_holes_middle)*aper_mir_nm
    ends_mir_aper_list = numpy.ones(num_mirror_holes_end)*aper_mir_nm
    nmax = (num_cavity_holes - 1) / 2.0
    aper_list_idx = numpy.linspace(0 - nmax, nmax, num_cavity_holes)
    aper_a = (aper_mir_nm - aper_cav_nm) / (nmax * nmax - 0.25)
    aper_c = (aper_mir_nm - 4 * aper_cav_nm * nmax * nmax) / (1 - 4 * nmax * nmax)
    cav_aper_list = [aper_a * x * x + aper_c for x in aper_list_idx]

    #there are PhCparams['num_PhC_per_wg'] cavities to define

    #add left mirror aper list - decide this will be smallest scaling
    aper_list.extend(ends_mir_aper_list)
    ellipse_dims_to_add = numpy.ones(len(ends_mir_aper_list))*param_sweep_scaling[0]
    ellipse_dims_x.extend(ellipse_dims_to_add*baseline_hx)
    ellipse_dims_y.extend(ellipse_dims_to_add*baseline_hy)

    #add cavity set of holes
    aper_list.extend(cav_aper_list)

    #this will be the first set of cavity holes and thus the smallest param scaling
    ellipse_dims_to_add = numpy.ones(len(cav_aper_list)) * param_sweep_scaling[0]
    ellipse_dims_x.extend(ellipse_dims_to_add * baseline_hx)
    ellipse_dims_y.extend(ellipse_dims_to_add * baseline_hy)

    #for however many more than 1 number of cavities per beam,
    for num_more_1_cav in range(1,num_PhC_per_beam):

        # add middle mirror set
        aper_list.extend(middle_mir_aper_list)
        ellipse_dims_to_add = numpy.ones(len(middle_mir_aper_list)) * param_sweep_scaling[num_more_1_cav]
        ellipse_dims_x.extend(ellipse_dims_to_add * baseline_hx)
        ellipse_dims_y.extend(ellipse_dims_to_add * baseline_hy)

        #add cavity set (last in for loop will be rightmost cavity on this beam)
        aper_list.extend(cav_aper_list)
        ellipse_dims_to_add = numpy.ones(len(cav_aper_list)) * param_sweep_scaling[num_more_1_cav]
        ellipse_dims_x.extend(ellipse_dims_to_add * baseline_hx)
        ellipse_dims_y.extend(ellipse_dims_to_add * baseline_hy)

    #add final set of end mirrors to aper list
    aper_list.extend(ends_mir_aper_list)
    ellipse_dims_to_add = numpy.ones(len(ends_mir_aper_list)) * param_sweep_scaling[-1]
    ellipse_dims_x.extend(ellipse_dims_to_add * baseline_hx)
    ellipse_dims_y.extend(ellipse_dims_to_add * baseline_hy)

    return [aper_list,ellipse_dims_x,ellipse_dims_y]

def generate_single_beam_set(aper_list,beam_ellipse_dims_x,beam_ellipse_dims_y,hole_center_x):
    ellipse_holes_single_beam = Device()
    for ellipse_index in range(len(aper_list)):
        # hole_center_x = aper_list[ellipse_index]/ 2.0  # SRP: aper_list is used to set center_x of ellipse
        PhC_ellipse = pg.ellipse(radii=(beam_ellipse_dims_x[ellipse_index]/2,beam_ellipse_dims_y[ellipse_index]/2),angle_resolution=2.5,layer=5)
        ellipse_holes_single_beam.add_ref(PhC_ellipse).movex(hole_center_x)
        # hole_center_y_new += aper # SRP: aper_list is used to set center_y of ellipse
        hole_center_x += aper_list[ellipse_index]
    return ellipse_holes_single_beam

def generate_PhC_holes(
    PhCparams,
    this_PhC_set_scale_min,
    this_PhC_set_scale_max,
    num_cavity_holes,
    num_mirror_holes_middle,
    num_mirror_holes_end
):
    scale_list_phc = numpy.linspace(this_PhC_set_scale_min, this_PhC_set_scale_max, PhCparams['num_PhC_per_GC'])  # param2['num_phc'] #SRP: this is where PhC scale list gets defined by sweeping over n total cavities - n/2 on the bottom waveguide and n/2 on the top waveguide
    bottom_and_top_sets = numpy.split(scale_list_phc,2)
    top_beam_scaling_phc = bottom_and_top_sets[0]
    bottom_beam_scaling_phc = bottom_and_top_sets[1]
    #get aper_list and ellipse scalings for the top beam
    [aper_list,top_beam_ellipse_dims_x,top_beam_ellipse_dims_y] = generate_aper_list(num_cavity_holes,num_mirror_holes_middle,num_mirror_holes_end,PhCparams['aper_mir'],PhCparams['aper_cav'],PhCparams['num_PhC_per_wg'],top_beam_scaling_phc,PhCparams['PhC_hx'],PhCparams['PhC_hy'])
    [aper_list_bot,bot_beam_ellipse_dims_x,bot_beam_ellipse_dims_y] = generate_aper_list(num_cavity_holes,num_mirror_holes_middle,num_mirror_holes_end,PhCparams['aper_mir'],PhCparams['aper_cav'],PhCparams['num_PhC_per_wg'],bottom_beam_scaling_phc,PhCparams['PhC_hx'],PhCparams['PhC_hy'])
    # top_beam_ellipse_dims = generate_ellipse_dims(num_cavity_holes,num_mirror_holes_middle,num_mirror_holes_end,top_beam_scaling_phc)

    hole_center_y_new = 0
    hole_center_x = 0
    ellipse_holes_top_beam = generate_single_beam_set(aper_list, top_beam_ellipse_dims_x, top_beam_ellipse_dims_y, hole_center_x)
    ellipse_holes_bottom_beam = generate_single_beam_set(aper_list_bot, bot_beam_ellipse_dims_x, bot_beam_ellipse_dims_y, hole_center_x)

    #get aper list and hole scalings for the bottom beam

    # for i in range(0, len(param['rad_list_mat'][:, 1])):
    #     hole_center_x = cell_edge_x + param['aper_list'][i] / 2.0  # SRP: aper_list is used to set center_x of ellipse
    #     hole_center_y_new += param['aper_list'][i]  # SRP: aper_list is used to set center_y of ellipse
    phc_beams = Device()
    phc_beams.add_ref(ellipse_holes_top_beam)
    phc_beams.add_ref(ellipse_holes_bottom_beam)

    phc_beams_group = Group([ellipse_holes_top_beam,ellipse_holes_bottom_beam])
    # phc_beams.write_gds('twobeams.gds', unit=1e-9, precision=1e-12)

    print('top_beam_ellipse_dims' + str(top_beam_ellipse_dims_x))
    print('bottom beam ellipse dims x' + str(bot_beam_ellipse_dims_x))
    return [phc_beams,aper_list]

def generate_PhC_skeleton(
    PhCparams,
    beam_len_x
):
    phc_skeleton_layer = 6
    phc_skeleton_silicon = Device()
    beam_path = pp.straight(length=beam_len_x)
    # P2 = pp.straight(length=test_PhC_bus_len)
    phc_beam_ref = beam_path.extrude(PhCparams['PhC_wy'],layer=phc_skeleton_layer)
    phc_top_beam = phc_skeleton_silicon << phc_beam_ref
    phc_bot_beam = phc_skeleton_silicon << phc_beam_ref
    phc_bus_ref = beam_path.extrude(PhCparams['bus_wg_width'],layer=phc_skeleton_layer)
    phc_bus_wg = phc_skeleton_silicon << phc_bus_ref

    #now align the elements with respect to each other
    phc_bus_wg.ymax = phc_top_beam.ymin-PhCparams['bus_wg_to_phc_wg_spacing']
    phc_bot_beam.ymax = phc_bus_wg.ymin - PhCparams['bus_wg_to_phc_wg_spacing']
    # pad_ref = P1.extrude(GCparams['grating_pad_width'], layer=2)
    # PhC_WG_ref = P2.extrude(GCparams['PhC_wy'], layer=2)
    return phc_skeleton_silicon