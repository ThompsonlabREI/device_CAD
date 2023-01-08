import numpy
import gdspy
from phidl import Device, Path, CrossSection
from phidl import Group
import phidl.geometry as pg
import phidl.path as pp
import pickle

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
    # print('param sweep scaling' + str(param_sweep_scaling))


    #for however many more than 1 number of cavities per beam,
    for num_more_1_cav in range(1,num_PhC_per_beam):

        # add middle mirror set
        aper_list.extend(middle_mir_aper_list)

        #split the number of middle mirror holes in 2 - left mirror holes match param scaling from left cavity and right mirror holes match param scaling from right cavity
        #left mirror holes get left scaling

        ellipse_dims_to_add_left = numpy.ones(len(middle_mir_aper_list[:len(middle_mir_aper_list)//2])) * param_sweep_scaling[num_more_1_cav-1]
        ellipse_dims_x.extend(ellipse_dims_to_add_left * baseline_hx)
        ellipse_dims_y.extend(ellipse_dims_to_add_left * baseline_hy)

        ellipse_dims_to_add_right = numpy.ones(len(middle_mir_aper_list[len(middle_mir_aper_list)//2:])) * param_sweep_scaling[num_more_1_cav]
        ellipse_dims_x.extend(ellipse_dims_to_add_right * baseline_hx)
        ellipse_dims_y.extend(ellipse_dims_to_add_right * baseline_hy)

        #add cavity set (last in for loop will be rightmost cavity on this bxeam)
        aper_list.extend(cav_aper_list)
        ellipse_dims_to_add = numpy.ones(len(cav_aper_list)) * param_sweep_scaling[num_more_1_cav]
        ellipse_dims_x.extend(ellipse_dims_to_add * baseline_hx)
        ellipse_dims_y.extend(ellipse_dims_to_add * baseline_hy)

    #add final set of end mirrors to aper list
    aper_list.extend(ends_mir_aper_list)
    ellipse_dims_to_add = numpy.ones(len(ends_mir_aper_list)) * param_sweep_scaling[-1]
    ellipse_dims_x.extend(ellipse_dims_to_add * baseline_hx)
    ellipse_dims_y.extend(ellipse_dims_to_add * baseline_hy)
    # print("aper list" + str(aper_list))

    return [aper_list,ellipse_dims_x,ellipse_dims_y]

def generate_single_beam_set(aper_list,beam_ellipse_dims_x,beam_ellipse_dims_y,hole_center_x,num_ellipse_points):
    ellipse_holes_single_beam = Device()
    phi_res = round((360.0/num_ellipse_points),1)
    print("angle res for ellipses" + str(phi_res))
    for ellipse_index in range(len(aper_list)):
        # hole_center_x = aper_list[ellipse_index]/ 2.0  # SRP: aper_list is used to set center_x of ellipse
        PhC_ellipse = pg.ellipse(radii=(beam_ellipse_dims_x[ellipse_index]/2,beam_ellipse_dims_y[ellipse_index]/2),angle_resolution=phi_res,layer=5)
        ellipse_holes_single_beam.add_ref(PhC_ellipse).movex(hole_center_x)
        # hole_center_y_new += aper # SRP: aper_list is used to set center_y of ellipse
        hole_center_x += aper_list[ellipse_index]
    return ellipse_holes_single_beam

def generate_PhC_holes(
    PhCparams,
    this_PhC_set_scale_min,
    this_PhC_set_scale_max
):
    scale_list_phc = numpy.linspace(this_PhC_set_scale_min, this_PhC_set_scale_max, PhCparams['num_PhC_per_GC'])  # param2['num_phc'] #SRP: this is where PhC scale list gets defined by sweeping over n total cavities - n/2 on the bottom waveguide and n/2 on the top waveguide
    bottom_and_top_sets = numpy.split(scale_list_phc,2)
    top_beam_scaling_phc = bottom_and_top_sets[0]
    bottom_beam_scaling_phc = bottom_and_top_sets[1]
    #get aper_list and ellipse scalings for the top beam
    [aper_list,top_beam_ellipse_dims_x,top_beam_ellipse_dims_y] = generate_aper_list(PhCparams['num_cavity_holes'],PhCparams['num_mirror_holes_middle'],PhCparams['num_mirror_holes_end'] ,PhCparams['aper_mir'],PhCparams['aper_cav'],PhCparams['num_PhC_per_wg'],top_beam_scaling_phc,PhCparams['PhC_hx'],PhCparams['PhC_hy'])
    [aper_list_bot,bot_beam_ellipse_dims_x,bot_beam_ellipse_dims_y] = generate_aper_list(PhCparams['num_cavity_holes'],PhCparams['num_mirror_holes_middle'],PhCparams['num_mirror_holes_end'] ,PhCparams['aper_mir'],PhCparams['aper_cav'],PhCparams['num_PhC_per_wg'],bottom_beam_scaling_phc,PhCparams['PhC_hx'],PhCparams['PhC_hy'])
    # top_beam_ellipse_dims = generate_ellipse_dims(num_cavity_holes,num_mirror_holes_middle,num_mirror_holes_end,top_beam_scaling_phc)

    hole_center_y_new = 0
    hole_center_x = 0
    ellipse_dims_check_savename_top = 'ellipse_dims_top_x_min_scale_' + str(round(this_PhC_set_scale_min,4)) + '.pickle'
    ellipse_dims_check_savename_bottom = 'ellipse_dims_bottom_xmin_scale_' + str(round(this_PhC_set_scale_min,4)) + '.pickle'
    # print("top beam ellipse dims y" + str(top_beam_ellipse_dims_y))

    with open(ellipse_dims_check_savename_top, 'wb') as handle_top:
        pickle.dump(top_beam_ellipse_dims_x, handle_top, protocol=pickle.HIGHEST_PROTOCOL)

    with open(ellipse_dims_check_savename_bottom, 'wb') as handle_bot:
        pickle.dump(bot_beam_ellipse_dims_x, handle_bot, protocol=pickle.HIGHEST_PROTOCOL)
    # with open(ellipse_dims_check_savename, 'rb') as handle:
    #     top_beam_ellipse_dims_x_loaded = pickle.load(handle)
    #     print("top beam list loaded" + str(top_beam_ellipse_dims_x_loaded))
    #     print("top beam list type " + str(type(top_beam_ellipse_dims_x_loaded)))
    ellipse_holes_top_beam = generate_single_beam_set(aper_list, top_beam_ellipse_dims_x, top_beam_ellipse_dims_y, hole_center_x, PhCparams['num_ellipse_points'])
    ellipse_holes_bottom_beam = generate_single_beam_set(aper_list_bot, bot_beam_ellipse_dims_x, bot_beam_ellipse_dims_y, hole_center_x, PhCparams['num_ellipse_points'])

    #get aper list and hole scalings for the bottom beam

    # for i in range(0, len(param['rad_list_mat'][:, 1])):
    #     hole_center_x = cell_edge_x + param['aper_list'][i] / 2.0  # SRP: aper_list is used to set center_x of ellipse
    #     hole_center_y_new += param['aper_list'][i]  # SRP: aper_list is used to set center_y of ellipse
    phc_beams = Device()
    phc_beams.add_ref(ellipse_holes_top_beam)
    phc_beams.add_ref(ellipse_holes_bottom_beam)

    phc_beams_group = Group([ellipse_holes_top_beam,ellipse_holes_bottom_beam])
    # phc_beams.write_gds('twobeams.gds', unit=1e-9, precision=1e-12)
    aper_list_reflector = numpy.ones(PhCparams['num_bus_reflector_mirrors'])*PhCparams['aper_mir']


    beam_ellipse_x_avg = numpy.average([numpy.average(top_beam_ellipse_dims_x),numpy.average(bot_beam_ellipse_dims_x)])
    # print("beam ellipse avg" +str(beam_ellipse_x_avg))
    beam_ellipse_y_avg = numpy.average([numpy.average(top_beam_ellipse_dims_y),numpy.average(bot_beam_ellipse_dims_y)])
    # print("beam ellipse avg" +str(beam_ellipse_y_avg))
    reflector_ellipse_dims_x = numpy.ones(PhCparams['num_bus_reflector_mirrors'])*beam_ellipse_x_avg
    reflector_ellipse_dims_y = numpy.ones(PhCparams['num_bus_reflector_mirrors']) * beam_ellipse_y_avg
    reflector_set = generate_single_beam_set(aper_list_reflector,reflector_ellipse_dims_x,reflector_ellipse_dims_y,hole_center_x, PhCparams['num_ellipse_points'])

    # print('top_beam_ellipse_dims' + str(top_beam_ellipse_dims_x))
    # print('bottom beam ellipse dims x' + str(bot_beam_ellipse_dims_x))

    return [ellipse_holes_top_beam,ellipse_holes_bottom_beam,aper_list,reflector_set]

def generate_PhC_skeleton(
    PhCparams,
    beam_len_x,
    reflector_len_x
):
    phc_skeleton_layer = 6
    phc_skeleton_subtractions_layer = 8
    phc_skeleton_silicon = Device()

    #calculate additional length needed

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

    #define y tether
    ytether_path = pp.straight(length=PhCparams['beam_tether_x'])
    ytether_len = phc_skeleton_silicon.ysize+2*PhCparams['outer_cutout_phc_wg']
    ytether = ytether_path.extrude(ytether_len,layer=phc_skeleton_subtractions_layer)
    ytether.center = phc_skeleton_silicon.center

    #define x taper
    phc_skeleton_silicon_tethered_layer = 9

    bus_center_taper_path = pp.straight(length=PhCparams['bus_reflector_taper_len'])
    bus_left_taper = bus_center_taper_path.extrude([PhCparams['bus_wg_width'],PhCparams['bus_taper_chonky_y']],layer=phc_skeleton_silicon_tethered_layer)
    bus_left_taper.y = ytether.y
    bus_left_taper.xmax = ytether.xmin

    bus_right_taper = pg.copy(bus_left_taper)
    bus_right_taper.rotate(180)
    bus_right_taper.xmin = ytether.xmax
    bus_right_taper.y=ytether.y

    #subtract y tether and x taper from 3 beams
    phc_skeleton_silicon_tethered = pg.boolean(A=phc_skeleton_silicon,B=ytether,operation='or',precision=1e-9,num_divisions=[1,1],layer=phc_skeleton_silicon_tethered_layer)
    phc_skeleton_silicon_tethered.write_gds('phcbones_tethered.gds', unit=1e-9, precision=1e-12)
    # phc_skeleton_silicon_tethered.add_ref(bus_left_taper)
    # phc_skeleton_silicon_tethered.add_ref(bus_right_taper)
    phc_skeleton_silicon_tether_tapered_layer = 10
    phc_skeleton_silicon_tether_tapered = pg.boolean(A=phc_skeleton_silicon_tethered,B=[bus_left_taper,bus_right_taper],operation='or',precision=1e-9,num_divisions=[1,1],layer=phc_skeleton_silicon_tether_tapered_layer)

    phc_skeleton_silicon_tether_tapered.write_gds('phcbones.gds', unit=1e-9, precision=1e-12)

    #define outer layer for cutout
    outer_layer_rect_ref = pg.rectangle(size=(phc_skeleton_silicon_tether_tapered.xsize,ytether_len),layer=7)
    phc_outer_box = Device()
    out_layer_rect = phc_outer_box << outer_layer_rect_ref
    out_layer_rect.y = phc_skeleton_silicon_tether_tapered.y

    phc_cutout_only = pg.boolean(A=phc_outer_box,B=phc_skeleton_silicon_tether_tapered,operation='not',precision=1e-9,num_divisions=[1,1],layer=0)

    #generate reflector skeleton
    #taper part
    #reflector beam
    reflector_beam_path = pp.straight(length=reflector_len_x+PhCparams['aper_mir']) #add the aper on for buffer between end of ellipse cutouts and start of non-suspended region
    reflector_taper_path = pp.straight(length=PhCparams['bus_reflect_taper_len_x'])
    reflector_beam_ref = reflector_beam_path.extrude(PhCparams['PhC_wy'],layer=phc_skeleton_layer)
    reflector_taper_ref = reflector_taper_path.extrude([PhCparams['bus_wg_width'],PhCparams['PhC_wy']],layer=phc_skeleton_layer)

    reflector_skeleton = Device()
    phc_reflector_beam = reflector_skeleton << reflector_beam_ref
    phc_reflector_taper = reflector_skeleton << reflector_taper_ref
    phc_reflector_taper.xmax=phc_reflector_beam.xmin
    reflector_skeleton.write_gds('bus_reflector_skeleton.gds', unit=1e-9, precision=1e-12)
    reflector_outer_rect_ref = pg.rectangle(size=(reflector_skeleton.xsize,PhCparams['bus_wg_width']+2*PhCparams['bus_wg_to_phc_wg_spacing']),layer=7)
    reflector_outer_rect_ref.center = reflector_skeleton.center
    reflector_cutout = pg.boolean(A=reflector_outer_rect_ref,B=reflector_skeleton,operation='not',precision=1e-9,num_divisions=[1,1],layer=0)
    reflector_cutout.write_gds('bus_reflector_skeleton_cutout.gds', unit=1e-9, precision=1e-12)
    return [phc_cutout_only,reflector_cutout]