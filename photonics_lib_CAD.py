######################################################################
#                                                                    #
#  Copyright 2009 Lucas Heitzmann Gabrielli.                         #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################


import numpy
import gdspy
from grating_def import *
from GC_def import *
import phidl.geometry as pg
from phidl import quickplot as qp
from multiplexed_cavity_def import *
from define_params import *
import pickle

def generate_single_device(GC_scale,phc_scale_min,phc_scale_max,num_tether_along_taper,num_cavities_per_wg,gc_taper_len_x):
    PhCparams['num_PhC_per_wg'] = num_cavities_per_wg
    PhCparams['num_PhC_per_GC'] = 2 * PhCparams['num_PhC_per_wg']
    [ellipse_holes_top_beam_ref,ellipse_holes_bottom_beam_ref,aper_list,reflector_hole_set_ref] = generate_PhC_holes(PhCparams,phc_scale_min,phc_scale_max)
    cavity_check_savename = 'matching_cavities_OR_code' + '.pickle'
    with open(cavity_check_savename, 'wb') as handle:
        pickle.dump(aper_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(cavity_check_savename, 'rb') as handle:
        aper_list_loaded = pickle.load(handle)
        print("aper list loaded" + str(aper_list_loaded))
        print("aper list type " + str(type(aper_list_loaded)))
    # print(aper_list)
    # reflector_hole_set.write_gds('reflector_hole_set.gds',unit=1e-9,precision=1e-12)
    # reflector_hole_set = generate_single_beam_set(aper_list,beam_ellipse_dims_x,beam_ellipse_dims_y,hole_center_x)
    #beam hole set will give length of top and bottom beams (adjust later for mirror at the end)

    reflector_len_x = reflector_hole_set_ref.xmax-reflector_hole_set_ref.xmin
    PhC_beam_len = ellipse_holes_top_beam_ref.xmax-ellipse_holes_top_beam_ref.xmin + PhCparams['phc_beam_buffer_x']
    # PhC_beam_len += PhCparams['phc_beam_buffer_x']
    [PhC_beam_skeleton_ref,reflector_skeleton_ref] = generate_PhC_skeleton(PhCparams,PhC_beam_len,reflector_len_x)
    PhC_holes_and_skeleton = Device()

    PhC_beam_skeleton = PhC_holes_and_skeleton << PhC_beam_skeleton_ref
    # PhC_holes_and_skeleton.add_ref(PhC_beam_skeleton)
    ellipse_holes_top_beam = PhC_holes_and_skeleton << ellipse_holes_top_beam_ref
    ellipse_holes_top_beam.x = PhC_beam_skeleton.x
    # PhC_holes_and_skeleton.add_ref(ellipse_holes_top_beam)

    reflector_skeleton = PhC_holes_and_skeleton << reflector_skeleton_ref
    reflector_skeleton.xmin = PhC_beam_skeleton.xmax
    reflector_skeleton.y=PhC_beam_skeleton.y
    reflector_hole_set = PhC_holes_and_skeleton << reflector_hole_set_ref
    reflector_hole_set.y=reflector_skeleton.y
    reflector_hole_set.xmax=reflector_skeleton.xmax

    ellipse_holes_bottom_beam = PhC_holes_and_skeleton << ellipse_holes_bottom_beam_ref
    ellipse_holes_bottom_beam.mirror(p1=(PhC_beam_skeleton.xmin,PhC_beam_skeleton.y),p2=(PhC_beam_skeleton.xmax,PhC_beam_skeleton.y))
    # PhC_holes_and_skeleton.add_ref(ellipse_holes_bottom_beam_ref).mirror(p1=(PhC_beam_skeleton.xmin,PhC_beam_skeleton.y),p2=(PhC_beam_skeleton.xmax,PhC_beam_skeleton.y))
    ellipse_holes_bottom_beam.x = PhC_beam_skeleton.x

    #combine a single GC and PhC
    single_device = Device()
    GCparams['taper_length'] = gc_taper_len_x
    grating_coupler_single = subwavelength_grating(air_hole_diameter_list_base, GC_scale, GCparams, grating_pad_center,
                                            GC_hole_layer, num_tether_along_taper, GC_tether_x)

    single_device.add_ref(grating_coupler_single)
    single_device_phc_half = single_device << PhC_holes_and_skeleton
    single_device_phc_half.xmin=grating_coupler_single.xmax
    single_device_phc_half.y=grating_coupler_single.y
    return single_device

if __name__ == "__main__":
    lib = gdspy.GdsLibrary(unit=1e-09)

    #get params
    [PhCparams, GCparams] = generate_photonics_params()

    grating_x_num = range(GCparams['num_grating_periods_x'])  # range starts 0,1,2...num_grating_period_x-1
    air_hole_diameter_list_base = [-0.09594 * GCparams['a_2DPhC'] * (GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) ** 4 + 0.9546 * GCparams['a_2DPhC'] * (GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) ** 3 - 3.586 * GCparams['a_2DPhC'] * (
            GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) ** 2 + 5.542 * GCparams['a_2DPhC'] * (
                                           GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) - 1.931 * GCparams['a_2DPhC'] for x in grating_x_num]

    #relevant params for sweep
    num_GC_scalings = 2
    GC_scale_min = 0.95
    GC_scale_max = 1.05
    GC_scales = numpy.linspace(GC_scale_min,GC_scale_max,num=num_GC_scalings)
    # print("gc scales" + str(GC_scales))

    xspacing = 50000
    yspacing = 50000

    # GC_scales = [0.97,1.01,1.05]
    # grating_coupler_list = []
    # silicon_skeleton_list = []
    grating_pad_center = [0,0]
    GC_hole_layer = 0
    GC_tether_x = 1000 #could try sweeping
    # num_tether_along_taper = 3 #could try sweeping

    # for GC_scale in GC_scales:
    #     grating_coupler = subwavelength_grating(air_hole_diameter_list_base,GC_scale,GCparams,grating_pad_center,GC_hole_layer, num_tether_along_taper, GC_tether_x)
    #     grating_coupler_list.append(grating_coupler)
    #     # silicon_skeleton_list.append(silicon_skeleton)
    # # grid_of_holes = pg.grid(grating_coupler_list, spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    # # grid_of_silicon = pg.grid(silicon_skeleton_list, spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    # grid_of_GCs = pg.grid(grating_coupler_list,spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')

    #create PhC part of cavity
    # num_cavity_holes = 12
    # num_mirror_holes_middle = 2
    # num_mirror_holes_end = 11
    # phc_scale_min = 0.87815
    # phc_scale_max = 1.0036

    # num_phc_in_sweep = [2,4, 6, 8, 10]
    num_phc_in_sweep=[4]
    # num_PhC_sweep = 4
    # phc_scale_min = 0.8
    # phc_scale_max = 1.1
    phc_scale_min = 1.0
    phc_scale_max = 1.0
    #overlap
    gc_phc_param_sweep_devices = []
    # num_tethers_sweep = [1,2,4,10]
    # min_num_tethers = 1
    # max_num_tethers = 10
    min_num_tethers = 3
    max_num_tethers = 5
    num_tethers_sweep = numpy.arange(min_num_tethers,max_num_tethers,2)
    print('num tethers sweep' + str(num_tethers_sweep))

    # GC_taper_min_x = 20000
    GC_taper_min_x = 185000
    GC_taper_max_x = 185000
    num_taper_len_sweep = 1
    gc_taper_lens = numpy.round(numpy.linspace(GC_taper_min_x,GC_taper_max_x,num=num_taper_len_sweep),0)
    print("gc taper lengths" + str(gc_taper_lens))
    # phc_freq_scalings = numpy.linspace(phc_scale_min,phc_scale_max,num=num_PhC_sweep+1)

    # print(GC_scales)
    for num_phc_sweep in num_phc_in_sweep:
        print('num phc on this wg' + str(num_phc_sweep))
        # phc_freq_scalings = numpy.linspace(phc_scale_min, phc_scale_max, num=num_phc_in_sweep[num_phc_sweep_index] + 1)
        # for phc_freq_scale_index in range(len(phc_freq_scalings)-1):
        for num_tethers in num_tethers_sweep:
            print('num tethers' + str(num_tethers))
            # GCparams['']
            for gc_taper_len in gc_taper_lens:
                print('gc taper len' + str(gc_taper_len))
                single_device = generate_single_device(GC_scales[1],phc_scale_min,phc_scale_max,num_tethers,num_phc_sweep,gc_taper_len)
                gc_phc_param_sweep_devices.append(single_device)

    # for num_phc_sweep_index in range(len(num_phc_in_sweep)):
    #     phc_freq_scalings = numpy.linspace(phc_scale_min, phc_scale_max, num=num_phc_in_sweep[num_phc_sweep_index] + 1)
    # for GC_scale_index in range(len(GC_scales)):
    #     for phc_freq_scale_index in range(len(phc_freq_scalings)-1):
    #         for num_tethers_index in range(min_num_tethers,max_num_tethers):
    #             single_device = generate_single_device(GC_scales[GC_scale_index],phc_freq_scalings[phc_freq_scale_index],phc_freq_scalings[phc_freq_scale_index+1])
    #             gc_phc_param_sweep_devices.append(single_device)

    # single_device_check = generate_single_device(GC_scales[0],phc_scale_min,phc_scale_max)
    grid_of_devices = pg.grid(gc_phc_param_sweep_devices,spacing=(xspacing,yspacing),separation=True,shape=(len(num_phc_in_sweep),len(num_tethers_sweep)*len(gc_taper_lens)),align_x='x',align_y='y',edge_x='x',edge_y='ymax')

    #get holes and GC together
    # pg.gridsweep()

    grid_of_devices.write_gds('checking cell spacing cavity section.gds',unit=1e-9,precision=1e-12)
    # single_device_check.write_gds('checking combined device.gds',unit=1e-9,precision=1e-12)

    # gdspy.LayoutViewer(lib)
