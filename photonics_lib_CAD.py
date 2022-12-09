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

if __name__ == "__main__":
    lib = gdspy.GdsLibrary(unit=1e-09)

    #SWG_GC params
    #copied from old CAD code
    GCparams = {}
    GCparams['a_2DPhC'] = 232.61  # looks like this is the same as the
    GCparams['num_grating_periods_x'] = 20
    GCparams['num_grating_periods_y'] = 31  # SC total grating y width 12.6 um, maximum overlap with fiber
    GCparams['grating_period_x_start'] = 700
    GCparams['grating_pad_width'] = 12600  # assuming this is the same as the taperendwidth parameter in the lumerical geometry file
    GCparams['grating_first_index'] = 2.55
    GCparams['grating_delta_index'] = 0.02
    GCparams['phaseFactor'] = 0.582
    GCparams['a_2DPhC'] = 232.61  # looks like this is the same as the
    GCparams['grating_pad_length'] = 17000
    GCparams['grating_pad_buffer'] = (GCparams['grating_pad_width'] - (GCparams['num_grating_periods_y'] - 1) * GCparams['a_2DPhC'] * numpy.sqrt(3)) / 2.0
    GCparams['grating_pad_offset'] = 2000
    GCparams['grating_pad_spacing'] = 3000
    GCparams['resonance']=1544
    GCparams['n_circle_points_GC'] = 110
    GCparams['PhC_wy'] = 600
    GCparams['taper_length']=185000
    grating_x_num = range(GCparams['num_grating_periods_x'])  # range starts 0,1,2...num_grating_period_x-1
    air_hole_diameter_list_base = [-0.09594 * GCparams['a_2DPhC'] * (GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) ** 4 + 0.9546 * GCparams['a_2DPhC'] * (GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) ** 3 - 3.586 * GCparams['a_2DPhC'] * (
            GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) ** 2 + 5.542 * GCparams['a_2DPhC'] * (
                                           GCparams['grating_first_index'] - GCparams['grating_delta_index'] * x) - 1.931 * GCparams['a_2DPhC'] for x in grating_x_num]

    #relevant params for sweep

    xspacing = 135000
    yspacing = 250000

    GC_scales = [0.97,1.01,1.05]
    grating_coupler_list = []
    # silicon_skeleton_list = []
    grating_pad_center = [0,0]
    GC_hole_layer = 0
    GC_tether_x = 1000 #could try sweeping
    num_tether_along_taper = 3 #could try sweeping

    for GC_scale in GC_scales:
        grating_coupler = subwavelength_grating(air_hole_diameter_list_base,GC_scale,GCparams,grating_pad_center,GC_hole_layer, num_tether_along_taper, GC_tether_x)
        grating_coupler_list.append(grating_coupler)
        # silicon_skeleton_list.append(silicon_skeleton)
    # grid_of_holes = pg.grid(grating_coupler_list, spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    # grid_of_silicon = pg.grid(silicon_skeleton_list, spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    grid_of_GCs = pg.grid(grating_coupler_list,spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')

    #create PhC part of cavity
    PhCparams = {}
    PhCparams['num_PhC_per_wg'] = 2
    PhCparams['num_PhC_per_GC'] = 2*PhCparams['num_PhC_per_wg']
    PhCparams['PhC_hx'] = 236.8
    PhCparams['PhC_hy'] = 407.03
    PhCparams['aper_cav'] = 349.1  # 298
    PhCparams['aper_mir'] = 435.2  # 343
    num_cavity_holes = 12
    num_mirror_holes_middle = 2
    num_mirror_holes_end = 11
    PhCparams['bus_wg_width']=392
    PhCparams['PhC_wy']=600
    PhCparams['bus_wg_to_phc_wg_spacing']=500
    PhCparams['outer_cutout_phc_wg']=1500
    PhCparams['num_bus_reflector_mirrors'] = 6
    PhCparams['bus_reflector_taper_len']=2000
    PhCparams['beam_tether_x'] = 500
    PhCparams['bus_taper_xlen'] = 500
    PhCparams['bus_taper_chonky_y'] = 1000
    PhCparams['bus_reflect_taper_len_x'] = PhCparams['bus_reflector_taper_len']

    phc_scale_min = 0.87815
    phc_scale_max = 1.0036

    [ellipse_holes_top_beam,ellipse_holes_bottom_beam,aper_list] = generate_PhC_holes(PhCparams,phc_scale_min,phc_scale_max, num_cavity_holes,num_mirror_holes_middle,num_mirror_holes_end)
    print(aper_list)

    #beam hole set will give length of top and bottom beams (adjust later for mirror at the end)
    PhC_beam_len = ellipse_holes_top_beam.xmax-ellipse_holes_top_beam.xmin
    PhC_beam_skeleton = generate_PhC_skeleton(PhCparams,PhC_beam_len)
    PhC_holes_and_skeleton = Device()
    PhC_holes_and_skeleton.add_ref(PhC_beam_skeleton)
    PhC_holes_and_skeleton.add_ref(ellipse_holes_top_beam)
    PhC_holes_and_skeleton.add_ref(ellipse_holes_bottom_beam).mirror(p1=(PhC_beam_skeleton.xmin,PhC_beam_skeleton.y),p2=(PhC_beam_skeleton.xmax,PhC_beam_skeleton.y))

    #combine a single GC and PhC
    single_device = Device()
    GC_scale_single_check = 0.97
    grating_coupler_single = subwavelength_grating(air_hole_diameter_list_base, GC_scale_single_check, GCparams, grating_pad_center,
                                            GC_hole_layer, num_tether_along_taper, GC_tether_x)

    single_device.add_ref(grating_coupler_single)
    single_device_phc_half = single_device << PhC_holes_and_skeleton
    single_device_phc_half.xmin=grating_coupler_single.xmax
    single_device_phc_half.y=grating_coupler_single.y

    #get holes and GC together

    PhC_holes_and_skeleton.write_gds('phcbeamholes.gds',unit=1e-9,precision=1e-12)
    single_device.write_gds('checking combined device.gds',unit=1e-9,precision=1e-12)
    # gdspy.LayoutViewer(lib)
