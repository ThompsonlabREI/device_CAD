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
    for GC_scale in GC_scales:
        grating_coupler = subwavelength_grating(air_hole_diameter_list_base,GC_scale,GCparams,grating_pad_center,GC_hole_layer)
        grating_coupler_list.append(grating_coupler)
        # silicon_skeleton_list.append(silicon_skeleton)
    # grid_of_holes = pg.grid(grating_coupler_list, spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    # grid_of_silicon = pg.grid(silicon_skeleton_list, spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    grid_of_GCs = pg.grid(grating_coupler_list,spacing=(xspacing,yspacing),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')

    # grid_of_holes_2 = pg.grid(grating_coupler_list, spacing=(5,1),separation=True,shape=(1,3),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    # big_boy = Device()
    # grid1 = big_boy << grid_of_holes
    # grid2 = big_boy << grid_of_holes
    # silicon_structure1 = big_boy << grid_of_silicon
    # silicon_structure2 = big_boy << grid_of_silicon
    # bounding_box_grid1 = grid1.get_bounding_box()
    # silicon_structure_box_grid1 = silicon_structure1.get_bounding_box()
    # grid2.move(origin=bounding_box_grid1[0],destination=bounding_box_grid1[1])
    # silicon_structure2.move(origin=silicon_structure_box_grid1[0],destination=silicon_structure_box_grid1[1])

    # big_boy.remap_layers({0:1})
    # D = pg.gridsweep(
    #     function=pg.circle,
    #     param_x={'radius': [10, 20, 30, 40, 50]},
    #     param_y={'layer': [0, 5, 9]},
    #     param_defaults={},
    #     param_override={},
    #     spacing=(30, 10),
    #     separation=True,
    #     align_x='x',
    #     align_y='y',
    #     edge_x='x',
    #     edge_y='ymax',
    #     label_layer=None)
    #
    # # qp(D)
    # D.write_gds('phidlcheck.gds')
    # grid_of_holes.write_gds('gridsweepcheck.gds',unit=1e-9,precision=1e-12)
    grid_of_GCs.write_gds('checking_deviceofdevice.gds',unit=1e-9,precision=1e-12)
    # gdspy.LayoutViewer(lib)
