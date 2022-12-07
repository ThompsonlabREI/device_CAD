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
    # Examples
    lib = gdspy.GdsLibrary(unit=1e-09)

    # Negative resist example
    width = 0.45
    bend_radius = 50.0
    ring_radius = 20.0
    taper_len = 50.0
    input_gap = 150.0
    io_gap = 500.0
    wg_gap = 20.0
    ring_gaps = [0.06 + 0.02 * i for i in range(8)]
    GC_hole_radius = 0.153

    ring = lib.new_cell("NRing")
    ring.add(
        gdspy.Round((ring_radius, 0), ring_radius, ring_radius - width, tolerance=0.001)
    )

    grat = lib.new_cell("NGrat")
    grat.add(
        grating(
            0.626,
            28,
            0.5,
            19,
            (0, 0),
            "+y",
            1.55,
            numpy.sin(numpy.pi * 8 / 180),
            21.5,
            width,
            tolerance=0.001,
        )
    )

    #SWG_GC params
    #copied from old CAD code
    param = {}
    param['a_2DPhC'] = 232.61  # looks like this is the same as the
    param['num_grating_periods_x'] = 20
    param['num_grating_periods_y'] = 31  # SC total grating y width 12.6 um, maximum overlap with fiber
    param['grating_period_x_start'] = 700
    param['grating_pad_width'] = 12600  # assuming this is the same as the taperendwidth parameter in the lumerical geometry file
    param['grating_first_index'] = 2.55
    param['grating_delta_index'] = 0.02
    param['phaseFactor'] = 0.582
    param['a_2DPhC'] = 232.61  # looks like this is the same as the
    param['grating_pad_length'] = 17000
    param['grating_pad_buffer'] = (param['grating_pad_width'] - (param['num_grating_periods_y'] - 1) * param['a_2DPhC'] * numpy.sqrt(3)) / 2.0
    param['grating_pad_offset'] = 2000
    param['grating_pad_spacing'] = 3000
    param['resonance']=1544
    param['n_circle_points_GC'] = 110
    grating_x_num = range(param['num_grating_periods_x'])  # range starts 0,1,2...num_grating_period_x-1
    air_hole_diameter_list_base = [-0.09594 * param['a_2DPhC'] * (param['grating_first_index'] - param['grating_delta_index'] * x) ** 4 + 0.9546 * param['a_2DPhC'] * (param['grating_first_index'] - param['grating_delta_index'] * x) ** 3 - 3.586 * param['a_2DPhC'] * (
                                               param['grating_first_index'] - param['grating_delta_index'] * x) ** 2 + 5.542 * param['a_2DPhC'] * (
                                               param['grating_first_index'] - param['grating_delta_index'] * x) - 1.931 * param['a_2DPhC'] for x in grating_x_num]

    constgrating_airholescale = 1.01
    SW_GC_cell = lib.new_cell("N_SW_GC")
    [SW_GC,disregard] = subwavelength_grating(air_hole_diameter_list_base,
        constgrating_airholescale,
        param['grating_period_x_start'],
        param['phaseFactor'],
        param['grating_delta_index'],
        param["resonance"],
        param['num_grating_periods_x'],
        param['num_grating_periods_y'],
        param['n_circle_points_GC'],
        param['a_2DPhC'],
        param['grating_pad_length'],
        param['grating_pad_width'],
        param['grating_pad_buffer'],
        0,
        0,
        SW_GC_cell)


    # def subwavelength_grating(
    #         air_hole_diameter_list_base,
    #         constgrating_airhole_scale_factor,
    #         grating_period_x_start,
    #         phase_factor,
    #         grating_delta_index,
    #         resonance_freq,
    #         num_grating_periods_x,
    #         num_grating_periods_y,
    #         num_points,
    #         a_2DPhC_nm,
    #         grating_pad_length,
    #         grating_pad_width,
    #         grating_pad_buffer,
    #         grating_pad_center_x,
    #         grating_pad_center_y,
    #         GC_holes_cell,
    # ):

    taper = lib.new_cell("NTaper")
    taper.add(gdspy.Path(0.12, (0, 0)).segment(taper_len, "+y", final_width=width))

    circle_array = lib.new_cell("N_GC_holes")

    circle_array.add(gdspy.Round((GC_hole_radius, 0), GC_hole_radius, tolerance=0.001))

    c = lib.new_cell("Negative")
    for i, gap in enumerate(ring_gaps):
        path = gdspy.FlexPath(
            [(input_gap * i, taper_len)],
            width=width,
            corners="circular bend",
            bend_radius=bend_radius,
            gdsii_path=True,
        )
        path.segment((0, 600 - wg_gap * i), relative=True)
        path.segment((io_gap, 0), relative=True)
        path.segment((0, 300 + wg_gap * i), relative=True)
        c.add(path)
        c.add(gdspy.CellReference(ring, (input_gap * i + width / 2 + gap, 300)))
    c.add(gdspy.CellArray(taper, len(ring_gaps), 1, (input_gap, 0), (0, 0)))
    c.add(
        gdspy.CellArray(
            grat, len(ring_gaps), 1, (input_gap, 0), (io_gap, 900 + taper_len)
        )
    )

    c.add(gdspy.CellArray(circle_array, len(ring_gaps), 1, (input_gap, 0), (0, 0)))
    c.add(SW_GC)

    # Save to a gds file and check out the output
    lib.write_gds("photonics_1layersaved.gds",cells=[c])

    GC_scales = [0.97,1.01,1.05]
    grating_coupler_list = []
    for GC_scale in GC_scales:
        [GC_holes,device_format_holes] = subwavelength_grating(air_hole_diameter_list_base,
        GC_scale,
        param['grating_period_x_start'],
        param['phaseFactor'],
        param['grating_delta_index'],
        param["resonance"],
        param['num_grating_periods_x'],
        param['num_grating_periods_y'],
        param['n_circle_points_GC'],
        param['a_2DPhC'],
        param['grating_pad_length'],
        param['grating_pad_width'],
        param['grating_pad_buffer'],
        0,
        0,
        SW_GC_cell)
        grating_coupler_list.append(device_format_holes)
    grid_of_holes = pg.grid(grating_coupler_list, spacing=(5,1),separation=True,shape=(3,1),align_x='x',align_y='y',edge_x='x',edge_y='ymax')
    D = pg.gridsweep(
        function=pg.circle,
        param_x={'radius': [10, 20, 30, 40, 50]},
        param_y={'layer': [0, 5, 9]},
        param_defaults={},
        param_override={},
        spacing=(30, 10),
        separation=True,
        align_x='x',
        align_y='y',
        edge_x='x',
        edge_y='ymax',
        label_layer=None)

    # qp(D)
    D.write_gds('phidlcheck.gds')
    grid_of_holes.write_gds('gridsweepcheck.gds')
    # gdspy.LayoutViewer(lib)
