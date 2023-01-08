import numpy
import gdspy
from phidl import Device, Path, CrossSection
from phidl import Group
import phidl.geometry as pg
import phidl.path as pp

def generate_grating_holes(air_hole_diameter_list_base,
    constgrating_airhole_scale_factor,
    GCparams,
    grating_pad_origin,
    GC_hole_layer,
    num_grating_airhole_points):
    SWG_holes = Device()

    air_hole_diameter_list = numpy.asarray(
        air_hole_diameter_list_base) * constgrating_airhole_scale_factor  # SRP: varies from 0.85 to 1.05 for first run
    # print("scaled airhole diameter list for GC" + str(air_hole_diameter_list))
    # air_hole_diameter_list = numpy.asarray(air_hole_diameter_list_base) + constgrating_airholescale_list[iB] # take off the 1.13 scaling factor so that GC hole sizes match lumerical
    grating_start_x = grating_pad_origin[0] - GCparams['grating_pad_length'] / 2.0 + GCparams['grating_pad_buffer']
    grating_start_y = grating_pad_origin[1] - GCparams['grating_pad_width'] / 2.0 + GCparams['grating_pad_buffer']
    grating_periods_x_list = numpy.zeros(len(air_hole_diameter_list)) + GCparams['grating_period_x_start']
    grating_phi_res = round((360.0/num_grating_airhole_points),1)
    print("angle res for grating" + str(grating_phi_res))

    for iG in range(1, len(air_hole_diameter_list), 1):
        grating_periods_x_list[iG] = grating_periods_x_list[iG - 1] + GCparams['phaseFactor'] * GCparams['grating_delta_index'] * \
                                     grating_periods_x_list[iG - 1] * grating_periods_x_list[iG - 1] / GCparams['resonance']

    for iGx in range(GCparams['num_grating_periods_x']):
        for iGy in range(GCparams['num_grating_periods_y']):

            grating_column_start_x = grating_start_x + air_hole_diameter_list[iGx] / 2.0 + numpy.sum(
                grating_periods_x_list[:iGx])
            grating_column_start_y = grating_start_y + numpy.sqrt(3) * GCparams['a_2DPhC'] * iGy
            hole = pg.circle(radius=air_hole_diameter_list[iGx] / 2.0, angle_resolution=grating_phi_res, layer=GC_hole_layer)
            hole1 = SWG_holes << hole
            hole2 = SWG_holes << hole
            if iGy < (GCparams['num_grating_periods_y'] - 1):
                hole3 = SWG_holes << hole
                hole3.move(origin=[0, 0], destination=[grating_column_start_x + GCparams['a_2DPhC'] / 2.0,
                                                       grating_column_start_y + GCparams['a_2DPhC'] / 2.0 * numpy.sqrt(3)])

            hole1.move(origin=[0, 0], destination=[grating_column_start_x, grating_column_start_y])
            hole2.move(origin=[0, 0], destination=[grating_column_start_x + GCparams['a_2DPhC'], grating_column_start_y])
            # SWG_holes.add([hole1,hole2,hole3])
            # GC_holes_cell.add(gdspy.Round([grating_column_start_x, grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
            # GC_holes_cell.add(gdspy.Round([grating_column_start_x + a_2DPhC_nm, grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
            # if iGy < (num_grating_periods_y - 1):
            #     GC_holes_cell.add(gdspy.Round( [grating_column_start_x + a_2DPhC_nm / 2.0,grating_column_start_y + a_2DPhC_nm / 2.0 * numpy.sqrt(3)],air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
    return SWG_holes


def generate_silicon_skeleton(
    GCparams,
    num_tether_along_taper,
    GC_tether_x_nm
):
    test_PhC_bus_len = 14000
    P1 = pp.straight(length=GCparams['grating_pad_length'])
    P2 = pp.straight(length=test_PhC_bus_len)
    pad_ref = P1.extrude(GCparams['grating_pad_width'], layer=2)
    PhC_WG_ref = P2.extrude(GCparams['bus_wg_width'], layer=2)

    # define the taper
    P3 = pp.straight(length=GCparams['taper_length'])
    # Use the transitional CrossSection to create a Device
    WG_trans = P3.extrude([GCparams['grating_pad_width'], GCparams['bus_wg_width']], layer=2)

    # add them all to a device together
    silicon_skeleton = Device()
    pad = silicon_skeleton << pad_ref
    PhC_WG = silicon_skeleton << PhC_WG_ref
    PhC_WG.movex(GCparams['taper_length'] + GCparams['grating_pad_length'])
    wg3 = silicon_skeleton << WG_trans
    wg3.movex(GCparams['grating_pad_length'])

    silicon_skeleton_cutout_layer = 1

    pad_and_taper = Device()
    taper_wg = pad_and_taper << WG_trans
    taper_wg.movex(GCparams['grating_pad_length'])
    # pad_and_taper.add_ref(WG_trans).movex(GCparams['grating_pad_length'])
    pad_and_taper.add_ref(pad_ref)
    # outline = pg.outline(pad_and_taper,distance=1000,precision=1e-9,max_points=8000,layer=silicon_skeleton_cutout_layer)
    inversion = pg.invert(pad_and_taper, border=GCparams['cutout_around_GC_taper'], precision=1e-9, layer=silicon_skeleton_cutout_layer)
    skeleton_and_outline = Device()
    skeleton_and_outline.add_ref(inversion)
    skeleton_and_outline.add_ref(pad_and_taper)

    outline_taper = pg.outline(wg3,distance=GCparams['cutout_around_GC_taper'],precision=1e-9,layer=3)
    skeleton_and_outline.add_ref(outline_taper)
    silicon_cutout_only = pg.boolean(A=outline_taper,B=pad_and_taper,operation='not',precision=1e-9,num_divisions=[1,1],layer=4)
    # skeleton_and_outline.add_ref(silicon_cutout_only)
    # skeleton_larger = pad_and_taper.copy(name='newcell',scale=1.1)
    # skeleton_larger.remap_layers({2:3})
    #get the outline itself that should remain
    tether_spacing = GCparams['taper_length']/(1+num_tether_along_taper)
    tethers_collection = Device()
    # pad_start_x = GCparams['grating_pad_length']
    for tether_index in range(0,num_tether_along_taper+2):
        #calculate tether x
        tether_path = pp.straight(length=GC_tether_x_nm)
        tether_ref = tether_path.extrude(silicon_cutout_only.ysize,layer=2)
        if tether_index == 0:
            tether_ref.xmin = silicon_cutout_only.xmin
        else:
            tetherx = tether_index * tether_spacing
            tether_ref.movex(tetherx+GCparams['grating_pad_length'])
        tethers_collection.add_ref(tether_ref)
        # tether_ref.move(origin=[taper_wg.xmin,pad_ref.y],destination=[taper_wg.xmin+tetherx,pad_ref.y])
        # tether_ref.movex(tetherx)

    cutout_including_tethers = pg.boolean(A=silicon_cutout_only,B=tethers_collection,operation='not',precision=1e-9,num_divisions=[1,1],layer=0)
    skeleton_and_outline.add_ref(cutout_including_tethers)
    return [skeleton_and_outline,pad.center]

def subwavelength_grating(
    air_hole_diameter_list_base,
    constgrating_airhole_scale_factor,
    GCparams,
    grating_pad_center,
    GC_hole_layer,
    num_tether_along_taper,
    GC_tether_x_nm
):
    grating_holes = generate_grating_holes(air_hole_diameter_list_base,
    constgrating_airhole_scale_factor,
    GCparams,
    grating_pad_center,
    GC_hole_layer,
    GCparams['num_GC_circle_points'])

    [silicon_stuff_all_layers,pad_center] = generate_silicon_skeleton(GCparams,num_tether_along_taper,GC_tether_x_nm)
    # silicon_stuff_all_layers.write_gds('checking_cutout_gc.gds',unit=1e-9,precision=1e-12)

    #align the silicon skeleton and holes
    grating_holes.center=pad_center
    # grating_holes.rotate(180)
    grating_holes.mirror(p1=(grating_holes.x,grating_holes.ymin),p2=(grating_holes.x,grating_holes.ymax))
    # silicon_stuff.y = grating_holes.y

    #subtract the outline to get just the silicon remaining
    grating_coupler=Device()
    silicon_stuff = pg.extract(silicon_stuff_all_layers,layers=[0])
    # silicon_stuff.remove_layers()
    grating_coupler.add_ref(grating_holes)
    grating_coupler.add_ref(silicon_stuff)
    # grating_coupler.add_ref(skeleton_outline)

    #generate the taper to connect both

    GC_cutout_taper_path = pp.straight(length=GCparams['cutout_taper_length'])
    connector_taper_y_start = GCparams['bus_wg_width'] + 2*GCparams['cutout_around_GC_taper']
    connector_taper_y_end = GCparams['bus_wg_width'] + 2 * GCparams['bus_wg_to_phc_wg_spacing']
    gc_to_bus_taper_outline_ref = GC_cutout_taper_path.extrude([connector_taper_y_start, connector_taper_y_end], layer=10)
    bus_wg_connector_ref = GC_cutout_taper_path.extrude(GCparams['bus_wg_width'],layer=10)
    gc_to_bus_taper_cutout_ref = pg.boolean(A=gc_to_bus_taper_outline_ref,B=bus_wg_connector_ref,operation='not',precision=1e-9,num_divisions=[1,1],layer=0)
    gc_to_bus_taper_cutout = grating_coupler << gc_to_bus_taper_cutout_ref
    gc_to_bus_taper_cutout.y=silicon_stuff.y
    gc_to_bus_taper_cutout.xmin=silicon_stuff.xmax

    grating_coupler.write_gds('grating_coupler_check.gds', unit=1e-9, precision=1e-12)

    #define and add holes air cutout
    GC_pad_side_cutout_path = pp.straight(length=grating_holes.xsize)
    GC_pad_side_cutout_ref = GC_pad_side_cutout_path.extrude(GCparams['grating_pad_spacing'], layer=0)
    GC_pad_top_cutout = grating_coupler << GC_pad_side_cutout_ref
    GC_pad_top_cutout.x = grating_holes.x
    GC_pad_top_cutout.ymin = grating_holes.ymax + GCparams['gc_holes_tether_y']

    GC_pad_bot_cutout = grating_coupler << GC_pad_side_cutout_ref
    GC_pad_bot_cutout.x = grating_holes.x
    GC_pad_bot_cutout.ymax = grating_holes.ymin - GCparams['gc_holes_tether_y']

    return grating_coupler

