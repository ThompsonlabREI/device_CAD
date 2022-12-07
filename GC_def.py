import numpy
import gdspy
from phidl import Device, Path, CrossSection
from phidl import Group
import phidl.geometry as pg
import phidl.path as pp

def generate_grating_holes(air_hole_diameter_list_base,
    constgrating_airhole_scale_factor,
    grating_period_x_start,
    phase_factor,
    grating_delta_index,
    resonance_freq,
    num_grating_periods_x,
    num_grating_periods_y,
    num_points,
    a_2DPhC_nm,
    grating_pad_length,
    grating_pad_width,
    grating_pad_buffer,
    grating_pad_center_x,
    grating_pad_center_y,
    layer_num):
    SWG_holes = Device()

    air_hole_diameter_list = numpy.asarray(
        air_hole_diameter_list_base) * constgrating_airhole_scale_factor  # SRP: varies from 0.85 to 1.05 for first run
    # air_hole_diameter_list = numpy.asarray(air_hole_diameter_list_base) + constgrating_airholescale_list[iB] # take off the 1.13 scaling factor so that GC hole sizes match lumerical
    grating_start_x = grating_pad_center_x - grating_pad_length / 2.0 + grating_pad_buffer
    grating_start_y = grating_pad_center_y - grating_pad_width / 2.0 + grating_pad_buffer
    grating_periods_x_list = numpy.zeros(len(air_hole_diameter_list)) + grating_period_x_start
    for iG in range(1, len(air_hole_diameter_list), 1):
        grating_periods_x_list[iG] = grating_periods_x_list[iG - 1] + phase_factor * grating_delta_index * \
                                     grating_periods_x_list[iG - 1] * grating_periods_x_list[iG - 1] / resonance_freq

    for iGx in range(num_grating_periods_x):
        for iGy in range(num_grating_periods_y):

            grating_column_start_x = grating_start_x + air_hole_diameter_list[iGx] / 2.0 + numpy.sum(
                grating_periods_x_list[:iGx])
            grating_column_start_y = grating_start_y + numpy.sqrt(3) * a_2DPhC_nm * iGy
            hole = pg.circle(radius=air_hole_diameter_list[iGx] / 2.0, angle_resolution=2.5, layer=layer_num)
            # hole2 = pg.circle(radius = air_hole_diameter_list[iGx] / 2.0, angle_resolution=2.5,layer=0)
            # hole3 = pg.circle(radius = air_hole_diameter_list[iGx] / 2.0, angle_resolution=2.5,layer=0)
            hole1 = SWG_holes << hole
            hole2 = SWG_holes << hole
            if iGy < (num_grating_periods_y - 1):
                hole3 = SWG_holes << hole
                hole3.move(origin=[0, 0], destination=[grating_column_start_x + a_2DPhC_nm / 2.0,
                                                       grating_column_start_y + a_2DPhC_nm / 2.0 * numpy.sqrt(3)])

            hole1.move(origin=[0, 0], destination=[grating_column_start_x, grating_column_start_y])
            hole2.move(origin=[0, 0], destination=[grating_column_start_x + a_2DPhC_nm, grating_column_start_y])
            # SWG_holes.add([hole1,hole2,hole3])
            # GC_holes_cell.add(gdspy.Round([grating_column_start_x, grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
            # GC_holes_cell.add(gdspy.Round([grating_column_start_x + a_2DPhC_nm, grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
            # if iGy < (num_grating_periods_y - 1):
            #     GC_holes_cell.add(gdspy.Round( [grating_column_start_x + a_2DPhC_nm / 2.0,grating_column_start_y + a_2DPhC_nm / 2.0 * numpy.sqrt(3)],air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
    return SWG_holes

def subwavelength_grating(
    air_hole_diameter_list_base,
    constgrating_airhole_scale_factor,
    grating_period_x_start,
    phase_factor,
    grating_delta_index,
    resonance_freq,
    num_grating_periods_x,
    num_grating_periods_y,
    num_points,
    a_2DPhC_nm,
    grating_pad_length,
    grating_pad_width,
    grating_pad_buffer,
    grating_pad_center_x,
    grating_pad_center_y,
    layer_num,
    PhC_wy
):
    grating_holes = generate_grating_holes(air_hole_diameter_list_base,
    constgrating_airhole_scale_factor,
    grating_period_x_start,
    phase_factor,
    grating_delta_index,
    resonance_freq,
    num_grating_periods_x,
    num_grating_periods_y,
    num_points,
    a_2DPhC_nm,
    grating_pad_length,
    grating_pad_width,
    grating_pad_buffer,
    grating_pad_center_x,
    grating_pad_center_y,
    layer_num)

    #define taper

    # cross_section_width = CrossSection()

    #define pad


    #define bit at the end

    # X1 = CrossSection()
    # X1.add(width=1.2, offset=0, layer=2, name='wg')
    # X1.add(width=2.2, offset=0, layer=3, name='etch')
    # X1.add(width=1.1, offset=3, layer=1, name='wg2')

    # Create the second CrossSection that we want to transition to
    # X2 = CrossSection()
    # X2.add(width=1, offset=0, layer=2, name='wg')
    # X2.add(width=3.5, offset=0, layer=3, name='etch')
    # X2.add(width=3, offset=5, layer=1, name='wg2')

    # To show the cross-sections, let's create two Paths and
    # create Devices by extruding them
    P1 = pp.straight(length=grating_pad_length)
    P2 = pp.straight(length=14000)
    WG1 = P1.extrude(grating_pad_width,layer=2)
    WG2 = P2.extrude(PhC_wy,layer=2)

    # Place both cross-section Devices and quickplot them

    # Create the transitional CrossSection
    # Xtrans = pp.transition(cross_section1=X1,
    #                        cross_section2=X2,
    #                        width_type='sine')
    # # Create a Path for the transitional CrossSection to follow
    P3 = pp.straight(length=180000)
    # Use the transitional CrossSection to create a Device
    WG_trans = P3.extrude([grating_pad_width,PhC_wy],layer=2)

    silicon_skeleton = Device()
    wg1 = silicon_skeleton << WG1
    wg2 = silicon_skeleton << WG2
    wg2.movex(180000)
    wg3 = silicon_skeleton << WG_trans

    return [grating_holes,silicon_skeleton]

