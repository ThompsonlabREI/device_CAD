import numpy
import gdspy

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
    GC_holes_cell,
):

    air_hole_diameter_list = numpy.asarray(air_hole_diameter_list_base) * constgrating_airhole_scale_factor  # SRP: varies from 0.85 to 1.05 for first run
# air_hole_diameter_list = numpy.asarray(air_hole_diameter_list_base) + constgrating_airholescale_list[iB] # take off the 1.13 scaling factor so that GC hole sizes match lumerical
    grating_start_x = grating_pad_center_x - grating_pad_length / 2.0 + grating_pad_buffer
    grating_start_y = grating_pad_center_y - grating_pad_width / 2.0 + grating_pad_buffer
    grating_periods_x_list = numpy.zeros(len(air_hole_diameter_list)) + grating_period_x_start
    for iG in range(1, len(air_hole_diameter_list), 1):
        grating_periods_x_list[iG] = grating_periods_x_list[iG - 1] + phase_factor * grating_delta_index * \
                                 grating_periods_x_list[iG - 1] * grating_periods_x_list[iG - 1] / resonance_freq

    for iGx in range(num_grating_periods_x):
        for iGy in range(num_grating_periods_y):

            grating_column_start_x = grating_start_x + air_hole_diameter_list[iGx] / 2.0 + numpy.sum(grating_periods_x_list[:iGx])
            grating_column_start_y = grating_start_y + numpy.sqrt(3) * a_2DPhC_nm * iGy

            GC_holes_cell.add(gdspy.Round([grating_column_start_x, grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
            GC_holes_cell.add(gdspy.Round([grating_column_start_x + a_2DPhC_nm, grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
            if iGy < (num_grating_periods_y - 1):
                GC_holes_cell.add(gdspy.Round( [grating_column_start_x + a_2DPhC_nm / 2.0,grating_column_start_y + a_2DPhC_nm / 2.0 * numpy.sqrt(3)],air_hole_diameter_list[iGx] / 2.0, number_of_points=num_points))
    return GC_holes_cell