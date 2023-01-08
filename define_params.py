import numpy
import gdspy
from phidl import Device, Path, CrossSection
from phidl import Group
import phidl.geometry as pg
import phidl.path as pp

def generate_photonics_params():
    PhCparams = {}
    PhCparams['num_PhC_per_wg'] = 2
    PhCparams['num_PhC_per_GC'] = 2 * PhCparams['num_PhC_per_wg']
    PhCparams['PhC_hx'] = 236.8
    PhCparams['PhC_hy'] = 407.03
    PhCparams['aper_cav'] = 349.1  # 298
    PhCparams['aper_mir'] = 435.2  # 343
    PhCparams['bus_wg_width'] = 392
    PhCparams['PhC_wy'] = 600
    PhCparams['bus_wg_to_phc_wg_spacing'] = 500
    PhCparams['outer_cutout_phc_wg'] = 1500
    PhCparams['num_bus_reflector_mirrors'] = 10
    PhCparams['bus_reflector_taper_len'] = 2000
    PhCparams['beam_tether_x'] = 500
    PhCparams['bus_taper_xlen'] = 500
    PhCparams['bus_taper_chonky_y'] = 1000
    PhCparams['bus_reflect_taper_len_x'] = PhCparams['bus_reflector_taper_len']
    PhCparams['num_cavity_holes']=12
    PhCparams['num_mirror_holes_middle'] = 2
    PhCparams['num_mirror_holes_end'] = 8
    PhCparams['phc_beam_buffer_x'] = 2*PhCparams['beam_tether_x']
    PhCparams['num_ellipse_points']=199

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
    GCparams['grating_pad_buffer'] = (GCparams['grating_pad_width'] - (GCparams['num_grating_periods_y'] - 1) *
                                      GCparams['a_2DPhC'] * numpy.sqrt(3)) / 2.0
    GCparams['grating_pad_offset'] = 2000
    GCparams['grating_pad_spacing'] = 3000
    GCparams['resonance'] = 1544
    GCparams['n_circle_points_GC'] = 110
    GCparams['PhC_wy'] = 600
    GCparams['taper_length'] = 185000
    GCparams['cutout_taper_length'] = 2000
    GCparams['bus_wg_width']=PhCparams['bus_wg_width']
    GCparams['cutout_around_GC_taper'] = 1000
    GCparams['bus_wg_to_phc_wg_spacing']=PhCparams['bus_wg_to_phc_wg_spacing']
    GCparams['gc_holes_tether_y'] = 175
    GCparams['num_GC_circle_points']=144

    return [PhCparams,GCparams]