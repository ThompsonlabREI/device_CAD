#!/usr/bin/python

import numpy
from copy import *
import math
from itertools import chain
import gdspy
# from gdsii.library import Library
# from gdsii.structure import *
# from gdsii.elements import *

# layers
# 1: waveguide
# 2: outerbox
# 3: notches
# 4: holes
# 5: guard
lGuard = 5

def write_beams(cell, param):
	""" write holes with the specified parameters, appending to the structure cell """
	#define holes as a boundary with 15 points at these angles
	# since 0 = 2*pi, the boundary will be closed

	global iX,iY

	philist=numpy.linspace(0,2*numpy.pi,param['n_circle_points'])

	if param['straight_q'] is True:
		# outerbox_x_max = param['array_orig_x'] + param['beam_len']
		# outerbox_x_min = param['array_orig_x']# - param['beam_len']/2.0

		outerbox_x_max = param['array_orig_x']

		if param['meander'] is True:
			outerbox_x_min = param['array_orig_x'] - param['beam_len'] - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2 + param2['vert_linker_width_left'] + 1.5 * param['beam_width']
		else:
			outerbox_x_min = param['array_orig_x'] - param['beam_len']

	else:

		outerbox_x_max = param['array_orig_x'] + param['beam_len'] + numpy.cos(param['bend_angle'])*(param['stick_len']) + numpy.sin(param['bend_angle'])*param['bend_rad'] + 500
		outerbox_x_min = param['array_orig_x']# - param['beam_len']/2.0

	# asymmetric shape if using vertical flag beams (vflagbeam) for alignment along the edge of the YSO
	if param['vflagbeam_q'] is True:
		outerbox_y_min = param['array_orig_y'] - param['beam_spacing']
		outerbox_y_max = param['array_orig_y'] + (param['num_beams']+3)*param['beam_spacing']
		outerbox_y_max_left = param['array_orig_y'] + (param['num_beams']+1)*param['beam_spacing']
		outerbox_x_min_left = param['array_orig_x'] - param['vert_align_offset'] - (param['end_width_taper_length'] + param['end_width_length'])
 		outerbox_pts = gdspy.Polygon(2,[
		(outerbox_x_min, outerbox_y_min),
		(outerbox_x_min, outerbox_y_max_left),
		(outerbox_x_min_left, outerbox_y_max_left),
		(outerbox_x_min_left, outerbox_y_max),
		(outerbox_x_max, outerbox_y_max),
		(outerbox_x_max, outerbox_y_min),
		(outerbox_x_min, outerbox_y_min)
		])
		outerbox_pts.fillet(2000) # AD Use fillet command to round corners by 2 um radius, relieves stress on structures
		cell.add(outerbox_pts)

	# normal rounded rectangle if not using vertical flag beams

	else:

		if param['meander'] is True:
			if param["2-axes"]==False:
				outerbox_y_min = param['array_orig_y'] -param['length_total_meander']- param['supporting_bar_width']- param['spacing_phc-vertical_buffer']- param['box_buffer']-param['beam_width'] / 2  #-param['width_pad_above_meander'] # SC make top and bottom outerbox symmetric
				outerbox_y_max = param['array_orig_y'] +param['length_total_meander']+ param['supporting_bar_width']+ param['spacing_phc-vertical_buffer']+ param['box_buffer']+param['beam_width'] / 2 +param['beam_spacing'] # +param['width_pad_above_meander']
				# + (param['num_beams']) * param['beam_spacing'] + param['box_buffer'] + param['length_total_meander']
			if param["2-axes"]==True:
				outerbox_y_min = param['array_orig_y'] -param['length_total_meander']- param['supporting_bar_width']- param['box_buffer']-param['beam_width'] / 2 -param['y_spacing_between_wabguides'] - param['Len_bus_waveguide_vertical'] #-param['width_pad_above_meander'] # SC make top and bottom outerbox symmetric
				outerbox_y_max = param['array_orig_y'] +param['length_total_meander']+ param['supporting_bar_width']+ param['spacing_phc-vertical_buffer']+ param['box_buffer']+param['beam_width'] / 2 +param['beam_spacing'] # +param['width_pad_above_meander']

		else:
			# outerbox_y_min = param['array_orig_y'] - param['beam_spacing'] + param['grating_pad_spacing']
			outerbox_y_min = param['array_orig_y'] - param['beam_spacing'] - param['box_buffer']  # SC make top and bottom outerbox symmetric
			outerbox_y_max = param['array_orig_y'] + (param['num_beams']) * param['beam_spacing'] + param['box_buffer']

		outerbox_pts = gdspy.Polygon(2,[
		(outerbox_x_min,outerbox_y_min),
		(outerbox_x_min,outerbox_y_max),
		(outerbox_x_max,outerbox_y_max),
		(outerbox_x_max,outerbox_y_min),
		(outerbox_x_min,outerbox_y_min)
		])

		outerbox_pts.fillet(2000) # AD Use fillet command to round corners by 2 um radius, relieves stress on structures

		cell.add(outerbox_pts)

	if iX==0 and iY==0:

		first_WF_x=(outerbox_x_min + outerbox_x_max)/2.0
		first_WF_y=(outerbox_y_min + outerbox_y_max)/2.0

		first_WF_pts=gdspy.Polygon(8,[
			(first_WF_x, first_WF_y - 1.0e3),
			(first_WF_x + 1.0e3, first_WF_y),
			(first_WF_x, first_WF_y + 1.0e3),
			(first_WF_x - 1.0e3, first_WF_y)])

		#cell.add(first_WF_pts)

	if param['support_tether_q'] is True: # SC adding horizantal supporting tether
		if param['meander'] is True:
			beam_center_x = param['array_orig_x'] - param['beam_len'] / 2.0
			circle_center_x = beam_center_x - param['beam_len'] / 2.0 + 2 * param2['vert_linker_width_left'] + param2['vert_linker_offset'] + 1.5 * param['beam_width']
			cell.add(gdspy.Polygon(1, [  # bottom block holder
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
				(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
				(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2)
			]))
			cell.add(gdspy.Polygon(1, [  # top block holder
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
				(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
				(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2)
			]))

		else:
			cell.add(gdspy.Polygon(1, [  # bottom block holder
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
				(outerbox_x_min + param['vert_linker_offset'] + 2 * param['vert_linker_width_left'] + param['vert_linker_gap'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
				(outerbox_x_min + param['vert_linker_offset'] + 2 * param['vert_linker_width_left'] + param['vert_linker_gap'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2)
			]))
			cell.add(gdspy.Polygon(1, [  # top block holder
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
				(outerbox_x_min + param['vert_linker_offset'] + 2 * param['vert_linker_width_left'] + param['vert_linker_gap'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
				(outerbox_x_min + param['vert_linker_offset'] + 2 * param['vert_linker_width_left'] + param['vert_linker_gap'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
				(outerbox_x_max - param['grating_pad_offset'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2)
			]))

		cell.add(gdspy.Polygon(1, [# SC add a vertical supporting beam holding left end of grating pad
			(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
			(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
			(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
			(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - param['support_tether_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
			(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - param['support_tether_width'] / 2,outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width'])
		]))

		for iT in range(param['num_tether_taper']):
			# cell.add(gdspy.Polygon(1, [
			# 	(outerbox_x_max - param['grating_pad_offset'] - (iT + 1) * (param['grating_taper_length'] + param['grating_pad_length']) / (param['num_tether_taper'] + 1) - param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
			# 	(outerbox_x_max - param['grating_pad_offset'] - (iT + 1) * (param['grating_taper_length'] + param['grating_pad_length']) / (param['num_tether_taper'] + 1) + param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
			# 	(outerbox_x_max - param['grating_pad_offset'] - (iT + 1) * (param['grating_taper_length'] + param['grating_pad_length']) / (param['num_tether_taper'] + 1) + param['support_tether_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
			# 	(outerbox_x_max - param['grating_pad_offset'] - (iT + 1) * (param['grating_taper_length'] + param['grating_pad_length']) / (param['num_tether_taper'] + 1) - param['support_tether_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
			# 	(outerbox_x_max - param['grating_pad_offset'] - (iT + 1) * (param['grating_taper_length'] + param['grating_pad_length']) / (param['num_tether_taper'] + 1) - param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width'])
			# ]))
			cell.add(gdspy.Polygon(1, [#SC adding vertical supporting tether beam
				(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - (iT + 1) * (param['grating_taper_length']) / (param['num_tether_taper'] + 1) - param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
				(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - (iT + 1) * (param['grating_taper_length']) / (param['num_tether_taper'] + 1) + param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
				(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - (iT + 1) * (param['grating_taper_length']) / (param['num_tether_taper'] + 1) + param['support_tether_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
				(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - (iT + 1) * (param['grating_taper_length']) / (param['num_tether_taper'] + 1) - param['support_tether_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
				(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] - (iT + 1) * (param['grating_taper_length']) / (param['num_tether_taper'] + 1) - param['support_tether_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width'])
			]))

		if param['meander'] is True:
			for iTT in range(param['num_tether_device']):  # SC adding vertical block pinch points
				cell.add(gdspy.Polygon(1, [  # bottom block linkers
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_min + param['beam_width'] / 2 - param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_min + param['beam_width'] / 2 - param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_min + param['beam_width'] / 2 - param['beam_width'] / 2)
				]))
				cell.add(gdspy.Polygon(3, [  # bottom block pinch points
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] - param['notch_end_width'] / 2 + param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] + param['notch_end_width'] / 2 + param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2 - param['notch_end_depth'], outerbox_y_min + param['notch_end_offset'] + param['beam_width'] / 2)]))

				cell.add(gdspy.Polygon(3, [  # bottom block pinch points
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] - param['notch_end_width'] / 2 + param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] + param['notch_end_width'] / 2 + param['beam_width'] / 2),
					(outerbox_x_min+ param['beam_width'] / 2+param['grating_pad_offset'] + (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2 + param['notch_end_depth'], outerbox_y_min + param['notch_end_offset'] + param['beam_width'] / 2)]))

				cell.add(gdspy.Polygon(1, [  # top block linkers
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_max - param['beam_width'] / 2 + param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_max - param['beam_width'] / 2 + param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_max - param['beam_width'] / 2 + param['beam_width'] / 2)
				]))
				cell.add(gdspy.Polygon(3, [  # top block pinch points
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] - param['notch_end_width'] / 2 - param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] + param['notch_end_width'] / 2 - param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) + param['beam_width'] / 2 - param['notch_end_depth'], outerbox_y_max - param['notch_end_offset'] - param['beam_width'] / 2)]))
				cell.add(gdspy.Polygon(3, [  # top block pinch points
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] - param['notch_end_width'] / 2 - param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] + param['notch_end_width'] / 2 - param['beam_width'] / 2),
					(outerbox_x_max- param['beam_width'] / 2-param['grating_pad_offset'] - (iTT ) * param['beam_len'] / (param['num_tether_device'] ) - param['beam_width'] / 2 + param['notch_end_depth'], outerbox_y_max - param['notch_end_offset'] - param['beam_width'] / 2)]))

		else:
			for iTT in range(param['num_tether_device']):  # SC adding vertical block pinch points
				cell.add(gdspy.Polygon(1, [  # bottom block linkers
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_min),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_min),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_min + param['box_buffer'] - param['beam_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_min + param['box_buffer'] - param['beam_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_min)
				]))
				cell.add(gdspy.Polygon(3, [  # bottom block pinch points
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] - param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] + param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2 - param['notch_end_depth'], outerbox_y_min + param['notch_end_offset'])]))
				cell.add(gdspy.Polygon(3, [  # bottom block pinch points
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] - param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_min + param['notch_end_offset'] + param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2 + param['notch_end_depth'], outerbox_y_min + param['notch_end_offset'])]))

				cell.add(gdspy.Polygon(1, [  # top block linkers
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_max),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_max),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_max - param['box_buffer'] + param['beam_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_max - param['box_buffer'] + param['beam_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_max)
				]))
				cell.add(gdspy.Polygon(3, [  # top block pinch points
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] - param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] + param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) + param['beam_width'] / 2 - param['notch_end_depth'], outerbox_y_max - param['notch_end_offset'])]))
				cell.add(gdspy.Polygon(3, [  # top block pinch points
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] - param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2, outerbox_y_max - param['notch_end_offset'] + param['notch_end_width'] / 2),
					(outerbox_x_max - (iTT + 1) * param['beam_len'] / (param['num_tether_device'] + 1) - param['beam_width'] / 2 + param['notch_end_depth'], outerbox_y_max - param['notch_end_offset'])]))

	for iB in range(-1,param['num_beams'] + 1): # Use range(0,param['num_beams'] + 1) with vertical linker
		#SC, starting from -1 to add bottom horizontal linker
		cell_edge_x = param['array_orig_x'] + param['beam_len']/2.0
		hole_center_y = param['array_orig_y'] + iB*param['beam_spacing']
		beam_center_x = param['array_orig_x'] - param['beam_len']/2.0

		# write the straight part of the beam
		circle_center_x = beam_center_x - param['beam_len'] / 2.0 + 2 * param2['vert_linker_width_left'] + param2['vert_linker_offset'] + 1.5 * param['beam_width']
		if param['meander'] is True:
			circle_center_x = beam_center_x - param['beam_len'] / 2.0 + 2 * param2['vert_linker_width_left'] + param2['vert_linker_offset'] + 1.5 * param['beam_width']

			if iB == 0 or iB == 1:
				cell.add(gdspy.Polygon(1, [
					(circle_center_x, hole_center_y - param['width_taper'] / 2.0),
					(circle_center_x, hole_center_y + param['width_taper'] / 2.0),
					(beam_center_x + param['beam_len'] / 2.0 - param2['vert_linker_offset'], hole_center_y + param['width_taper'] / 2.0),
					(beam_center_x + param['beam_len'] / 2.0 - param2['vert_linker_offset'], hole_center_y - param['width_taper'] / 2.0),
					(circle_center_x, hole_center_y - param['width_taper'] / 2.0)
				]))
		else:
			cell.add(gdspy.Polygon(1, [
				(beam_center_x - param['beam_len'] / 2.0, hole_center_y - param['beam_width'] / 2.0),
				(beam_center_x - param['beam_len'] / 2.0, hole_center_y + param['beam_width'] / 2.0),
				(beam_center_x + param['beam_len'] / 2.0, hole_center_y + param['beam_width'] / 2.0),
				(beam_center_x + param['beam_len'] / 2.0, hole_center_y - param['beam_width'] / 2.0),
				(beam_center_x - param['beam_len'] / 2.0, hole_center_y - param['beam_width'] / 2.0)
			]))

		if param['guard_q'] is True:
			cell.add(gdspy.Polygon(lGuard,[
				(beam_center_x - param['beam_len']/2.0, hole_center_y - param['beam_width']/2.0 - param['guard_width']),
				(beam_center_x - param['beam_len']/2.0, hole_center_y + param['beam_width']/2.0 + param['guard_width']),
				(beam_center_x + param['beam_len']/2.0, hole_center_y + param['beam_width']/2.0 + param['guard_width']),
				(beam_center_x + param['beam_len']/2.0, hole_center_y - param['beam_width']/2.0 - param['guard_width']),
				(beam_center_x - param['beam_len']/2.0, hole_center_y - param['beam_width']/2.0 - param['guard_width'])
				]))

		# write the grating pad
		if param['grating_pad_q'] is True and iB < param['num_beams'] and iB > -1:

			grating_pad_center_x = beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset'] - param['grating_pad_length']/2.0
			grating_pad_center_y = hole_center_y

			grating_start_x = grating_pad_center_x - param['grating_pad_length']/2.0 + param['grating_pad_buffer']
			grating_start_y = grating_pad_center_y - param['grating_pad_width']/2.0 + param['grating_pad_buffer']


			# add rectangular pad for grating coupler
			grating_pad_pts = gdspy.Polygon(1,[
				(grating_pad_center_x - param['grating_pad_length']/2.0, grating_pad_center_y - param['grating_pad_width']/2.0),
				(grating_pad_center_x - param['grating_pad_length']/2.0, grating_pad_center_y + param['grating_pad_width']/2.0),
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y + param['grating_pad_width']/2.0),
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y - param['grating_pad_width']/2.0),
				(grating_pad_center_x - param['grating_pad_length']/2.0, grating_pad_center_y - param['grating_pad_width']/2.0)
				])

			cell.add(grating_pad_pts)

			# add top linker beam for grating coupler
			grating_beam_top_pts=gdspy.Polygon(1,[
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y + param['grating_pad_width']/2.0),
				(grating_pad_center_x + param['grating_pad_length']/2.0 + param['grating_pad_offset'], grating_pad_center_y + param['grating_pad_width']/2.0),
				(grating_pad_center_x + param['grating_pad_length']/2.0 + param['grating_pad_offset'], grating_pad_center_y + param['grating_pad_width']/2.0 - param['beam_width']),
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y + param['grating_pad_width']/2.0 - param['beam_width']),
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y + param['grating_pad_width']/2.0)
				])

			cell.add(grating_beam_top_pts)

			# add bottom linker beam for grating coupler
			grating_beam_bot_pts=gdspy.Polygon(1,[
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y - param['grating_pad_width']/2.0),
				(grating_pad_center_x + param['grating_pad_length']/2.0 + param['grating_pad_offset'], grating_pad_center_y - param['grating_pad_width']/2.0),
				(grating_pad_center_x + param['grating_pad_length']/2.0 + param['grating_pad_offset'], grating_pad_center_y - param['grating_pad_width']/2.0 + param['beam_width']),
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y - param['grating_pad_width']/2.0 + param['beam_width']),
				(grating_pad_center_x + param['grating_pad_length']/2.0, grating_pad_center_y - param['grating_pad_width']/2.0)
				])

			cell.add(grating_beam_bot_pts)

			# Now write the 2D triangular air holes for the coupling grating
			# air_hole_diameter_list = numpy.asarray(air_hole_diameter_list_base) + grating_airholescale_list[iB] / (param['num_grating_periods_x'] - 1.0) * numpy.asarray(range(param['num_grating_periods_x']))#SC AirGapwidth linear Offset
			# air_hole_diameter_list = (numpy.asarray(air_hole_diameter_list_base) + constgrating_airholescale_list[iB])/1.13 # SO 2021
			air_hole_diameter_list = numpy.asarray(air_hole_diameter_list_base) + constgrating_airholescale_list[iB] # take off the 1.13 scaling factor so that GC hole sizes match lumerical

			grating_periods_x_list = numpy.zeros(len(air_hole_diameter_list)) + param['grating_period_x_start']
			for iG in range(1, len(air_hole_diameter_list), 1):
				grating_periods_x_list[iG] = grating_periods_x_list[iG - 1] + param['phaseFactor'] * param['grating_delta_index'] * grating_periods_x_list[iG - 1] * grating_periods_x_list[iG - 1] / param2["resonance"]


			for iGx in range(param['num_grating_periods_x']):
				for iGy in range(param['num_grating_periods_y']):

					grating_column_start_x = grating_start_x + air_hole_diameter_list[iGx] / 2.0 + numpy.sum(grating_periods_x_list[:iGx])
					grating_column_start_y = grating_start_y + numpy.sqrt(3) * param['a_2DPhC'] *iGy

					cell.add(gdspy.Round(3, [grating_column_start_x, grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points = param['n_circle_points_GC']))
					cell.add(gdspy.Round(3, [grating_column_start_x + param['a_2DPhC'], grating_column_start_y], air_hole_diameter_list[iGx] / 2.0, number_of_points = param['n_circle_points_GC']))
					if iGy < (param['num_grating_periods_y'] -1):
						cell.add(gdspy.Round(3, [grating_column_start_x + param['a_2DPhC'] / 2.0, grating_column_start_y + param['a_2DPhC'] / 2.0 *numpy.sqrt(3)], air_hole_diameter_list[iGx] / 2.0, number_of_points = param['n_circle_points_GC']))

		if iB < param['num_beams']:

			# original beam narrowing taper (NO GRATINGS)
			if param['straight_q'] is True and param['grating_pad_q'] is False and param['end_width'] < param['beam_width'] and iB > -1:
				cell.add(gdspy.Polygon(3,[
					(beam_center_x + param['beam_len']/2.0, hole_center_y + param['end_width']/2.0),
					(beam_center_x + param['beam_len']/2.0, hole_center_y + param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['end_width_length'] - param['end_width_taper_length'], hole_center_y + param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['end_width_length'], hole_center_y + param['end_width']/2.0),
					(beam_center_x + param['beam_len']/2.0, hole_center_y + param['end_width']/2.0)
					]))

				cell.add(gdspy.Polygon(3,[
					(beam_center_x + param['beam_len']/2.0, hole_center_y - param['end_width']/2.0),
					(beam_center_x + param['beam_len']/2.0, hole_center_y - param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['end_width_length'] - param['end_width_taper_length'], hole_center_y - param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['end_width_length'], hole_center_y - param['end_width']/2.0),
					(beam_center_x + param['beam_len']/2.0, hole_center_y - param['end_width']/2.0)
					]))

			# SO meander Semi-circles
			if param['straight_q'] is True and param['meander'] is True and iB > -1:

				if iB==1: #face down
					S=1
				else:
					S=-1 #face up

				circle_center_x = beam_center_x - param['beam_len'] / 2.0 + 2 * param2['vert_linker_width_left'] + param2['vert_linker_offset'] + 1.5 * param['beam_width']
				circle_center_y = hole_center_y + S* param['bus_bend_radius']

				outer_radius350 = param['bus_bend_radius'] + param['ww'] / 2.0 #radius of circle of width ww
				inner_radius350 = param['bus_bend_radius'] - param['ww'] / 2.0

				outer_radius650 = param['bus_bend_radius'] + param['beam_width'] / 2.0  # radius of circle of width beam width (650 or 680nm)
				inner_radius650 = param['bus_bend_radius'] - param['beam_width'] / 2.0

				outer_radius = param['bus_bend_radius'] + param['width_taper'] / 2.0 #radius of circle of width equal to the width of the taper
				inner_radius = param['bus_bend_radius'] - param['width_taper'] / 2.0

				#add semi-circle rotation
				if param['waveguide_with_end_mirror'] % 2 == 1:
					q = 2
				else:
					q = 1
				for iM in range(int(param['waveguide_with_end_mirror']+q)/2):
					if iM != param2['waveguide_with_end_mirror']-1 or iM-1<1:
						cell.add(  # left semi-circle facing up
							gdspy.Round(1, (circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides'])),
										outer_radius, inner_radius,
										S * numpy.pi, S * 3 * numpy.pi / 2.0)
						)

						if iB==1:
							cell.add(  # left semi-circle facing down
								gdspy.Round(1, (circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + param['y_spacing_between_wabguides'] + 2 * iM * param['y_spacing_between_wabguides'])),
											outer_radius, inner_radius,
											S * numpy.pi, S * numpy.pi / 2.0)
							)


				for iM in range(int(param['waveguide_with_end_mirror']) / 2):
					cell.add( #right semi-circle facing up
						gdspy.Round(1, (circle_center_x + param['Len_bus_waveguide'],
										circle_center_y+S*(param['y_spacing_between_wabguides']+2*param['bus_bend_radius'] + 4* iM * param['bus_bend_radius']+ 2*iM*param['y_spacing_between_wabguides'])),
									outer_radius,
									inner_radius,
									-S*numpy.pi / 2.0, 0)
					)

					cell.add( #right semi-circle facing down
						gdspy.Round(1, (circle_center_x + param['Len_bus_waveguide'],
										circle_center_y+ S*(2*param['y_spacing_between_wabguides']+2*param['bus_bend_radius'] + 4* iM * param['bus_bend_radius']+ 2*iM*param['y_spacing_between_wabguides'])),
									outer_radius,
									inner_radius,
									S*numpy.pi / 2, 0)
					)

					#separation waveguide on right
					cell.add(gdspy.Polygon(1, [
						(circle_center_x+ param['Len_bus_waveguide'] - param['width_taper']/2+param['bus_bend_radius'], circle_center_y+  S*(param['y_spacing_between_wabguides']+ 4* iM * param['bus_bend_radius']+ 2*iM*param['y_spacing_between_wabguides']+ param['y_spacing_between_wabguides']+2*param['bus_bend_radius'])),
						(circle_center_x+ param['Len_bus_waveguide'] - param['width_taper']/2+param['bus_bend_radius'], circle_center_y+  S*(4* iM * param['bus_bend_radius']+ 2*iM*param['y_spacing_between_wabguides']+ param['y_spacing_between_wabguides']+2*param['bus_bend_radius'])),
						(circle_center_x+ param['Len_bus_waveguide'] + param['width_taper']/2+param['bus_bend_radius'], circle_center_y+  S*(4* iM * param['bus_bend_radius']+ 2*iM*param['y_spacing_between_wabguides']+ param['y_spacing_between_wabguides']+2*param['bus_bend_radius'])),
						(circle_center_x+ param['Len_bus_waveguide'] + param['width_taper']/2+param['bus_bend_radius'], circle_center_y+  S*(param['y_spacing_between_wabguides']+ 4* iM * param['bus_bend_radius']+ 2*iM*param['y_spacing_between_wabguides']+ param['y_spacing_between_wabguides']+2*param['bus_bend_radius'])),
					]))

				#support blocks at the left end of the devices:
					for iM in range(int(param['waveguide_with_end_mirror'])):
						# right side of the meander
						if iM % 2 == 0 and iM > 0:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param['Len_bus_waveguide'], circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * (iM - 1) * param['bus_bend_radius'] + (1 + (iM - 1)) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * (iM - 1) * param['bus_bend_radius'] + (1 + (iM - 1)) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'], circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))
						# left side of the meander
						if iM % 2 == 1:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * (iM - 1) * param['bus_bend_radius'] + (1 + (iM - 1)) * param['y_spacing_between_wabguides'])),
								(circle_center_x, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * (iM - 1) * param['bus_bend_radius'] + (1 + (iM - 1)) * param['y_spacing_between_wabguides'])),
								(circle_center_x, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))

				# taper waveguides
				for iM in range(0,int(param['waveguide_with_end_mirror']) / 2+1):

					#tapering the bus waveguide for tethering
					if iM < int(param['waveguide_with_end_mirror']+2) / 2 :
						for s1 in range(0, 2):
							if iM==0 and S==1: #for the first waveguide on top
								s1=0
							if iM==0 and S==-1: #for the first waveguide on bottom
								s1=1
							if iM==int(param['waveguide_with_end_mirror']+2) / 2-1and S==1and param['waveguide_with_end_mirror'] % 2 == 0: # for last waveguide on top, even
								s1=1
							if iM==int(param['waveguide_with_end_mirror']+2) / 2-1and S==-1and param['waveguide_with_end_mirror'] % 2 == 0: # for last waveguide on bottom, even
								s1=0

							if iB==1:
								C=0

								#right
								cell.add(gdspy.Polygon(1, [
									(circle_center_x+ param2['Len_bus_waveguide']/2+param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S *(4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] + param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x+ param2['Len_bus_waveguide']/2+param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S *(4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] - param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x+ param2['Len_bus_waveguide']/2 + param2['bus_taper_len']+param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S *(4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x+ param2['Len_bus_waveguide']/2 + param2['bus_taper_len']+param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S *(4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								]))

								#left
								cell.add(gdspy.Polygon(1, [
									(circle_center_x + param2['Len_bus_waveguide']/2-param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S * (4* iM *  param['bus_bend_radius'] + 2* iM *  param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] - param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']/2-param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S * (4* iM *  param['bus_bend_radius'] + 2* iM *  param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] + param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']/2-param2['length_wider_bus_wavguide_part']/2 - param2['bus_taper_len']+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']/2-param2['length_wider_bus_wavguide_part']/2 - param2['bus_taper_len']+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM *  param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								]))

								#center
								cell.add(gdspy.Polygon(1, [
									(circle_center_x + param2['Len_bus_waveguide']/2+param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] + param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']/2-param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] + param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']/2-param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] - param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']/2+param2['length_wider_bus_wavguide_part']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + param['bus_bend_radius'] - param['width_taper_middle'] / 2.0 - 2 * s1 * param['bus_bend_radius']),
								]))

								# tether for waveguide
								cell.add(gdspy.Polygon(1, [
									(circle_center_x + param2['Len_bus_waveguide'] / 2 - param2['width_tether_bus_waveguide'] / 2+2*param['aper_cav'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius'] + param2['wg'] - (outer_radius - inner_radius - param2['ww']) / 2),
									(circle_center_x + param2['Len_bus_waveguide'] / 2 + param2['width_tether_bus_waveguide'] / 2+2*param['aper_cav'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius'] + param2['wg'] - (outer_radius - inner_radius - param2['ww']) / 2),
									(circle_center_x + param2['Len_bus_waveguide'] / 2 + param2['width_tether_bus_waveguide'] / 2+2*param['aper_cav'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius'] - param2['wg'] + (outer_radius - inner_radius - param2['ww']) / 2),
									(circle_center_x + param2['Len_bus_waveguide'] / 2 - param2['width_tether_bus_waveguide'] / 2+2*param['aper_cav'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius'] - param2['wg'] + (outer_radius - inner_radius - param2['ww']) / 2),
								]))

							elif param["2-axes"]==True:
								C = 1
								#tether for waveguide
								cell.add(gdspy.Polygon(1, [
									(circle_center_x -param2['width_tether_bus_waveguide']/2+param['y_spacing_between_wabguides']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + outer_radius - 2 * s1 * param['bus_bend_radius']+param2['wg']-(outer_radius-inner_radius-param2['ww'])/2),
									(circle_center_x +param2['width_tether_bus_waveguide']/2+param['y_spacing_between_wabguides']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + outer_radius - 2 * s1 * param['bus_bend_radius']+param2['wg']-(outer_radius-inner_radius-param2['ww'])/2),
									(circle_center_x +param2['width_tether_bus_waveguide']/2+param['y_spacing_between_wabguides']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + inner_radius - 2 * s1 * param['bus_bend_radius']-param2['wg']+(outer_radius-inner_radius-param2['ww'])/2),
									(circle_center_x -param2['width_tether_bus_waveguide']/2+param['y_spacing_between_wabguides']/2+2*param['aper_cav'], circle_center_y + S * (4* iM * param['bus_bend_radius'] + 2* iM * param['y_spacing_between_wabguides']+ C*(param['Len_bus_waveguide_vertical']+2*param['y_spacing_between_wabguides'])) + inner_radius - 2 * s1 * param['bus_bend_radius']-param2['wg']+(outer_radius-inner_radius-param2['ww'])/2),
								]))


					if iM < int(param['waveguide_with_end_mirror']) / 2 and iM>0 and param['waveguide_with_end_mirror'] % 2 == 0 :
						# rest of taper waveguides on the left, even
						for s1 in range(0, 2):
							cell.add(gdspy.Polygon(1, [
								(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x+ param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x+ param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
							]))

						# rest of taper waveguides on the right
						for s1 in range(0, 2):
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
							]))

					elif iM < int(param['waveguide_with_end_mirror']) / 2 and iM != int(param['waveguide_with_end_mirror']) / 2 and iM>0 and param['waveguide_with_end_mirror'] % 2 == 1 :
						# rest of taper waveguides on the left, odd
						for s1 in range(0, 2):
							cell.add(gdspy.Polygon(1, [
								(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
							]))

							if iM == int(param['waveguide_with_end_mirror']) / 2-1 :

								cell.add(gdspy.Polygon(1, [
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
								]))
								cell.add(gdspy.Polygon(1, [
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
								]))

						# rest of taper waveguides on the right
						for s1 in range(0, 2):
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
							]))

							if iM == int(param['waveguide_with_end_mirror']) / 2 - 1:

								cell.add(gdspy.Polygon(1, [
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
								]))

								cell.add(gdspy.Polygon(1, [
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius'] - 2 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius'] - 2 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius'] - 2 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius'] - 2 * param['bus_bend_radius']),
								]))
					if iB==1:
						if iM == int(param['waveguide_with_end_mirror']) / 2 and param['waveguide_with_end_mirror'] % 2 == 1 and S==-1:
							# last for odd number of waveguides
								# right
								cell.add(gdspy.Polygon(1, [
									(circle_center_x+ param2['Len_bus_waveguide']-param2['bus_taper_len']+4*param['aper_mir'], circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 *  param['bus_bend_radius']),
									(circle_center_x+ param2['Len_bus_waveguide']-param2['bus_taper_len']+4*param['aper_mir'], circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 *  param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']+4*param['aper_mir'], circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + outer_radius650 - 2 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']+4*param['aper_mir'], circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + inner_radius650 - 2 * param['bus_bend_radius']),
								]))
								# left
								cell.add(gdspy.Polygon(1, [
									(circle_center_x +param2['bus_taper_len'], circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * param['bus_bend_radius']),
									(circle_center_x +param2['bus_taper_len'], circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * param['bus_bend_radius']),
									(circle_center_x , circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + inner_radius - 2 * param['bus_bend_radius']),
									(circle_center_x , circle_center_y + S * (4 * (iM) * param['bus_bend_radius'] + 2 * (iM) * param['y_spacing_between_wabguides']) + outer_radius - 2 * param['bus_bend_radius']) ,
								]))

						if iM == int(param['waveguide_with_end_mirror']) / 2 and param['waveguide_with_end_mirror'] % 2 == 1 and S==1: # for odd number of waveguides
								# right
								cell.add(gdspy.Polygon(1, [
									(circle_center_x+ param2['Len_bus_waveguide']-param2['bus_taper_len']+4*param['aper_mir'], circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 *  param['bus_bend_radius']),
									(circle_center_x+ param2['Len_bus_waveguide']-param2['bus_taper_len']+4*param['aper_mir'], circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 *  param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']+4*param['aper_mir'], circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + outer_radius650 - 2 * param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide']+4*param['aper_mir'], circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + inner_radius650 - 2 * param['bus_bend_radius']),
								]))
								# left
								cell.add(gdspy.Polygon(1, [
									(circle_center_x +param2['bus_taper_len'], circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * param['bus_bend_radius']),
									(circle_center_x +param2['bus_taper_len'], circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * param['bus_bend_radius']),
									(circle_center_x , circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + inner_radius - 2 * param['bus_bend_radius']),
									(circle_center_x , circle_center_y + S * (4 * (iM+0.5) * param['bus_bend_radius'] + 2 * (iM+0.5) * param['y_spacing_between_wabguides']) + outer_radius - 2 * param['bus_bend_radius']) ,
								]))

					elif param2["2-axes"]==True:
						#taper the half circle to waveguide after songtao's large taper
						cell.add(gdspy.Polygon(1, [
							(circle_center_x- (inner_radius+outer_radius)/2- param2['ww'] / 2, circle_center_y+ S*(param['y_spacing_between_wabguides']+ param2['bus_taper_len'])),
							(circle_center_x- (inner_radius+outer_radius)/2+ param2['ww'] / 2, circle_center_y+ S*(param['y_spacing_between_wabguides']+ param2['bus_taper_len'])),
							(circle_center_x- inner_radius, circle_center_y+ S*param['y_spacing_between_wabguides']),
							(circle_center_x- outer_radius, circle_center_y+ S*param['y_spacing_between_wabguides']),
						]))

						cell.add(gdspy.Polygon(1, [
							(circle_center_x- (inner_radius+outer_radius)/2- param2['ww'] / 2, circle_center_y+ S*(param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']-param["beam_width"])),
							(circle_center_x- (inner_radius+outer_radius)/2+ param2['ww'] / 2, circle_center_y+ S*(param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']-param["beam_width"])),
							(circle_center_x- (inner_radius+outer_radius)/2+param["beam_width"]/2, circle_center_y+ S*(param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param2['bus_taper_len']-param["beam_width"])),
							(circle_center_x- (inner_radius+outer_radius)/2-param["beam_width"]/2, circle_center_y+ S*(param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param2['bus_taper_len']-param["beam_width"])),
						]))


					if param['waveguide_with_end_mirror'] != 1:
						#taper at first waveguide
						if iM==0:
							if S==1:
								s1=0
							elif S==-1:
								s1=1

							cell.add(gdspy.Polygon(1, [
								(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
							]))

							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']),
							]))

							if param['waveguide_with_end_mirror']==3:
								cell.add(gdspy.Polygon(1, [
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
								]))

								cell.add(gdspy.Polygon(1, [
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']+2*param['bus_bend_radius']),
								]))

								cell.add(gdspy.Polygon(1, [
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x, circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x + param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
								]))

								cell.add(gdspy.Polygon(1, [
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
									(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * s1 * param['bus_bend_radius']-2*param['bus_bend_radius']),
								]))

					# taper at the last bus waveguide at the top
					if iM == int(param['waveguide_with_end_mirror']) / 2 and S==1 :
						if param['waveguide_with_end_mirror'] % 2 == 0: # for even number of waveguides
							# left
							cell.add(gdspy.Polygon(1, [
								(circle_center_x-4*param['aper_mir'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius650 - 2 * param['bus_bend_radius']),
								(circle_center_x-4*param['aper_mir'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius650 - 2 * param['bus_bend_radius']),
								(circle_center_x + param2['bus_taper_len']-4*param['aper_mir'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 *  param['bus_bend_radius']),
								(circle_center_x + param2['bus_taper_len']-4*param['aper_mir'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 *  param['bus_bend_radius']),
							]))
							# right
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius - 2 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * param['bus_bend_radius']),
								(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * param['bus_bend_radius']),
							]))

					# taper at the last bus waveguide in bottom, even
					if iM == int(param['waveguide_with_end_mirror']) / 2 and param['waveguide_with_end_mirror'] % 2 == 0 and S==-1 and iM>0 :
						# left
						cell.add(gdspy.Polygon(1, [
							(circle_center_x-4*param['aper_mir'], circle_center_y + S * (4 * (iM - 1) * param['bus_bend_radius'] + 2 * (iM - 1) * param['y_spacing_between_wabguides']) + outer_radius650 - 2 * 2 * param['bus_bend_radius']),
							(circle_center_x-4*param['aper_mir'], circle_center_y + S * (4 * (iM - 1) * param['bus_bend_radius'] + 2 * (iM - 1) * param['y_spacing_between_wabguides']) + inner_radius650 - 2 * 2 * param['bus_bend_radius']),
							(circle_center_x + param2['bus_taper_len']-4*param['aper_mir'], circle_center_y + S * (4 * (iM - 1) * param['bus_bend_radius'] + 2 * (iM - 1) * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * 2 * param['bus_bend_radius']),
							(circle_center_x + param2['bus_taper_len']-4*param['aper_mir'], circle_center_y + S * (4 * (iM - 1) * param['bus_bend_radius'] + 2 * (iM - 1) * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * 2 * param['bus_bend_radius']),
						]))
						# right
						cell.add(gdspy.Polygon(1, [
							(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * (iM-1) * param['bus_bend_radius'] + 2 * (iM-1) * param['y_spacing_between_wabguides']) + inner_radius - 2 * 2 * param['bus_bend_radius']),
							(circle_center_x + param2['Len_bus_waveguide'], circle_center_y + S * (4 * (iM-1) * param['bus_bend_radius'] + 2 * (iM-1) * param['y_spacing_between_wabguides']) + outer_radius - 2 * 2 * param['bus_bend_radius']),
							(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * (iM-1) * param['bus_bend_radius'] + 2 * (iM-1) * param['y_spacing_between_wabguides']) + outer_radius - (param['width_taper'] - param2['ww']) / 2 - 2 * 2 * param['bus_bend_radius']),
							(circle_center_x + param2['Len_bus_waveguide'] - param2['bus_taper_len'], circle_center_y + S * (4 * (iM-1) * param['bus_bend_radius'] + 2 * (iM-1) * param['y_spacing_between_wabguides']) + inner_radius + (param['width_taper'] - param2['ww']) / 2 - 2 * 2 * param['bus_bend_radius']),
						]))

				# add Bus waveguides and PhCs
				for iM in range(int(param['waveguide_with_end_mirror'])):
					#Bus waveguide
					if iB==1:
						cell.add(gdspy.Polygon(1, [
							(circle_center_x, circle_center_y+  S*(inner_radius350+ 2*iM*param['bus_bend_radius'] +(1+iM)*param['y_spacing_between_wabguides'])),
							(circle_center_x, circle_center_y+  S*(outer_radius350+ 2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides'])),
							(circle_center_x + param['Len_bus_waveguide'], circle_center_y+  S*(outer_radius350+ 2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides'])),
							(circle_center_x + param['Len_bus_waveguide'], circle_center_y+  S*(inner_radius350+ 2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides'])),
						]))
					elif param["2-axes"]==True and iB!=1:

						#veritical waveguide
						cell.add(gdspy.Polygon(1, [
							(circle_center_x-outer_radius350, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical'])),
							(circle_center_x-inner_radius350, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical'])),
							(circle_center_x-inner_radius350, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides'])),
							(circle_center_x-outer_radius350, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides'])),
						]))
						#phc waveguide paralle to the verical meander
						#left
						cell.add(gdspy.Polygon(1, [
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]-param['beam_width'] / 2- param['ww'] / 2.0-param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])),
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]-param['beam_width'] / 2- param['ww'] / 2.0+param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])),
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]-param['beam_width'] / 2- param['ww'] / 2.0+param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides'])),
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]-param['beam_width'] / 2- param['ww'] / 2.0-param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(1+iM)*param['y_spacing_between_wabguides'])),
						]))
						#right
						cell.add(gdspy.Polygon(1, [
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2, circle_center_y + S * (2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'] + param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])),
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 + param['beam_width'] / 2, circle_center_y + S * (2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'] + param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])),
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 + param['beam_width'] / 2, circle_center_y + S * (2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2, circle_center_y + S * (2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
						]))

						#thether support left vertical phc waveguide
						cell.add(gdspy.Polygon(1, [ #bottom
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]- param['ww'] / 2.0,                                circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])+param2['width_tether_phc']),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])+param2['width_tether_phc']),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]- param['ww'] / 2.0,                                 outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
						]))
						cell.add(gdspy.Polygon(1, [#top
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]- param['ww'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])-param2['width_tether_phc']),
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]- param['ww'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])-param2['width_tether_phc']),
						]))

						cell.add(gdspy.Polygon(1, [#middle
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]- param['ww'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)+param2['width_tether_phc']/2-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)+param2['width_tether_phc']/2-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)-param2['width_tether_phc']/2-2*param['aper_cav']),
							(circle_center_x-(outer_radius350+inner_radius350)/2-param2["wg"]- param['ww'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)-param2['width_tether_phc']/2-2*param['aper_cav']),
						]))

						#thether support right vertical phc waveguide
						cell.add(gdspy.Polygon(1, [#bottom
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])+param2['width_tether_phc']),
							(circle_center_x - param['bus_bend_radius'] / 2,                                                                                                   circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']+param['y_spacing_between_wabguides'])+param2['width_tether_phc']),
							(circle_center_x - param['bus_bend_radius'] / 2,                                                                                                    outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2,  outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
						]))
						cell.add(gdspy.Polygon(1, [#top
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])),
							(circle_center_x - param['bus_bend_radius'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])),
							(circle_center_x - param['bus_bend_radius'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])-param2['width_tether_phc']),
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical'])-param2['width_tether_phc']),
						]))
						cell.add(gdspy.Polygon(1, [#middle
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)+param2['width_tether_phc']/2 -2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)+param2['width_tether_phc']/2-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)-param2['width_tether_phc']/2-2*param['aper_cav']),
							(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 - param['beam_width'] / 2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2)-param2['width_tether_phc']/2-2*param['aper_cav']),
						]))

						#tapering for tethering vertical waveguide
						# bottom
						cell.add(gdspy.Polygon(1, [
							(circle_center_x - param['bus_bend_radius']- param['ww'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2+ param2['length_wider_bus_wavguide_part'] / 2+param2['bus_taper_len'])-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']+ param['ww'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2+ param2['length_wider_bus_wavguide_part'] / 2+param2['bus_taper_len'])-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']+ param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2+ param2['length_wider_bus_wavguide_part'] / 2)-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']- param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2+ param2['length_wider_bus_wavguide_part'] / 2)-2*param['aper_cav']),
						]))

						# top
						cell.add(gdspy.Polygon(1, [
							(circle_center_x - param['bus_bend_radius']- param['ww'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2- 1.5*param2['length_wider_bus_wavguide_part'] +param2['length_wider_bus_wavguide_part']-param2['bus_taper_len'])-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']+ param['ww'] / 2.0,  circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2- 1.5*param2['length_wider_bus_wavguide_part'] +param2['length_wider_bus_wavguide_part']-param2['bus_taper_len'])-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']+ param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2- 1.5*param2['length_wider_bus_wavguide_part'] +param2['length_wider_bus_wavguide_part'])-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']- param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2- 1.5*param2['length_wider_bus_wavguide_part'] +param2['length_wider_bus_wavguide_part'])-2*param['aper_cav']),
						]))

						# center
						cell.add(gdspy.Polygon(1, [
							(circle_center_x - param['bus_bend_radius']  - param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2+ param2['length_wider_bus_wavguide_part'] / 2)-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']  + param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2+ param2['length_wider_bus_wavguide_part'] / 2)-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']  + param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2- param2['length_wider_bus_wavguide_part'] / 2)-2*param['aper_cav']),
							(circle_center_x - param['bus_bend_radius']  - param['width_taper_middle'] / 2.0, circle_center_y+ S*(2*iM*param['bus_bend_radius']+(iM)*param['y_spacing_between_wabguides']+ param['Len_bus_waveguide_vertical']- param['Len_bus_waveguide_vertical']/2- param2['length_wider_bus_wavguide_part'] / 2)-2*param['aper_cav']),
						]))

					# add the extra waveguide for mirror at the end of the meander for even number of waveguide
					if (iM + 1) == param['waveguide_with_end_mirror'] and param['waveguide_with_end_mirror'] % 2 == 0:
						cell.add(gdspy.Polygon(1, [
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (inner_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (outer_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							(circle_center_x -4*param['aper_mir'], circle_center_y + S * (outer_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							(circle_center_x -4*param['aper_mir'], circle_center_y + S * (inner_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
						]))

						# add the extra waveguide for mirror at the end of the meander for odd number of waveguide
					if iB==1:
						if (iM + 1) == param['waveguide_with_end_mirror'] and param['waveguide_with_end_mirror'] % 2 == 1:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, circle_center_y + S * (inner_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, circle_center_y + S * (outer_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] +4*param['aper_mir'], circle_center_y + S * (outer_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] +4*param['aper_mir'], circle_center_y + S * (inner_radius650 + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))
					elif param2["2-axes"] == True:
						#vertical waveguide end mirror
						cell.add(gdspy.Polygon(1, [
							(circle_center_x-(outer_radius350+inner_radius350)/2-param['beam_width']/2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+param['Len_bus_waveguide_vertical']+(1+iM)*param['y_spacing_between_wabguides']+param2['bus_taper_len']-param["beam_width"])),
							(circle_center_x-(outer_radius350+inner_radius350)/2+param['beam_width']/2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+param['Len_bus_waveguide_vertical']+(1+iM)*param['y_spacing_between_wabguides']+param2['bus_taper_len']-param["beam_width"])),
							(circle_center_x-(outer_radius350+inner_radius350)/2+param['beam_width']/2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+param['Len_bus_waveguide_vertical']+(1+iM)*param['y_spacing_between_wabguides']+18 * param['aper_mir']+param2['bus_taper_len']-param["beam_width"])),
							(circle_center_x-(outer_radius350+inner_radius350)/2-param['beam_width']/2, circle_center_y+ S*(2*iM*param['bus_bend_radius']+param['Len_bus_waveguide_vertical']+(1+iM)*param['y_spacing_between_wabguides']+18 * param['aper_mir']+param2['bus_taper_len']-param["beam_width"])),
						]))

					#PhC waveguide top
					if iB==1:
						if iM % 2==0:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'], circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'], circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) *param['y_spacing_between_wabguides'])),
							]))

						elif iM==param['waveguide_with_end_mirror']-1:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) *param['y_spacing_between_wabguides'])),
							]))
							y_top=circle_center_y + (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])
							y_bottom = circle_center_y - (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])

						else:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x, circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) *param['y_spacing_between_wabguides'])),
								(circle_center_x, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, circle_center_y + S * (outer_radius350 + param['wg'] + param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) *param['y_spacing_between_wabguides'])),
							]))


						# PhC waveguide bottom
						if iM % 2==1:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (inner_radius350 - param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'], circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'], circle_center_y + S * (inner_radius350 - param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))

						else:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x, circle_center_y + S * (inner_radius350 - param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, circle_center_y + S * (inner_radius350 - param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))

					#support pad tethers connected to the left frame:
					if iM % 2 == 1:

						middle=( circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])  + circle_center_y + S * (inner_radius350 - param['wg'] + 2 * (iM-1) * param['bus_bend_radius'] + ( iM) * param['y_spacing_between_wabguides']))/2
						cell.add(gdspy.Polygon(1, [
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, middle + param2['width_support_tather']/2),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, middle - param2['width_support_tather']/2),
							(circle_center_x + param['Len_bus_waveguide']+param['bus_bend_radius']/2, middle - param2['width_support_tather']/2),
							(circle_center_x + param['Len_bus_waveguide']+param['bus_bend_radius']/2, middle + param2['width_support_tather']/2),
						]))

						#tether connecting bus waveguie to the left frame connecting the bent circle
						cell.add(gdspy.Polygon(1, [
							(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, middle + param2['width_tether_bus_waveguide'] / 2), #+ param['buffer_siliconpad_bendwaveguide']
							(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, middle - param2['width_tether_bus_waveguide'] / 2),
							(circle_center_x + param['Len_bus_waveguide']+param['bus_bend_radius']/2, middle - param2['width_tether_bus_waveguide'] / 2),
							(circle_center_x + param['Len_bus_waveguide']+param['bus_bend_radius']/2, middle + param2['width_tether_bus_waveguide'] / 2),
						]))

						if iB==1:
							#tethers connecting phc to the left frame at the end point connecting the bent circle
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param['Len_bus_waveguide']- param2['width_tether_phc'], circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide']- param2['width_tether_phc'], middle + S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide'], middle + S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide'], circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param['Len_bus_waveguide']- param2['width_tether_phc'], middle + S * (-param2['width_support_tather'] / 2 - inner_radius350 + param2['width_support_tather'] / 2 + param['wg'] + param['beam_width'])),
								(circle_center_x + param['Len_bus_waveguide']- param2['width_tether_phc'], middle - S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide'], middle - S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide'], middle + S*(-param2['width_support_tather'] / 2- inner_radius350 + param2['width_support_tather']/2+ param['wg']+ param['beam_width'])),
							]))

							#tethers connecting phc to the left frame at the middle point
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param['Len_bus_waveguide']/2- param2['width_tether_phc']/2, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param['Len_bus_waveguide']/2- param2['width_tether_phc']/2, middle + S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide']/2+param2['width_tether_phc']/2, middle + S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide']/2+param2['width_tether_phc']/2, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param['Len_bus_waveguide']/2- param2['width_tether_phc']/2, middle + S * (-param2['width_support_tather'] / 2 - inner_radius350 + param2['width_support_tather'] / 2 + param['wg'] + param['beam_width'])),
								(circle_center_x + param['Len_bus_waveguide']/2- param2['width_tether_phc']/2, middle - S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide']/2+param2['width_tether_phc']/2, middle - S*param2['width_support_tather'] / 2 ),
								(circle_center_x + param['Len_bus_waveguide']/2+param2['width_tether_phc']/2, middle + S*(-param2['width_support_tather'] / 2- inner_radius350 + param2['width_support_tather']/2+ param['wg']+ param['beam_width'])),
							]))

					#support pad tethers connected to the pad:
					if iM % 2 == 0:

						if iM==0 and S==1:
							middle = (circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides']) + circle_center_y + S * (inner_radius350 - param['wg'] + 2 * (iM - 1) * param['bus_bend_radius'] + (iM) * param['y_spacing_between_wabguides'])) / 2

							if iB==1:
								cell.add(gdspy.Polygon(1, [
									(circle_center_x - param['bus_bend_radius'] / 2, middle + param2['width_support_tather'] / 2),
									(circle_center_x - param['bus_bend_radius'] / 2, param['extra_silicon_pad_gap'] * 3.5 + hole_center_y),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, param['extra_silicon_pad_gap'] * 3.5 + hole_center_y),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, middle + param2['width_support_tather'] / 2),
								]))

							elif param["2-axes"]==True:
								cell.add(gdspy.Polygon(1, [
									(circle_center_x - param['bus_bend_radius'] / 2, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] - param2['spacing_phc-vertical_buffer'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'] + param['y_spacing_between_wabguides'] + param['Len_bus_waveguide_vertical'])),
									(circle_center_x - param['bus_bend_radius'] / 2, param['extra_silicon_pad_gap'] * 3.5 + hole_center_y),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, param['extra_silicon_pad_gap'] * 3.5 + hole_center_y),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, circle_center_y + S * (inner_radius350 - param['wg'] - param2['spacing_phc-vertical_buffer'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'] + param['y_spacing_between_wabguides'] + param['Len_bus_waveguide_vertical'])),
								]))

						elif iM==0 and S==-1:
							middle = (circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides']) + circle_center_y + S * (inner_radius350 - param['wg'] + 2 * (iM - 1) * param['bus_bend_radius'] + (iM) * param['y_spacing_between_wabguides'])) / 2

							if iB==1:
								cell.add(gdspy.Polygon(1, [
									(circle_center_x - param['bus_bend_radius'] / 2, param['extra_silicon_pad_gap'] * (-3.5) + hole_center_y),
									(circle_center_x - param['bus_bend_radius'] / 2, middle - param2['width_support_tather'] / 2),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, middle - param2['width_support_tather'] / 2),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, param['extra_silicon_pad_gap'] * (-3.5) + hole_center_y),
								]))
							elif param["2-axes"]==True:
								cell.add(gdspy.Polygon(1, [
									(circle_center_x - param['bus_bend_radius'] / 2, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
									(circle_center_x - param['bus_bend_radius'] / 2, param['extra_silicon_pad_gap'] * (-3.5) + hole_center_y),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, param['extra_silicon_pad_gap'] * (-3.5) + hole_center_y),
									(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide'] + outer_radius, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
								]))

						else:
							middle=( circle_center_y + S * (outer_radius350 + param['wg'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides']) + circle_center_y + S * (inner_radius350 - param['wg'] + 2 * (iM-1) * param['bus_bend_radius'] + (iM) * param['y_spacing_between_wabguides']))/2
							cell.add(gdspy.Polygon(1, [
								(circle_center_x-param['bus_bend_radius']/2, middle + param2['width_support_tather']/2),
								(circle_center_x-param['bus_bend_radius']/2, middle - param2['width_support_tather']/2),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, middle - param2['width_support_tather']/2),
								(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, middle + param2['width_support_tather']/2),
							]))

						if iB==1:
							#tether connecting bus waveguie to the right pad
							cell.add(gdspy.Polygon(1, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, middle + param2['width_tether_bus_waveguide'] / 2), #+ param['buffer_siliconpad_bendwaveguide']
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, middle - param2['width_tether_bus_waveguide'] / 2),
								(circle_center_x-param['bus_bend_radius']/2, middle - param2['width_tether_bus_waveguide'] / 2),
								(circle_center_x-param['bus_bend_radius']/2, middle + param2['width_tether_bus_waveguide'] / 2),
							]))
					 	else:
							#tether connecting bus waveguie to the right pad
							cell.add(gdspy.Polygon(1, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, middle + param2['width_tether_bus_waveguide'] / 2-param2['Len_bus_waveguide_vertical']/2-2*param['aper_cav']), #+ param['buffer_siliconpad_bendwaveguide']
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, middle - param2['width_tether_bus_waveguide'] / 2-param2['Len_bus_waveguide_vertical']/2-2*param['aper_cav']),
								(circle_center_x-param['bus_bend_radius']/2, middle - param2['width_tether_bus_waveguide'] / 2-param2['Len_bus_waveguide_vertical']/2-2*param['aper_cav']),
								(circle_center_x-param['bus_bend_radius']/2, middle + param2['width_tether_bus_waveguide'] / 2-param2['Len_bus_waveguide_vertical']/2-2*param['aper_cav']),
							]))


						if iM>0:
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param2['width_tether_phc'], middle + S * (-param2['width_support_tather'] / 2 - inner_radius350 + param2['width_support_tather'] / 2 + param['wg'] + param['beam_width'])),
								(circle_center_x + param2['width_tether_phc'], middle - S * param2['width_support_tather'] / 2),
								(circle_center_x, middle - S * param2['width_support_tather'] / 2),
								(circle_center_x, middle + S * (-param2['width_support_tather'] / 2 - inner_radius350 + param2['width_support_tather'] / 2 + param['wg'] + param['beam_width'])),
							]))

						if iB==1:
							# tethers connecting phc to the pad at the end point
							cell.add(gdspy.Polygon(1, [
								(circle_center_x + param2['width_tether_phc'], circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x + param2['width_tether_phc'], middle + S * param2['width_support_tather'] / 2),
								(circle_center_x, middle + S * param2['width_support_tather'] / 2),
								(circle_center_x, circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))
							# tethers connecting phc to the pad at the middle point
							cell.add(gdspy.Polygon(1, [
								(circle_center_x- param2['width_tether_phc']/2+ param['Len_bus_waveguide']/2+2*param['aper_cav'], circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
								(circle_center_x- param2['width_tether_phc']/2+ param['Len_bus_waveguide']/2+2*param['aper_cav'], middle + S*param2['width_support_tather'] / 2 ),
								(circle_center_x+ param2['width_tether_phc']/2+ param['Len_bus_waveguide']/2+2*param['aper_cav'], middle + S*param2['width_support_tather'] / 2 ),
								(circle_center_x+ param2['width_tether_phc']/2+ param['Len_bus_waveguide']/2+2*param['aper_cav'], circle_center_y + S * (inner_radius350 - param['wg'] - param['beam_width'] + 2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])),
							]))
							if iM>0:
								cell.add(gdspy.Polygon(1, [
									(circle_center_x - param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], middle + S * (-param2['width_support_tather'] / 2 - inner_radius350 + param2['width_support_tather'] / 2 + param['wg'] + param['beam_width'])),
									(circle_center_x - param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], middle - S * param2['width_support_tather'] / 2),
									(circle_center_x + param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], middle - S * param2['width_support_tather'] / 2),
									(circle_center_x + param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], middle + S * (-param2['width_support_tather'] / 2 - inner_radius350 + param2['width_support_tather'] / 2 + param['wg'] + param['beam_width'])),
								]))


				# add cuts inside top pad for HF:
				if iB==1:
					aa=14
					bb=int((param['waveguide_with_end_mirror']) / 2) + 1
				else:
					aa=19
					bb=6

				for iM in range(bb+1):

					if param2["2-axes"] == True and iB!=1:
						for z in range(0, aa):
							if S == 1:
								s1 = 0
							if S == -1:
								s1 = 1
							cell.add(gdspy.Round(3,
								 (outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2 - param2['width_silicon_between_HF_cuts'] / 2 - z * (param2['width_silicon_between_HF_cuts'] + 2 * param2['width_HF_cuts']),
								  circle_center_y + S * (4 * iM * param['bus_bend_radius']/3 ) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								 param2['width_HF_cuts']))
							cell.add(gdspy.Polygon(3, [  # rectangular cuts on the far right
								(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2 + param2['width_silicon_between_HF_cuts'] / 2,
								 circle_center_y + S * (4 * iM * param['bus_bend_radius']/3 ) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2 + param2['width_silicon_between_HF_cuts'] / 2,
								 circle_center_y + S * (4 * iM * param['bus_bend_radius']/3 ) + outer_radius - 2 * s1 * param['bus_bend_radius'] + S * param2['width_HF_cuts']),
								(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'] - param2['width_silicon_between_HF_cuts'] / 2,
								 circle_center_y + S * (4 * iM * param['bus_bend_radius']/3 ) + outer_radius - 2 * s1 * param['bus_bend_radius'] + S * param2['width_HF_cuts']),
								(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'] - param2['width_silicon_between_HF_cuts'] / 2,
								 circle_center_y + S * (4 * iM * param['bus_bend_radius']/3 ) + outer_radius - 2 * s1 * param['bus_bend_radius']),
							]))

					if iB==1:
						if iM == 0:  # the circles closest to the big taper
							for z in range(0, aa-1):
								if S == 1:
									s1 = 0
								if S == -1:
									s1 = 1
								cell.add(gdspy.Round(3,
													 (outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2 - param2['width_silicon_between_HF_cuts'] / 2 - z * (param2['width_silicon_between_HF_cuts'] + 2 * param2['width_HF_cuts']),
													  circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
													 param2['width_HF_cuts']))
								cell.add(gdspy.Polygon(3, [  # rectangular cuts on the far right
									(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2 + param2['width_silicon_between_HF_cuts'] / 2,
									 circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
									(outerbox_x_max - param['grating_pad_offset'] - param['grating_pad_length'] + param['support_tether_width'] / 2 + param2['width_silicon_between_HF_cuts'] / 2,
									 circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius'] + S * param2['width_HF_cuts']),
									(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'] - param2['width_silicon_between_HF_cuts'] / 2,
									 circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius'] + S * param2['width_HF_cuts']),
									(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'] - param2['width_silicon_between_HF_cuts'] / 2,
									 circle_center_y + S * (4 * iM * param['bus_bend_radius'] + 2 * iM * param['y_spacing_between_wabguides']) + outer_radius - 2 * s1 * param['bus_bend_radius']),
								]))

			#add block and tethers for the topmost and bottom most PhCs:
			if param['meander'] is True:
				for i in range(0, 2):
					if i == 0:  # for bottom one
						y = outerbox_y_min
						s2 = -1.0
					if i == 1:  # for top one
						y = outerbox_y_max
						s2 = 1.0
					# start with the tether at the center of the PhCs
				if iB==1:
					cell.add(gdspy.Polygon(1, [
						(circle_center_x - param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], y - s2 * param['box_buffer'] - s2 * param['beam_width']/2.003 - s2 * param['supporting_bar_width']),
						(circle_center_x - param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], y - s2 * param['box_buffer'] - s2 * param['beam_width']/2.003 - s2 * param['supporting_bar_width'] - s2 * param2['spacing_phc-vertical_buffer']-30),
						(circle_center_x + param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], y - s2 * param['box_buffer'] - s2 * param['beam_width']/2.003 - s2 * param['supporting_bar_width'] - s2 * param2['spacing_phc-vertical_buffer']-30),
						(circle_center_x + param2['width_tether_phc'] / 2 + param['Len_bus_waveguide'] / 2+2*param['aper_cav'], y - s2 * param['box_buffer'] - s2 * param['beam_width']/2.003 - s2 * param['supporting_bar_width'])
					]))

					# now add a block at the right and left ends of the PhCs

				if iB==1:
					cell.add(gdspy.Polygon(1, [  # right above phc
						(circle_center_x + param['Len_bus_waveguide'], y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * (param['supporting_bar_width']+2*param["wg"]+param["beam_width"]+param2['spacing_phc-vertical_buffer']+param["ww"])),
						(circle_center_x + param['Len_bus_waveguide'], y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * (param['supporting_bar_width']+2*param["wg"]+param["beam_width"]+param2['spacing_phc-vertical_buffer']+param["ww"]) - param['bus_bend_radius']),
						(circle_center_x + param['Len_bus_waveguide'] + param['bus_bend_radius'] + param2['buffer_siliconpad_bendwaveguide'] + param2['width_taper'] / 2, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * (param['supporting_bar_width']+2*param["wg"]+param["beam_width"]+param2['spacing_phc-vertical_buffer']+param["ww"]) - param['bus_bend_radius']),
						(circle_center_x + param['Len_bus_waveguide'] + param['bus_bend_radius'] + param2['buffer_siliconpad_bendwaveguide'] + param2['width_taper'] / 2, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * (param['supporting_bar_width']+2*param["wg"]+param["beam_width"]+param2['spacing_phc-vertical_buffer']+param["ww"])),
						]))
					cell.add(gdspy.Polygon(1, [  # right
						(circle_center_x + param['Len_bus_waveguide'], y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width']),
						(circle_center_x + param['Len_bus_waveguide'], y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width'] - s2 * param2['spacing_phc-vertical_buffer'] - s2 * param['beam_width']-30),
						(circle_center_x + param['Len_bus_waveguide'] + param['bus_bend_radius'] + param2['buffer_siliconpad_bendwaveguide'] + param2['width_taper'] / 2, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width'] - s2 * param['beam_width'] - s2 * param2['spacing_phc-vertical_buffer']-30),
						(circle_center_x + param['Len_bus_waveguide'] + param['bus_bend_radius'] + param2['buffer_siliconpad_bendwaveguide'] + param2['width_taper'] / 2, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width']),
					]))

					cell.add(gdspy.Polygon(1, [  # left
						(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width']),
						(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width'] - s2 * param2['spacing_phc-vertical_buffer']-30),
						(circle_center_x, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width'] - s2 * param2['spacing_phc-vertical_buffer']-30),
						(circle_center_x, y - s2 * param['box_buffer'] - s2 * param['beam_width'] / 2.003 - s2 * param['supporting_bar_width'])
					]))


			# modified beam widening taper (GRATINGS)
			# SC incoporate the nonlinear waveguide taper design, adding one more variable "taper_nonlinear_order"
			# SC this part now also include the extra silicon pad to decrease the exposure time
			if param['straight_q'] is True and param['grating_pad_q'] is True and param['end_width'] > param['beam_width'] and iB > -1:

				xposlist = numpy.linspace(0, param['grating_taper_length'], 195)

				if param['meander'] is True:
					pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x,
									 hole_center_y + param['end_width'] / 2.0 + (param['width_taper'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order']) for x in xposlist]

					pts_up_1 = numpy.append([
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y + param['end_width'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y + param['end_width'] / 2.0)],
						# (beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - param['grating_taper_length'], hole_center_y + param['beam_width']/2.0)],
						pts_up_taper, axis=0)

					pts_up = numpy.append(pts_up_1,[
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'],hole_center_y + param['width_taper'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y + param['end_width'] / 2.0)], axis = 0)

					pts_down_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x,
									   hole_center_y - param['end_width'] / 2.0 - (param['width_taper'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param2['taper_nonlinear_order']) for x in xposlist]

					pts_down_1 = numpy.append([
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y - param['end_width'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y - param['end_width'] / 2.0)],
						# (beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - param['grating_taper_length'], hole_center_y - param['beam_width']/2.0),
						pts_down_taper, axis=0)

					pts_down = numpy.append(pts_down_1, [
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y - param['width_taper'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y - param['end_width'] / 2.0)], axis=0)

				else:
					xposlist = numpy.linspace(0, param['grating_taper_length'], 195)
					pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x, hole_center_y + param['end_width'] / 2.0 + (param['beam_width'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order']) for x in xposlist]

					pts_up_1 = numpy.append([
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y + param['end_width'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y + param['end_width'] / 2.0)],
						# (beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - param['grating_taper_length'], hole_center_y + param['beam_width']/2.0)],
						pts_up_taper, axis=0)

					pts_up = numpy.append(pts_up_1, [
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y + param['beam_width'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y + param['end_width'] / 2.0)], axis=0)

					pts_down_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x,
									   hole_center_y - param['end_width'] / 2.0 - (param['beam_width'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param2['taper_nonlinear_order']) for x in xposlist]

					pts_down_1 = numpy.append([
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y - param['end_width'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y - param['end_width'] / 2.0)],
						# (beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - param['grating_taper_length'], hole_center_y - param['beam_width']/2.0),
						pts_down_taper, axis=0)

					pts_down = numpy.append(pts_down_1, [
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y - param['beam_width'] / 2.0),
						(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'], hole_center_y - param['end_width'] / 2.0)], axis=0)

				cell.add(gdspy.Polygon(1, pts_up))
				cell.add(gdspy.Polygon(1,pts_down))

				# extra silicon pad
				if param['extra_silicon_pad_q'] is True and iB==param['num_beams']/2-1:
					# first half of middle the silicon pad
					xposlist = numpy.linspace(param['grating_taper_length'], 0, 195)

					if param['meander'] is True:
						pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x,
										 param['extra_silicon_pad_gap'] + hole_center_y + param['end_width'] / 2.0 + (param['width_taper'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order']) for x in xposlist]

						pts_up_1 = numpy.append(pts_up_taper, [
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0)],
												axis=0)

						pts_up = numpy.append(pts_up_1, [
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, hole_center_y + param['beam_width'] / 2.0 + param['extra_silicon_pad_gap'])], axis=0)

						cell.add(gdspy.Polygon(1, pts_up))

						# SC adding center cut inside silicon pad for HF undercut
						if iB < (param['num_beams'] - 1):
							cell.add(gdspy.Polygon(3, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 - param['centercut_silicon_pad'] / 2.0),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 + param['centercut_silicon_pad'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - param['grating_taper_length'] / 2.0, hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 + param['centercut_silicon_pad'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param[
									'grating_pad_length'] - param['grating_taper_length'] / 2.0, hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 - param['centercut_silicon_pad'] / 2.0),
							]))

					#bottom silicon pad
					xposlist = numpy.linspace(param['grating_taper_length']- 2* param['bus_bend_radius']-param['Len_bus_waveguide'] -param['buffer_siliconpad_bendwaveguide'], 0, 195)
					if param['meander'] is True:
						pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x,
										 - param['extra_silicon_pad_gap'] + hole_center_y - param['end_width'] / 2.0 - (param['width_taper'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order'])
										for x in xposlist]

						pts_up_1 = numpy.append(pts_up_taper, [
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
							(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width'])
						], axis=0)
						pts_up = numpy.append(pts_up_1, [
							(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius,
							 param['extra_silicon_pad_gap']*(-3.5) + hole_center_y)], axis=0)

						cell.add(gdspy.Polygon(1, pts_up))

				#second half of middle silicon pad
				if param['extra_silicon_pad_q'] is True and iB==param['num_beams']/2:
					xposlist = numpy.linspace(param['grating_taper_length'], 0, 195)
					if param['meander'] is True:
						pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x,
										 - param['extra_silicon_pad_gap'] + hole_center_y - param['end_width'] / 2.0 - (param['width_taper'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order']) for x in xposlist]
						pts_up_1 = numpy.append(pts_up_taper, [
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y - (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0),
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, hole_center_y - (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0)], axis=0)
						pts_up = numpy.append(pts_up_1, [
							(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param2['width_taper'] / 2, hole_center_y - param['beam_width'] / 2.0 - param['extra_silicon_pad_gap'])], axis=0)
						cell.add(gdspy.Polygon(1, pts_up))

						#top silicon pad
						xposlist1 = numpy.linspace(param['grating_taper_length']- 2* param['bus_bend_radius']-param['Len_bus_waveguide'] -param['buffer_siliconpad_bendwaveguide'], 0, 195)

						pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x,
										 param['extra_silicon_pad_gap'] + hole_center_y + param['end_width'] / 2.0 + (param['width_taper'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order']) for x in xposlist1]

						pts_up_1 = numpy.append(pts_up_taper, [
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
							(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width'])
						], axis=0)
						pts_up = numpy.append(pts_up_1, [
							(circle_center_x + param['Len_bus_waveguide'] + param['buffer_siliconpad_bendwaveguide']+outer_radius, param['extra_silicon_pad_gap']*3.5+ hole_center_y)], axis=0)

						cell.add(gdspy.Polygon(1, pts_up))

				if param['meander'] is False:
					if param['extra_silicon_pad_q'] is True:
						# upper part of the silicon pad
						xposlist = numpy.linspace(param['grating_taper_length'], 0, 195)

						pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x, param['extra_silicon_pad_gap']+ hole_center_y + param['end_width'] / 2.0 + (param['beam_width'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order']) for x in xposlist]

						if iB == (param['num_beams'] - 1):  # for the top most extra silicon pad, make the connection to right tether also 1um wide
							pts_up_1 = numpy.append(pts_up_taper, [
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 + (param['grating_pad_spacing'] - 2.0 * param['extra_silicon_pad_gap']) / 2.0),
								(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'], hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 + (param['grating_pad_spacing'] - 2.0 * param['extra_silicon_pad_gap']) / 2.0)],
													axis=0)
						else:
							pts_up_1 = numpy.append(pts_up_taper, [
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0),
								(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'], hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0)],
													axis=0)

						pts_up = numpy.append(pts_up_1, [
							(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'],
							 hole_center_y + param['beam_width'] / 2.0 + param['extra_silicon_pad_gap'])], axis=0)

						cell.add(gdspy.Polygon(1, pts_up))

							# SC adding center cut inside silicon pad for HF undercut
						if iB < (param['num_beams'] - 1):
							cell.add(gdspy.Polygon(3, [
								(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'], hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 - param['centercut_silicon_pad'] / 2.0),
								(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'], hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 + param['centercut_silicon_pad'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - param['grating_taper_length'] / 2.0, hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 + param['centercut_silicon_pad'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - param['grating_taper_length'] / 2.0, hole_center_y + (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 - param['centercut_silicon_pad'] / 2.0),
								]))

							# lower part of the silicon pad
						xposlist = numpy.linspace(param['grating_taper_length'], 0, 195)

						pts_up_taper = [(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] - x, - param['extra_silicon_pad_gap']+ hole_center_y - param['end_width'] / 2.0 - (param['beam_width'] - param['end_width']) / 2.0 * (x / param['grating_taper_length']) ** param['taper_nonlinear_order']) for x in xposlist]
						if iB == 0:  # for the top most extra silicon pad, make the connection to right tether also 1um wide
							pts_up_1 = numpy.append(pts_up_taper, [
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y - (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 - (param['grating_pad_spacing'] - 2.0 * param['extra_silicon_pad_gap']) / 2.0),
								(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'], hole_center_y - (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0 - (param['grating_pad_spacing'] - 2.0 * param['extra_silicon_pad_gap']) / 2.0)],
													axis=0)
						else:
							pts_up_1 = numpy.append(pts_up_taper, [
								(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y - (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0),
								(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'], hole_center_y - (param['grating_pad_width'] + param['grating_pad_spacing']) / 2.0)],
													axis=0)

						pts_up = numpy.append(pts_up_1, [
							(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2.0 * param['supporting_bar_width'],
							 hole_center_y - param['beam_width'] / 2.0 - param['extra_silicon_pad_gap'])], axis=0)

						cell.add(gdspy.Polygon(1, pts_up))

			if param['meander'] is False:
				# write the notches pinch points: left side
				cell.add(gdspy.Polygon(3, [
					(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] - param['notch_end_width'] / 2.0, hole_center_y - param['beam_width'] / 2.0),
					(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'], hole_center_y - param['beam_width'] / 2.0 + param['notch_end_depth']),
					(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] + param['notch_end_width'] / 2.0, hole_center_y - param['beam_width'] / 2.0)]))

				cell.add(gdspy.Polygon(3, [
					(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] - param['notch_end_width'] / 2.0, hole_center_y + param['beam_width'] / 2.0),
					(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'], hole_center_y + param['beam_width'] / 2.0 - param['notch_end_depth']),
					(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] + param['notch_end_width'] / 2.0, hole_center_y + param['beam_width'] / 2.0)]))
				# SC write the notches: left side, adding two more notches evenly in "beam_spacing"

				if param['extra_left_notches'] is True:
					cell.add(gdspy.Polygon(3, [
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] - param['notch_end_width'] / 2.0, hole_center_y - param['beam_width'] / 2.0 + param['beam_spacing'] / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'], hole_center_y - param['beam_width'] / 2.0 + param['notch_end_depth'] + param['beam_spacing'] / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] + param['notch_end_width'] / 2.0, hole_center_y - param['beam_width'] / 2.0 + param['beam_spacing'] / 3.0)]))

					cell.add(gdspy.Polygon(3, [
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] - param['notch_end_width'] / 2.0, hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'], hole_center_y + param['beam_width'] / 2.0 - param['notch_end_depth'] + param['beam_spacing'] / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] + param['notch_end_width'] / 2.0, hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] / 3.0)]))

					cell.add(gdspy.Polygon(3, [
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] - param['notch_end_width'] / 2.0, hole_center_y - param['beam_width'] / 2.0 + param['beam_spacing'] * 2.0 / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'], hole_center_y - param['beam_width'] / 2.0 + param['notch_end_depth'] + param['beam_spacing'] * 2.0 / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] + param['notch_end_width'] / 2.0, hole_center_y - param['beam_width'] / 2.0 + param['beam_spacing'] * 2.0 / 3.0)]))

					cell.add(gdspy.Polygon(3, [
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] - param['notch_end_width'] / 2.0, hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] * 2.0 / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'], hole_center_y + param['beam_width'] / 2.0 - param['notch_end_depth'] + param['beam_spacing'] * 2.0 / 3.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['notch_end_offset'] + param['notch_end_width'] / 2.0, hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] * 2.0 / 3.0)]))

			if param['straight_q'] is True and param['grating_pad_q'] is False and iB > -1:
				# right side, has been modified to allow different notch parameters for left and right sides
				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] - param['notch_taper_end_width']/2.0, hole_center_y - param['end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'], hole_center_y - param['end_width']/2.0 + param['notch_taper_end_depth']),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] + param['notch_taper_end_width']/2.0, hole_center_y - param['end_width']/2.0)] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] - param['notch_taper_end_width']/2.0, hole_center_y + param['end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'], hole_center_y + param['end_width']/2.0 - param['notch_taper_end_depth']),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] + param['notch_taper_end_width']/2.0, hole_center_y + param['end_width']/2.0)] ) )

			if param['straight_q'] is True and param['grating_pad_q'] is True and iB > -1:
				# right side, has been modified to allow different notch parameters for left and right sides
				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] - param['notch_taper_end_width']/2.0, hole_center_y - param['beam_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'], hole_center_y - param['beam_width']/2.0 + param['notch_taper_end_depth']),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] + param['notch_taper_end_width']/2.0, hole_center_y - param['beam_width']/2.0)] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] - param['notch_taper_end_width']/2.0, hole_center_y + param['beam_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'], hole_center_y + param['beam_width']/2.0 - param['notch_taper_end_depth']),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] + param['notch_taper_end_width']/2.0, hole_center_y + param['beam_width']/2.0)] ) )

				# notches in extra top linker beams for grating pad
				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 - param['notch_taper_end_width']/2.0, hole_center_y + param['grating_pad_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 + param['notch_taper_end_width']/2.0, hole_center_y + param['grating_pad_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0, hole_center_y + param['grating_pad_width']/2.0 - param['notch_taper_end_depth'])] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 - param['notch_taper_end_width']/2.0, hole_center_y + param['grating_pad_width']/2.0 - param['beam_width']),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 + param['notch_taper_end_width']/2.0, hole_center_y + param['grating_pad_width']/2.0 - param['beam_width']),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0, hole_center_y + param['grating_pad_width']/2.0 - param['beam_width'] + param['notch_taper_end_depth'])] ) )

				# notches in extra bottom linker beams for grating pad
				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 - param['notch_taper_end_width']/2.0, hole_center_y - param['grating_pad_width']/2.0 + param['beam_width']),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 + param['notch_taper_end_width']/2.0, hole_center_y - param['grating_pad_width']/2.0 + param['beam_width']),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0, hole_center_y - param['grating_pad_width']/2.0 + param['beam_width'] - param['notch_taper_end_depth'])] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 - param['notch_taper_end_width']/2.0, hole_center_y - param['grating_pad_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0 + param['notch_taper_end_width']/2.0, hole_center_y - param['grating_pad_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['grating_pad_offset']/2.0, hole_center_y - param['grating_pad_width']/2.0 + param['notch_taper_end_depth'])] ) )

		if param['meander'] is False:

			if iB == param['num_beams'] or iB == -1: # Special exception used horizontal dummy beam

					# Four different notches on the left and right sides of horizontal dummy beam w/ no holes in it
				cell.add(gdspy.Polygon(3, [
						(beam_center_x - param['beam_len']/2.0 + param['notch_end_offset'] - param['notch_end_width']/2.0, hole_center_y - param['beam_width']/2.0),
						(beam_center_x - param['beam_len']/2.0 + param['notch_end_offset'], hole_center_y - param['beam_width']/2.0 + param['notch_end_depth']),
						(beam_center_x - param['beam_len']/2.0 + param['notch_end_offset'] + param['notch_end_width']/2.0, hole_center_y - param['beam_width']/2.0)] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x - param['beam_len']/2.0 + param['notch_end_offset'] - param['notch_end_width']/2.0, hole_center_y + param['beam_width']/2.0),
						(beam_center_x - param['beam_len']/2.0 + param['notch_end_offset'], hole_center_y + param['beam_width']/2.0 - param['notch_end_depth']),
						(beam_center_x - param['beam_len']/2.0 + param['notch_end_offset'] + param['notch_end_width']/2.0, hole_center_y + param['beam_width']/2.0)] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] - param['notch_end_width']/2.0, hole_center_y - param['beam_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'], hole_center_y - param['beam_width']/2.0 + param['notch_end_depth']),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] + param['notch_end_width']/2.0, hole_center_y - param['beam_width']/2.0)] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] - param['notch_end_width']/2.0, hole_center_y + param['beam_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'], hole_center_y + param['beam_width']/2.0 - param['notch_end_depth']),
						(beam_center_x + param['beam_len']/2.0 - param['notch_end_offset'] + param['notch_end_width']/2.0, hole_center_y + param['beam_width']/2.0)] ) )

		if iB == param['num_beams'] and param['vflagbeam_q'] is True: # Special exception used for vertical "flag" beams to help with alignment to edge of YSO
			# Vertical "flag" beam to help with alignment
			cell.add(gdspy.Polygon(1, [
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'], hole_center_y + param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['beam_width'], hole_center_y + param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['beam_width'], hole_center_y + 3*param['beam_spacing']),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'], hole_center_y + param['beam_width']/2.0 + 3*param['beam_spacing'] - param['beam_width']/2.0)]))

			# 2nd vertical "flag" beam to help with alignment
			cell.add(gdspy.Polygon(1, [
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'], hole_center_y + param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'] + param['beam_width'], hole_center_y + param['beam_width']/2.0),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'] + param['beam_width'], hole_center_y + 3*param['beam_spacing']),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'], hole_center_y + param['beam_width']/2.0 + 3*param['beam_spacing'] - param['beam_width']/2.0)]))

			# Horizontal between "flag" beams
			cell.add(gdspy.Polygon(1, [
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['beam_width'], hole_center_y - param['beam_width']/2.0 + 2*param['beam_spacing']),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'], hole_center_y - param['beam_width']/2.0 + 2*param['beam_spacing']),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'], hole_center_y + param['beam_width']/2.0 + 2*param['beam_spacing']),
					(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['beam_width'], hole_center_y + param['beam_width']/2.0 + 2*param['beam_spacing'])]))

			if param['meander'] is False:
				# Notches in vertical beam
				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['beam_width'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] + param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['beam_width'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] - param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['beam_width'] - param['notch_end_depth'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset']) ] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] + param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] - param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset'] + param['notch_end_depth'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset']) ] ) )

				# Notches in 2nd vertical beam
				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'] + param['beam_width'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] + param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'] + param['beam_width'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] - param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'] + param['beam_width'] - param['notch_end_depth'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset']) ] ) )

				cell.add(gdspy.Polygon(3, [
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] + param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset'] - param['notch_end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_align_offset2'] + param['notch_end_depth'], hole_center_y + 3*param['beam_spacing'] - param['notch_end_offset']) ] ) )

		#cell_edge_x = beam_center_x - param['extra_len_far']
		cell_edge_x = beam_center_x - param['beam_len']/2.0 + param['extra_len_far']
		cell_edge_y = hole_center_y + param['extra_len_far']+ param['extra_len_far']

		if param['meander'] is True:
			scale_list_phc = numpy.linspace(param2['start_sweep_meander'], param2['end_sweep_meander'], param2['num_phc_total'])  # param2['num_phc']
			hole_scale_list_phc = numpy.zeros(2 * param['waveguide_with_end_mirror'] * (param2['num_phc'] * param2['cavity_len'] + 4 * param['num_mirror_holes'] + (param2['num_phc'] - 2) * param2['middle_mirror_len']))  # param2['num_phc']

			holes_in_waveguide = param2['num_phc'] * param2['cavity_len'] + 4 * param['num_mirror_holes'] + (param2['num_phc'] - 2) * param2['middle_mirror_len']
			holes_in_half_waveguide = param2['num_phc'] / 2 * param2['cavity_len'] + 2 * param['num_mirror_holes'] + (param2['num_phc'] / 2 - 1) * param2['middle_mirror_len']

			for y in range(4 * param['waveguide_with_end_mirror']):
				for iHx in range(param['num_mirror_holes'] + param2['cavity_len'] + param2['middle_mirror_len'] / 2):
					hole_scale_list_phc[iHx + y * holes_in_half_waveguide] = scale_list_phc[0 + y * param2['num_phc'] / 2]

				for i in range(param2['num_phc'] / 2 - 2):
					for iHx in range(param2['cavity_len'] + 2):
						hole_scale_list_phc[i * (param2['cavity_len'] + 2) + iHx + param['num_mirror_holes'] + param2['cavity_len'] + param2['middle_mirror_len'] / 2 + y * holes_in_half_waveguide] = scale_list_phc[1 + i + y * param2['num_phc'] / 2]

				for iHx in range(param['num_mirror_holes'] + param2['cavity_len'] + param2['middle_mirror_len'] / 2):
					hole_scale_list_phc[(param2['num_phc'] / 2 - 2) * (param2['cavity_len'] + 2) + iHx + param['num_mirror_holes'] + param2['cavity_len'] + param2['middle_mirror_len'] / 2 + y * holes_in_half_waveguide] = scale_list_phc[param2['num_phc'] / 2 - 2 + 1 + y * param2['num_phc'] / 2]

			# print(scale_list_phc)
			# print(hole_scale_list_phc)

		hole_center_y_new = 0#cell_edge_y

		if iB < param['num_beams'] and iB > -1:

			if param['holes_q'] is True:

				for i in range(0, len(param['rad_list_mat'][:, 1])):

					hole_center_x = cell_edge_x + param['aper_list'][i] / 2.0
					hole_center_y_new += param['aper_list'][i]

					rad = param['rad_list_mat'][i, iB]
					rad2 = param['rad2_list_mat'][i, iB]  # this is the half hole size perp. to propagation axis

					if param['meander'] is True:

						if iB == 1:  # face down
							S = 1
						else:
							S = -1  # face up

						for iM in range(int(param['waveguide_with_end_mirror'])):

							if param['sweep_hole_size_phc_meander'] is True:
								# Phcs on waveguides with sweep
								rad_top = param['rad_list_mat'][i, iB] * hole_scale_list_phc[i + (2 * iM + 1) * holes_in_waveguide]  # for phc above meander waveguide
								print("rad top" + str(rad_top))
								rad2_top = param['rad2_list_mat'][i, iB] * hole_scale_list_phc[i + (2 * iM + 1) * holes_in_waveguide]  # for phc above meander waveguide
								print("rad2 top" + str(rad2_top))
								rad_bottom = param['rad_list_mat'][i, iB] * hole_scale_list_phc[i + 2 * iM * holes_in_waveguide]  # for phc bollow meander waveguide
								print("rad bottom" + str(rad_bottom))
								rad2_bottom = param['rad2_list_mat'][i, iB] * hole_scale_list_phc[i + 2 * iM * holes_in_waveguide]  # for phc bollow meander waveguide
								print("rad2 bottom" + str(rad2_bottom))

								if iB == 1:
									pts = [(hole_center_x + param2['bus_taper_len'] / 2 + rad_top * numpy.cos(x)+3* param['aper_mir'], hole_center_y + rad2_top * numpy.sin(x) + S * (2 * (iM + 1) * param['bus_bend_radius'] + param['wg'] + param['ww'] / 2 + param['beam_width'] / 2 + (1 + iM) * param['y_spacing_between_wabguides'])) for x in philist]
									pts2 = [(hole_center_x + param2['bus_taper_len'] / 2 + rad_bottom * numpy.cos(x)+3* param['aper_mir'], hole_center_y + rad2_bottom * numpy.sin(x) + S * (2 * (iM + 1) * param['bus_bend_radius'] - param['wg'] - param['ww'] / 2 - param['beam_width'] / 2 + (1 + iM) * param['y_spacing_between_wabguides'])) for x in philist]

									cell.add(gdspy.Polygon(3, pts))
									cell.add(gdspy.Polygon(3, pts2))
								if param2["2-axes"] == True and iB!=1:
									pts = [(circle_center_x - (outer_radius350 + inner_radius350) / 2 - param2["wg"] - param['beam_width'] / 2 - param['ww'] / 2.0 + rad2_top * numpy.sin(x), circle_center_y + S * (2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'] + hole_center_y_new + rad_top * numpy.cos(x) + 5 * param['aper_mir'])) for x in philist]
									pts2 = [(circle_center_x - (outer_radius350 + inner_radius350) / 2 + param2["wg"] + param['beam_width'] / 2 + param['ww'] / 2.0 + rad2_bottom * numpy.sin(x), circle_center_y + S * (2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'] + hole_center_y_new + rad_bottom * numpy.cos(x) + 5 * param['aper_mir'])) for x in philist]
									cell.add(gdspy.Polygon(3, pts))
									cell.add(gdspy.Polygon(3, pts2))

							# add the mirror at the end of the meander for even number of waveguide
							if iB==1:
								for y in range(0, int(param['num_mirror_holes_end_meander'] / 6)):
									if (iM + 1) == param['waveguide_with_end_mirror'] and i < 6 and param['waveguide_with_end_mirror'] % 2 == 0:  # param['waveguide_with_end_mirror']
										pts3 = [(hole_center_x + rad * (hole_scale_list_phc[-1] + hole_scale_list_phc[0]) / 2 * numpy.cos(x) - 1536 / 4 - param2['bus_taper_len'] - 6 * y * param['aper_mir'],
												 hole_center_y + rad2 * (hole_scale_list_phc[-1] + hole_scale_list_phc[0]) / 2 * numpy.sin(x) + S * (2 * (iM + 1) * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])) for x in philist]
										cell.add(gdspy.Polygon(3, pts3))

									# add the mirror at the end of the meander for odd number of waveguide
									elif (iM + 1) == param['waveguide_with_end_mirror'] and i < 6 and param['waveguide_with_end_mirror'] % 2 == 1:  # > (len(param['rad_list_mat'][:, 1]) - 6)
										pts3 = [(hole_center_x + rad * (hole_scale_list_phc[-1] + hole_scale_list_phc[0]) / 2 * numpy.cos(x) + param2['bus_taper_len'] + param['Len_bus_waveguide'] + 6 * y * param['aper_mir']-param["beam_width"],
												 hole_center_y + rad2 * (hole_scale_list_phc[-1] + hole_scale_list_phc[0]) / 2 * numpy.sin(x) + S * (2 * (iM + 1) * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'])) for x in philist]
										cell.add(gdspy.Polygon(3, pts3))
							else:
								for y in range(0, int(param['num_mirror_holes_end_meander'] / 6)):
									if (iM + 1) == param['waveguide_with_end_mirror'] and i < 6 and param['waveguide_with_end_mirror'] % 2 == 1:  # > (len(param['rad_list_mat'][:, 1]) - 6)
										pts3 = [(circle_center_x - (outer_radius350 + inner_radius350) / 2  + rad2 * (hole_scale_list_phc[-1] + hole_scale_list_phc[0]) / 2 * numpy.sin(x),
												 circle_center_y + S * (2 * iM * param['bus_bend_radius'] + (1 + iM) * param['y_spacing_between_wabguides'] + hole_center_y_new + rad * (hole_scale_list_phc[-1] + hole_scale_list_phc[0]) / 2 * numpy.cos(x) + param['Len_bus_waveguide_vertical']-param["beam_width"]+param2['bus_taper_len']+ 6 * y * param['aper_mir'])) for x in philist]
										cell.add(gdspy.Polygon(3, pts3))
					else:
						rad = param['rad_list_mat'][i, iB]
						rad2 = param['rad2_list_mat'][i, iB]  # this is the half hole size perp. to propagation axis
						pts = [(hole_center_x + rad * numpy.cos(x), hole_center_y + rad2 * numpy.sin(x)) for x in philist]
						cell.add(gdspy.Polygon(3, pts))

					# now modify the beam width by adding to notch layer

					if i == 0:
						last_beam_width_delta = 0
						next_beam_width_delta = param['beam_width_delta_list'][i + 1]
					elif i == len(param['rad_list']) - 1:
						last_beam_width_delta = param['beam_width_delta_list'][i - 1]
						next_beam_width_delta = 0
					else:
						last_beam_width_delta = param['beam_width_delta_list'][i - 1]
						next_beam_width_delta = param['beam_width_delta_list'][i + 1]

					dy_pts = [
						(cell_edge_x, hole_center_y + param['beam_width'] / 2.0 + (last_beam_width_delta + param['beam_width_delta_list'][i]) / 4.0),
						(cell_edge_x + param['aper_list'][i] / 2.0, hole_center_y + param['beam_width'] / 2.0 + param['beam_width_delta_list'][i] / 2.0),
						(cell_edge_x + param['aper_list'][i], hole_center_y + param['beam_width'] / 2.0 + (next_beam_width_delta + param['beam_width_delta_list'][i]) / 4.0),
						(cell_edge_x + param['aper_list'][i], hole_center_y + param['beam_width'] / 2.0),
						(cell_edge_x, hole_center_y + param['beam_width'] / 2.0),
						(cell_edge_x, hole_center_y + param['beam_width'] / 2.0 + (last_beam_width_delta + param['beam_width_delta_list'][i]) / 4.0)
					]

					dy_pts_bot = [
						(cell_edge_x, hole_center_y - param['beam_width'] / 2.0 - (last_beam_width_delta + param['beam_width_delta_list'][i]) / 4.0),
						(cell_edge_x + param['aper_list'][i] / 2.0, hole_center_y - param['beam_width'] / 2.0 - param['beam_width_delta_list'][i] / 2.0),
						(cell_edge_x + param['aper_list'][i], hole_center_y - param['beam_width'] / 2.0 - (next_beam_width_delta + param['beam_width_delta_list'][i]) / 4.0),
						(cell_edge_x + param['aper_list'][i], hole_center_y - param['beam_width'] / 2.0),
						(cell_edge_x, hole_center_y - param['beam_width'] / 2.0),
						(cell_edge_x, hole_center_y - param['beam_width'] / 2.0 - (last_beam_width_delta + param['beam_width_delta_list'][i]) / 4.0)
					]

					cell.add(gdspy.Polygon(3, dy_pts))
					cell.add(gdspy.Polygon(3, dy_pts_bot))

					cell_edge_x = cell_edge_x + param['aper_list'][i]
		# AD Vertical linkers to facilitate transfer of all beams at the same time
		# SC two more horizontal linkers to accommodate two more notches
		if param['vert_linker_q'] is True and iB < param['num_beams']:
				if param['meander'] is False:
					# Left side first
					cell.add(gdspy.Polygon(1, [
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], hole_center_y + param['beam_width'] / 2.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + param['vert_linker_width_left'], hole_center_y + param['beam_width'] / 2.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + param['vert_linker_width_left'], hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] - param['beam_width']),
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] - param['beam_width'])]))

					#  Add 2nd vertical linker on left side for stability and improved adhesion
					cell.add(gdspy.Polygon(1, [
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + param['vert_linker_width_left'] + param['vert_linker_gap'], hole_center_y + param['beam_width'] / 2.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2 * param['vert_linker_width_left'] + param['vert_linker_gap'], hole_center_y + param['beam_width'] / 2.0),
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + 2 * param['vert_linker_width_left'] + param['vert_linker_gap'], hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] - param['beam_width']),
						(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + param['vert_linker_width_left'] + param['vert_linker_gap'], hole_center_y + param['beam_width'] / 2.0 + param['beam_spacing'] - param['beam_width'])]))

				# SO left side tether:
				if iB==0 and param['meander'] is True:

					cell.add(gdspy.Polygon(1, [
						(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide']-param2['width_taper']/2, outerbox_y_max- param['box_buffer'] - param['beam_width'] / 2),
						(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2, outerbox_y_max- param['box_buffer'] - param['beam_width'] / 2),
						(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2, outerbox_y_min+ param['box_buffer'] + param['beam_width'] / 2),
						(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide']-param2['width_taper']/2,outerbox_y_min+ param['box_buffer'] + param['beam_width'] / 2)]))
					# cell.add(gdspy.Polygon(1, [
					# 	(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
					# 	(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + param['vert_linker_width_left'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
					# 	(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'] + param['vert_linker_width_left'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
					# 	(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2)]))

					#right side
					cell.add(gdspy.Polygon(1, [
						(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
						(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2),
						(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2),
						(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2)]))

				if param['extra_left_notches'] is True:

					#  Add linkers and notches pinches in left and right frame:
					if param['meander'] is True:
						a=int(outerbox_y_max - outerbox_y_min- param['box_buffer'] - param['beam_width'])/(1+param2['separation_notches_meander_sides']+param['beam_width'])

						for i in range(a): #3+param['waveguide_with_end_mirror']*3
							# left
							cell.add(gdspy.Polygon(1, [
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left'] - param2['width_taper'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) - param['beam_width'] / 2.0),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left'] - param2['width_taper'] / 2, outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) + param['beam_width'] / 2.0),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left'] - param2['width_taper'] / 2 - param['vert_linker_offset'], outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) + param['beam_width'] / 2.0),
								(circle_center_x - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left'] - param2['width_taper'] / 2 - param['vert_linker_offset'], outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) - param['beam_width'] / 2.0)
							]))
							cell.add(gdspy.Polygon(3, [
								(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2 - param['notch_end_offset'] - param['notch_end_width'] / 2.0, outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) - param['beam_width'] / 2.0),
								(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2 - param['notch_end_offset'], outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) + param['notch_end_depth'] - param['beam_width'] / 2.0),
								(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2 - param['notch_end_offset'] + param['notch_end_width'] / 2.0, outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) - param['beam_width'] / 2.0 )]))

							cell.add(gdspy.Polygon(3, [
								(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2 - param['notch_end_offset'] - param['notch_end_width'] / 2.0, outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) + param['beam_width'] / 2.0 ),
								(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2- param['notch_end_offset'], outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) - param['notch_end_depth'] + param['beam_width'] / 2.0),
								(circle_center_x  - param['bus_bend_radius'] - param2['buffer_siliconpad_bendwaveguide'] - param['vert_linker_width_left']-param2['width_taper']/2 - param['notch_end_offset'] + param['notch_end_width'] / 2.0,outerbox_y_max - param['box_buffer'] - param['beam_width'] - (i) * (param2['separation_notches_meander_sides']+param['beam_width']) + param['beam_width'] / 2.0)]))
						#right
						for i in chain(range(a/2+1), range(a/2+7, a)):
							cell.add(gdspy.Polygon(1, [
								(beam_center_x + param['beam_len'] / 2.0, outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) - param['beam_width'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0, outerbox_y_min  + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) + param['beam_width'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) + param['beam_width'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) - param['beam_width'] / 2.0)
							]))

							cell.add(gdspy.Polygon(3, [
								(beam_center_x + param['beam_len'] / 2.0 - param['notch_end_offset'] - param['notch_end_width'] / 2.0, outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) - param['beam_width'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['notch_end_offset'], outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) + param['notch_end_depth'] - param['beam_width'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['notch_end_offset'] + param['notch_end_width'] / 2.0, outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) - param['beam_width'] / 2.0)]))

							cell.add(gdspy.Polygon(3, [
								(beam_center_x + param['beam_len'] / 2.0 - param['notch_end_offset'] - param['notch_end_width'] / 2.0, outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) + param['beam_width'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['notch_end_offset'], outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) - param['notch_end_depth'] + param['beam_width'] / 2.0),
								(beam_center_x + param['beam_len'] / 2.0 - param['notch_end_offset'] + param['notch_end_width'] / 2.0, outerbox_y_min + param['box_buffer'] + param['beam_width'] + (i) * (param2['separation_notches_meander_sides'] + param['beam_width']) + param['beam_width'] / 2.0)]))

							#extra in right:
							# param['array_orig_y'] + iB * param['beam_spacing']
							for ix in range(2):
								cell.add(gdspy.Polygon(1, [
									(beam_center_x + param['beam_len'] / 2.0, param['array_orig_y'] + ix * param['beam_spacing'] - param['beam_width'] / 2.0),
									(beam_center_x + param['beam_len'] / 2.0, param['array_orig_y'] + ix * param['beam_spacing'] + param['beam_width'] / 2.0),
									(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], param['array_orig_y'] + ix * param['beam_spacing'] + param['beam_width'] / 2.0),
									(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], param['array_orig_y'] + ix * param['beam_spacing'] - param['beam_width'] / 2.0)
								]))

					else:
						#  Add first horizontal linker for notch
						cell.add(gdspy.Polygon(1, [
							(beam_center_x - param['beam_len'] / 2.0, hole_center_y - param['beam_width'] / 2.0 + param2['beam_spacing'] / 3.0),
							(beam_center_x - param['beam_len'] / 2.0, hole_center_y + param['beam_width'] / 2.0 + param2['beam_spacing'] / 3.0),
							(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], hole_center_y + param['beam_width'] / 2.0 + param2['beam_spacing'] / 3.0),
							(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], hole_center_y - param['beam_width'] / 2.0 + param2['beam_spacing'] / 3.0)
						]))

						#  Add second horizontal linker for notch
						cell.add(gdspy.Polygon(1, [
							(beam_center_x - param['beam_len'] / 2.0, hole_center_y - param['beam_width'] / 2.0 + param2['beam_spacing'] * 2.0 / 3.0),
							(beam_center_x - param['beam_len'] / 2.0, hole_center_y + param['beam_width'] / 2.0 + param2['beam_spacing'] * 2.0 / 3.0),
							(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], hole_center_y + param['beam_width'] / 2.0 + param2['beam_spacing'] * 2.0 / 3.0),
							(beam_center_x - param['beam_len'] / 2.0 + param['vert_linker_offset'], hole_center_y - param['beam_width'] / 2.0 + param2['beam_spacing'] * 2.0 / 3.0)
						]))

				if param['meander'] is False:
					#  Add number of mirror holes as text between vertical linkers to better identify intended Q after transfer
					cell.add(gdspy.Text(3,
										str(param['num_mirror_holes']), 2750,  # int(2*param['vert_linker_offset']),
										(beam_center_x - param['beam_len'] / 2.0 + 1.25 * param['vert_linker_offset'], hole_center_y + param['beam_spacing'] / 5)))

					#  Add hole group # as text between vertical linkers to better identify intended resonance wavelength after transfer
					cell.add(gdspy.Text(3,
										str(param['hole_group']), 1500,  # int(2*param['vert_linker_offset']),
										(beam_center_x - param['beam_len'] / 2.0 + 1.4 * param['vert_linker_offset'] + param['vert_linker_width_left'] + param['vert_linker_gap'], hole_center_y + param['beam_spacing'] / 10.0)))

					#  Add "S" or "A" as text between vertical linkers to better identify symmetric vs asymmetric cavities
					if param['sym_hole_taper_q'] is False:
						cell.add(gdspy.Text(3, "A", 1500,  # int(2*param['vert_linker_offset']),
											(beam_center_x - param['beam_len'] / 2.0 + 1.4 * param['vert_linker_offset'] + param['vert_linker_width_left'] + param['vert_linker_gap'], hole_center_y + param['beam_spacing'] / 2.0)))

					else:
						cell.add(gdspy.Text(3, "S", 1500,  # int(2*param['vert_linker_offset']),
											(beam_center_x - param['beam_len'] / 2.0 + 1.4 * param['vert_linker_offset'] + param['vert_linker_width_left'] + param['vert_linker_gap'], hole_center_y + param['beam_spacing'] / 2.0)))

				#  Now right side linkers
				if param['grating_pad_q'] is True and iB == -1: # SC: bottom vertical linker

					cell.add(gdspy.Polygon(1, [
						(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'], hole_center_y + param['beam_spacing'] - param['end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['beam_spacing'] - param['end_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['beam_width']/2.0),
						(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'], hole_center_y + param['beam_width']/2.0)]))

				if param['grating_pad_q'] is True and iB < (param['num_beams']-1) and iB > -1:
					cell.add(gdspy.Polygon(1, [
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'], hole_center_y + param['end_width']/2.0),
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['end_width']/2.0),
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['end_width']/2.0 + param['beam_spacing'] - param['end_width']),
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'], hole_center_y + param['end_width']/2.0 + param['beam_spacing'] - param['end_width'])]))

				if param['grating_pad_q'] is True and iB == (param['num_beams']-1):
					if param['meander'] is True:
						# cell.add(gdspy.Polygon(1, [
						# 	(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], hole_center_y + param['end_width'] / 2.0),
						# 	(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['end_width'] / 2.0),
						# 	(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['beam_spacing'] - param['beam_width'] ),
						# 	(beam_center_x + param['beam_len'] / 2.0 - param['vert_linker_offset'], hole_center_y + param['beam_spacing'] - param['beam_width'] / 2.0)]))

						# top right linker
						cell.add(gdspy.Polygon(1, [
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] + param2['grating_pad_length'] - param['vert_linker_width_right'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], outerbox_y_max - param['box_buffer'] - param['beam_width'] / 2 - param['supporting_bar_width']),
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y + param['beam_spacing'] - param2['grating_pad_width']/2 ),
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] + param2['grating_pad_length'] - param['vert_linker_width_right'], hole_center_y + param['beam_spacing'] - param2['grating_pad_width']/2 )]))

						# bottom right linker
						cell.add(gdspy.Polygon(1, [
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] + param2['grating_pad_length'] - param['vert_linker_width_right'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], outerbox_y_min + param['box_buffer'] + param['beam_width'] / 2 + param['supporting_bar_width']),
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'], hole_center_y - 2 * param['grating_pad_spacing'] -param2['grating_pad_width']*3/2 ),
							(beam_center_x + param['beam_len'] / 2.0 - param['grating_pad_offset'] - param['grating_pad_length'] + param2['grating_pad_length'] - param['vert_linker_width_right'], hole_center_y - 2 * param['grating_pad_spacing'] -param2['grating_pad_width']*3/2)]))


					else:
						cell.add(gdspy.Polygon(1, [
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'], hole_center_y + param['end_width']/2.0),
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['end_width']/2.0),
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'] - param['vert_linker_width_right'], hole_center_y + param['beam_spacing'] - param['beam_width']/2.0),
							(beam_center_x + param['beam_len']/2.0 - param['vert_linker_offset'], hole_center_y + param['beam_spacing'] - param['beam_width']/2.0)]))

def write_arrow(cell,center_x,center_y,angle,scale=1):

	shaft_width = scale*20000
	shaft_len = scale*60000
	head_size = scale*40000

	pts = [ (-shaft_width/2.0, -shaft_len/2.0),
			(-shaft_width/2.0, shaft_len/2.0),
			(-head_size/2.0, shaft_len/2.0),
			(0,shaft_len/2.0 + head_size),
			(head_size/2.0, shaft_len/2.0),
			(shaft_width/2.0, shaft_len/2.0),
			(shaft_width/2.0, -shaft_len/2.0),
			(-shaft_width/2.0, -shaft_len/2.0) ]

	pts_rot = [ (center_x + x[0]*math.cos(angle*3.142/180) + x[1]*math.sin(angle*3.142/180), center_y - x[0]*math.sin(angle*3.142/180) + x[1]*math.cos(angle*3.142/180)) for x in pts]

	cell.add(gdspy.Polygon(2,pts_rot))

def write_alignment(cell,center_x,center_y,scale=1,angle=0,bounding=None,layer=2,aspect=1):

	shaft_width=2e4*scale*aspect
	shaft_len=2e4*scale
	bounding_width = 3e5

	marker_eps = 10

	pts = [ (0, 0),
			(0 + shaft_len, 0),
			(0 + shaft_len, 0 + shaft_width),
			(0 + shaft_width, 0 + shaft_width),
			(0 + shaft_width, 0 + shaft_len),
			(0, 0 + shaft_len),
			(0, 0) ]


	pts_rot = [ (center_x + marker_eps + x[0]*math.cos(angle*3.142/180) + x[1]*math.sin(angle*3.142/180), center_y + marker_eps - x[0]*math.sin(angle*3.142/180) + x[1]*math.cos(angle*3.142/180)) for x in pts]
	cell.add(gdspy.Polygon(layer,pts_rot))

	pts1 = [ (0, 0),
			(0 - shaft_len, 0),
			(0 - shaft_len, 0 - shaft_width),
			(0 - shaft_width, 0 - shaft_width),
			(0 - shaft_width, 0 - shaft_len),
			(0, 0 - shaft_len),
			(0, 0) ]

	pts1_rot = [ (center_x - marker_eps + x[0]*math.cos(angle*3.142/180) + x[1]*math.sin(angle*3.142/180), center_y - marker_eps - x[0]*math.sin(angle*3.142/180) + x[1]*math.cos(angle*3.142/180)) for x in pts1]
	cell.add(gdspy.Polygon(layer,pts1_rot))

	if bounding is not None:
		bounding.add(gdspy.Polygon(2,[
			(center_x + bounding_width,center_y + bounding_width),
			(center_x - bounding_width,center_y + bounding_width),
			(center_x - bounding_width,center_y - bounding_width),
			(center_x + bounding_width,center_y - bounding_width),
			(center_x + bounding_width,center_y + bounding_width)
			]))

def write_alignment_small(cell,center_x,center_y,scale=1,angle=0,bounding=None,layer=2):

	shaft_width=30e3*scale
	shaft_len=74e3*scale
	bounding_width = 3e5

	marker_eps = 10

	pts = [ (0, 0),
			(0 + shaft_width, 0),
			(0 + shaft_width, 0 + shaft_len),
			(0, 0 + shaft_len),
			(0, 0) ]

	pts_rot = [ (center_x + marker_eps + x[0]*math.cos(angle*3.142/180) + x[1]*math.sin(angle*3.142/180), center_y + marker_eps - x[0]*math.sin(angle*3.142/180) + x[1]*math.cos(angle*3.142/180)) for x in pts]
	cell.add(gdspy.Polygon(layer,pts_rot))

	pts1 = [ (0, 0),
			(0 - shaft_width, 0),
			(0 - shaft_width, 0 - shaft_len),
			(0, 0 - shaft_len),
			(0, 0) ]

	pts1_rot = [ (center_x - marker_eps + x[0]*math.cos(angle*3.142/180) + x[1]*math.sin(angle*3.142/180), center_y - marker_eps - x[0]*math.sin(angle*3.142/180) + x[1]*math.cos(angle*3.142/180)) for x in pts1]
	cell.add(gdspy.Polygon(layer,pts1_rot))

	if bounding is not None:
		cell.add(gdspy.Polygon(9,[
			(center_x + bounding_width,center_y + bounding_width),
			(center_x - bounding_width,center_y + bounding_width),
			(center_x - bounding_width,center_y - bounding_width),
			(center_x + bounding_width,center_y - bounding_width),
			(center_x + bounding_width,center_y + bounding_width)
			]))
def write_heater_line(cell, center_x, center_y, height, width,layer=9):

	pts = [ (center_x, center_y),
			(center_x, center_y + height),
			(center_x + width, center_y + height),
			(center_x + width, center_y),
			(center_x, center_y)
			 ]

	cell.add(gdspy.Polygon(layer,pts)) #Previously used for single vertical lines at the bottom of the pattern

def write_line(cell,center_x,center_y,height,width,vbeamwidth,layer=1):

	########### AD This has been heavily modified to create a vertical beam + breakpoints to aid in the alignment of transferred beams to trenches in YSO

	# pts = [ (center_x-width/2.0,center_y-len/2.0),
	# 		(center_x+width/2.0,center_y-len/2.0),
	# 		(center_x+width/2.0,center_y+len/2.0),
	# 		(center_x-width/2.0,center_y+len/2.0),
	# 		(center_x-width/2.0,center_y-len/2.0) ]

	# cell.add(gdspy.Polygon(layer,pts)) Previously used for single vertical lines at the bottom of the pattern

	pts_outer = gdspy.Polygon(2,[ (center_x-width/2.0,center_y-height/2.0),
			(center_x+width/2.0,center_y-height/2.0),
			(center_x+width/2.0,center_y+height/2.0),
			(center_x-width/2.0,center_y+height/2.0),
			(center_x-width/2.0,center_y-height/2.0) ])

	pts_outer.fillet(2000) # AD Use fillet command to round corners by 2 um radius

	width2 = vbeamwidth/2.0
	height2 = height/2.0

	vnotch_depth=vbeamwidth/2.3
	vnotch_height=numpy.sqrt(2.0)*vnotch_depth

	pts_inner = [ (center_x-width2,center_y-height2),
			(center_x+width2,center_y-height2),
			(center_x+width2,center_y+height2),
			(center_x-width2,center_y+height2),
			(center_x-width2,center_y-height2) ]

	left_arrow = [(center_x-width2,center_y),(center_x-width2-1500,center_y),(center_x-width2,center_y+height2-2*vnotch_height),(center_x-width2,center_y)]
	right_arrow = [(center_x+width2,center_y),(center_x+width2+1500,center_y),(center_x+width2,center_y+height2-2*vnotch_height),(center_x+width2,center_y)]
	hole_arrow = [(center_x-width2/2.0,center_y+height/16.0),
			(center_x+width2/2.0,center_y+height/16.0),
			(center_x+width2/2.0,center_y+height/4.0),
			(center_x-width2/2.0,center_y+height/4.0),
			(center_x-width2/2.0,center_y-height/16.0)]

	cell.add(pts_outer)
	cell.add(gdspy.Polygon(1,pts_inner))
	cell.add(gdspy.Polygon(1,left_arrow))
	cell.add(gdspy.Polygon(1,right_arrow))
	cell.add(gdspy.Polygon(3,hole_arrow))


	vnotch_leftdown=[ (center_x-width2,center_y-height2),
			(center_x-width2+vnotch_depth,center_y-height2+vnotch_height),
			(center_x-width2,center_y-height2+2.0*vnotch_height),(center_x-width2,center_y-height2) ]
	vnotch_rightdown=[ (center_x+width2,center_y-height2),
			(center_x+width2-vnotch_depth,center_y-height2+vnotch_height),
			(center_x+width2,center_y-height2+2.0*vnotch_height),(center_x+width2,center_y-height2) ]
	vnotch_leftup=[ (center_x-width2,center_y+height2),
			(center_x-width2+vnotch_depth,center_y+height2-vnotch_height),
			(center_x-width2,center_y+height2-2.0*vnotch_height),(center_x-width2,center_y+height2) ]
	vnotch_rightup=[ (center_x+width2,center_y+height2),
			(center_x+width2-vnotch_depth,center_y+height2-vnotch_height),
			(center_x+width2,center_y+height2-2.0*vnotch_height),(center_x+width2,center_y+height2) ]
	cell.add(gdspy.Polygon(3,vnotch_leftdown))
	cell.add(gdspy.Polygon(3,vnotch_rightdown))
	cell.add(gdspy.Polygon(3,vnotch_leftup))
	cell.add(gdspy.Polygon(3,vnotch_rightup))

# Used for writing trenches
# def write_trench(cell,center_x,center_y,trench_length,trench_width,layer=1):

# 	pts_trench = gdspy.Polygon(2,[ (center_x-trench_width/2.0,center_y-trench_length/2.0),
# 			(center_x+trench_width/2.0,center_y-trench_length/2.0),
# 			(center_x+trench_width/2.0,center_y+trench_length/2.0),
# 			(center_x-trench_width/2.0,center_y+trench_length/2.0),
# 			(center_x-trench_width/2.0,center_y-trench_length/2.0) ])

# 	pts_trench.fillet(2000) # AD Use fillet command to round corners by 2 um radius, can remove later when switching to photoresist

# 	cell.add(pts_trench)

def f_eta(x,eta):
	# helper function for taper_type == chan below
	if x<=0.5:
		return 0.5*(2*x)**eta
	else:
		return 1-0.5*(2*(1-x))**eta

def make_cavity_params_tm_refl(param):

	param['aper_list'] = []
	param['rad_list'] = []
	param['rad_list2'] = []
	param['beam_width_delta_list'] = []

	if param['mirror_len'] >= 0:
		# and set the other parameters

		if param['end_taper_n'] > 0:
			taper_scale_list = numpy.linspace(0,1.0,param['end_taper_n']+2)
			taper_scale_list = taper_scale_list[1:-1]
			#rad_list_mir_taper = numpy.linspace(0,param['hole_rad2'],param['end_taper_n']+2)
			#rad_list_mir_taper = rad_list_mir_taper[1:-1] # drop first and last element
			#print rad_list_mir_taper
			taper_scale_list_trunc = numpy.flipud(taper_scale_list)
			taper_scale_list_left = numpy.flipud(taper_scale_list_trunc)

		nmax = (param['cavity_len'] - 1)/2.0
		aper_list_idx = numpy.linspace(0-nmax,nmax,param['cavity_len'])
		slot_list_idx = numpy.linspace(0-nmax,nmax,param['cavity_len'])

		if param['taper_type'] == 'parabola':

			aper_a = (param['aper_mir'] - param['aper_cav'])/(nmax*nmax - 0.25)
			aper_c = (param['aper_mir'] - 4*param['aper_cav']*nmax*nmax)/(1-4*nmax*nmax)

			aper_list_cav = [aper_a*x*x + aper_c for x in aper_list_idx]

		elif param['taper_type'] == 'chan':
			# function from Jasper Chan's PhD thesis which has smooth derivatives
			# 2x**3 - 3x**2 + 1
			aper_list_x = [f_eta(abs(x/(nmax+1)),param['eta']) for x in aper_list_idx]
			aper_list_cav = [param['aper_mir'] - (param['aper_mir']-param['aper_cav'])*(2*x**3 - 3*x**2 + 1) for x in aper_list_x]

		# Parameters for the slot
		slot_width_list_cav = [0 for x in slot_list_idx]
		slot_width_list_cav[(param['cavity_len']-1)/2] = param['cavity_slot_width']

		if param['cavity_slot_q'] is True and param['cavity_len']%2==1:
			aper_list_cav[(param['cavity_len']-1)/2] = param['cavity_slot_length']

		aper_list_mir = [param['aper_mir'] for x in range(param['mirror_len'] + param['end_taper_n'])]
		aper_list_mir_middle = [param['aper_mir'] for x in range(param['middle_mirror_len'])] #mirror between PhCs (i.e not the end mirrors)

		slot_width_mir = [0 for x in aper_list_mir]
		slot_width_mir_middle = [0 for x in aper_list_mir_middle]

		if param['meander'] is True:
			param['slot_width_list'] = copy(slot_width_mir)
			param['slot_width_list'].extend(slot_width_list_cav)
			param['slot_width_list'].extend(slot_width_mir_middle)
			for ip in range(param['num_phc']-2):
				param['slot_width_list'].extend(slot_width_list_cav)
				param['slot_width_list'].extend(slot_width_mir_middle)
			param['slot_width_list'].extend(slot_width_list_cav)
			param['slot_width_list'].extend(slot_width_mir)

		else:
			param['slot_width_list'] = copy(slot_width_mir)
			param['slot_width_list'].extend(slot_width_list_cav)
			param['slot_width_list'].extend(slot_width_mir)

		if param['taper_q'] is False:
			aper_list_cav = [param['aper_mir'] for x in aper_list_idx]

		if param['sym_hole_taper_q'] is True:
			# print "Inside sym taper q"
			param['aper_list'].extend(aper_list_mir)
			param['aper_list'].extend(aper_list_cav)


			if param['meander'] is True:
				param['aper_list'].extend(aper_list_mir_middle)
				for ip in range(param['num_phc'] - 2):
					param['aper_list'].extend(aper_list_cav)
					param['aper_list'].extend(aper_list_mir_middle)
				param['aper_list'].extend(aper_list_cav)

			param['aper_list'].extend(aper_list_mir)

		else:
			param['aper_list'] = copy(aper_list_cav)
			param['aper_list'].extend(aper_list_mir)

		param['rad_list'] = [param['hole_rad'] for x in param['aper_list']]
		param['rad_list2'] = [param['hole_rad2'] for x in param['aper_list']]

		if param['cavity_slot_q'] is True and param['cavity_len']%2==1:
			param['rad_list'][(param['cavity_len']-1)/2] = 0
			param['rad_list2'][(param['cavity_len']-1)/2] = 0

		if param['end_taper_n'] > 0 and param['sym_hole_taper_q'] is False:
			param['rad_list2'][-param['end_taper_n']:] = taper_scale_list[-1::-1] * param['hole_rad2']
			param['rad_list'][-param['end_taper_n']:] = taper_scale_list[-1::-1] * param['hole_rad']

		if param['end_taper_n'] > 0 and param['sym_hole_taper_q'] is True:

			# Swap first end_taper_n elements with taper_scale_list_left
			param['rad_list'][0:param['end_taper_n']] = taper_scale_list_left * param['hole_rad']
			param['rad_list2'][0:param['end_taper_n']] = taper_scale_list_left * param['hole_rad2']

			# Swap last end_taper_n elements with taper_scale_list_left
			param['rad_list2'][-param['end_taper_n']:] = taper_scale_list_trunc * param['hole_rad2']
			param['rad_list'][-param['end_taper_n']:] = taper_scale_list_trunc * param['hole_rad']

		param['beam_width_delta_list'] = numpy.linspace(0,0,len(param['aper_list']))

	#now build list of things that should go on the left side
	param['aper_list_tm'] = []
	param['hx_list_tm'] = []
	param['hy_list_tm'] = []
	param['beam_width_delta_list_tm'] = []

	ap_list_tm = numpy.linspace(param['aper_mir'],param['aper_mir'],param['extra_te_mir'] + param['mirror_len'])
	wy_list_tm = numpy.linspace(param['beam_width'],param['beam_width'],len(ap_list_tm))
	dum_list = numpy.linspace(1,1,len(ap_list_tm))
	hx_list_tm = dum_list*param['hole_rad']*2.0
	hy_list_tm = dum_list*param['hole_rad2']*2.0
	delta_list_tm = (wy_list_tm - param['beam_width'])

	param['aper_list_tm'].extend(ap_list_tm)
	param['hx_list_tm'].extend(hx_list_tm/2.0)
	param['hy_list_tm'].extend(hy_list_tm/2.0)
	param['beam_width_delta_list_tm'].extend(delta_list_tm)

	#param['aper_list_tm'].extend(numpy.linspace(param['aper_mir'],param['aper_mir'],param['extra_te_mir']))

	#hx

	if param['extra_te_mir'] > 0 and param['sym_hole_taper_q'] is False:

		# print "Inside last if statement"
		param['aper_list_tm'].extend(param['aper_list'])
		param['hx_list_tm'].extend(param['rad_list'])
		param['hy_list_tm'].extend(param['rad_list2'])
		param['beam_width_delta_list_tm'].extend(param['beam_width_delta_list'])


		param['aper_list'] = param['aper_list_tm']
		param['rad_list'] = param['hx_list_tm']
		param['rad_list2'] = param['hy_list_tm']
		param['beam_width_delta_list'] = param['beam_width_delta_list_tm']

	# print repr(param['aper_list'])
	if param['meander'] is True:
		param['beam_len'] = param['extra_len_near'] + param['extra_len_far']
		param['beam_center_x'] = (param['extra_len_far'] - param['extra_len_near']) + numpy.sum(param['aper_list']) / 2.0

	else:
		param['beam_len'] = numpy.sum(param['aper_list']) + param['extra_len_near'] + param['extra_len_far']
		param['beam_center_x'] = (param['extra_len_far'] - param['extra_len_near']) + numpy.sum(param['aper_list']) / 2.0

	# Create new matrices that will contain all rad list data for different hole sizes
	param['rad_list_mat'] = numpy.zeros((len(param['aper_list']), param['num_beams']))
	param['rad2_list_mat'] = numpy.zeros((len(param['aper_list']), param['num_beams']))

	for iHH in range(0,param['num_beams']):

		param['rad_list_mat'][:,iHH]=numpy.array(param['rad_list'])*hole_scale_list[iHH]
		param['rad2_list_mat'][:,iHH]=numpy.array(param['rad_list2'])*hole_scale_list[iHH]

# don't know what the options are... presumably 1e-9 is the unit (in meters)
beams = gdspy.Cell('beams')

param = {}

param['vary_number_meander_waveguide']=False
param['sweep_hole_size_phc_meander']=True

###YSO PhC hole scale
param['start_sweep_GC10'] = 0.95
param['end_sweep_GC10'] = 0.95
#meander sweep
xin = numpy.linspace(-0.015, 0.025, 6)
startsweep = 1 + xin
endsweep = 1.03+xin

#CaWO4 PhC hole scale

num_rows = 2 #repetition in y
num_cols = 2 #reptition in x

device_y_um = 108
device_y_nm = device_y_um * 1000

device_x_and_spacing_x_um = 347
device_x_and_spacing_x_nm = device_x_and_spacing_x_um * 1000

param['start_sweep_GC10_CaWO4'] = 0.960
param['end_sweep_GC10_CaWO4'] = 1.2

#meander sweep
#PhC resonance sweep parameters
#sweep goes from left to right
#xin  used to calculate the scaling factor of the PhC hole sizes hx,hy; hx = startsweep*hx_0; hy = startsweep*hy_0
#xin is used for central PhC in meander, if we use 4 columns then  xin will be 4 element array
#xin = numpy.linspace(-0.1 -0.004*10 , 0, 4)
# xin = numpy.linspace(-0.005*8, 0.005*4, 4) #expected shift of centr. wav. -8 nm -4 nm 0 4nm
xin = numpy.linspace(-0.005*8, 0.005*4, num_cols)
#startsweep will be the scaling factor for the smallest holes in the meander
startsweep =(1 +xin)*0.93
#endsweep will be the scaling factor for the biggest holes in the meander
endsweep =(1 +xin)*0.99
# for 0.93 - 0.99 I expect around 12 nm sweep of PhC wav. res. for single meander

startsweep_CaWO4=startsweep
endsweep_CaWO4=endsweep

# print ("YSOstart" + str(startsweep))
# print ("YSOend" + str(endsweep))
print ("calciumstart" + str(startsweep_CaWO4))
print ("calciumend" + str(endsweep_CaWO4))

target_GC_center_lambda_nm = 1544

mirror_num = 7
mirror_list = mirror_num*numpy.ones(num_rows, dtype=int) #number of mirrors is held constant for all the PhCs

#assume cavity list is the same for now
cavity_num = 12
cavity_list = cavity_num*numpy.ones(num_rows, dtype=int)

grating_taper_length_nm = 185000
grating_taper_length_list = grating_taper_length_nm*numpy.ones(num_cols, dtype=int)

param['hole_size_scale_meander_no_sweep']=1.088 # scale in case of no sweep for meander

param['vert_linker_q'] = True
param['vert_linker_offset'] = 3000
param['vert_linker_width_left'] = 630 + 50 # Added 50 nm for sample on YSO
param['vert_linker_width_right'] = 2500
param['vert_linker_gap'] = 0
param['vert_align_offset'] = 5000
param['vert_align_offset2'] = 2000

param['supporting_bar_width'] = 5000	#SC supporting bar top and bottom

param['beam_spacing'] = 4500
param['beam_width'] = 630+35
param['num_mirror_holes'] = 6

param['meander'] = False
param['hoste']=1 #1==YSO, 2==Tio2
param['line_defect']=False
param['num_hole_size'] = 10	#4 SC trying to decrease number of devices in each block
param['num_hole_cluster'] = 1	#3 SC trying to decrease number of devices in each block
param['num_hole_clusterfuck'] = 2
param['num_beams'] = param['num_hole_size'] * param['num_hole_cluster']
param['hole_group'] = 1

# parameters that apply to all things
param['box_buffer'] = 2000	#SC
param['support_tether_width'] = 1000	#SC changes from 650 to 1000 after adding those extra silicon pads
param['support_tether_q'] = False	#SC
param['num_tether_taper'] = 1	#SC tether line within taper
# param['num_tether_device'] =10	#SC tether line to hold the device

param['holes_q'] = True
param['square_holes_q'] = False
param['taper_q'] = True
param['stepped_width_q'] = False

param['end_taper_n'] = 0
param['extra_len_near'] = 2000
param['extra_len_far'] = 2000

# parameters for grating
param['grating_pad_q'] = True
param['num_grating_periods_x'] = 20
param['num_grating_periods_y'] = 33	  #SC total grating y width 12.6 um, maximum overlap with fiber
param['grating_period_x_start'] = 620
param['grating_pad_width'] = 12600 #along y
param['a_2DPhC'] =213.54  #214.97 ########
param['grating_pad_buffer'] = (param['grating_pad_width'] - (param['num_grating_periods_y']-1) * param['a_2DPhC'] * numpy.sqrt(3)) / 2.0
param['grating_pad_offset'] = 2000
param['grating_pad_length'] = 18000
param['grating_pad_spacing'] =3000 #(150000- 2* param['supporting_bar_width']-2*param['box_buffer']-2.0* param['vert_linker_width_left'] - param['num_beams']*param['grating_pad_width'])/(param['num_beams']+1)	#SC changing from 1000 to 3000
param['phaseFactor'] = 0.588
param['vflagbeam_q'] = False

param['extra_te_mir'] =0

param['extra_silicon_pad_q'] = True
param['centercut_silicon_pad'] = 1000 #SC a center line cut for HF undercut
param['extra_silicon_pad_gap'] = 1000 #SC extra silicon pad to decrease the exposure time

# parameters for the break-off notch left side
param['notch_end_offset'] = 1000 #offset from the end of the beam
param['notch_end_depth'] = 75 #depth of the notch
param['notch_end_width'] = param['notch_end_depth']*numpy.sqrt(2.0) # width of the notch base
# param['notch_end_offset'] = param['notch_end_width']/2.0
param['notch_end_offset'] = 1000

# parameters for the break-off notch right (tapered) side
param['notch_taper_end_depth'] = 95
param['notch_taper_end_width'] = param['notch_taper_end_depth']*numpy.sqrt(2.0)
param['notch_end_depth_blade'] = 0 #depth of the notch
param['notch_end_width_blade'] = param['notch_end_depth_blade']*numpy.sqrt(2.0) # width of the notch base
param['notch_end_offset_blade'] = param['notch_end_width_blade']/2.0

param['end_width'] = 230
param['end_width_length'] = 11000
param['end_width_taper_length'] = 3000

param['notch_q'] = True	# whether to make a notch

# parameters for end cut
param['end_cut_q'] = False
param['end_cut_width'] = 1000
param['end_cut_offset'] = 1000

# parameters for the orientation flag
param['flag_end_offset'] = 2000 # offset from the far end of the beam (increasing x value)
param['flag_width'] = 350 # size of flag along the axis of the beam
param['flag_length'] = 500 #distance that flag protrudes from the side
param['flag_q'] = False

param['frame_thickness'] = 150

param['frame_q'] = False
param['window_q'] = True
param['frame_holes_q'] = False

param['grating_q'] = False
param['grating_aper_list'] = [700,700,700]
param['grating_rad_list'] = [120,120,120]
param['grating_offset'] = 1000

param['bend_rad'] = 3000
param['bend_angle'] = numpy.pi/4.0
param['stick_len'] = 7100
param['blade_width'] = 50 # width of the coupling region of the blade
param['blade_taper_len'] = 7000 #length of the taper transition to coupling width
param['normal_side_len'] = 100

param['left_taper_q'] = False
param['left_taper_min_width'] = 200
param['left_taper_len'] = 0

param['num_disks'] = 1
param['disk_spacing'] = 20000
param['disk_rad'] = 3000
param['wg_offset'] = 500
param['wg_width'] = 250

param['sym_hole_taper_q'] = True

param['tm_len'] = 5
param['tm_a_taper_len'] = 3
param['tm_w_taper_len'] = 4
param['tm_ff'] = 0.2
param['tm_aspect'] = 1.0
param['tm_aper'] = 350
param['tm_hx_offset'] = 25
param['tm_hy_offset'] = 25
param['tm_wy'] = 700 + 25

param['pad_end_offset'] = 3000
param['pad_rad'] = 1200
param['pad_q'] = False

param['guard_q'] = False
param['guard_width'] = 500

param['taper_type'] = 'parabola'
param['eta'] = 1.3

# param['n_circle_points'] = 61
param['n_circle_points'] = 199	#SC from 121 to 199, to make the ellipses as smooth as possible
param['n_circle_points_GC'] = 110

ichip = 0

# Now we want to have multiple hole sizes in each cluster, right now this is defined as 4 different hole sizes, but 3 beams of each hole size cluster together

scale_list1 = numpy.linspace(param['start_sweep_GC10'],param['end_sweep_GC10'], param['num_hole_size']) # SC make the sweeping range narrow 1.07-1.126
scale_list2 = numpy.linspace(param['start_sweep_GC10'],param['end_sweep_GC10'], param['num_hole_size']) # for sweeping the size of the holes
scale_list = numpy.append(scale_list1, scale_list2)

scale_list_CaWO4_1 = numpy.linspace(param['start_sweep_GC10_CaWO4'],param['end_sweep_GC10_CaWO4'], param['num_hole_size']) # SC make the sweeping range narrow 1.07-1.126
scale_list_CaWO4_2 = numpy.linspace(param['start_sweep_GC10_CaWO4'],param['end_sweep_GC10_CaWO4'], param['num_hole_size']) # for sweeping the size of the holes
scale_list_CaWO4 = numpy.append(scale_list_CaWO4_1, scale_list_CaWO4_2)

if param['sweep_hole_size_phc_meander'] is True:
	scale_meander = 1
if param['sweep_hole_size_phc_meander'] is False:
	scale_meander = param['hole_size_scale_meander_no_sweep']


hole_scale_list1 = numpy.zeros(param['num_beams'])
hole_scale_list2 = numpy.zeros(param['num_beams'])
hole_scale_list1CaWO4 = numpy.zeros(param['num_beams'])
hole_scale_list2CaWO4 = numpy.zeros(param['num_beams'])
hole_scale_list = numpy.zeros(param['num_beams'])
hole_scale_list_meander= numpy.zeros(param['num_beams'])

# Parameters for 540 nm wide beam w/ center slot in air
param['cavity_slot_q'] = False
param['cavity_slot_width'] = 36 # In nm
param['cavity_slot_length'] = 378 # This number is doubled because cavity length is supposed to be an even number

# parameters for the gash type cavity down in the middle of a narrow beam
param['cavity_gash_q'] = False
param['cavity_gash_width'] = 50
param['cavity_gash_length'] = 14 # This number should be even

for iHx in range(param['num_hole_cluster']):
	for iHy in range(param['num_hole_size']):
		hole_scale_list1[(iHx + param['num_hole_cluster']*iHy)] = scale_list1[iHy]
		hole_scale_list2[(iHx + param['num_hole_cluster']*iHy)] = scale_list2[iHy]

		hole_scale_list1CaWO4[(iHx + param['num_hole_cluster'] * iHy)] = scale_list_CaWO4_1[iHy]
		hole_scale_list2CaWO4[(iHx + param['num_hole_cluster'] * iHy)] = scale_list_CaWO4_2[iHy]

for iHx in range(param['num_hole_clusterfuck']):
	for iHy in range(param['num_hole_size']):
		hole_scale_list_meander[iHx] = scale_meander

for iXM in range(1):
	for iYM in range(1):

		ichip = ichip+1

		chipCenterX = 0#5e6+10e6*(1+iXM)
		chipCenterY = 0#1.2e6*(iYM)

		n_tm_list = [0]

		cavity_list_CaWO4 =cavity_list

		grating_airholescale_list1 = numpy.linspace(-50, 0, param['num_beams'])	# SC air hole diameter sweep
		grating_airholescale_list2 = numpy.linspace(-50, 0, param['num_beams'])  # SC air hole diameter sweep
		constgrating_airholescale_list = [5,5]   # air hole diameter sweep constant

		grating_airholescale_list1CaWO4 = numpy.linspace(-50, 0, param['num_beams'])	# SC air hole diameter sweep linear
		grating_airholescale_list2CaWO4 = numpy.linspace(-50, 0, param['num_beams'])  # SC air hole diameter sweep linear
		# constgrating_airholescale_list1CaWO4 = [0,0]  # SC air hole diameter sweep constant

		# grating_taper_length_list=[10000,25000,50000,100000,150000]
		# notch_depth_list = [250,260,270,280,290,250,260,270,280,290]

		blade_width_list = [200]
		for iX in range(num_cols): #4
			if iX==0:
				y_start=0 #0
				# off_in_end=4#4
			else:
				y_start=0
				# off_in_end = 4#4

			# for iY in range(y_start, num_rows-off_in_end):
			for iY in range(y_start, num_rows):
				# now make a specific design
				param2 = copy(param)

				param2['meander'] = True
				if iX==0:
					param2["meander"]=True
				if iX==1:
					param2["meander"]=True

				param2['grating_pad_length'] = 18000  # 30000 #SC 2D PhC requires less length
				param2['grating_first_index'] = 2.95
				param2['grating_delta_index'] = 0.057
				param2['taper_nonlinear_order'] = 1.3

				param2['y_spacing_between_wabguides'] = 0

				param2['hoste'] = 3

				mirror_listCaWO4 = mirror_list

				if param2['meander']== True:
					mirror_listCaWO4 = mirror_list

					if param2['hoste'] == 1:
						param2['start_sweep_meander'] = startsweep[iX]
						param2['end_sweep_meander'] = endsweep[iX]
					if param2['hoste'] == 3:
						param2['start_sweep_meander'] = startsweep_CaWO4[iX]
						param2['end_sweep_meander'] = endsweep_CaWO4[iX]

				# param2['mirror_len'] = 9
				param2['support_tether_q'] = True
				param2['cavity_len'] = 10

				if param2['hoste'] == 1:  # 1==YSO
					param2['wg'] = 485
					param2["2-axes"] = True
					param2['cavity_len'] = cavity_list[iY]
					param2['mirror_len'] = mirror_list[iY]
					param2['num_mirror_holes'] = param2['mirror_len']
					param2['aper_cav'] = 295.4
					param2['aper_mir'] = 340.4
					param2['beam_width'] = 685 -12 +13 +5
					param2['hole_rad'] = (143.1 - 20.0) / 2.0  # SC changes -40 to -20 for cold developing Decided to add offset to correct for bulk broadening from exposure and/or etching
					param2['hole_rad2'] = (314.9 - 10.0) / 2.0  # SC changes -30 to -10 for cold developing Decided to add offset to correct for bulk broadening from exposure and/or etching
					param2['vflagbeam_q'] = False
					param2['y_spacing_between_wabguides'] = 00
					param2["resonance"] = 1536
					if iY % 2 == 0:
						grating_airholescale_list = grating_airholescale_list1

					if iY % 2 != 0:
						grating_airholescale_list = grating_airholescale_list2

					if param2['meander'] == True:
						# xs=grating_airholescale_list1[iY]
						# grating_airholescale_list=[xs,xs]

						if iY == 0:
							constgrating_airholescale_list = [-15, -15]
						if iY == 1:
							constgrating_airholescale_list = [-17, -17]
						if iY == 2:
							constgrating_airholescale_list = [-23, -23]
						if iY == 3:
							constgrating_airholescale_list = [-25, -25]
						if iY == 4:
							constgrating_airholescale_list = [-27, -27]
						if iY == 5:
							constgrating_airholescale_list = [-29, -29]

				if param2['hoste'] == 3:  # 3==CaWO4
					# pass
					param2['wg'] = 500 #figure out what this waveguide thing does
					param2["resonance"] = target_GC_center_lambda_nm
					param2["2-axes"] = True
					param2['cavity_len'] = cavity_list_CaWO4[iY]
					param2['mirror_len'] = mirror_listCaWO4[iY]
					param2['num_mirror_holes'] = param2['mirror_len']
					param2['aper_cav'] = 349.1 #298
					param2['aper_mir'] = 435.2 #343
					param2['beam_width'] =  600 # SRP: I think this defines the PhC wy
					param2['hole_rad'] = 236.8 / 2.0 # (145.6-20) / 2.0  # SC changes -40 to -20 for cold developing Decided to add offset to correct for bulk broadening from exposure and/or etching
					param2['hole_rad2'] = 407.03 / 2.0  #(307.8-10) / 2.0  # SC changes -30 to -10 for cold developing Decided to add offset to correct for bulk broadening from exposure and/or etching
					param2['vflagbeam_q'] = False
					param2['y_spacing_between_wabguides'] = 00

					param2['num_grating_periods_x'] = 20
					param2['num_grating_periods_y'] = 31  # SC total grating y width 12.6 um, maximum overlap with fiber
					param2['grating_period_x_start'] = 700
					param2['grating_pad_width'] = 12600 #assuming this is the same as the taperendwidth parameter in the lumerical geometry file
					param2['grating_first_index'] = 2.55
					param2['grating_delta_index'] = 0.02
					param2['phaseFactor'] = 0.582
					param2['a_2DPhC'] = 232.61 #looks like this is the same as the
					param2['grating_pad_length'] = 17000
					param2['grating_pad_buffer'] = (param2['grating_pad_width'] - (param2['num_grating_periods_y'] - 1) * param2['a_2DPhC'] * numpy.sqrt(3)) / 2.0
					param2['grating_pad_offset'] = 2000
					param2['grating_pad_spacing'] = 3000
					param2['vflagbeam_q'] = False

					#GC sweep parameters
					#iY is counted from bottom to top
					if param2['meander'] == True:
						constgrating_airholescale_list = [0, 0]

						# if iY == 0:
						# 	constgrating_airholescale_list = [10, 10]
						# if iY == 1:
						# 	constgrating_airholescale_list = [0, 0]
						# if iY == 2:
						# 	constgrating_airholescale_list = [-10, -10]
						# if iY == 3:
						# 	constgrating_airholescale_list = [-20, -20]
						# if iY == 4:
						# 	constgrating_airholescale_list = [-30, -30]
						# if iY == 5:
						# 	constgrating_airholescale_list = [-40, -40]

				grating_x_num = range(param2['num_grating_periods_x']) #range starts 0,1,2...num_grating_period_x-1
				# air_hole_diameter_list_base = [-0.09594*param2['a_2DPhC']*(param2['grating_first_index']-param2['grating_delta_index']*x)**4 + 0.9546*param2['a_2DPhC']*(param2['grating_first_index']-param2['grating_delta_index']*x)**3 - 3.586*param2['a_2DPhC']*(param2['grating_first_index']-param2['grating_delta_index']*x)**2 + 5.542*param2['a_2DPhC']*(param2['grating_first_index']-param2['grating_delta_index']*x) - 1.931*param2['a_2DPhC'] for x in grating_x_num]
				air_hole_diameter_list_base = [-0.09594 * param2['a_2DPhC'] * (param2['grating_first_index'] - param2['grating_delta_index'] * x) ** 4 + 0.9546 * param2['a_2DPhC'] * (param2['grating_first_index'] - param2['grating_delta_index'] * x) ** 3 - 3.586 * param2['a_2DPhC'] * (param2['grating_first_index'] - param2['grating_delta_index'] * x) ** 2 + 5.542 * param2['a_2DPhC'] * (param2['grating_first_index'] - param2['grating_delta_index'] * x) - 1.931 * param2['a_2DPhC'] for x in grating_x_num]
				#print air_hole_diameter_list_base


				param2['grating_taper_length'] = grating_taper_length_list[iX]
				param2['end_width'] = param2['grating_pad_width']
				param2['beam_spacing'] = param2['grating_pad_width'] + param2['grating_pad_spacing']
				param2['extra_len_near'] = param2['grating_pad_offset'] + param2['grating_taper_length'] + param2['grating_pad_length']
				param2['end_width_taper_length'] = 0
				param2['end_width_length'] = 0
				# param2['end_width_length'] = param2['extra_len_near']-param2['end_width_taper_length'] % this generally used to make width taper start directly at end of PhC holes
				param2['extra_len_far'] = 13000 # SC was originally 8000; now add more to accomondate wider vertical linker. AD, was originally 9000 w/ heater pad
				param2['taper_q'] = True
				#param2['beam_width'] = 630 + 50 # Going to vary beam width slightly to compensate for over exposure
				#param2['beam_width'] = beam_width_list[iY]

				param2['middle_mirror_len'] = 2  # number of mirror hole between Phcs
				if param2['meander'] is True:
					param2['box_buffer'] = param['box_buffer']

					param2['num_mirror_holes'] = param2['mirror_len']
				else:
					param2['box_buffer'] = param['box_buffer'] + param2['beam_width'] / 2

				# param2['num_tether_device'] = 30

				param2['square_holes_q'] = False

				param2['tm_ff'] = 0.2*scale_list[0]**2
				param2['tm_len'] = 0
				param2['tm_a_taper_len'] = 0
				param2['tm_w_taper_len'] = 0

				param2['notch_q'] = True
				param2['extra_left_notches'] = True
				# param2['notch_end_depth'] = 285 + 25# The +25 is only added for 650 nm beam width
				param2['notch_end_depth'] = 285 + 25 - 30
				param2['notch_end_width'] = numpy.sqrt(2.0)*param2['notch_end_depth']
				# param2['notch_end_depth'] = notch_depth_list[iX] # used for sweeping notch depth with notch_depth_list
				param2['notch_taper_end_depth'] = 285 + 25 - 30	#SC all notch points have same geometry 285+25-20
				param2['notch_taper_end_width'] = numpy.sqrt(2.0)*param2['notch_taper_end_depth']
				param2['end_taper_n'] = 0
				param2['straight_q'] = True

				param2['vert_linker_offset'] = 2000
				param2['vert_linker_width_left'] = 5000 #SC changing from 2500 to 5000
				param2['vert_linker_width_right'] = 2500	#SC enabling different linker width
				param2['vert_align_offset'] = param2['extra_len_near']
				param2['vert_align_offset2'] = param2['end_width_taper_length'] + param2['end_width_length']

				param2['taper_type'] = 'parabola'
				#SO added stuff
				param2['num_hole_size'] = 1
				param2['num_hole_cluster'] = 2
				param2['num_beams'] = param2['num_hole_size'] * param2['num_hole_cluster']

				param2['num_phc'] = 8  # number of PhCs per waveguide slub
				param2['bus_taper_len'] = 2000
				param2['bus_bend_radius'] = 6000
				param2['Len_bus_waveguide'] = (2*param2['mirror_len']+ (param2['num_phc'])*cavity_list_CaWO4[1]+(param2['num_phc']-1)*param2['middle_mirror_len']+2)*param2['aper_mir']
				param2['Len_bus_waveguide_vertical'] =(2*param2['mirror_len']+ param2['num_phc']*cavity_list_CaWO4[1]+(param2['num_phc']-1)*param2['middle_mirror_len']+2)*param2['aper_mir']
				param2['num_bus_waveguide'] = 1
				param2['ww'] = 350+35+7  # width bus waveguide
				param2['buffer_siliconpad_bendwaveguide'] = 2000
				param2['waveguide_with_end_mirror'] = param2['num_bus_waveguide']  # where to end the meander and put the mirror
				param2['num_phc_total'] = param2['num_phc'] * 2 * param2['waveguide_with_end_mirror']
				if param['vary_number_meander_waveguide'] is True:
					#left column has all even numbers of waveguides
					for i in range(6):
						if iY == 6-i and iX==0:
							param2['waveguide_with_end_mirror'] =1# i
							param2['num_phc_total'] = param2['num_phc'] * 2 * param2['waveguide_with_end_mirror']
					# right column has all odd numbers of waveguides
					for i in range(6):
						if iY == i and iX==1:
							param2['waveguide_with_end_mirror'] =1# i+1
							param2['num_phc_total'] = param2['num_phc'] * 2 * param2['waveguide_with_end_mirror']

				param2['length_total_meander'] = param2['waveguide_with_end_mirror'] * (param2['y_spacing_between_wabguides'] + 2 * param2['bus_bend_radius']) + param2['ww']/2+ param2['wg'] +param['beam_width']
				param2['spacing_phc-vertical_buffer'] = 1600
				param2['width_support_tather']=3500 #support tether added between phc
				param2['width_tether_bus_waveguide']=200
				param2['width_taper'] = 1000+40
				param2['length_wider_bus_wavguide_part'] = 600
				param2['width_taper_middle'] = 1000
				param2['radius_line_defect_circles']=100
				param2['aspect_line_defect_circles']=1.1379
				param2['a_defetc']=164.1*2
				param2['a_defetc_y'] = 124/2 +param2['radius_line_defect_circles']
				param2['width_tether_phc'] = 500

				param2['width_silicon_between_HF_cuts']=8000
				param2['width_HF_cuts']=1000
				param2['num_mirror_holes_end_meander']=18 # has to be a multiple of 6
				param2['separation_notches_meander_sides']=4500
				param2['width_wider_bus_wavguide_tether'] = 200

				if param2['meander'] is True:
					param2['num_tether_device'] = 20  # SC tether line to hold the device
					param2['extra_te_mir'] = 0
					param2['sym_hole_taper_q'] = True

				if param2['meander'] is False:
					param2['end_taper_n'] = 4
					param2['num_hole_size'] = 5  # 4 SC trying to decrease number of devices in each block
					param2['num_hole_cluster'] = 2  # 3 SC trying to decrease number of devices in each block
					param2['num_beams'] = param2['num_hole_size'] * param2['num_hole_cluster']
					param2['num_tether_device'] = 10  # SC tether line to hold the device
					param2['extra_te_mir'] = 10
					param2['sym_hole_taper_q'] = False
					if iY % 2 == 0:
						if param2['hoste'] == 1:
							hole_scale_list = hole_scale_list1
							param2['hole_group'] = 1
						if param2['hoste'] == 3:
							hole_scale_list = hole_scale_list1CaWO4
							param2['hole_group'] = 1

					if iY % 2 !=0:
						if param2['hoste'] == 1:
							hole_scale_list = hole_scale_list2
							param2['hole_group'] = 2
						if param2['hoste'] == 3:
							hole_scale_list = hole_scale_list2CaWO4
							param2['hole_group'] = 2

				if param2['meander'] is True:
					if iY % 2 == 0:
						hole_scale_list = hole_scale_list_meander
						param2['hole_group'] = 1

					if iY % 2 != 0:
						hole_scale_list = hole_scale_list_meander
						param2['hole_group'] = 2

				make_cavity_params_tm_refl(param2)

				param2['array_orig_x'] = chipCenterX + device_x_and_spacing_x_nm * iX

				y_offset = iY*num_cols*device_y_nm
				# param2['array_orig_y'] = chipCenterY + 190e3 * iY  #chipCenterY + 350e3 * iY -50e3*iY
				param2['array_orig_y'] = chipCenterY + iX*device_y_nm + y_offset
				# if iX==0:
				# 	param2['array_orig_y'] = chipCenterY + 350e3 * iY -50e3*iY
				print("writing sweep combo x" + str(param2['array_orig_x']) + "y" + str(param2['array_orig_y']))
				write_beams(beams, param2)
				param3 = copy(param2)

gdspy.gds_print('mirror_num_cavity_num_in_code_match_sim.gds', unit=1.0e-9, precision=1.0e-10)
	