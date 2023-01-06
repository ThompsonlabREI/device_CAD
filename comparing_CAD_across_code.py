import pickle
import matplotlib.pyplot as plt
import numpy as np

num_phc_per_GC_check = 8
num_PhC_center_sweep = 5
phc_scale_min_param_sweep = 0.8
phc_scale_max_param_sweep = 1.1
phc_freq_scalings = np.linspace(phc_scale_min_param_sweep, phc_scale_max_param_sweep, num=num_PhC_center_sweep + 2)

OR_code_filename = 'matching_cavities_OR_code' + '.pickle'
old_code_filename = 'matching_cavities_jankycode' + '.pickle'

with open(OR_code_filename, 'rb') as handle:
    aper_list_OR_design = pickle.load(handle)
    print("aper list loaded from object oriented" + str(aper_list_OR_design))

with open(old_code_filename, 'rb') as f:
    aper_list_old_design= pickle.load(f,encoding='latin1')
    print("aper list loaded confirmed" + str(aper_list_old_design))

[fig, ax] = plt.subplots()
x = np.arange(len(aper_list_OR_design))

ax.scatter(x, aper_list_OR_design, color='b', label='object oriented')
    # ax.plot(piezo_voltages[peaks], interferometer_signal[peaks], "o", mfc='coral', mec='black', mew=0.5, color='C0',
    #         linewidth=1, markersize=4, label='find peaks result')
x2 = np.arange(len(aper_list_old_design))
ax.scatter(x2, aper_list_old_design, marker='^', color='r', label='old design')
scalings_label_top = 'freq multiplex, top ellipse dims x'
scalings_label_bot = 'freq multiplex, bot ellipse dims x'
axright = ax.twinx()

for center_freq_scaling_cav_index in range(num_PhC_center_sweep):
    ellipse_dims_check_savename_top = 'ellipse_dims_top_x_min_scale_' + str(round(phc_freq_scalings[center_freq_scaling_cav_index],4))  +'.pickle'
    ellipse_dims_check_savename_bottom = 'ellipse_dims_bottom_xmin_scale_' + str(round(phc_freq_scalings[center_freq_scaling_cav_index],4)) +'.pickle'

    with open(ellipse_dims_check_savename_top, 'rb') as hole_scalings_file_top:
        hole_scalings_x_top = pickle.load(hole_scalings_file_top)
        print("hole scalings x top" + str(hole_scalings_x_top))

    with open(ellipse_dims_check_savename_bottom, 'rb') as hole_scalings_file_bot:
        hole_scalings_x_bot = pickle.load(hole_scalings_file_bot)
        print("hole scalings x bot" + str(hole_scalings_x_bot))

    axright.scatter(x2, hole_scalings_x_top, marker='*', color='g', label=scalings_label_top)
    axright.scatter(x2, hole_scalings_x_bot, marker='*', color='black', label=scalings_label_bot)

axright.set_ylabel('hole scaling [nm]', fontsize=12)
    # ax.scatter(piezo_voltages[peaks], interferometer_signal[peaks], s=0.5, color='orange')
ax.legend(fontsize='x-small', loc='lower right')
axright.legend(fontsize='x-small', loc='center right')
    # ax.set_xlim([68, 72])
ax.set_xlabel('aper cell index', fontsize=12)
ax.set_ylabel('aper [nm]', fontsize=12)
comparison_title_apers = 'aper match check for ' + str(num_phc_per_GC_check) + ' cavities per GC, scale = ' + str(phc_scale_min_param_sweep) + ' to' + str(phc_scale_max_param_sweep)
ax.set_title(comparison_title_apers, fontsize=7)
savename_before_fits = np.array([comparison_title_apers + '.svg'])[0]
plt.savefig(savename_before_fits)
