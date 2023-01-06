import pickle
import matplotlib.pyplot as plt
import numpy as np

num_phc_per_GC_check = 8
OR_code_filename = 'matching_cavities_OR_code' + '.pickle'
old_code_filename = 'matching_cavities_jankycode' + '.pickle'
ellipse_dims_check_savename = 'ellipse_dims_top_x' + '.pickle'
with open(OR_code_filename, 'rb') as handle:
    aper_list_OR_design = pickle.load(handle)
    print("aper list loaded from object oriented" + str(aper_list_OR_design))

with open(old_code_filename, 'rb') as f:
    aper_list_old_design= pickle.load(f,encoding='latin1')
    print("aper list loaded confirmed" + str(aper_list_old_design))

with open(ellipse_dims_check_savename, 'rb') as hole_scalings_file:
    hole_scalings_x_top = pickle.load(hole_scalings_file)
    print("hole scalings x top" + str(hole_scalings_x_top))

[fig, ax] = plt.subplots()
x = np.arange(len(aper_list_OR_design))

ax.scatter(x,aper_list_OR_design, color='b',label='object oriented')
    # ax.plot(piezo_voltages[peaks], interferometer_signal[peaks], "o", mfc='coral', mec='black', mew=0.5, color='C0',
    #         linewidth=1, markersize=4, label='find peaks result')
x2 = np.arange(len(aper_list_old_design))
ax.scatter(x2,aper_list_old_design,marker='^', color='r', label='old design')

scalings_label = 'freq multiplex, top ellipse dims x'
axright = ax.twinx()
axright.scatter(x2,hole_scalings_x_top,marker='*', color='g', label=scalings_label)
axright.set_ylabel('hole scaling [nm]', fontsize=12)
    # ax.scatter(piezo_voltages[peaks], interferometer_signal[peaks], s=0.5, color='orange')
ax.legend(fontsize='x-small', loc='lower right')
axright.legend(fontsize='x-small', loc='center right')
    # ax.set_xlim([68, 72])
comparison_title_apers = 'aper match check for ' + str(num_phc_per_GC_check) + ' cavities per GC'
ax.set_title(comparison_title_apers, fontsize=7)
ax.set_xlabel('aper cell index', fontsize=12)
ax.set_ylabel('aper [nm]', fontsize=12)
savename_before_fits = np.array([comparison_title_apers + '.svg'])[0]
# ax.set_xlim([68, 72])
plt.savefig(savename_before_fits)