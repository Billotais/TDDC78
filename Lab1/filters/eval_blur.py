import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from subprocess import check_output

# List of images
images = {"im1", "im2", "im3", "im4"}

num_jobs = range(0, 7)  #7
num_jobs = [2**i for i in num_jobs]

filter_rad = range(1,8)  #8
filter_rad = [2**i + 1 for i in filter_rad]

for img in images:
	time_all = np.zeros((len(num_jobs), len(filter_rad)))

	for n_id, n in enumerate(num_jobs, 0):
		for r_id, r in enumerate(filter_rad, 0):			
			cmd_to_run = "mpirun -n " + str(n) + " ./blur " + str(r) + " ../" + img + ".ppm ../" + img + "_out_blur.ppm"

			program_out = check_output(cmd_to_run, shell=True)

			# Find the time
			tmp = program_out.split("Average filtering time: ")[1]
			run_time = tmp.split(" secs\n")[0]
			run_time = float(run_time)

			time_all[n_id, r_id] = run_time

	# plot
	for r_id, r in enumerate(filter_rad, 0):
		plot_label = "Radius = " + str(r)
		plt.plot(num_jobs, time_all[:, r_id], linewidth=2.0, label=plot_label)

	plt.ylabel('Time in secs')
	plt.xlabel('Num workers')
	plt.xscale('log', basex=2)
	
	xi = [2**n for n in range(0, len(num_jobs))]
	plt.xticks(num_jobs, xi)
	plt.title('Time vs Num workers ('+img+')')

	plt.legend(loc='best')

	plt.savefig('plots_blur/'+img+'_plot.png', bbox_inches='tight')
	plt.close()

	# Plot speedup
	speedup = time_all / time_all[0,:]
	speedup = 1 / speedup

	for r_id, r in enumerate(filter_rad, 0):
		plot_label = "Radius = " + str(r)
		plt.plot(num_jobs, speedup[:, r_id], linewidth=2.0, label=plot_label)

	plt.ylabel('Speedup')
	plt.xlabel('Num workers')
	plt.xscale('log', basex=2)
	plt.yscale('log', basey=2)
	
	xi = [2**n for n in range(0, len(num_jobs))]
	plt.xticks(num_jobs, xi)

	plt.yticks(num_jobs, xi)

	plt.title('Speedup vs Num workers ('+img+')')

	plt.legend(loc='best')

	plt.savefig('plots_blur/'+img+'_plot_speedup.png', bbox_inches='tight')
	plt.close()



