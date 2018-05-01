import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from subprocess import check_output

num_jobs = range(0, 5)  
num_jobs = [2**i for i in num_jobs]

grid_size = [500, 1000, 2000, 4000, 8000, 16000]  

time_all = np.zeros((len(num_jobs), len(grid_size)))

for n_id, n in enumerate(num_jobs, 0):
	for g_id, g in enumerate(grid_size, 0):			
		cmd_to_run = "export OMP_NUM_THREADS=" + str(n) + "; echo \""+str(g)+"\" | ./laplsolv "

		program_out = check_output(cmd_to_run, shell=True)

		# Find the time
		tmp = program_out.split("Time: ")[1]
		run_time = tmp.split(" ::")[0]
		run_time = float(run_time.replace("\x00", ""))
		print("Runtime is:  " + str(run_time) + "  n jobs:  " + str(n) + "n grid:  " + str(g))

		time_all[n_id, g_id] = run_time

# plot
for g_id, g in enumerate(grid_size, 0):
	plot_label = "Grid Size = " + str(g)
	plt.plot(num_jobs, time_all[:, g_id], linewidth=2.0, label=plot_label)

plt.ylabel('Time in secs')
plt.xlabel('Num workers')
plt.xscale('log', basex=2)
#plt.yscale('log', basey=2)

xi = [2**n for n in range(0, len(num_jobs))]
plt.xticks(num_jobs, xi)
plt.title('Time vs Num workers')

plt.legend(loc='best')

plt.savefig('plots/plot.png', bbox_inches='tight')
plt.close()

# Plot speedup
speedup = time_all / (time_all[0,:] + 1e-9 )
speedup = 1 / speedup

for g_id, g in enumerate(grid_size, 0):
	plot_label = "Grid Size = " + str(g)
	plt.plot(num_jobs, speedup[:, g_id], linewidth=2.0, label=plot_label)

plt.ylabel('Speedup')
plt.xlabel('Num workers')
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)

xi = [2**n for n in range(0, len(num_jobs))]
plt.xticks(num_jobs, xi)

plt.yticks(num_jobs, xi)
plt.title('Speedup vs Num workers')

plt.legend(loc='best')

plt.savefig('plots/plot_speedup.png', bbox_inches='tight')
plt.close()



