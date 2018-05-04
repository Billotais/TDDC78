import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from subprocess import check_output

num_jobs = range(0, 5) # 6    
num_jobs = [2**i for i in num_jobs]

num_particles = [1000, 2000, 5000, 7500, 10000]
box_sizes = [5000, 6000, 7500, 10000]

base_num_particles = 5000
base_box_size = 5000

sim_time = 100


# Vary num particles
time_all = np.zeros((len(num_jobs), len(num_particles)))

pressure = np.zeros(len(num_particles))

for n_id, n in enumerate(num_jobs, 0):
	for g_id, g in enumerate(num_particles, 0):			
		cmd_to_run = "mpirun -n " + str(n) + " ./parallel " + str(sim_time) + "  " + str(base_box_size) + "  " + str(g)
		print("\n\nCommand is:  " + cmd_to_run)

		program_out = check_output(cmd_to_run, shell=True)

		print("Output is:  " + program_out)
		# Find the time
		tmp = program_out.split("Average filtering time: ")[1]
		run_time = tmp.split(" secs\n")[0]
		run_time = float(run_time)
		print("Runtime is:  " + str(run_time) + "  n jobs:  " + str(n) + "  n particles:  " + str(g))

		time_all[n_id, g_id] = run_time

		# Find pressure
		tmp = program_out.split("Average pressure = ")[1]
		new_pressure = tmp.split(" units\n")[0]
		new_pressure = float(new_pressure)

		pressure[g_id] += new_pressure 


pressure = pressure / len(num_jobs)


# plot
plt.plot(num_particles, pressure, linewidth=2.0)

plt.ylabel('Pressure')
plt.xlabel('Num particles')
# plt.xscale('log', basex=2)
# plt.yscale('log', basey=2)

#xi = [2**n for n in range(0, len(num_jobs))]
#plt.xticks(num_jobs, xi)
plt.title('Pressure vs Num particles')

plt.legend(loc='best')

plt.savefig('plots_t2/plot_pressure_vs_num_particles.png', bbox_inches='tight')
plt.close()

# plot
for g_id, g in enumerate(num_particles, 0):
	plot_label = "Num Particles = " + str(g)
	plt.plot(num_jobs, time_all[:, g_id], linewidth=2.0, label=plot_label)

plt.ylabel('Time in secs')
plt.xlabel('Num workers')
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)

xi = [2**n for n in range(0, len(num_jobs))]
plt.xticks(num_jobs, xi)
plt.title('Time vs Num workers')

plt.legend(loc='best')

plt.savefig('plots_t2/plot_num_particles.png', bbox_inches='tight')
plt.close()

# Plot speedup
speedup = time_all / (time_all[0,:] + 1e-9 )
speedup = 1 / speedup

for g_id, g in enumerate(num_particles, 0):
	plot_label = "Num Particles = " + str(g)
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

plt.savefig('plots_t2/plot_num_particles_speedup.png', bbox_inches='tight')
plt.close()






# Do the same for box size
# Vary num particles
time_all = np.zeros((len(num_jobs), len(box_sizes)))
pressure = np.zeros(len(box_sizes))

for n_id, n in enumerate(num_jobs, 0):
	for g_id, g in enumerate(box_sizes, 0):			
		cmd_to_run = "mpirun -n " + str(n) + " ./parallel " + str(sim_time) + "  " + str(g) + "  " + str(base_num_particles)

		program_out = check_output(cmd_to_run, shell=True)

		# Find the time
		tmp = program_out.split("Average filtering time: ")[1]
		run_time = tmp.split(" secs\n")[0]
		run_time = float(run_time)
		print("Runtime is:  " + str(run_time) + "  n jobs:  " + str(n) + "n grid:  " + str(g))

		time_all[n_id, g_id] = run_time

		# Find pressure
		tmp = program_out.split("Average pressure = ")[1]
		new_pressure = tmp.split(" units\n")[0]
		new_pressure = float(new_pressure)

		pressure[g_id] += new_pressure 


pressure = pressure / len(num_jobs)


# plot
inv_box_vol = [b**-2 for b in box_sizes]
plt.plot(inv_box_vol, pressure, linewidth=2.0)

plt.ylabel('Pressure')
plt.xlabel('1 / Box Volume')
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)

#xi = [2**n for n in range(0, len(box_sizes))]
x_t = [ "1/"+str(b**2) for b in box_sizes]
plt.xticks(inv_box_vol, x_t)
plt.title('Pressure vs 1 /(Box Volume)')

plt.legend(loc='best')

plt.savefig('plots_t2/plot_pressure_vs_box_vol.png', bbox_inches='tight')
plt.close()



# plot
for g_id, g in enumerate(box_sizes, 0):
	plot_label = "Box size = " + str(g)
	plt.plot(num_jobs, time_all[:, g_id], linewidth=2.0, label=plot_label)

plt.ylabel('Time in secs')
plt.xlabel('Num workers')
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)

xi = [2**n for n in range(0, len(num_jobs))]
plt.xticks(num_jobs, xi)
plt.title('Time vs Num workers')

plt.legend(loc='best')

plt.savefig('plots_t2/plot_box_size.png', bbox_inches='tight')
plt.close()

# Plot speedup
speedup = time_all / (time_all[0,:] + 1e-9 )
speedup = 1 / speedup

for g_id, g in enumerate(box_sizes, 0):
	plot_label = "Box size = " + str(g)
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

plt.savefig('plots_t2/plot_box_size_speedup.png', bbox_inches='tight')
plt.close()



