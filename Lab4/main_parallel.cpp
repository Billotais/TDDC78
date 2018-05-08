#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "mpi.h"
#include <VT.h>


#include "coordinate.h"
#include "definitions.h"
#include "physics.h"

#include <list>
#include <vector>



using namespace std;


//Feel free to change this program to facilitate parallelization.

float rand1(){
	return (float)( rand()/(float) RAND_MAX );
}

void init_collisions(bool *collisions, unsigned int max){
	for(unsigned int i=0;i<max;++i)
		collisions[i]=0;
}


int main(int argc, char** argv){


	unsigned int time_stamp = 0, time_max, total_num_particles, box_size, wall_length;
	double stime, etime;
	float pressure = 0;


	// parse arguments
	if(argc != 4) {
		fprintf(stderr, "Usage: %s simulation_time box_size total_num_particles \n", argv[0]);
		fprintf(stderr, "For example: %s 10 1000 500\n", argv[0]);
		exit(1);
	}

	time_max = atoi(argv[1]);
	box_size = atoi(argv[2]);
	total_num_particles = atoi(argv[3]);

	wall_length = (2.0*box_size+2.0*box_size);

	// Init MPI
	MPI_Init(NULL, NULL);
    int myid, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    // pcord_t Datatype for MPI
    pcord_t item;
	MPI_Datatype pcord_t_mpi;
	int block_length[] = {1, 1, 1, 1};
	MPI_Datatype block_types [] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
	MPI_Aint start, displ[4];

	MPI_Address(&item, &start);
	MPI_Address(&item.x, &displ[0]);
	MPI_Address(&item.y, &displ[1]);
	MPI_Address(&item.vx, &displ[2]);
	MPI_Address(&item.vy, &displ[3]);

	displ[0] -= start;
	displ[1] -= start;
	displ[2] -= start;
	displ[3] -= start;

	MPI_Type_struct(4, block_length, displ, block_types, &pcord_t_mpi);
	MPI_Type_commit(&pcord_t_mpi);
	MPI_Status status;
	
	/* Initialize */
	// 1. set the walls
	cord_t wall;
	wall.y0 = wall.x0 = 0;
	wall.x1 = box_size;
	wall.y1 = box_size;
	
	int section_size = box_size / mpi_size;
	float up_limit = myid * section_size; // Included
	float down_limit = ((myid+1) * section_size) -1; // Included
	down_limit = down_limit >= box_size ? box_size - 1 : down_limit;
	
	
	///////////////////////
	int vt_alloc;
	VT_funcdef("Alloc", VT_NOCLASS, &vt_alloc);
    VT_begin(vt_alloc);
    
	// 2. allocate particle bufer and initialize the particles
	list<pcord_t> particles;
	bool *collisions=(bool *)malloc(total_num_particles*sizeof(bool) );

	srand(time(NULL) + 1234 );

	float r, a;
	for(int i=0; i<total_num_particles/mpi_size; i++){
		// initialize random position
		pcord_t p;
		p.x = wall.x0 + rand1()*box_size;
		p.y = up_limit + rand1()*section_size;

		// initialize random velocity
		r = rand1()*MAX_INITIAL_VELOCITY;
		a = rand1()*2*PI;
		p.vx = r*cos(a);
		p.vy = r*sin(a);
		particles.push_front(p);
	}
	VT_end(vt_alloc);


	unsigned int p, pp;

    stime = MPI_Wtime();
    
    
    int vt_mainloop;
	VT_funcdef("Main Loop", VT_NOCLASS, &vt_mainloop);
    
    VT_begin(vt_mainloop);
    
	/* Main loop */
	for (time_stamp=0; time_stamp<time_max; time_stamp++) { // for each time stamp
		
		///////////////////////
		int vt_simulate;
		VT_funcdef("Simulate", VT_NOCLASS, &vt_simulate);
		VT_begin(vt_simulate);
		
		init_collisions(collisions, total_num_particles);
		
		//for(p=0; p<total_num_particles; p++) { // for all particles
		int p_index = 0;
		for (auto it_p = particles.begin(); it_p != particles.end(); ++it_p) {
			if(collisions[p_index]) continue; // If this particle has already collided qith a previous particle, ignore

			/* check for collisions with all next particles*/
			int pp_index = 0;
			
			for(auto it_pp = next(it_p,1);it_pp != particles.end(); it_pp++){
				if(collisions[p_index + pp_index]) continue; // If the particle we want to colide with has already collided, ignore
				float t=collide(&(*it_p), &(*it_pp));
				if(t!=-1){ // collision
					collisions[p_index]=collisions[pp_index]=1;
					interact(&(*it_p), &(*it_pp), t);
					break; // only check collision of two particles
				}
				pp_index++;
			}
			p_index++;
		}
		p_index = 0;
		// move particles that has not collided with another
		for(auto it_p = particles.begin(); it_p != particles.end(); ++it_p)
			if(!collisions[p_index]){
				feuler(&(*it_p), 1);

				/* check for wall interaction and add the momentum */
				pressure += wall_collide(&(*it_p), wall);
			}
		
		
		VT_end(vt_simulate);
		/////////////////////	
		
		
		///////////////////////
		int vt_tosend;
		VT_funcdef("Find to send", VT_NOCLASS, &vt_tosend);
		VT_begin(vt_tosend);
		
		
		// We now find the particls we have to send to other threads
		vector<pcord_t> send_up;
		vector<pcord_t> send_down;
		
		auto it_p = particles.begin(); 
		while (it_p != particles.end())
		{
			if ((*it_p).y < up_limit && myid > 0) // If the particle goes up
			{
				auto old_it = it_p;
				it_p++;
				send_up.push_back(*old_it);
				particles.erase(old_it);
			}
			else if ((*it_p).y > up_limit && myid < mpi_size-1) // f the particle goes down
			{
				auto old_it = it_p;
				it_p++;
				send_down.push_back(*old_it);
				particles.erase(old_it);
			}
			else it_p++;
		}
		
		
		VT_end(vt_tosend);
		///////////////////
		int vt_send_size;
		VT_funcdef("Send sizes", VT_NOCLASS, &vt_send_size);
		VT_begin(vt_send_size);
		
		// Put in static arrays
		
		int up_size = send_up.size();
		int down_size = send_down.size();
		
		int up_receive_size=0;
		int down_receive_size=0;
		
		
		// Send the sizes
		if (myid > 0) MPI_Send(&up_size, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD);
		if (myid < mpi_size - 1) MPI_Send(&down_size, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD);
		
		if (myid > 0) MPI_Recv(&up_receive_size, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, &status);
		if (myid < mpi_size - 1) MPI_Recv(&down_receive_size, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD, &status);
		
		
		VT_end(vt_send_size);
		////////////////////
		
		
		
		pcord_t* receive_up = (pcord_t*)malloc(up_receive_size*sizeof(pcord_t));
		pcord_t* receive_down = (pcord_t*)malloc(down_receive_size*sizeof(pcord_t));
		

		///////////////////////
		int vt_send_part;
		VT_funcdef("Send particles", VT_NOCLASS, &vt_send_part);
		VT_begin(vt_send_part);
		
		// Send the data
		
		if(up_size > 0){ // Send it up
			pcord_t* up_array = &send_up[0];
			MPI_Send(up_array, up_size, pcord_t_mpi, myid-1, 1, MPI_COMM_WORLD);
		}
		
		if(down_size > 0) { 
			pcord_t* down_array = &send_down[0];
			MPI_Send(down_array, down_size, pcord_t_mpi, myid+1, 1, MPI_COMM_WORLD);
		}
		 
		// receive the data
		if(up_receive_size > 0){ // receive it up
			MPI_Recv(receive_up, up_receive_size, pcord_t_mpi, myid-1, 1, MPI_COMM_WORLD, &status);
		}
		
		if(down_receive_size > 0) { 
			MPI_Recv(receive_down, down_receive_size, pcord_t_mpi, myid+1, 1, MPI_COMM_WORLD, &status);
		}
		VT_end(vt_send_part);
		//////////////////////
		
		///////////////////////
		int vt_update_part;
		VT_funcdef("Send particles", VT_NOCLASS, &vt_update_part);
		VT_begin(vt_update_part);
		
		
		// Add to particles
		for (int particle_id=0; particle_id < up_receive_size; particle_id++) {
			particles.push_front(receive_up[particle_id]);
		}
		
		for (int particle_id=0; particle_id < down_receive_size; particle_id++) {
			particles.push_front(receive_down[particle_id]);
		}
		
		free(receive_up);
		free(receive_down);
		VT_end(vt_update_part);
		//////////////////////
	}
	
	VT_end(vt_mainloop);
	///////////////////////



	///////////////////////
	int vt_gather_pressure;
	VT_funcdef("Gather pressure", VT_NOCLASS, &vt_gather_pressure);
    VT_begin(vt_gather_pressure);
    
    
	// Gather all pressure
	float pressure_all;
	MPI_Reduce(&pressure, &pressure_all, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	VT_end(vt_gather_pressure);
	///////////////////////
	
	etime = MPI_Wtime();

	if (myid==0) {
		printf("Average pressure = %f units\n\n", pressure_all / (wall_length*time_max));
	}
	///////////////////////
	int vt_gather_time;
	VT_funcdef("Gather time", VT_NOCLASS, &vt_gather_time);
    VT_begin(vt_gather_time);
    
	double total_time = etime - stime;
	double total_time_all;

	// Get time from all workers
	MPI_Reduce(&total_time, &total_time_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	total_time_all = total_time_all / (float)mpi_size;
	
	VT_end(vt_gather_time);
	///////////////////////
	/* write result */
	if (myid == 0) {
		printf("Average simulation time: %g secs\n", total_time_all) ;
	}
	//free(particles);
	free(collisions);

	MPI_Finalize();
	
	return 0;

}

