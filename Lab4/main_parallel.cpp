#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "mpi.h"


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


	unsigned int time_stamp = 0, time_max;
	float pressure = 0;


	// parse arguments
	if(argc != 2) {
		fprintf(stderr, "Usage: %s simulation_time\n", argv[0]);
		fprintf(stderr, "For example: %s 10\n", argv[0]);
		exit(1);
	}

	time_max = atoi(argv[1]);



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
	
	
	/* Initialize */
	// 1. set the walls
	cord_t wall;
	wall.y0 = wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;
	
	int section_size = BOX_VERT_SIZE / mpi_size;
	float up_limit = myid * section_size; // Included
	float down_limit = ((myid+1) * section_size) -1; // Included
	down_limit = down_limit >= BOX_VERT_SIZE ? BOX_VERT_SIZE - 1 : down_limit;
	
	

	
    
	// 2. allocate particle bufer and initialize the particles
	list<pcord_t> particles;
	//pcord_t *particles = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
	//forward_list<bool> collisions();
	bool *collisions=(bool *)malloc(INIT_NO_PARTICLES*sizeof(bool) );

	srand( time(NULL) + 1234 );

	float r, a;
	for(int i=0; i<INIT_NO_PARTICLES/mpi_size; i++){
		// initialize random position
		pcord_t p;
		p.x = wall.x0 + rand1()*BOX_HORIZ_SIZE;
		p.y = up_limit + rand1()*section_size;

		// initialize random velocity
		r = rand1()*MAX_INITIAL_VELOCITY;
		a = rand1()*2*PI;
		p.vx = r*cos(a);
		p.vy = r*sin(a);
		particles.push_front(p);
	}


	unsigned int p, pp;

	/* Main loop */
	for (time_stamp=0; time_stamp<time_max; time_stamp++) { // for each time stamp

		
		init_collisions(collisions, INIT_NO_PARTICLES);
		
		//for(p=0; p<INIT_NO_PARTICLES; p++) { // for all particles
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
			
		vector<pcord_t> send_up;
		vector<pcord_t> send_down;
		
		auto it_p = particles.begin(); 
		while (it_p != particles.end())
		{
			if ((*it_p).y < up_limit && myid > 0) 
			{
				auto old_it = it_p;
				it_p++;
				send_up.push_back(*old_it);
				particles.erase(old_it);
			}
			else if ((*it_p).y > up_limit && myid < mpi_size-1) 
			{
				auto old_it = it_p;
				it_p++;
				send_down.push_back(*old_it);
				particles.erase(old_it);
			}
			else it_p++;
		}
		
		// Put in static arrays
		
		int up_size = send_up.size();
		int down_size = send_down.size();
		
		int up_receive_size=0;
		int down_receive_size=0;
		
		if (myid > 0) MPI_Send(&up_size, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD);
		if (myid < mpi_size - 1) MPI_Send(&down_size, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD);
		
		if (myid > 0) MPI_Recv(&up_receive_size, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD, NULL);
		if (myid < mpi_size - 1) MPI_Recv(&down_receive_size, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD, NULL);
		
		
		
		
		pcord_t* receive_up = (pcord_t*)malloc(up_receive_size*sizeof(pcord_t));
		pcord_t* receive_down = (pcord_t*)malloc(down_receive_size*sizeof(pcord_t));
		

		if(up_size > 0){ // Send it up
			pcord_t* up_array = &send_up[0];
			MPI_Send(up_array, up_size, pcord_t_mpi, myid-1, 1, MPI_COMM_WORLD);
		}
		
		if(down_size > 0) { 
			pcord_t* down_array = &send_down[0];
			MPI_Send(down_array, down_size, pcord_t_mpi, myid+1, 1, MPI_COMM_WORLD);
		 }
		 
		 
		 if(up_receive_size > 0){ // receive it up
			MPI_Recv(receive_up, up_receive_size, pcord_t_mpi, myid-1, 1, MPI_COMM_WORLD, NULL);
		}
		
		if(down_receive_size > 0) { 
			MPI_Recv(receive_down, down_receive_size, pcord_t_mpi, myid+1, 1, MPI_COMM_WORLD, NULL);
		 }
		
		
		// Add to particles
		for (int particle_id=0; particle_id < up_receive_size; particle_id++) {
			particles.push_front(receive_up[particle_id]);
		}
		
		for (int particle_id=0; particle_id < down_receive_size; particle_id++) {
			particles.push_front(receive_down[particle_id]);
		}
		
		free(receive_up);
		free(receive_down);
	}

	// Gather all pressure
	float pressure_all;
	MPI_Reduce(&pressure, &pressure_all, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (myid==0) {
		printf("Average pressure = %f\n", pressure_all / (WALL_LENGTH*time_max));
	}
	
	//free(particles);
	free(collisions);

	MPI_Finalize();
	
	return 0;

}

