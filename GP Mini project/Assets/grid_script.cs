using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using System.Linq;
using System.Threading;

public class grid_script : MonoBehaviour
{
    public GameObject grid_cell;
    public GameObject particle;
    public float timestep = 0.02f;

    public const int rho = 1;

    private float time;
    private Renderer[,,] renderers_grid;
    private Renderer[] renderers_particle;
    private Thread[] particlePool;

    private bool is_running;
    
    enum state {air, water, solid}

    enum face { x_face, y_face, z_face }

    //MAC Grid
    float[,,] pressure ;
    state[,,] cell_state;

    float[,,] ux;
    float[,,] uy;
    float[,,] uz;
    float[,,] ux_new;
    float[,,] uy_new;
    float[,,] uz_new;

    //marker particles
    Vector3[] particles ; //particles positions. we also need velocities

    int sim_step = 0;
    int num_particles = 0;

    List<Vector3> water_cells;


    //parameters to alter 5x5x5
    public int grid_dim_x = 5;
    public int grid_dim_y = 5;
    public int grid_dim_z = 5;

    //determines which cells the marker particles are initialy placed in 3x3x3
    public int waterfall_dim_x = 1;
    public int waterfall_dim_y = 1;
    public int waterfall_dim_z = 1;

   public int part_spawn = 1000; //for each cell 10/20
    static int jacobi_iterations = 30; //number of iterations for the Jacobi solver
    private const bool threading = false;
    private const bool clamp_part_pos = true; //when true, sets particle position to border when trying to escape >:D
    public const float p_atm = 0.0f; //0 is 'normal' and maybe best?


    // Use this for initialization
    void Start()
    {
        time = 0;
        particles = new Vector3[waterfall_dim_x * waterfall_dim_y * waterfall_dim_z * part_spawn];
        renderers_grid = new Renderer[grid_dim_x, grid_dim_y, grid_dim_z];
        renderers_particle = new Renderer[waterfall_dim_x * waterfall_dim_y * waterfall_dim_z * part_spawn];

        ux = new float[grid_dim_x + 1, grid_dim_y, grid_dim_z];
        uy = new float[grid_dim_x, grid_dim_y + 1, grid_dim_z];
        uz = new float[grid_dim_x, grid_dim_y, grid_dim_z + 1];
        pressure = new float[grid_dim_x, grid_dim_y, grid_dim_z];
        cell_state = new state[grid_dim_x, grid_dim_y, grid_dim_z];
        //grid visualisation
        for (int x = 0; x < grid_dim_x + 1; x++)
        {
            for (int y = 0; y < grid_dim_y + 1; y++)
            {
                for (int z = 0; z < grid_dim_z + 1; z++)
                {
                    if (x < grid_dim_x && y < grid_dim_y && z < grid_dim_z)
                    {
                        Vector3 pos = new Vector3(x, y, z);
                        GameObject gc = Instantiate(grid_cell, pos, new Quaternion());
                        Renderer r = gc.GetComponent<Renderer>();
                        renderers_grid[x, y, z] = r;
                    }
                }
            }
        }
    }

    void spawnParticles(int sim_step)
    {
        num_particles += waterfall_dim_x * waterfall_dim_y * waterfall_dim_z;
        for (int x = 0; x < grid_dim_x; x++)
        {
            for (int y = 0; y < grid_dim_y; y++)
            {
                for (int z = 0; z < grid_dim_z; z++)
                {
                    if (x < waterfall_dim_x && y >= (grid_dim_y - waterfall_dim_y) && y < grid_dim_y && z < waterfall_dim_z)
                    {
                        Vector3 pos = new Vector3(x, y, z);
                        GameObject p = Instantiate(particle, pos, new Quaternion());
                        Renderer pr = p.GetComponent<Renderer>();
                        int ind = (sim_step - 1) * waterfall_dim_x * waterfall_dim_y * waterfall_dim_z + x * waterfall_dim_y * waterfall_dim_z + (y - (grid_dim_y - waterfall_dim_y)) * waterfall_dim_z + z;
                        particles[ind] = pos;
                        renderers_particle[ind] = pr;
                    }
                }
            }
        }
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKeyDown("s"))
        {
            updateTime();
            simStep();
        }

        if (Input.GetKeyDown("space"))
        {
            if (is_running)
                is_running = false;
            else
                is_running = true;

        }

        if (is_running)
        {
            updateTime();            
            simStep();
        }
    }

    void updateTime()
    {
        setTimeStep();
        time +=timestep;
        print("time: " + time);
    }

    void setTimeStep()
    {
        float max = getMaxVelocity();
        timestep = Math.Max(0.01f,Math.Min(1/ max,0.1f));
    }

    float getMaxVelocity()
    {
        float max = 0.0f;
        for (int x = 0; x < grid_dim_x + 1; x++)
        {
            for (int y = 0; y < grid_dim_y + 1; y++)
            {
                for (int z = 0; z < grid_dim_z + 1; z++)
                {
                    if (y < grid_dim_y && z < grid_dim_z)
                        if (Math.Abs(ux[x, y, z]) > max)
                            max = Math.Abs(ux[x, y, z]);

                    if (x < grid_dim_x && z < grid_dim_z)
                        if (Math.Abs(uy[x, y, z]) > max)
                            max = Math.Abs(uy[x, y, z]);

                    if (x < grid_dim_x && y < grid_dim_y)
                        if (Math.Abs(uz[x, y, z]) > max)
                            max = Math.Abs(uz[x, y, z]);
                }
            }
        }
        return max;
    }
    void simStep()
    {
        //initial velocities?
        sim_step++;

        if (sim_step <= part_spawn)
            spawnParticles(sim_step);

        //marks each cell as either water or air
        markAirWaterCells();

        //3a advect velocities for each dimension -> behavior at edges is still a thing
        advect(); //forward euler advection
        //advectRK2(); //RK2 advection

        //3b body forces -> add gravity to velocity
        addGravity();
        float[,,] divergence2 = get3DDivergence();

        //if is cell has air, overwrite pressure value to 0
        //setAirPressure();

        //3g set boundary velocities to zero if they point into the walls
        setBoundary();

        //3d solve pressure -> get pressure for each cell by solving system of linear equations
        solvePressure();

        //3e adjust velocity field acording to pressure
        adjustVel();

        float[,,] divergence = get3DDivergence();

        //3f move air velocities
        advectAir();

        //3g set boundary velocities to zero if they point into the walls
        setBoundary();

        //4 move particles according to interpolated velocities
        moveParticles();

        adjustWaterColor();
    }

    void adjustWaterColor()
    {
        for (int x = 0; x < grid_dim_x; x++)
        {
            for (int y = 0; y < grid_dim_y; y++)
            {
                for (int z = 0; z < grid_dim_z; z++)
                {
                    if (cell_state[x, y, z] == state.water)
                    {
                        renderers_grid[x, y, z].material.color = new Color(0, 0.4f, 1, 0.5f);
                        renderers_grid[x, y, z].enabled = true;
                    }
                    else
                    {
                        renderers_grid[x, y, z].material.color = new Color(0, 1, 1, 0.1f);
                        renderers_grid[x, y, z].enabled = false;
                    }
                }
            }
        }
    }
    void setAirPressure()
    {
        for (int x = 0; x < grid_dim_x; x++)
        {
            for (int y = 0; y < grid_dim_y; y++)
            {
                for (int z = 0; z < grid_dim_z; z++)
                {
                    if(cell_state[x,y,z] == state.air)
                        pressure[x, y, z] = 0; //overwrite presure solve -> set value to 0
                }
            }
        }
    }

    //for all velocity components that border walls, set to 0 if it points into the wall
    void setBoundary()
    {
        for (int x = 0; x < grid_dim_x + 1; x++)
        {
            for (int y = 0; y < grid_dim_y + 1; y++)
            {
                for (int z = 0; z < grid_dim_z + 1; z++)
                {
                    if (z < grid_dim_z && y < grid_dim_y && x == 0 && ux[x, y, z] < 0) //left wall
                        ux[x, y, z] = 0;

                    if (z < grid_dim_z && y < grid_dim_y && x == grid_dim_x && ux[x, y, z] > 0) //right wall
                        ux[x, y, z] = 0;
        
                    if (z < grid_dim_z && x < grid_dim_x && y == 0 && uy[x, y, z] < 0) //down wall
                        uy[x, y, z] = 0;

                    if (z < grid_dim_z && x < grid_dim_x && y == grid_dim_y && uy[x, y, z] > 0) //up wall
                        uy[x, y, z] = 0;

                    if (x < grid_dim_x && y < grid_dim_y && z == 0 && uz[x, y, z] < 0) //front wall
                        uz[x, y, z] = 0;

                    if (x < grid_dim_x && y < grid_dim_y && z == grid_dim_z && uz[x, y, z] > 0) //back wall
                        uz[x, y, z] = 0;
                }
            }
        }
    }

    void markAirWaterCells()
    {
        cellStateSetZero();
        for (int p = 0; p < num_particles; p++)
        {
            Vector3 pos = particles[p];
            //print(pos);
            int ind_x = (int)Math.Floor(pos.x + 0.5);
            int ind_y = (int)Math.Floor(pos.y + 0.5);
            int ind_z = (int)Math.Floor(pos.z + 0.5);

            if (ind_x == grid_dim_x)
                ind_x = grid_dim_x - 1;

            if (ind_y == grid_dim_y)
                ind_y = grid_dim_y - 1;

            if (ind_z == grid_dim_z)
                ind_z = grid_dim_z - 1;

            //print(ind_x + " " + ind_y + " " + ind_z);
            cell_state[ind_x, ind_y, ind_z] = state.water;//water

        }
        //print(cell_state[grid_dim_x - 1,grid_dim_y - 1,grid_dim_z - 1].ToString());
    }

    //returns true if the given face has water on either side
    bool stateIsWater(int x, int y, int z, face f)
    {
        int ind_x = x;
        int ind_y = y;
        int ind_z = z;

        int ind_x2 = x - 1;
        int ind_y2 = y - 1;
        int ind_z2 = z - 1;

        if (ind_x > grid_dim_x - 1)
            ind_x = grid_dim_x - 1;

        if (ind_y > grid_dim_y - 1)
            ind_y = grid_dim_y - 1;

        if (ind_z > grid_dim_z - 1)
            ind_z = grid_dim_z - 1;

        if (ind_x2 < 0)
            ind_x2 = 0;

        if (ind_y2 < 0)
            ind_y2 = 0;

        if (ind_z2 < 0)
            ind_z2 = 0;

        state s1 = state.air, s2;
        s2 = cell_state[ind_x, ind_y, ind_z];
        //print("cell: " + s1 + ": " + s2)
        if (f == face.x_face)
        {
            s1 = cell_state[ind_x2, ind_y, ind_z];
        }

        else if (f == face.y_face)
        {
            s1 = cell_state[ind_x, ind_y2, ind_z];
        }
        else if (f == face.z_face)
        {
            s1 = cell_state[ind_x, ind_y, ind_z2];
        }
        else
            print("face enum error");

        if (s1 == state.water || s2 == state.water)
            return true;
        else
            return false;

    }

    //solve the pressure equations
    void solvePressure()
    {
        float[] divergence = getDivergence();
        float[,] AMatrix = getAMatrix();

        adjustDivergence(divergence);
        JacobiMethod(AMatrix, divergence, jacobi_iterations); //pressure solve
    }

    //adjusts the divergence for the pressure solve
    void adjustDivergence(float[] divergence)
    {
        for (int i = 0; i < divergence.Length; i++)
        {
            Vector3 pos = water_cells.ElementAt(i);
            List<Vector3> neighbors = getNeighbors((int)pos.x, (int)pos.y, (int)pos.z);
            int air_neighbors = 0;
            for(int j = 0; j < neighbors.Count; j++)
            {
                Vector3 pos2 = neighbors.ElementAt(j);
                if (cell_state[(int)pos2.x, (int)pos2.y, (int)pos2.z] == state.air)
                    air_neighbors++;
            }

            divergence[i] = (rho / timestep) * divergence[i] - air_neighbors * p_atm;
        }

    }

    //adds gravity to all y velocity components that border a fluid cell
    void addGravity()
    {
        for (int x = 0; x < grid_dim_x; x++)
        {
            for (int y = 0; y < grid_dim_y + 1; y++)
            {
                for (int z = 0; z < grid_dim_z; z++)
                {
                    if (stateIsWater(x, y, z, face.y_face))
                       uy[x, y, z] -= 9.81f * timestep;
                }
            }
        }
    }

    //adjusts the velocities based on the pressure
    void adjustVel()
    {
        for (int x = 0; x < grid_dim_x + 1; x++)
            for (int y = 0; y < grid_dim_y + 1; y++)
                for (int z = 0; z < grid_dim_z + 1; z++)
                {
                    if (z < grid_dim_z && y < grid_dim_y)
                    {
                        if (!(x == 0 || x == grid_dim_x)) //left wall and right wall
                        {
                            float pressure_x = pressure[x - 1, y, z];
                            float pressure_x2 = pressure[x, y, z];

                            if (stateIsWater(x, y, z, face.x_face))
                                ux[x, y, z] = ux[x, y, z] - (timestep / rho) * (pressure_x2 - pressure_x);
                        }
                    }

                    if (x < grid_dim_x && z < grid_dim_z)
                    {
                        if (!(y == 0 || y == grid_dim_y)) //down wall and up wall
                        {
                            float pressure_y = pressure[x, y - 1, z];
                            float pressure_y2 = pressure[x, y, z];

                            if (stateIsWater(x, y, z, face.y_face))
                                uy[x, y, z] = uy[x, y, z] - (timestep / rho) * (pressure_y2 - pressure_y);
                        }
                    }

                    if (x < grid_dim_x && y < grid_dim_y)
                    {
                        if (!(z == 0 || z == grid_dim_z)) //back and front wall
                        {
                            float pressure_z = pressure[x, y, z - 1];
                            float pressure_z2 = pressure[x, y, z];

                            if (stateIsWater(x, y, z, face.z_face))
                                uz[x, y, z] = uz[x, y, z] - (timestep / rho) * (pressure_z2 - pressure_z);
                        }
                    }
                }
    }

    void pressureSetZero()
    {
        for (int x = 0; x < grid_dim_x; x++)
            for (int y = 0; y < grid_dim_y; y++)
                for (int z = 0; z < grid_dim_z; z++)
                    pressure[x, y, z] = 0;
    }

    void cellStateSetZero()
    {
        for (int x = 0; x < grid_dim_x; x++)
            for (int y = 0; y < grid_dim_y; y++)
                for (int z = 0; z < grid_dim_z; z++)
                    cell_state[x, y, z] = 0;
    }

    //with the A matrix and adjusted divergence as input, solves for pressure
    void JacobiMethod(float[,] inputMatrix, float[] expectedOutcome, int iterations)
    {
        pressureSetZero();

        for (int p = 0; p < iterations; p++)
        {
            for (int i = 0; i < expectedOutcome.Length; i++)
            {
                int ix = (int) water_cells.ElementAt(i).x;
                int iy = (int) water_cells.ElementAt(i).y;
                int iz = (int) water_cells.ElementAt(i).z;

                float sigma = 0f;
                for (int j = 0; j < expectedOutcome.Length; j++)
                {
                    int jx = (int)water_cells.ElementAt(j).x;
                    int jy = (int)water_cells.ElementAt(j).y;
                    int jz = (int)water_cells.ElementAt(j).z;

                    if (j != i)
                    {
                        sigma += inputMatrix[i, j] * pressure[jx, jy, jz];
                    }
                    pressure[ix, iy, iz] = (expectedOutcome[i] - sigma) / inputMatrix[i, i];
                }
            }
        }
    }

    void AverageMethod(float[,,] divergence, int iterations)
    {
        float[,,] tempPressure = (float[,,])pressure.Clone();
        float xn = 0, xp = 0, yn = 0, yp = 0, zn = 0, zp = 0;

        for (int i = 0; i < iterations; i++)
        {
            for (int x = 0; x < grid_dim_x; x++)
            {
                for (int y = 0; y < grid_dim_y; y++)
                {
                    for (int z = 0; z < grid_dim_z; z++)
                    {

                        if (x == 0) //left wall?
                            xn = pressure[x, y, z]; 
                        else
                            xn = pressure[x - 1, y, z];

                        if (x == grid_dim_x)
                            xp = pressure[x - 1, y, z]; 
                        else
                            xp = pressure[x, y, z];

                        if (y == 0) //up wall
                            yn = pressure[x, y, z];
                        else
                            yn = pressure[x, y - 1, z];

                        if (y == grid_dim_y)
                            yp = pressure[x, y - 1, z];
                        else
                            yp = pressure[x, y, z];

                        if (z == 0) //back wall
                            zn = pressure[x, y, z];
                        else
                            zn = pressure[x, y, z - 1];

                        if (z == grid_dim_z) //front wall
                            zp = pressure[x, y, z - 1];
                        else
                            zp = pressure[x, y, z];

                        if (cell_state[x, y, z] == state.water)
                            tempPressure[x, y, z] = (xn + xp + yn + yp + zn + zp - divergence[x, y, z]) / 6;
                        else
                            tempPressure[x, y, z] = 0;
                    }
                }
            }
            pressure = (float[,,])tempPressure.Clone();
        }
    }


    //get the divergence values for all water cells
    float[] getDivergence()
    {
        water_cells = getWaterCells();
        float[] divergence = new float[water_cells.Count];

        for (int i = 0; i < water_cells.Count; i++)
        {
            int x = (int)water_cells.ElementAt(i).x;
            int y = (int)water_cells.ElementAt(i).y;
            int z = (int)water_cells.ElementAt(i).z;

            float ux_term = (ux[x + 1, y, z] - ux[x, y, z]);
            float uy_term = (uy[x, y + 1, z] - uy[x, y, z]);
            float uz_term = (uz[x, y, z + 1] - uz[x, y, z]);

            divergence[i] = (ux_term + uy_term + uz_term);
        }

        return divergence;
    }

    float[,,] get3DDivergence()
    {
        float[,,] divergence = new float[grid_dim_x, grid_dim_y, grid_dim_z];
        for (int x = 0; x < grid_dim_x; x++)
        {
            for (int y = 0; y < grid_dim_y; y++)
            {
                for (int z = 0; z < grid_dim_z; z++)
                {
                    float ux_term = (ux[x + 1, y, z] - ux[x, y, z]);
                    float uy_term = (uy[x, y + 1, z] - uy[x, y, z]);
                    float uz_term = (uz[x, y, z + 1] - uz[x, y, z]);

                    divergence[x, y, z] = (ux_term + uy_term + uz_term);
                }
            }
        }
        return divergence;
    }

    bool hasWaterNeigbor(int x, int y, int z)
    {
        List<Vector3> neighbors = getNeighbors(x, y, z);
        for (int i = 0; i < neighbors.Count; i++)
        {
            int x2 = (int)neighbors.ElementAt(i).x;
            int y2 = (int)neighbors.ElementAt(i).y;
            int z2 = (int)neighbors.ElementAt(i).z;
            if (cell_state[x2, y2, z2] == state.water)
                return true;
        }
        return false;
    }

    void advectAir()
    {
        for (int x = 0; x < grid_dim_x; x++)
        {
            for (int y = 0; y < grid_dim_y; y++)
            {
                for (int z = 0; z < grid_dim_z; z++)
                {
                    if (cell_state[x, y, z] == state.air) //not acounting for solid cells
                    {
                        if(hasWaterNeigbor(x, y, z))
                        {
                            if (!stateIsWater(x, y, z, face.x_face))
                                ux[x, y, z] = getAverageOfWaterNeighbors(x, y, z, face.x_face);
                            if (!stateIsWater(x + 1, y, z, face.x_face))
                                ux[x + 1, y, z] = getAverageOfWaterNeighbors(x + 1, y, z, face.x_face, true);

                            if (!stateIsWater(x, y, z, face.y_face))
                                uy[x, y, z] = getAverageOfWaterNeighbors(x, y, z, face.y_face);
                            if (!stateIsWater(x, y + 1, z, face.y_face))
                                uy[x, y + 1, z] = getAverageOfWaterNeighbors(x, y + 1, z, face.y_face, true);

                            if (!stateIsWater(x, y, z, face.z_face))
                                uz[x, y, z] = getAverageOfWaterNeighbors(x, y, z, face.z_face);
                            if (!stateIsWater(x + 1, y, z, face.z_face))
                                uz[x, y, z + 1] = getAverageOfWaterNeighbors(x, y, z + 1, face.z_face, true);
                        }
                    }
                }
            }
        }
    }

    //returns water and air neighbors
    List<Vector3> getNeighbors(int x, int y, int z)
    {
        List<Vector3> neighbors = new List<Vector3>();
        for (int x2 = x - 1; x2 <= x + 1; x2 += 2)
        {
            if (x2 >= 0 && x2 < grid_dim_x && y >= 0 && y < grid_dim_y && z >= 0 && z < grid_dim_z)
                neighbors.Add(new Vector3(x2, y, z));
        }

        for (int y2 = y - 1; y2 <= y + 1; y2 += 2)
        {
            if (x >= 0 && x < grid_dim_x && y2 >= 0 && y2 < grid_dim_y && z >= 0 && z < grid_dim_z)
                neighbors.Add(new Vector3(x, y2, z));
        }

        for (int z2 = z - 1; z2 <= z + 1; z2 += 2)
        {
            if (x >= 0 && x < grid_dim_x && y >= 0 && y < grid_dim_y && z2 >= 0 && z2 < grid_dim_z)
                neighbors.Add(new Vector3(x, y, z2));
        }
        return neighbors;
    }

    float getAverageOfWaterNeighbors(int x, int y, int z, face f, bool right_face = false)
    {
        int num_water_neighbors = 0;
        float sum = 0;
        List<Vector3> neighbors = new List<Vector3>();

        if (right_face)
        {
            if (f == face.x_face)
                neighbors = getNeighbors(x - 1, y, z);
            else if (f == face.y_face)
                neighbors = getNeighbors(x, y - 1, z);
            else if (f == face.z_face)
                neighbors = getNeighbors(x, y, z - 1);
        }
        else
            neighbors = getNeighbors(x, y, z);

        for(int i = 0; i < neighbors.Count; i++)
        {
            int x2 = (int)neighbors.ElementAt(i).x;
            int y2 = (int)neighbors.ElementAt(i).y;
            int z2 = (int)neighbors.ElementAt(i).z;

            if (cell_state[x2, y2, z2] == state.water)
            {
                num_water_neighbors++;
                if (f == face.x_face)
                {
                    if (right_face)
                        sum += ux[x2 + 1, y2, z2];
                    else
                        sum += ux[x2, y2, z2];
                }
                if (f == face.y_face)
                {
                    if (right_face)
                        sum += uy[x2, y2 + 1, z2];
                    else
                        sum += uy[x2, y2, z2];

                }
                if (f == face.z_face)
                {
                    if (right_face)
                        sum += uz[x2, y2, z2 + 1];
                    else
                        sum += uz[x2, y2, z2];
                }
            }
        }
        return sum / num_water_neighbors;
    }

        //update velocity values of MAC Grid?
    void advect()
    {
        ux_new = (float[,,]) ux.Clone();
        uy_new = (float[,,]) uy.Clone();
        uz_new = (float[,,]) uz.Clone();

        //for each grid cell
        for (int x = 0; x < grid_dim_x + 1; x++)
        {
            for (int y = 0; y < grid_dim_y + 1; y++)
            {
                for (int z = 0; z < grid_dim_z + 1; z++)
                {
                    //advection uses forward euler integrator for now TODO: upgrade to better integrator
                    if (y < grid_dim_y && z < grid_dim_z && stateIsWater(x, y, z, face.x_face))
                        advectUx2(x, y, z);
                    if (x < grid_dim_x && z < grid_dim_z && stateIsWater(x, y, z, face.y_face))
                        advectUy2(x, y, z);
                    if (x < grid_dim_x && y < grid_dim_y && stateIsWater(x, y, z, face.z_face))
                        advectUz2(x, y, z);
                }
            }
        }
        ux = ux_new;
        uy = uy_new;
        uz = uz_new;
    }

    void advectRK2()
    {
        ux_new = (float[,,])ux.Clone();
        uy_new = (float[,,])uy.Clone();
        uz_new = (float[,,])uz.Clone();

        //for each grid cell
        for (int x = 0; x < grid_dim_x; x++)
        {
            for (int y = 0; y < grid_dim_y; y++)
            {
                for (int z = 0; z < grid_dim_z; z++)
                {
                    Vector3 new_location = traceParticle(x, y, z);
                    Vector3 new_speed = getVelocity(new_location.x, new_location.y, new_location.z);
                    if(stateIsWater(x,y,z,face.x_face))
                        ux_new[x, y, z] = new_speed.x;
                    if(stateIsWater(x,y,z,face.y_face))
                        uy_new[x, y, z] = new_speed.y;
                    if(stateIsWater(x,y,z,face.z_face))
                        uz_new[x, y, z] = new_speed.z;
                }
            }
        }
        ux = ux_new;
        uy = uy_new;
        uz = uz_new;
    }

    Vector3 traceParticle(float x, float y, float z)
    {
        Vector3 V = getVelocity(x, y, z);
        V = getVelocity(x + 0.5f * timestep * V.x, y + 0.5f * timestep * V.y, z + 0.5f * timestep * V.z);
        Vector3 point = new Vector3(x, y, z);
        return point + timestep * V;
    }

    Vector3 getVelocity(float x, float y, float z)
    {
        Vector3 V;
        V.x = getInterpolatedXValue(x, y - 0.5f, z - 0.5f);
        V.y = getInterpolatedYValue(x - 0.5f, y, z - 0.5f);
        V.z = getInterpolatedZValue(x - 0.5f, y - 0.5f, z);
        return V;
    }

    float getInterpolatedXValue(float x, float y, float z)
    {
        int i = (int)Math.Floor(x);
        int j = (int)Math.Floor(y);
        int k = (int)Math.Floor(z);
        List<float> values = new List<float>();

        if (indexExists(i, j, k))
            values.Add((i + 1 - x) * (j + 1 - y) * (k + 1 - z) * ux[i, j, k]);

        if (indexExists(i + 1, j, k))
            values.Add((x - i) * (j + 1 - y) * (k + 1 - z) * ux[i + 1, j, k]);

        if (indexExists(i, j + 1, k))
            values.Add((i + 1 - x) * (y - j) * (k + 1 - z) * ux[i, j + 1, k]);

        if (indexExists(i + 1, j + 1, k))
            values.Add((x - i) * (y - j) * (k + 1 - z) * ux[i + 1, j + 1, k]);

        if (indexExists(i, j, k + 1))
            values.Add((i + 1 - x) * (j + 1 - y) * (z - k) * ux[i, j, k + 1]);

        if (indexExists(i + 1, j, k + 1))
            values.Add((x - i) * (j + 1 - y) * (z - k) * ux[i + 1, j, k + 1]);

        if (indexExists(i, j + 1, k + 1))
            values.Add((i + 1 - x) * (y - j) * (z - k) * ux[i, j + 1, k + 1]);

        if (indexExists(i + 1, j + 1, k + 1))
            values.Add((x - i) * (y - j) * (z - k) * ux[i + 1, j + 1, k + 1]);

        float sum = 0;
        foreach (float f in values)
        {
            sum += f;
        }
        return sum;
    }

    float getInterpolatedYValue(float x, float y, float z)
    {
        int i = (int)Math.Floor(x);
        int j = (int)Math.Floor(y);
        int k = (int)Math.Floor(z);
        List<float> values = new List<float>();

        if (indexExists(i, j, k))
            values.Add((i + 1 - x) * (j + 1 - y) * (k + 1 - z) * uy[i, j, k]);

        if (indexExists(i + 1, j, k))
            values.Add((x - i) * (j + 1 - y) * (k + 1 - z) * uy[i + 1, j, k]);

        if (indexExists(i, j + 1, k))
            values.Add((i + 1 - x) * (y - j) * (k + 1 - z) * uy[i, j + 1, k]);

        if (indexExists(i + 1, j + 1, k))
            values.Add((x - i) * (y - j) * (k + 1 - z) * uy[i + 1, j + 1, k]);

        if (indexExists(i, j, k + 1))
            values.Add((i + 1 - x) * (j + 1 - y) * (z - k) * uy[i, j, k + 1]);

        if (indexExists(i + 1, j, k + 1))
            values.Add((x - i) * (j + 1 - y) * (z - k) * uy[i + 1, j, k + 1]);

        if (indexExists(i, j + 1, k + 1))
            values.Add((i + 1 - x) * (y - j) * (z - k) * uy[i, j + 1, k + 1]);

        if (indexExists(i + 1, j + 1, k + 1))
            values.Add((x - i) * (y - j) * (z - k) * uy[i + 1, j + 1, k + 1]);

        float sum = 0;
        foreach (float f in values)
        {
            sum += f;
        }
        return sum;
    }

    bool indexExists(int x, int y, int z)
    {
        if (x < grid_dim_x && x > -1 && y < grid_dim_y && y > -1 && z < grid_dim_z && z > -1)
            return true;

        return false;
    }


    float getInterpolatedZValue(float x, float y, float z)
    {
        int i = (int)Math.Floor(x);
        int j = (int)Math.Floor(y);
        int k = (int)Math.Floor(z);
        List<float> values = new List<float>();

        if (indexExists(i, j, k))
            values.Add((i + 1 - x) * (j + 1 - y) * (k + 1 - z) * uz[i, j, k]);

        if (indexExists(i + 1, j, k))
            values.Add((x - i) * (j + 1 - y) * (k + 1 - z) * uz[i + 1, j, k]);

        if (indexExists(i, j + 1, k))
            values.Add((i + 1 - x) * (y - j) * (k + 1 - z) * uz[i, j + 1, k]);

        if (indexExists(i + 1, j + 1, k))
            values.Add((x - i) * (y - j) * (k + 1 - z) * uz[i + 1, j + 1, k]);

        if (indexExists(i, j, k + 1))
            values.Add((i + 1 - x) * (j + 1 - y) * (z - k) * uz[i, j, k + 1]);

        if (indexExists(i + 1, j, k + 1))
            values.Add((x - i) * (j + 1 - y) * (z - k) * uz[i + 1, j, k + 1]);

        if (indexExists(i, j + 1, k + 1))
            values.Add((i + 1 - x) * (y - j) * (z - k) * uz[i, j + 1, k + 1]);

        if (indexExists(i + 1, j + 1, k + 1))
            values.Add((x - i) * (y - j) * (z - k) * uz[i + 1, j + 1, k + 1]);

        float sum = 0;
        foreach (float f in values)
        {
            sum += f;
        }
        return sum;
    }

    void advectUx2(int x, int y, int z)
    {
        Vector3 vel = getVelocity(x, y - 0.5f, z);

        //calculate spatial position of Qi, j, k -> X
        float x_cur = x - 0.5f;

        //calculate Xprev
        float x_prev = x_cur - vel.x * timestep;

        //clamp index, better methods exist for this
        if (x_prev < -0.5)
            x_prev = -0.5f;
        if (x_prev > grid_dim_x - 0.5)
            x_prev = grid_dim_x - 0.5f;

        int ind_right = (int)Math.Ceiling(x_prev + 0.5f);
        int ind_left = (int)Math.Floor(x_prev + 0.5f);
        float dist = x_prev - (ind_left - 0.5f); //distance from left grid edge
        float vel_right = ux[ind_right, y, z];
        float vel_left = ux[ind_left, y, z];

        //calculate speed at x_prev by linear interpolation -> not in paper but in slides
       // float new_vel = getVelocity(x_prev, y, z).x;
        float new_vel = Mathf.Lerp(vel_left, vel_right, dist);
        ux_new[x, y, z] = new_vel;
    }

    void advectUy2(int x, int y, int z)
    {
        Vector3 vel = getVelocity(x - 0.5f, y, z);

        //calculate spatial position of Qi, j, k -> X
        float y_cur = y - 0.5f;

        //calculate Xprev
        float y_prev = y_cur - vel.y * timestep;

        //clamp index, better methods exist for this
        if (y_prev < -0.5)
            y_prev = -0.5f;
        if (y_prev > grid_dim_y - 0.5)
            y_prev = grid_dim_y - 0.5f;

        int ind_up = (int)Math.Ceiling(y_prev + 0.5f);
        int ind_down = (int)Math.Floor(y_prev + 0.5f);
        float dist = y_prev - (ind_down - 0.5f); //distance from left grid edge
        float vel_up = uy[x, ind_up, z];
        float vel_down = uy[x, ind_down, z];

        //calculate speed at x_prev by linear interpolation -> not in paper but in slides
        //float new_vel = getVelocity(x, y_prev,z).y;
        float new_vel = Mathf.Lerp(vel_down, vel_up, dist);
        uy_new[x, y, z] = new_vel;
    }

    void advectUz2(int x, int y, int z)
    {
        Vector3 vel = getVelocity(x, y, z - 0.5f);

        //calculate spatial position of Qi, j, k -> X
        float z_cur = z - 0.5f;

        //calculate Xprev
        float z_prev = z_cur - vel.z * timestep;

        //clamp index, TODO: better methods exist for this
        if (z_prev < -0.5)
            z_prev = -0.5f;
        if (z_prev > grid_dim_z - 0.5)
            z_prev = grid_dim_z - 0.5f;

        int ind_back = (int)Math.Ceiling(z_prev + 0.5f);
        int ind_front = (int)Math.Floor(z_prev + 0.5f);
        float dist = z_prev - (ind_front - 0.5f); //distance from left grid edge
        float vel_back = uz[x, y, ind_back];
        float vel_front = uz[x, y, ind_front];


        //calculate speed at x_prev by linear interpolation -> not in paper but in slides
        //  float new_vel = getVelocity(x, y, z_prev).z;
        float new_vel = Mathf.Lerp(vel_front, vel_back, dist);
        uz_new[x, y, z] = new_vel;
    }

    void advectUx(int x, int y, int z)
    {
        float x_vel;
        if (x == grid_dim_x)
            x_vel = ux[x, y, z];
        else
            x_vel = ux[x + 1, y, z];

        //calculate -dQ/dt using central differencing
        float dif_x = (x_vel - ux[x, y, z]) / (2 * timestep); //average value per timestep in this cell

        //calculate spatial position of Qi, j, k -> X
        float x_cur = x - 0.5f;

        //calculate Xprev
        float x_prev = x_cur - dif_x * timestep; 


        //clamp index, better methods exist for this
        if (x_prev < -0.5)
            x_prev = -0.5f;
        if (x_prev > grid_dim_x - 0.5)
            x_prev = grid_dim_x - 0.5f;

        int ind_right = (int)Math.Ceiling(x_prev + 0.5f);
        int ind_left = (int)Math.Floor(x_prev + 0.5f);
        float dist = x_prev - (ind_left - 0.5f); //distance from left grid edge
        float vel_right = ux[ind_right, y, z];
        float vel_left = ux[ind_left, y, z];

       //calculate speed at x_prev by linear interpolation -> not in paper but in slides
        float new_vel = Mathf.Lerp(vel_left, vel_right, dist);
        ux_new[x, y, z] = new_vel;
    }

    void advectUy(int x, int y, int z)
    {
        float y_vel;
        if (y == grid_dim_y)
            y_vel = uy[x, y, z];
        else
            y_vel = uy[x, y + 1, z];

        //calculate -dQ/dt using central differencing
        float dif_y = (y_vel - uy[x, y, z]) / (2 * timestep); //average value per timestep in this cell

        //calculate spatial position of Qi, j, k -> X
        float y_cur = y - 0.5f;

        //calculate Xprev
        float y_prev = y_cur - dif_y * timestep; //do stuff if x_prev falls exactly on an edge?

        //clamp index, better methods exist for this
        if (y_prev < -0.5)
            y_prev = -0.5f;
        if (y_prev > grid_dim_y - 0.5)
            y_prev = grid_dim_y - 0.5f;

        int ind_up = (int)Math.Ceiling(y_prev + 0.5f);
        int ind_down = (int)Math.Floor(y_prev + 0.5f);
        float dist = y_prev - (ind_down - 0.5f); //distance from left grid edge
        float vel_up = uy[x, ind_up, z];
        float vel_down = uy[x, ind_down, z];

        //calculate speed at x_prev by linear interpolation -> not in paper but in slides
        float new_vel = Mathf.Lerp(vel_down, vel_up, dist);
        uy_new[x, y, z] = new_vel;
    }

    void advectUz(int x, int y, int z)
    {
        float z_vel;
        if (z == grid_dim_z)
            z_vel = uz[x, y, z];
        else
            z_vel = uz[x, y, z + 1];

        //calculate -dQ/dt using central differencing
        float dif_z = (z_vel - uz[x, y, z]) / (2 * timestep); //average value per timestep in this cell

        //calculate spatial position of Qi, j, k -> X
        float z_cur = z - 0.5f;

        //calculate Xprev
        float z_prev = z_cur - dif_z * timestep; //do stuff if x_prev falls exactly on an edge?

        //clamp index, TODO: better methods exist for this
        if (z_prev < -0.5)
            z_prev = -0.5f;
        if (z_prev > grid_dim_z - 0.5)
            z_prev = grid_dim_z - 0.5f;

        int ind_back = (int)Math.Ceiling(z_prev + 0.5f);
        int ind_front = (int)Math.Floor(z_prev + 0.5f);
        float dist = z_prev - (ind_front - 0.5f); //distance from left grid edge
        float vel_back = uz[x, y, ind_back];
        float vel_front = uz[x, y, ind_front];


        //calculate speed at x_prev by linear interpolation -> not in paper but in slides
        float new_vel = Mathf.Lerp(vel_front, vel_back, dist);
        uz_new[x, y, z] = new_vel;
    }

    void moveParticles()
    {
        particlePool = new Thread[num_particles];
        for (int p = 0; p < num_particles; p++)
        {
            int index = p;
            if (threading) {
                particlePool[p] = new Thread(moveParticle);
                particlePool[p].Start(index);
            }
            else
                moveParticle(p);
        }
        for (int p = 0; p < num_particles; p++)
        {
            if(threading)
            {
                particlePool[p].Join();
            }            
            renderers_particle[p].transform.position = particles[p];
        }
    }

    void moveParticle(System.Object index)
    {
        //get particle position
        Vector3 pos = particles[(int)index];
        
        float new_x = moveX(pos);
        float new_y = moveY(pos);
        float new_z = moveZ(pos);

        pos.x = new_x;
        pos.y = new_y;
        pos.z = new_z;

        particles[(int)index] = pos;

    }

    float moveX(Vector3 pos)
    {
        //get particle cell indices (x, y, z)
        int ind_right = (int)Math.Ceiling(pos.x + 0.5);
        int ind_left = (int)Math.Floor(pos.x + 0.5);

        int ind_y = (int)Math.Floor(pos.y + 0.5);
        int ind_z = (int)Math.Floor(pos.z + 0.5);

        if (ind_y > grid_dim_y - 1)
            ind_y = grid_dim_y - 1;

        if (ind_z > grid_dim_z - 1)
            ind_z = grid_dim_z - 1;

        float dist = pos.x - (ind_left - 0.5f); //distance from left grid edge

        //interpolate new speed
        float vel_right = ux[ind_right, ind_y, ind_z];
        float vel_left = ux[ind_left, ind_y, ind_z];
        float new_vel_x = Mathf.Lerp(vel_left, vel_right, dist);
        float new_pos_x = pos.x + new_vel_x * timestep;

        //if particle wants to move out of grid, put it at grid edge
        if (clamp_part_pos)
        {
            if (new_pos_x > grid_dim_x - 0.5)
                new_pos_x = grid_dim_x - 0.5f;
            if (new_pos_x < -0.5)
                new_pos_x = -0.5f;
        }

        return new_pos_x;
    }

    float moveY(Vector3 pos)
    {
        //get particle cell indices (x, y, z)
        int ind_up = (int)Math.Ceiling(pos.y + 0.5);
        int ind_down = (int)Math.Floor(pos.y + 0.5);

        int ind_x = (int)Math.Floor(pos.x + 0.5);
        int ind_z = (int)Math.Floor(pos.z + 0.5);

        if (ind_x > grid_dim_x - 1)
            ind_x = grid_dim_x - 1;

        if (ind_z > grid_dim_z - 1)
            ind_z = grid_dim_z - 1;

        float dist = pos.y - (ind_down - 0.5f); //distance from down grid edge

        //interpolate new speed
        float vel_up = uy[ind_x, ind_up, ind_z];
        float vel_down = uy[ind_x, ind_down, ind_z];
        float new_vel_y = Mathf.Lerp(vel_down, vel_up, dist);
        float new_pos_y = pos.y + new_vel_y * timestep;

        //if particle wants to move out of grid, put it at grid edge (TODO: needs to be changed at some point)
        if (clamp_part_pos)
        {
            if (new_pos_y > grid_dim_y - 0.5)
                new_pos_y = grid_dim_y - 0.5f;
            if (new_pos_y < -0.5)
                new_pos_y = -0.5f;
        }

        return new_pos_y;
    }

    float moveZ(Vector3 pos)
    {
        //get particle cell indices (x, y, z)
        int ind_back = (int)Math.Ceiling(pos.z + 0.5);
        int ind_front = (int)Math.Floor(pos.z + 0.5);

        int ind_x = (int)Math.Floor(pos.x + 0.5);
        int ind_y = (int)Math.Floor(pos.y + 0.5);

        if (ind_x > grid_dim_x - 1)
            ind_x = grid_dim_x - 1;

        if (ind_y > grid_dim_y - 1)
            ind_y = grid_dim_y - 1;

        float dist = pos.z - (ind_front - 0.5f); //distance from back grid edge

        //interpolate new speed
        float vel_back = uz[ind_x, ind_y, ind_back];
        float vel_front = uz[ind_x, ind_y, ind_front];
        float new_vel_z = Mathf.Lerp(vel_front, vel_back, dist);
        float new_pos_z = pos.z + new_vel_z * timestep;

        //if particle wants to move out of grid, put it at grid edge (TODO: needs to be changed at some point)
        if (clamp_part_pos)
        {
            if (new_pos_z > grid_dim_z - 0.5)
                new_pos_z = grid_dim_z - 0.5f;
            if (new_pos_z < -0.5)
                new_pos_z = -0.5f;
        }

        return new_pos_z;
    }

    List<Vector3> getWaterCells()
    {
        List<Vector3> water_cells = new List<Vector3>();
        for (int x = 0; x < grid_dim_x; x++)
            for (int y = 0; y < grid_dim_y; y++)
                for (int z = 0; z < grid_dim_z; z++)
                {
                    if (cell_state[x, y, z] == state.water)
                        water_cells.Add(new Vector3(x, y, z));
                }
        return water_cells;
    }

    float[,] getAMatrix()
    {
        List<Vector3> water_cells = getWaterCells();
        float[,] AMatrix = new float[water_cells.Count, water_cells.Count];

        for (int i = 0; i < water_cells.Count; i++)
        {
            int cell_x = (int) water_cells.ElementAt(i).x;
            int cell_y = (int) water_cells.ElementAt(i).y;
            int cell_z = (int) water_cells.ElementAt(i).z;

            for (int j = 0; j < water_cells.Count; j++)
            {
                int cell_x2 = (int)water_cells.ElementAt(j).x;
                int cell_y2 = (int)water_cells.ElementAt(j).y;
                int cell_z2 = (int)water_cells.ElementAt(j).z;

                if (i == j)
                {
                    int solid_neighbours = 0;
                    if (cell_x == 0 || cell_x == grid_dim_x - 1)
                        solid_neighbours++;

                    if (cell_y == 0 || cell_y == grid_dim_y - 1)
                        solid_neighbours++;

                    if (cell_z == 0 || cell_z == grid_dim_z - 1)
                        solid_neighbours++;

                    AMatrix[i, j] = -(6 - solid_neighbours);
                }
                else
                {
                    //cell is a neighbor is 2 of (x, y, z) are the same and the last one is either +/-1
                    if (isNeighbor(cell_x, cell_y, cell_z, cell_x2, cell_y2, cell_z2))
                        AMatrix[i, j] = 1;
                }
            }
        }
        return AMatrix;
    }

    bool isNeighbor(int ix, int iy, int iz, int jx, int jy, int jz)
    {
        if (ix == jx)
        {
            if (iy == jy)
            {
                if ((iz == jz - 1) || (iz == jz + 1))
                    return true;
                else return false;
            }
            else if (iz == jz)
            {
                if ((iy == jy - 1) || (iy == jy + 1))
                    return true;
                else return false;
            }
            else return false;
        }

        else if (iy == jy)
        {
            if (iz == jz)
            {
                if ((ix == jx - 1) || (ix == jx + 1))
                    return true;
                else return false;
            }
            else return false;
        }
        else return false;
    }
}
