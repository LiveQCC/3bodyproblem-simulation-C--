#include "glad/glad.h" // Make sure this is included BEFORE glfw
#define CALLBACK __stdcall
#include <GL/glu.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>       // Needed for sin/cos
#include <vector>
#include <cstdlib>     // For rand()
#include <ctime>       // For time()
#include <fstream>     // For file output
#include "../dependencies/stb/stb_easy_font.h"

using namespace std;

// Constants for the simulation
const double G = 10.0;        // Gravitational constant (properly scaled)
const double DT = 0.01;       // Time step (reduced for stability)
const int NUM_BODIES = 3;
const int NUM_SIMULATIONS = 10; // Number of concurrent simulations

// Structure to represent a celestial body
struct Body {
    double x, y, z;     // Position
    double vx, vy, vz;  // Velocity
    double mass;        // Mass
    double r, g, b;     // Color
    vector<tuple<double, double, double>> trail; // Trail positions
};

// Global variables
vector<vector<Body>> simulations(NUM_SIMULATIONS);
vector<vector<Body>> initial_states(NUM_SIMULATIONS); // Store initial states
GLFWwindow* window;
GLFWwindow* triangle_window;

// Function to calculate gravitational force between two bodies
void calculateForce(Body& b1, Body& b2, double& fx, double& fy, double& fz) {
    double dx = b2.x - b1.x;
    double dy = b2.y - b1.y;
    double dz = b2.z - b1.z;
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    if (dist < 1e-10) return; // Avoid division by zero
    
    double force = G * b1.mass * b2.mass / (dist * dist);
    fx = force * dx / dist;
    fy = force * dy / dist;
    fz = force * dz / dist;
}

// Function to update positions and velocities using Euler integration for a single simulation
void updateSimulation(vector<Body>& bodies) {
    // Calculate forces for each body
    vector<double> fx(NUM_BODIES, 0.0);
    vector<double> fy(NUM_BODIES, 0.0);
    vector<double> fz(NUM_BODIES, 0.0);
    
    for (int i = 0; i < NUM_BODIES; ++i) {
        for (int j = 0; j < NUM_BODIES; ++j) {
            if (i != j) {
                double fxi, fyi, fzi;
                calculateForce(bodies[i], bodies[j], fxi, fyi, fzi);
                fx[i] += fxi;
                fy[i] += fyi;
                fz[i] += fzi;
            }
        }
    }
    
    // Update velocities and positions
    for (int i = 0; i < NUM_BODIES; ++i) {
        // Update velocities
        bodies[i].vx += fx[i] / bodies[i].mass * DT;
        bodies[i].vy += fy[i] / bodies[i].mass * DT;
        bodies[i].vz += fz[i] / bodies[i].mass * DT;
        
        // Update positions
        bodies[i].x += bodies[i].vx * DT;
        bodies[i].y += bodies[i].vy * DT;
        bodies[i].z += bodies[i].vz * DT;
        
        // Add current position to trail
        bodies[i].trail.emplace_back(bodies[i].x, bodies[i].y, bodies[i].z);
        // No limit - trails stay forever
    }
    
}



// Function to initialize the bodies for a single simulation with randomization
void initializeSimulation(vector<Body>& bodies, int sim_idx) {
    bodies.resize(NUM_BODIES);
    
    // Randomize masses (between 1 and 200 for more variation)
    for (int i = 0; i < NUM_BODIES; ++i) {
        bodies[i].mass = 1.0 + (rand() % 200);
    }
    
    // Randomize positions (within -5 to 5 for stronger interactions)
    for (int i = 0; i < NUM_BODIES; ++i) {
        bodies[i].x = (rand() % 10) - 5.0;
        bodies[i].y = (rand() % 10) - 5.0;
        bodies[i].z = (rand() % 10) - 5.0;
    }
    
    // Randomize velocities (orbital scale, approximately -15 to 15)
    for (int i = 0; i < NUM_BODIES; ++i) {
        bodies[i].vx = (rand() % 300 - 150) / 10.0;
        bodies[i].vy = (rand() % 300 - 150) / 10.0;
        bodies[i].vz = (rand() % 300 - 150) / 10.0;
    }
    
    // Assign colors based on simulation index
    // Different color schemes for each simulation
    float colorSchemes[10][9] = {
        {1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f}, // Sim 0: Yellow, Green, Blue
        {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f}, // Sim 1: Red, Cyan, Magenta
        {1.0f, 0.5f, 0.0f, 0.5f, 1.0f, 0.0f, 0.0f, 0.5f, 1.0f}, // Sim 2: Orange, Lime, Sky Blue
        {1.0f, 0.0f, 0.5f, 0.0f, 1.0f, 0.5f, 0.5f, 0.0f, 1.0f}, // Sim 3: Pink, Spring Green, Purple
        {0.5f, 0.5f, 0.0f, 0.0f, 0.5f, 0.5f, 0.5f, 0.0f, 0.5f}, // Sim 4: Olive, Teal, Maroon
        {1.0f, 0.75f, 0.5f, 0.5f, 1.0f, 0.75f, 0.75f, 0.5f, 1.0f}, // Sim 5: Peach, Mint, Lavender
        {0.75f, 1.0f, 0.0f, 0.0f, 0.75f, 1.0f, 1.0f, 0.0f, 0.75f}, // Sim 6: Chartreuse, Azure, Rose
        {0.5f, 0.25f, 0.0f, 0.0f, 0.5f, 0.25f, 0.25f, 0.0f, 0.5f}, // Sim 7: Brown, Sea Green, Indigo
        {1.0f, 0.5f, 0.75f, 0.75f, 1.0f, 0.5f, 0.5f, 0.75f, 1.0f}, // Sim 8: Light Pink, Light Green, Light Blue
        {0.25f, 0.5f, 0.75f, 0.75f, 0.25f, 0.5f, 0.5f, 0.75f, 0.25f}  // Sim 9: Steel Blue, Plum, Lime Green
    };
    
    for (int i = 0; i < NUM_BODIES; ++i) {
        bodies[i].r = colorSchemes[sim_idx][i*3];
        bodies[i].g = colorSchemes[sim_idx][i*3 + 1];
        bodies[i].b = colorSchemes[sim_idx][i*3 + 2];
    }
}

// Function to initialize all simulations
void initializeBodies() {
    srand(time(NULL)); // Seed random number generator
    for (int sim = 0; sim < NUM_SIMULATIONS; ++sim) {
        initializeSimulation(simulations[sim], sim);
        initial_states[sim] = simulations[sim]; // Save initial state
    }
}

// Function to print a string at a given position in screen coordinates
void print_string(float x, float y, char *text, float r, float g, float b) {
    static char buffer[99999]; // ~500 chars
    int num_quads;

    num_quads = stb_easy_font_print(x, y, text, NULL, buffer, sizeof(buffer));

    glColor3f(r, g, b);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 16, buffer);
    glDrawArrays(GL_QUADS, 0, num_quads * 4);
    glDisableClientState(GL_VERTEX_ARRAY);
}

// OpenGL rendering function
void display() {
    // Ensure main window context is current
    glfwMakeContextCurrent(window);
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_SCISSOR_TEST);
    
    // Get current window size for dynamic grid
    int window_width, window_height;
    glfwGetWindowSize(window, &window_width, &window_height);
    
    const int GRID_COLS = 5;
    const int GRID_ROWS = 2;
    const int VIEWPORT_WIDTH = window_width / GRID_COLS;
    const int VIEWPORT_HEIGHT = window_height / GRID_ROWS;
    
    for (int sim_idx = 0; sim_idx < NUM_SIMULATIONS; ++sim_idx) {
        const auto& sim = simulations[sim_idx];
        
        // Compute center of this simulation
        double cx = 0.0, cy = 0.0, cz = 0.0;
        for (const auto& body : sim) {
            cx += body.x;
            cy += body.y;
            cz += body.z;
        }
        cx /= NUM_BODIES;
        cy /= NUM_BODIES;
        cz /= NUM_BODIES;
        
        // Compute max distance from center for scaling
        double max_dist = 0.0;
        for (const auto& body : sim) {
            double dist = sqrt((body.x - cx)*(body.x - cx) + 
                              (body.y - cy)*(body.y - cy) + 
                              (body.z - cz)*(body.z - cz));
            if (dist > max_dist) max_dist = dist;
        }
        if (max_dist < 1.0) max_dist = 1.0; // Minimum distance to avoid too close views
        
        // Set viewport and scissor for this simulation
        int row = sim_idx / GRID_COLS;
        int col = sim_idx % GRID_COLS;
        int vx = col * VIEWPORT_WIDTH;
        int vy = row * VIEWPORT_HEIGHT;
        glViewport(vx, vy, VIEWPORT_WIDTH, VIEWPORT_HEIGHT);
        glScissor(vx, vy, VIEWPORT_WIDTH, VIEWPORT_HEIGHT);
        
        // Clear this viewport
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Set projection for this viewport
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(60.0, (GLfloat)VIEWPORT_WIDTH / VIEWPORT_HEIGHT, 0.1, max_dist * 8.0);
        
        // Set up camera for this simulation, scaled to fit
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(cx, cy, cz + max_dist * 2.5,  // Eye position scaled to fit (increased multiplier)
                  cx, cy, cz,                   // Look at center
                  0.0, 1.0, 0.0);              // Up vector
        
        // Collect vertex and color data for this simulation
        vector<GLfloat> vertices;
        vector<GLfloat> colors;
        
        for (const auto& body : sim) {
            vertices.push_back(body.x);
            vertices.push_back(body.y);
            vertices.push_back(body.z);
            colors.push_back(body.r);
            colors.push_back(body.g);
            colors.push_back(body.b);
        }
        
        // Batch render points for this simulation
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, vertices.data());
        glColorPointer(3, GL_FLOAT, 0, colors.data());
        glPointSize(3.0f);  // Smaller points for grid view
        glDrawArrays(GL_POINTS, 0, vertices.size() / 3);
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        
        // Draw trails for this simulation
        for (const auto& body : sim) {
            if (body.trail.size() > 1) {
                glColor3f(body.r, body.g, body.b); // Full color for trails
                glBegin(GL_LINE_STRIP);
                for (const auto& pos : body.trail) {
                    glVertex3f(get<0>(pos), get<1>(pos), get<2>(pos));
                }
                glEnd();
            }
        }
        
        // Calculate force matrix
        double F[3][3] = {0};
        for (int i = 0; i < NUM_BODIES; ++i) {
            for (int j = i + 1; j < NUM_BODIES; ++j) {
                double dx = sim[j].x - sim[i].x;
                double dy = sim[j].y - sim[i].y;
                double dz = sim[j].z - sim[i].z;
                double dist = sqrt(dx * dx + dy * dy + dz * dz);
                if (dist > 0.0) {
                    double force = G * sim[i].mass * sim[j].mass / (dist * dist);
                    F[i][j] = F[j][i] = force;
                }
            }
        }
        
        // Draw opaque triangle connecting the three bodies
        glDisable(GL_DEPTH_TEST); // Disable depth test for solid triangle
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        glColor4f(1.0f, 1.0f, 0.0f, 0.7f); // Yellow with transparency
        glBegin(GL_TRIANGLES);
        glVertex3d( sim[0].x, sim[0].y, sim[0].z);
        glVertex3d(sim[1].x, sim[1].y, sim[1].z);
        glVertex3d(sim[2].x, sim[2].y, sim[2].z);
        glEnd();
        
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
    }
    
    glDisable(GL_SCISSOR_TEST);
    
    // Draw grid lines
    glViewport(0, 0, window_width, window_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, window_width, 0, window_height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glDisable(GL_DEPTH_TEST); // For 2D overlay
    glColor3f(0.5f, 0.5f, 0.5f); // Gray lines
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    
    // Vertical lines
    for (int i = 1; i < GRID_COLS; ++i) {
        int x = i * VIEWPORT_WIDTH;
        glVertex2i(x, 0);
        glVertex2i(x, window_height);
    }
    
    // Horizontal lines
    for (int i = 1; i < GRID_ROWS; ++i) {
        int y = i * VIEWPORT_HEIGHT;
        glVertex2i(0, y);
        glVertex2i(window_width, y);
    }
    
    glEnd();
    
    // Draw simulation numbers
    for (int sim_idx = 0; sim_idx < NUM_SIMULATIONS; ++sim_idx) {
        int row = sim_idx / GRID_COLS;
        int col = sim_idx % GRID_COLS;
        int vx = col * VIEWPORT_WIDTH;
        int vy = row * VIEWPORT_HEIGHT;
        
        char buf[16];
        sprintf(buf, "%d", sim_idx);
        print_string(vx + 5, vy + 5, buf, 1.0f, 1.0f, 1.0f);
    }
    
    glEnable(GL_DEPTH_TEST);
    
    glfwSwapBuffers(window);
}

// Function to display triangles in 2D grid
void displayTriangles() {
    // Make triangle window current
    glfwMakeContextCurrent(triangle_window);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Get triangle window size
    int window_width, window_height;
    glfwGetWindowSize(triangle_window, &window_width, &window_height);

    // Calculate grid dimensions
    int n_sims = NUM_SIMULATIONS;
    int n_cols = (int)ceil(sqrt(n_sims));
    int n_rows = (int)ceil((double)n_sims / n_cols);

    int cell_width = window_width / n_cols;
    int cell_height = window_height / n_rows;

    // Set up 2D orthographic projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, window_width, 0, window_height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glDisable(GL_DEPTH_TEST);

    for (int sim_idx = 0; sim_idx < n_sims; ++sim_idx) {
        const auto& sim = simulations[sim_idx];

        // Calculate cell position
        int row = sim_idx / n_cols;
        int col = sim_idx % n_cols;

        int cell_x = col * cell_width;
        int cell_y = row * cell_height;

        // Set scissor for this cell
        glScissor(cell_x, cell_y, cell_width, cell_height);
        glEnable(GL_SCISSOR_TEST);

        // Clear cell background
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Get body positions (2D projection)
        double positions[3][2];
        for (int i = 0; i < 3; ++i) {
            positions[i][0] = sim[i].vx;
            positions[i][1] = sim[i].vy;
        }

        // Find bounding box
        double min_x = positions[0][0], max_x = positions[0][0];
        double min_y = positions[0][1], max_y = positions[0][1];

        for (int i = 1; i < 3; ++i) {
            min_x = min(min_x, positions[i][0]);
            max_x = max(max_x, positions[i][0]);
            min_y = min(min_y, positions[i][1]);
            max_y = max(max_y, positions[i][1]);
        }

        double width = max_x - min_x;
        double height = max_y - min_y;

        // Scale to fit in cell with margin
        double margin = 0.1;
        double scale_x = (cell_width * (1 - 2*margin)) / width;
        double scale_y = (cell_height * (1 - 2*margin)) / height;
        double scale = min(scale_x, scale_y);

        if (width == 0 || height == 0) scale = 1.0;

        // Center in cell
        double offset_x = cell_x + cell_width/2.0 - (min_x + width/2.0) * scale;
        double offset_y = cell_y + cell_height/2.0 - (min_y + height/2.0) * scale;

        // Draw triangle
        glColor3f(1.0f, 1.0f, 0.0f); // Yellow triangle
        glBegin(GL_TRIANGLES);
        for (int i = 0; i < 3; ++i) {
            double x = offset_x + positions[i][0] * scale;
            double y = offset_y + positions[i][1] * scale;
            glVertex2d(x, y);
        }
        glEnd();

        // Draw body positions
        glPointSize(5.0f);
        glBegin(GL_POINTS);
        for (int i = 0; i < 3; ++i) {
            glColor3f(sim[i].r, sim[i].g, sim[i].b);
            double x = offset_x + positions[i][0] * scale;
            double y = offset_y + positions[i][1] * scale;
            glVertex2d(x, y);
        }
        glEnd();

        // Draw simulation index
        char buf[16];
        sprintf(buf, "%d", sim_idx);
        print_string(cell_x + 5, cell_y + 5, buf, 1.0f, 1.0f, 1.0f);
    }

    glDisable(GL_SCISSOR_TEST);
    glEnable(GL_DEPTH_TEST);
    glfwSwapBuffers(triangle_window);
    
    // Restore main window context
    glfwMakeContextCurrent(window);
}
// Function to save escaped simulation initial state
void saveEscapedSimulation(int sim_idx, const vector<Body>& initial_bodies) {
    ofstream file("escaped_simulations.txt", ios::app);
    if (file.is_open()) {
        file << "Simulation " << sim_idx << " escaped" << endl;
        for (int i = 0; i < NUM_BODIES; ++i) {
            const Body& b = initial_bodies[i];
            file << "Body " << i << ": pos(" << b.x << "," << b.y << "," << b.z 
                 << ") vel(" << b.vx << "," << b.vy << "," << b.vz 
                 << ") mass(" << b.mass << ")" << endl;
        }
        file << endl;
        file.close();
    }
}

// GLFW error callback
void error_callback(int error, const char* description) {
    fprintf(stderr, "Error: %s\n", description);
}

// GLFW key callback
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
}

// Main function
int main() {
    // Initialize GLFW
    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }
    
    // Create main window
    window = glfwCreateWindow(800, 600, "3-Body Problem Simulation", NULL, NULL);
    if (!window) {
        glfwTerminate();
        cout << "Failed to create GLFW window" << endl;
        exit(EXIT_FAILURE);
    }
    
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    
    // Create triangle visualization window (don't share context)
    triangle_window = glfwCreateWindow(800, 600, "Triangle Visualization", NULL, NULL);
    if (!triangle_window) {
        glfwTerminate();
        cout << "Failed to create triangle GLFW window" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Initialize GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    
    cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
    cout << "GLFW window created successfully" << endl;
    
    // Set up OpenGL
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    
    // Initialize bodies
    initializeBodies();
    
    // Main loop
    int step = 0;
    while (!glfwWindowShouldClose(window) && !glfwWindowShouldClose(triangle_window)) {  // Increased to 5000 steps for longer simulation
        // Update all simulations
        for (auto& sim : simulations) {
            updateSimulation(sim);
        }
        
        // Check for escaped simulations
        for (int sim_idx = 0; sim_idx < NUM_SIMULATIONS; ++sim_idx) {
            const auto& sim = simulations[sim_idx];
            double cx = 0.0, cy = 0.0, cz = 0.0;
            for (const auto& body : sim) {
                cx += body.x;
                cy += body.y;
                cz += body.z;
            }
            cx /= NUM_BODIES;
            cy /= NUM_BODIES;
            cz /= NUM_BODIES;
            
            double max_dist = 0.0;
            for (const auto& body : sim) {
                double dist = sqrt((body.x - cx)*(body.x - cx) + 
                                  (body.y - cy)*(body.y - cy) + 
                                  (body.z - cz)*(body.z - cz));
                if (dist > max_dist) max_dist = dist;
            }
            
            if (max_dist > 80.0) {  // Threshold for "flung away"
                //saveEscapedSimulation(sim_idx, initial_states[sim_idx]);
                // Clear old trails
                for (auto& body : simulations[sim_idx]) {
                    body.trail.clear();
                }
                initializeSimulation(simulations[sim_idx], sim_idx);
                initial_states[sim_idx] = simulations[sim_idx]; // Save new initial state
            }
        }
        
        display();
        displayTriangles();
        glfwPollEvents();
        
        // Print progress every 500 steps
        if (step % 500 == 0) {
            cout << "Step " << step << ": Simulations updated" << endl;
        }
        step++;
    }
    
    glfwDestroyWindow(window);
    glfwDestroyWindow(triangle_window);
    glfwTerminate();
    return 0;
}


