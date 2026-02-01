import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import re

class Body:
    def __init__(self, x, y, z, vx, vy, vz, mass):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.vx = float(vx)
        self.vy = float(vy)
        self.vz = float(vz)
        self.mass = float(mass)

def parse_escaped_simulations(filename):
    """Parse the escaped simulations file and return list of simulations"""
    simulations = []

    with open(filename, 'r') as f:
        content = f.read()

    # Split by simulation blocks
    sim_blocks = re.split(r'Simulation \d+ escaped', content)[1:]  # Skip empty first element

    for block in sim_blocks:
        bodies = []
        lines = block.strip().split('\n')

        for line in lines:
            if line.startswith('Body'):
                # Parse: Body X: pos(x,y,z) vel(vx,vy,vz) mass(mass)
                match = re.search(r'pos\(([^)]+)\) vel\(([^)]+)\) mass\(([^)]+)\)', line)
                if match:
                    pos_str = match.group(1)
                    vel_str = match.group(2)
                    mass_str = match.group(3)

                    pos_vals = [float(x) for x in pos_str.split(',')]
                    vel_vals = [float(x) for x in vel_str.split(',')]
                    mass = float(mass_str)

                    body = Body(pos_vals[0], pos_vals[1], pos_vals[2],
                              vel_vals[0], vel_vals[1], vel_vals[2], mass)
                    bodies.append(body)

        if len(bodies) == 3:  # Only include complete 3-body simulations
            simulations.append(bodies)

    return simulations

def calculate_triangle_features(bodies):
    """Calculate features that describe the triangle shape for color assignment"""
    positions = np.array([[b.x, b.y] for b in bodies])

    # Calculate side lengths
    sides = []
    for i in range(3):
        j = (i + 1) % 3
        side = np.linalg.norm(positions[j] - positions[i])
        sides.append(side)

    # Sort sides for consistent ordering
    sides = sorted(sides)

    # Calculate side ratios (normalize by largest side)
    max_side = max(sides)
    if max_side > 0:
        side_ratios = [s / max_side for s in sides]
    else:
        side_ratios = [1.0, 1.0, 1.0]

    # Calculate area using cross product
    v1 = positions[1] - positions[0]
    v2 = positions[2] - positions[0]
    area = abs(v1[0] * v2[1] - v1[1] * v2[0]) / 2.0

    # Calculate angles
    angles = []
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3

        # Vectors from vertex i to j and i to k
        v_ij = positions[j] - positions[i]
        v_ik = positions[k] - positions[i]

        # Cosine of angle at vertex i
        cos_angle = np.dot(v_ij, v_ik) / (np.linalg.norm(v_ij) * np.linalg.norm(v_ik))
        cos_angle = np.clip(cos_angle, -1, 1)  # Ensure valid range
        angle = np.arccos(cos_angle)
        angles.append(angle)

    # Sort angles for consistency
    angles = sorted(angles)

    # Return feature vector
    features = side_ratios + [area] + angles
    return np.array(features)

def calculate_forces(bodies):
    """Calculate gravitational forces between bodies"""
    G = 6.67430e-11  # Gravitational constant (scaled for simulation)
    forces = []

    for i in range(3):
        total_force = np.array([0.0, 0.0, 0.0])

        for j in range(3):
            if i != j:
                # Vector from body i to body j
                dx = bodies[j].x - bodies[i].x
                dy = bodies[j].y - bodies[i].y
                dz = bodies[j].z - bodies[i].z

                dist = np.sqrt(dx*dx + dy*dy + dz*dz)

                if dist > 0:
                    # Gravitational force magnitude
                    force_mag = G * bodies[i].mass * bodies[j].mass / (dist * dist)

                    # Force direction (normalized vector)
                    force_dir = np.array([dx, dy, dz]) / dist

                    # Add force vector
                    total_force += force_mag * force_dir

        # Store magnitude of total force on this body
        forces.append(np.linalg.norm(total_force))

    return forces

def draw_force_triangle(bodies, forces, ax, sim_idx, triangle_color):
    """Draw triangle with corners at body positions and side lengths proportional to forces"""

    # Body positions (2D projection - use x,y coordinates only)
    positions = np.array([[b.x, b.y] for b in bodies])

    # Draw the triangle connecting the three bodies with assigned color
    triangle = patches.Polygon(positions, closed=True, alpha=0.3, color=triangle_color,
                              label=f'Simulation {sim_idx}')
    ax.add_patch(triangle)

    # Plot body positions with triangle color
    for i, (body, force) in enumerate(zip(bodies, forces)):
        ax.scatter(body.x, body.y,
                  c=[triangle_color], s=50,
                  alpha=0.8, edgecolors='black')

    # Draw triangle edges
    for i in range(3):
        j = (i + 1) % 3
        ax.plot([positions[i, 0], positions[j, 0]],
                [positions[i, 1], positions[j, 1]],
                color='black', linewidth=2, alpha=0.7)

    # Add simulation index as text
    centroid = np.mean(positions, axis=0)
    ax.text(centroid[0], centroid[1], f'{sim_idx}', ha='center', va='center',
            fontsize=12, fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

def main():
    # Parse escaped simulations
    simulations = parse_escaped_simulations('../build/escaped_simulations.txt')

    if not simulations:
        print("No escaped simulations found in escaped_simulations.txt")
        return

    print(f"Found {len(simulations)} escaped simulations")

    # Calculate features for all triangles
    all_features = []
    for bodies in simulations:
        features = calculate_triangle_features(bodies)
        all_features.append(features)

    all_features = np.array(all_features)

    # Assign colors based on triangle similarity using distance-based approach
    if len(simulations) > 1:
        # Calculate similarity scores based on distance to mean triangle
        mean_features = np.mean(all_features, axis=0)
        distances = np.linalg.norm(all_features - mean_features, axis=1)

        # Normalize distances to [0, 1] for color mapping
        if np.max(distances) > 0:
            similarity_scores = 1 - (distances / np.max(distances))
        else:
            similarity_scores = np.ones(len(simulations)) * 0.5

        # Map similarity to colors (high similarity = warm colors, low = cool colors)
        triangle_colors = plt.cm.plasma(similarity_scores)

        # For display purposes, assign cluster-like labels
        cluster_labels = np.digitize(similarity_scores, np.linspace(0, 1, 6))  # 5 clusters
    else:
        triangle_colors = [plt.cm.plasma(0.5)] * len(simulations)
        cluster_labels = [1] * len(simulations)

    # Create grid layout
    n_sims = len(simulations)
    n_cols = int(np.ceil(np.sqrt(n_sims)))
    n_rows = int(np.ceil(n_sims / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()

    # Plot each simulation in its own subplot
    for sim_idx, (bodies, color) in enumerate(zip(simulations, triangle_colors)):
        ax = axes[sim_idx]
        forces = calculate_forces(bodies)
        draw_force_triangle(bodies, forces, ax, sim_idx, color)

        print(f"Simulation {sim_idx}: Forces = {[f'{f:.2e}' for f in forces]}")

        # Set equal aspect ratio for each subplot
        ax.set_aspect('equal')

        # Set title with cluster info
        cluster_id = cluster_labels[sim_idx] if len(simulations) > 1 else 1
        ax.set_title(f'Sim {sim_idx} (Cluster {cluster_id})', fontsize=10)

        # Remove axis labels for cleaner look
        ax.set_xlabel('')
        ax.set_ylabel('')

    # Hide empty subplots
    for i in range(n_sims, len(axes)):
        axes[i].set_visible(False)

    # Set overall title
    fig.suptitle('Escaped 3-Body Simulations - Triangle Shapes by Similarity', fontsize=16)

    plt.tight_layout()
    plt.savefig('triangle_similarity_grid.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
