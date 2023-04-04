from matplotlib import pyplot as plt
from numpy.random import random, seed
from numpy import array
from numpy.linalg import norm
from math import pi, cos, sin, sqrt
from time import time

from kdtree import Sphere, KDTree

seed(1234)

def main():
    start = time()
    domain = [(0, 1), (0, 1)]
    rejectedRadii = []
    sphere_list, packing_front = initializePacking(domain)
    for s in packing_front:
        s.setColor('g')
    # print("Initial Spheres")
    # print("\n".join([str(n) for n in sphere_list]))
    seed_sphere = sphere_list[0]
    # plotSpheres(sphere_list, domain)

    # tree = KDTree(sphere_list)
    # print(tree)
    # return

    # Generate new sphere
    while len(packing_front) > 0:
        newR = sampleRadius(rejected=rejectedRadii)
        # print(f"New R: {newR}")
        box_r = packing_front[0].r + 2.0*newR
        box = [ packing_front[0].x - box_r, packing_front[0].x + box_r, packing_front[0].y - box_r, packing_front[0].y + box_r ]
        neighbors = getNeighbors(box, sphere_list, packing_front[0])
        # print("Neighbors")
        # print("\n".join([str(n) for n in neighbors]))
        candidate_locations = generate_candidate_centers(neighbors, packing_front[0], newR)

        # filter intersections with neighbors and boundary
        eps = 1.0e-8
        intersects_neighbor = lambda x: any([norm(x - other.c) < (newR + other.r - eps) for other in neighbors])
        inside_domain = lambda s: s[0] > domain[0][0] + newR - eps and s[0] < domain[0][1] - newR + eps and s[1] > domain[1][0] + newR - eps and s[1] < domain[1][1] - newR + eps
        candidate_locations = [c for c in candidate_locations if not intersects_neighbor(c) and inside_domain(c)]
        # print("Candidates")
        # print(candidate_locations)

        if len(candidate_locations) > 0:
            sphere_list.append( Sphere(center=candidate_locations[0], radius=newR) )
            packing_front.append(sphere_list[-1])
            packing_front.sort(key=lambda x: norm(x.c - seed_sphere.c))
            for s in packing_front:
                s.setColor('g')
        else:
            rejectedRadii.append(newR)

        # Decide if we're done with the top of packing_front
        if len(candidate_locations) < 1:
            packing_front[0].setColor('b')
            packing_front.pop(0)
        
        # plotSpheres(sphere_list, domain)

    # Plot Spheres
    sphere_list[0].setColor('r')
    print(f"Packed {len(sphere_list)} spheres in {time() - start:.2f}s")
    plotSpheres(sphere_list, domain)


def sampleRadius(rMin=0.01, rMax=0.01, rejected=[]):
    return rejected.pop(0) if len(rejected) > 0 else random() * (rMax-rMin) + rMin

def initializePacking(domain):
    sphere_list = []
    packing_front = []

    # Initial sphere
    seed_center = array( [(domain[0][1]-domain[0][0])/2.0, (domain[1][1]-domain[1][0])/2.0] )
    seed_radius = sampleRadius()
    sphere_list.append( Sphere(center=seed_center, radius=seed_radius, color='r') )
    packing_front.append( sphere_list[-1] )
    # Second sphere
    # angle = random() * 2.0 * pi
    angle = pi/2.0
    newRadius = sampleRadius()
    halo_distance = sphere_list[-1].r + newRadius
    newOffset = array( [cos(angle)*halo_distance, sin(angle)*halo_distance] )
    sphere_list.append( Sphere( center=seed_center+newOffset, radius=newRadius) )
    packing_front.append( sphere_list[-1] )
    packing_front.sort(key=lambda x: norm(x.c - seed_center))

    return sphere_list, packing_front

def getNeighbors(box, sphere_list, exclude):
    insideBox = lambda s: s.x > box[0] - s.r and s.x < box[1] + s.r and s.y > box[2] - s.r and s.y < box[3] + s.r and s is not exclude
    return [sphere for sphere in sphere_list if insideBox(sphere)]
    
def generate_candidate_centers(neighbors, root, newR):
    intersection = []
    root_r = root.r + newR
    for sphere in neighbors:
        sphere_r = sphere.r + newR
        distance = norm(sphere.c - root.c)
        x = (distance**2 - sphere_r**2 + root_r**2) / (2.0 * distance)
        if x > root_r: continue

        unit_d = ( sphere.c - root.c ) / distance
        ortho_d = array( [-unit_d[1], unit_d[0]] )

        y = sqrt( sphere_r**2 - ((root_r**2 - sphere_r**2 - distance**2)/(2.0*distance))**2 )
        intersection.append( root.c + x*unit_d + y*ortho_d )
        intersection.append( root.c + x*unit_d - y*ortho_d )

    return intersection

def plotSpheres(sphere_list, domain):
    plt.cla(); plt.close()
    plt.figure(figsize=(8,8))
    for sphere in sphere_list:
        plt.gca().add_patch( sphere.draw() )
    plt.xlim([*domain[0]])
    plt.ylim([*domain[1]])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
    # import cProfile, pstats
    # with cProfile.Profile() as pr:
    #     main()
    #     pstats.Stats(pr).sort_stats("cumulative").print_stats()