from turtle import Vec2D
import pygame
from pygame.locals import *
from pygame.color import *
import pymunk as pm
import numpy as np
import sys
# import matplotlib
from matplotlib import cm
from simulation import Simulation


auto = True # automatically rigidify tiles

DISPLAY_SIZE = (500,500)
w,h = DISPLAY_SIZE

rows = 15
cols = 15
s = 15

learn_locations = True
select_tiles_only = True
explosive_k = 25
selection = np.argmin
can_select_implicitly_rigid = True

# do not set theta to 0 or the rigidity calculation will not work
theta = np.pi/6
fps = 15


def draw_shapes(screen, sim):
    # PYGAME COORDS: X,Y from upper left corner
    # draw tiles
    for i in range(rows):
        for j in range(cols):
            vertices = sim.get_tile_vertices(i,j)
            pygame.draw.polygon(screen, 
                                sim.get_tile_color(sim.tiles_rigid[i][j]),
                                vertices.astype(int))

    # draw rifts
    for i in range(rows - 1):
        for j in range(cols - 1):
            vertices = sim.get_rift_vertices(i,j)
            pygame.draw.polygon(screen, 
                                sim.get_rift_color(sim.rifts_rigid[i][j]),
                                vertices.astype(int))


# polygon: list of points
# return true if point in polygon
def is_within_polygon(polygon, point):
    A = []
    B = []
    C = []  
    for i in range(len(polygon)):
        p1 = polygon[i]
        p2 = polygon[(i + 1) % len(polygon)]
        
        # calculate A, B and C
        a = -(p2[1] - p1[1])
        b = p2[0] - p1[0]
        c = -(a * p1[0] + b * p1[1])

        A.append(a)
        B.append(b)
        C.append(c)

    D = []
    for i in range(len(A)):
        d = A[i] * point[0] + B[i] * point[1] + C[i]
        D.append(d)

    t1 = all(d >= 0 for d in D)
    t2 = all(d <= 0 for d in D)
    return t1 or t2



def main():
    sim = Simulation(rows, cols, s, theta, display_size = DISPLAY_SIZE, learn_locations=learn_locations)
    pygame.init()
    screen = pygame.display.set_mode(DISPLAY_SIZE, 0)
    font = pygame.font.Font(None, 16) 

    # run game
    clock = pygame.time.Clock()
    running = True
    while running:
        for event in pygame.event.get():
            LEFT = 1
            RIGHT = 3
            if event.type == QUIT:
                running = False

            elif event.type == KEYDOWN and event.key == K_r:
                sim.reset()
                
            elif event.type == MOUSEBUTTONDOWN and event.button == LEFT:
                mpos = pygame.mouse.get_pos()
                for i in range(rows):
                    for j in range(cols):
                        if is_within_polygon(sim.get_tile_vertices(i,j), mpos):
                            # print("tile", i,j)
                            sim.tile_onclick(i,j)
                            # print(sim.tiles_rigid)
                        if i < rows - 1 and j < cols - 1:
                            if is_within_polygon(sim.get_rift_vertices(i,j), mpos):
                                # print("rift", i,j)
                                sim.rift_onclick(i,j)
                                # print(sim.rifts_rigid)

        screen.fill(THECOLORS["white"]) # clear the screen
        if select_tiles_only:
            pct_rigidified = round(np.count_nonzero(sim.tiles_rigid_explicit) / sim.tiles_rigid_explicit.size, 3)
            screen.blit(font.render("r: " + str(pct_rigidified), 1, THECOLORS["darkgrey"]), (5,5))

        if auto:
            sim.rigidify_random_floppy_tile(k = explosive_k, selection = selection, select_tiles_only=select_tiles_only, can_select_implicitly_rigid=can_select_implicitly_rigid)

        draw_shapes(screen, sim)

        clock.tick(fps)

        pygame.display.flip()



if __name__ == '__main__':
    sys.exit(main())