from hashlib import new
from turtle import Vec2D
import pygame
from pygame.locals import *
from pygame.color import *
import pymunk as pm
import numpy as np
import sys
# import matplotlib
from matplotlib import cm

class Simulation(object):
    def __init__(self, rows, cols, s, theta, display_size = (500, 500), learn_locations = True):
        self.rows = rows
        self.cols = cols
        self.s = s
        self.theta = theta
        self.learn_locations = learn_locations
        w, h = display_size

        self.rounds = 0

        self.r1 = .5 * s * np.sqrt(2) * np.sqrt(1 - np.cos(theta))
        self.r2 = .5 * s * np.sqrt(2) * np.sqrt(1 + np.cos(theta))

        self.x_offset = (w - rows * (self.r1 + self.r2)) / 2 # r1 + r2 = height/width of the theta-tilted square's bounding box
        self.y_offset = (h - cols * (self.r1 + self.r2)) / 2

        # 0 = flexible, integer value k = became rigid in round k
        # these two arrays store both explicitly and implicitly rigidified tiles and holes
        self.tiles_rigid = np.zeros([rows, cols])
        self.rifts_rigid = np.zeros([rows - 1, cols - 1])

        # store only tiles/holes explicitly made rigid
        self.tiles_rigid_explicit = np.zeros([rows, cols])
        self.rifts_rigid_explicit = np.zeros([rows - 1, cols - 1])

        self.all_tile_vertices = []
        for i in range(rows):
            self.all_tile_vertices.append([])
            for j in range(cols):
                self.all_tile_vertices[-1].append(self.get_tile_vertices(i,j))
        self.all_tile_vertices = np.array(self.all_tile_vertices)

        self.all_rift_vertices = []
        for i in range(rows - 1):
            self.all_rift_vertices.append([])
            for j in range(cols - 1):
                self.all_rift_vertices[-1].append(self.get_rift_vertices(i,j))
        self.all_rift_vertices = np.array(self.all_rift_vertices)

        # print(self.all_rift_vertices[1][1])

        self.all_vertices = {tuple(x): [] for x in np.unique(self.all_tile_vertices.reshape(rows * cols * 4, 2), axis = 0)}
        for i in range(rows):
            for j in range(cols):
                for coords in self.all_tile_vertices[i][j]:
                    self.all_vertices[tuple(coords)].append(["tile", i,j])
                if i < rows - 1 and j < cols - 1:
                    for coords in self.all_rift_vertices[i][j]:
                        self.all_vertices[tuple(coords)].append(["rift", i,j])

        # print(self.all_vertices)


    def reset(self):
        self.tiles_rigid = np.zeros([self.rows, self.cols])
        self.rifts_rigid = np.zeros([self.rows - 1, self.cols - 1])
        self.tiles_rigid_explicit = np.zeros([self.rows, self.cols])
        self.rifts_rigid_explicit = np.zeros([self.rows - 1, self.cols - 1])
        self.rounds = 0

    def tile_onclick(self, i,j):
        if self.tiles_rigid_explicit[i][j] == 0:
            self.rounds += 1
            self.tiles_rigid_explicit[i][j] = self.rounds
            self.tiles_rigid[i][j] = self.rounds
            self.update_rigidity()

    def rift_onclick(self, i,j):
        if self.rifts_rigid_explicit[i][j] == 0:
            self.rounds += 1
            self.rifts_rigid_explicit[i][j] = self.rounds
            self.rifts_rigid[i][j] = self.rounds
            self.update_rigidity()

    def rigidify_random_floppy_tile(self, k = 1, selection = np.argmin, select_tiles_only = True, can_select_implicitly_rigid = True):

        if select_tiles_only:
            if can_select_implicitly_rigid: # do not select from tiles that have been explicitly rigidified
                floppy_tiles = np.array(np.where(self.tiles_rigid_explicit == 0)).T
            else: # do not select from tiles that have been explicitly or implicitly rigidified
                floppy_tiles = np.array(np.where(self.tiles_rigid == 0)).T

            kk = min(k, len(floppy_tiles))
            if kk != 0:
                samples = floppy_tiles[np.random.randint(len(floppy_tiles), size = kk)]
                    
                samples_results = []
                for sample in samples:
                    s_tiles = np.copy(self.tiles_rigid)
                    s_tiles[sample[0]][sample[1]] = self.rounds
                    s_rifts = np.copy(self.rifts_rigid)
                    sample_tiles_outcome = self.update_rigidity_helper(s_tiles, s_rifts)
                    samples_results.append(np.count_nonzero(sample_tiles_outcome[0]) + np.count_nonzero(sample_tiles_outcome[1]))
                
                sample = samples[selection(samples_results)]

                self.rounds += 1
                self.tiles_rigid_explicit[sample[0]][sample[1]] = self.rounds
                self.tiles_rigid[sample[0]][sample[1]] = self.rounds
                self.update_rigidity()

        else:
            if can_select_implicitly_rigid: # do not select from tiles/rifts that have been explicitly rigidified
                floppy_tiles = np.array(np.where(self.tiles_rigid_explicit == 0)).T
                floppy_rifts = np.array(np.where(self.rifts_rigid_explicit == 0)).T
            else: # do not select from tiles that have been explicitly or implicitly rigidified
                floppy_tiles = np.array(np.where(self.tiles_rigid == 0)).T
                floppy_rifts = np.array(np.where(self.rifts_rigid == 0)).T

            kk = min(k, len(floppy_tiles) + len(floppy_rifts))
            if kk != 0:
                sample_ids = np.random.randint(len(floppy_tiles) + len(floppy_rifts), size = kk)
                samples_results = []

                for sample_id in sample_ids:
                    s_tiles = np.copy(self.tiles_rigid)
                    s_rifts = np.copy(self.rifts_rigid)
                    if sample_id < len(floppy_tiles): # take from tiles
                        sample = floppy_tiles[sample_id] # convert sample id into coordinate of tile to rigidify
                        s_tiles[sample[0]][sample[1]] = self.rounds
                        sample_tiles_outcome = self.update_rigidity_helper(s_tiles, s_rifts)
                        samples_results.append(np.count_nonzero(sample_tiles_outcome[0]) + np.count_nonzero(sample_tiles_outcome[1]))

                    else: # sample from floppy rifts
                        sample = floppy_rifts[sample_id - len(floppy_tiles)]
                        s_rifts[sample[0]][sample[1]] = self.rounds
                        sample_rifts_outcome = self.update_rigidity_helper(s_tiles, s_rifts)
                        samples_results.append(np.count_nonzero(sample_rifts_outcome[0]) + np.count_nonzero(sample_rifts_outcome[1]))
            # print(samples_results)
            sample_id = sample_ids[selection(samples_results)] # these were already sampled in a random order using randint
            # print(sample_id) 
            self.rounds += 1
            if sample_id < len(floppy_tiles):
                sample = floppy_tiles[sample_id]
                self.tiles_rigid_explicit[sample[0]][sample[1]] = self.rounds
                self.tiles_rigid[sample[0]][sample[1]] = self.rounds
                self.update_rigidity()

            else:
                sample = floppy_rifts[sample_id - len(floppy_tiles)]
                self.rifts_rigid_explicit[sample[0]][sample[1]] = self.rounds
                self.rifts_rigid[sample[0]][sample[1]] = self.rounds
                self.update_rigidity()

    def update_rigidity_helper(self, rigid_tiles_array, rigid_rifts_array):
        new_tiles = np.copy(rigid_tiles_array)
        new_rifts = np.copy(rigid_rifts_array)
        changed = True
        while changed: 
            changed = False
            # iterate through vertices; if any vertex with 4 angles meeting now has 3 known, the 4th becomes rigid too
            if not self.learn_locations:
                for coord, vertices in self.all_vertices.items():
                    count_rigid = 0
                    for v in vertices:
                        if v[0] == "tile" and new_tiles[v[1]][v[2]] > 0:
                            count_rigid += 1
                        elif v[0] == "rift" and new_rifts[v[1]][v[2]] > 0:
                            count_rigid += 1

                    if len(vertices) == 4 and count_rigid == 3: # account for boundary nodes where there may not be 4 complete rifts/tiles meeting
                        changed = True
                        for v in vertices:
                            if v[0] == "tile" and new_tiles[v[1]][v[2]] == 0:
                                new_tiles[v[1]][v[2]] = self.rounds
                            elif v[0] == "rift" and new_rifts[v[1]][v[2]] == 0:
                                new_rifts[v[1]][v[2]] = self.rounds

            # if we know the locations of colored tiles, not only if they're rigid or not
            # now, if any component X has 3+ vertices which are fixed (because another component with a vertex in the same position is fixed), X becomes fixed
            if self.learn_locations:
                memoize_coord_rigidity = {} # save and look up whether a vertex is rigid
                for i in range(self.rows):
                    for j in range(self.cols):
                        # for tiles
                        if new_tiles[i][j] == 0:
                            vertices = self.all_tile_vertices[i][j]
                            if sum([self.is_vertex_fixed(tuple(v), new_tiles, new_rifts, memoize_coord_rigidity) for v in vertices]) >= 3:
                                changed = True
                                new_tiles[i][j] = self.rounds

                        # for rifts
                        if i < self.rows - 1 and j < self.cols - 1 and new_rifts[i][j] == 0:
                            vertices = self.all_rift_vertices[i][j]
                            if sum([self.is_vertex_fixed(tuple(v), new_tiles, new_rifts, memoize_coord_rigidity) for v in vertices]) >= 3:
                                changed = True
                                new_rifts[i][j] = self.rounds

        return new_tiles, new_rifts

    # see if any components with a vertex at these coordinates are rigid
    def is_vertex_fixed(self, coords, new_tiles, new_rifts, memoize_coord_rigidity):
        if coords in memoize_coord_rigidity:
            return memoize_coord_rigidity[coords]
        else: 
            vertices_at_coords = self.all_vertices[coords]
            for v in vertices_at_coords:
                if v[0] == "tile" and new_tiles[v[1]][v[2]] > 0:
                    memoize_coord_rigidity[coords] = True
                    return True
                elif v[0] == "rift" and new_rifts[v[1]][v[2]] > 0:
                    memoize_coord_rigidity[coords] = True
                    return True
            memoize_coord_rigidity[coords] = False
            return False




    def update_rigidity(self):
        new_tiles, new_rifts = self.update_rigidity_helper(self.tiles_rigid, self.rifts_rigid)
        self.tiles_rigid = new_tiles
        self.rifts_rigid = new_rifts



    # r = 0 for floppy, integers in order of which tiles were clicked
    def get_tile_color(self, r):
        if r > 0: 
            # return THECOLORS["red"]
            col = cm.autumn(1 - (1.0 * r) / self.rounds)
            # return pygame.Color(255 * 1. *  np.array(col))
            return pygame.Color(255 * .9 *  np.array(col))

        else:
            return THECOLORS["darkseagreen2"]


    def get_rift_color(self, r):
        if r > 0:
            # return THECOLORS["lightpink"]
            col = cm.autumn(1 - (1.0 * r) / self.rounds)
            # return pygame.Color(255 * 1. *  np.array(col) + [-30,0,30,0])
            return pygame.Color(255 * 1. *  np.array(col))

        else:
            return THECOLORS["honeydew1"]

    # get pygame vertices for the i,jth tile from the upper left
    def get_tile_vertices(self, i,j):
        r1 = self.r1
        r2 = self.r2
        if (i + j) % 2 == 0:
            vertices_base = [(r1, r1 + r2), (0, r1), (r2, 0), (r1 + r2, r2)] # rotated ccw
        else:
            vertices_base = [(0, r2), (r1, 0), (r1 + r2, r1), (r2, r1 + r2)] # rotated cw

        vertices = np.array(vertices_base) + np.array([j * (r1 + r2) + self.x_offset, i * (r1 + r2) + self.y_offset])
        # vertices = list(map(lambda x: (x[0] + j * (r1 + r2), x[1] + i * (r1 + r2)), vertices_base))
        return np.round(vertices, 6)

    
    # get pygame vertices for the i,jth rift from the upper left
    def get_rift_vertices(self, i,j):
        r1 = self.r1
        r2 = self.r2
        if (i + j) % 2 == 0:
            vertices_base = [(r1, r1 + r2), (r1 + r2, r2), (r1 + 2 * r2, r1 + r2), (r1 + r2, 2 * r1 + r2)] # horizontal
        else: 
            vertices_base = [(r2, r1 + r2), (r1 + r2, r1), (2 * r1 + r2, r1 + r2), (r1 + r2, r1 + 2 * r2)] #vertical

        vertices = np.array(vertices_base) + np.array([j * (r1 + r2) + self.x_offset, i * (r1 + r2) + self.y_offset])
        return np.round(vertices, 6)