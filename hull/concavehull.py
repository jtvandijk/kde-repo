#JTVD - 2020

#based on: https://github.com/joaofig/uk-accidents/blob/master/geomath/hulls.py

#import libraries
import numpy as np
from shapely.geometry import asPoint
from shapely.geometry import asLineString
from shapely.geometry import asPolygon
from scipy.spatial import distance

class ConcaveHull(object):

    def __init__(self,points,prime_ix=0):

        #data
        self.data_set = np.array(points)

        #index
        self.indices = np.ones(self.data_set.shape[0],dtype=bool)
        self.prime_k = np.array([4,5,6,7,8])
        self.prime_ix = prime_ix

    def get_next_k(self):
        if self.prime_ix < len(self.prime_k):
            return self.prime_k[self.prime_ix]
        else:
            return -1

    def distance(self, loc_ini, loc_end):
        rep = np.array([loc_ini,]*loc_end.shape[0])
        dst = map(distance.euclidean,rep,loc_end)
        return(dst)

    @staticmethod
    def get_lowest_y_index(points):
        indices = np.argsort(points[:,1])
        return indices[0]

    def get_k_nearest(self, ix, k):
        ixs = self.indices
        base_indices = np.arange(len(ixs))[ixs]
        distances = self.distance(self.data_set[ix, :], self.data_set[ixs, :])
        sorted_indices = np.argsort(distances)
        kk = min(k, len(sorted_indices))
        k_nearest = sorted_indices[range(kk)]
        return base_indices[k_nearest]

    def calculate_headings(self, ix, ixs, ref_heading=0.0):
        r_ix = self.data_set[ix, :]
        r_ixs = self.data_set[ixs, :]
        y = r_ixs[:, 0] - r_ix[0]
        x = r_ixs[:, 1] - r_ix[1]
        bearings = (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0 - ref_heading
        bearings[bearings < 0.0] += 360.0
        return bearings

    def recurse_calculate(self):
        recurse = ConcaveHull(self.data_set, self.prime_ix + 1)
        next_k = recurse.get_next_k()
        if next_k == -1:
            return None
        return recurse.calculate(next_k)

    def calculate(self, k=3):
                
        #check validity of input
        if self.data_set.shape[0] < 3:
            return None

        if self.data_set.shape[0] == 3:
            return self.data_set

        #starting point
        kk = min(k, self.data_set.shape[0])
        first_point = self.get_lowest_y_index(self.data_set)
        current_point = first_point

        #hull
        hull = np.reshape(np.array(self.data_set[first_point,:]),(1, 2))
        test_hull = hull

        #mask the first point
        self.indices[first_point] = False

        #loop
        prev_angle = 270
        step = 2
        stop = 2 + kk

        while ((current_point != first_point) or (step == 2)) and len(self.indices[self.indices]) > 0:
            if step == stop:
                self.indices[first_point] = True

            #indices nearest neighbours
            knn = self.get_k_nearest(current_point,kk)

            #headings between first_point and the knn points, return angles
            angles = self.calculate_headings(current_point, knn, prev_angle)

            #calculate the candidate indexes (largest angles first)
            candidates = np.argsort(-angles)

            i = 0
            invalid_hull = True

            while invalid_hull and i < len(candidates):
                candidate = candidates[i]

                #check for self-intersections
                next_point = np.reshape(self.data_set[knn[candidate]], (1,2))
                test_hull = np.append(hull, next_point, axis=0)
                line = asLineString(test_hull)
                invalid_hull = not line.is_simple
                i += 1

            if invalid_hull:
                return self.recurse_calculate()

            prev_angle = self.calculate_headings(knn[candidate], np.array([current_point]))
            current_point = knn[candidate]
            hull = test_hull

            self.indices[current_point] = False
            step += 1

        poly = asPolygon(hull)

        #check if all points within shape
        count = 0
        total = self.data_set.shape[0]
        for ix in range(total):
            pt = asPoint(self.data_set[ix, :])
            if poly.intersects(pt) or pt.within(poly):
                count += 1
            else:
                d = poly.distance(pt)
                if d < 100:
                    count += 1

        #if all points within
        if count == total:
            return hull
        #if not next iteration
        else:
            return self.recurse_calculate()