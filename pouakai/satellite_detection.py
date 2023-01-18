import numpy as np
import cv2
from copy import deepcopy
from astropy.stats import sigma_clipped_stats
class sat_streaks():
    def __init__(self,image,thickness=17,sigma=15,angle_tol=5,run=True) -> None:
        self.image = image - np.nanmedian(image)
        self._set_threshold(sigma)
        self.thickness = thickness
        self.angle_tol = angle_tol

        if run:
            self._dilate()
            self._edges()
            self._lines()
            self._consolidate_lines()
            self.make_mask()
            self._detected()


    def _set_threshold(self,sigma):
        mean, med, std = sigma_clipped_stats(self.image)
        self.threshold = mean + sigma*std

    def _detected(self):
        if len(self.consolidated_lines) > 0:
            self.satellite = True
        else:
            self.satellite = False
        self.sat_num = len(self.consolidated_lines)

    def _dilate(self):
        # set all values below the threshold to zero
        arr = deepcopy(self.image)
        arr[arr < self.threshold] = 0

        # create a structuring element for morphological dilation
        kernel = np.ones((9, 9))

        # dilate the array
        dilated = cv2.dilate(arr, kernel, iterations=1)
        dilated[dilated<10]  = 0
        # set all non-zero values in the dilated array that are not connected to other non-zero values to zero
        arr[(arr != 0) & (dilated == 0)] = 0

        d = (dilated > 0) * 1
        self.gray = (d*255/np.max(d)).astype('uint8')
        
    def _edges(self):
        low_threshold = 50
        high_threshold = 150
        self.edges = cv2.Canny(self.gray, low_threshold, high_threshold)

    def _lines(self):
        lines = cv2.HoughLinesP(
                                self.edges, # Input edge image
                                1, # Distance resolution in pixels
                                np.pi/180, # Angle resolution in radians
                                threshold=100, # Min number of votes for valid line
                                minLineLength=200, # Min allowed length of line
                                maxLineGap=50 # Max allowed gap between line for joining them
                                )
        if lines is not None:
            good = []
            for i in range(len(lines)):
                line = lines[i]
                x1, y1, x2, y2 = line[0]
                # calculate the angle of the line
                angle = np.arctan2(y2 - y1, x2 - x1) * 180 / np.pi
                if abs(angle) < 85:
                    good += [i]
            good = np.array(good,dtype=int)

            self.lines = lines[good]
        else:
            self.lines = []


    def _consolidate_lines(self):
        angle_tolerance = self.angle_tol
        # create an empty list to store the consolidated lines
        consolidated_lines = []

        # loop through the lines detected by HoughLinesP
        for line in self.lines:
            x1, y1, x2, y2 = line[0]
            # calculate the angle of the line
            angle = np.arctan2(y2 - y1, x2 - x1) * 180 / np.pi
            # round the angle to the nearest multiple of the angle tolerance
            angle = round(angle / angle_tolerance) * angle_tolerance
            # check if a line with a similar angle has already been consolidated
            found = False
            for consolidated_line in consolidated_lines:
                if abs(consolidated_line[0] - angle) < angle_tolerance:
                    # if a similar line has been found, add the current line to the consolidated line
                    consolidated_line[1] += [(x1, y1), (x2, y2)]
                    found = True
                    break
            if not found:
                # if no similar line has been found, create a new consolidated line
                consolidated_lines.append([angle, [(x1, y1), (x2, y2)]])
        self.consolidated_lines = consolidated_lines
    

    def make_mask(self,thickness=None):
        if thickness is None:
            thickness = self.thickness
        # create a black image with the same size as the input image
        mask = np.zeros_like(self.image)

        # loop through the consolidated lines
        sat_mask = []
        for consolidated_line in self.consolidated_lines:
            angle = consolidated_line[0]
            points = np.array(consolidated_line[1])
            x = points[:,0]
            y = points[:,1]
            coefs = np.polyfit(x, y, 1)
            slope = coefs[0]
            intercept = coefs[1]
            x1 = np.min(x)
            xx = np.arange(np.min(x),np.max(x))
            y1 = int(slope * x1 + intercept)
            x2 = np.max(x)
            y2 = int(slope * x2 + intercept)
            yy = (slope * xx + intercept).astype(int)
            # Draw the line on the black image using cv2.line
            tmp = cv2.line(mask, (x1, y1), (x2, y2), (255, 255, 255), thickness=thickness)
            tmp[tmp > 0] = 1
            sat_mask += [tmp]

        # The pixels that are covered by the line will have a value of 8 in the mask
        # You can use this mask to extract the pixels that belong to the lines
        sat_mask = np.array(sat_mask).astype(int)
        if len(sat_mask) == 0:
            sat_mask = mask
        
        self.mask = sat_mask

        tmp = np.nansum(self.mask,axis=0)
        tmp[tmp > 0] = 1
        self.total_mask = tmp.astype(int)