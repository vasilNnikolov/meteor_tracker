# meteor_tracker

A hobby meteor tracker which records meteors with a camera and analyses and stores the data while observing. The camera can locate itself by looking at the stars, so no need to calibrate position.

## Modules

## 1. take_picture

This module which handles taking pictures with a camera. It gets as input camera type, ISO, exposure and aperture and returns (non-blocking) a picture (.png, .tiff, ...), wrapped in a struct. Optionally receives feedback for camera settings.

## 2. coordinates

This module receives as input the Picture struct, parses it into star coordinates filtering everything else, and determines how the camera is facing. On the initial setup it solves the lost-in-space problem, and then it recieves feedback on its past position so it does not need to solve the whole problem again, just adjust the facing parameters a bit. If this fails(the camera has been moved), it solves the lost-in-space problem again. It returns the rotation matrix which converts from geocentric to camera coordinate system. If it cannot detect enough stars or cannot solve for the matrix it pauses execution and logs it in the output file. 

### 2.1 Coordinate system

#### Geocentric coordinate system

The $x$ axis is from the center of the Earth to the point of spring equinox $\gamma$, the $z$ axis is from the center of the Earth to the North Pole and the $y$ axis is such that the coordinate system is orthogonal and right-handed $(\hat{x} \times \hat{y} = \hat{z})$

#### Camera coordinate system 

The $x$ axis is from the center of the picture to the right, the $z$ axis is from the center of the picture going up, and the $y$ axis is perpendicular to the picture, going in the direction of shooting (to the sky).

### 2.1 Algorithm 

From a database of coordinates and 


## 3. Meteor detection

This module recieves the original Picture, transforms its features to 3D using the result of the coordinates module and then applies Hough transform??? to find all the lines. If internet connection is present it tries to filter out planes, satellites and such, if not it logs that this filtering needs to be done when there is enough info.
