# meteor_tracker

A hobby meteor tracker which records meteors with a camera and analyses and stores the data while observing. The camera can locate itself by looking at the stars, so no need to calibrate position.

## Modules

## take_picture

This module which handles taking pictures with a camera. It gets as input camera type, ISO, exposure and aperture and returns (non-blocking) a picture (.png, .tiff, ...), wrapped in a struct. Optionally receives feedback for camera settings.

## parse_stars

This module parses the picture taken to stars, and returns a vector of stars, in the coordinate system of the camera. The brightness of the stars is relative, maybe the same scale as magnitudes but the brightest star in the image is magnitude $0^m$.

## coordinates

This module receives as input the vector of stars and their relative brightnesses, and tries to determine how the camera is facing. On the initial setup it solves the lost-in-space problem, and then it receives feedback on its past position so it does not need to solve the whole problem again, just adjust the facing parameters a bit. If this fails(the camera has been moved), it solves the lost-in-space problem again. It returns the rotation matrix which converts from geocentric to camera coordinate system. If it cannot detect enough stars or cannot solve for the matrix it pauses execution and logs it in the output file. 

### Coordinate system

#### Geocentric coordinate system

The $x$ axis is from the center of the Earth to the point of spring equinox $\gamma$, the $z$ axis is from the center of the Earth to the North Pole and the $y$ axis is such that the coordinate system is orthogonal and right-handed $(\hat{x} \times \hat{y} = \hat{z})$

#### Sky coordinate system

This is the system which specifies how to generate the flower pattern for the catalogue stars. The $x$ axis will always be parallel to the $xy$ plane of the geocentric coordinate system, and when looked from the North pole should point clockwise. The $y$ axis is radial from the center of the celestial sphere. The $z$ axis is such that this coordinate system is right-handed. Specifying this coordinate system is important because the $x$ axis in particular is used to calculate the ordering of the $k$ brightest stars in the FOV around a central star, as described in the [flower pattern](#flower-pattern) section.


#### Camera coordinate system 

The $x$ axis is from the center of the picture to the right, the $z$ axis is from the center of the picture going up, and the $y$ axis is perpendicular to the picture, going in the direction of shooting (to the sky).

### Algorithm 

From the [HYG Star database](https://github.com/astronexus/HYG-Database) we extract right ascension, declination and visual magnitude and store them in a file, sorted by brightness in descending order, removing the stars dimmer than the limiting magnitude $M_{lim} \approx 4^{m}$. We expect to have around $N \approx 300$ stars. For each one of them we generate a flower pattern as described below, which has some invariances under rotation. When we take a picture we take a central star (or multiple central stars for robustness), generate a flower pattern, and then search the database for the best match.

#### Flower pattern

A central star is picked. For it the $k \approx 10$ brightest stars in the FOV are found, and we take note of the distance between the central star and each of these stars, as well as the angle between adjacent lines connecting the central star with the outer ones. We number the stars in the flower pattern with indexes $[1..k]$, where the star of index $1$ in the flower pattern is the first one, rotating anti-clockwise, after the $x$ axis of the "camera", as described in the [sky coordinate system](#sky-coordinate-system) section. Let the distances to the $i$-th brightest star be $r_i, i \in [1..k]$, and let the angular separation between stars $i$ and $i+1$ be $\delta_i$. The star with index $1$ from the outer ones is picked as the first past the $x$ axis, rotating counter-clockwise. Then in the database of stars we will store the fourier series of $r(i) = r_i$ and $\delta(i) = \delta_i$. These are invariant under translation, and under rotation the functions $r'(i) = r(i - \tau)$, and $\delta'(i) = \delta(i - \tau)$. The variable \tau can be found using phase correlation for both $r'(i)$ and $\delta'(i)$. The best match is selected, and from it the rotation matrix is computed. 

### Database

The database is fixed for a given FOV of the camera and a given number of stars to detect $k$. The database is a .json file of $N$ entries. The key is the star id according to HYG db, value are two arrays of complex numbers, which are the DFT coeffs for $r(i)$ and $\delta(i)$. The filename should encode the number $k$, as well as FOV, in a meaningful way. In the main program, on startup the file is loaded or generated if it does not exist, and is then used in the matching process. 

## Meteor detection

This module receives the original Picture, transforms its features to 3D using the result of the coordinates module and then applies Hough transform??? to find all the lines. If internet connection is present it tries to filter out planes, satellites and such, if not it logs that this filtering needs to be done when there is enough info.
