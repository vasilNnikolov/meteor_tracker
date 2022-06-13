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

### 2.2 Algorithm 

From the [HYG Star database](https://github.com/astronexus/HYG-Database) we extract right ascension, declination and visual magnitude and store them in a file, sorted by brightness in descending order, removing the stars dimmer than the limiting magnitude $M_{lim} \approx 4^{m}$. We expect to have around $N \approx 300$ stars. For each one of them we generate a flower pattern as described below, which has some invariances under rotation. When we take a picture we take a central star (or multiple central stars for robustness), generate a flower pattern, and then search the database for the best match.

#### Flower pattern

A central star is picked. For it the $k \approx 10$ brightest stars in the FOV are found, and we take note of the distance between the central star and each of these stars, as well as the angle between adjacent lines connecting the central star with the outer ones. Let the distances to the $i$-th brightest star be $r_i, i \in [1..k]$, and let the angular separation between stars $i$ and $i+1$ be $\delta_i$. The star with index $1$ from the outer ones is picked as the first past the $x$ axis, rotating counter-clockwise. Then the pattern for a given central star is 
\begin{equation}
    P = \begin{pmatrix}
        r_1 & r_2 & ... & r_k \\
        \delta_1 & \delta_2 & ... & \delta_k
    \end{pmatrix}
\end{equation}

If we pick a given star for the central star from our image and construct its star pattern in the same way, it will be the same as in the database, but column shifted by some number $\tau$. This will be due to the difference in orientation of the $x$ axis. The search step consists of going one by one through the $N$ stars in the database, and through $k$ different options for $\tau$. Since for each option it takes $k$ comparisons to check if it is close enough, this search algorithm has time complexity $O(Nk^2)$. Not great, not terrible. A better search algorithm storing the DFT of the rows of the flower pattern, computing the DFT of the observed pattern, and finding $\tau$ using phase correlation may be implemented. 

## 3. Meteor detection

This module receives the original Picture, transforms its features to 3D using the result of the coordinates module and then applies Hough transform??? to find all the lines. If internet connection is present it tries to filter out planes, satellites and such, if not it logs that this filtering needs to be done when there is enough info.
