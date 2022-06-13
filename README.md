# meteor_tracker
A hobby meteor tracker which records meteors with a camera and analyses and stores the data while observing. The camera can locate itself by looking at the stars, so no need to calibrate position.

## Modules

### take_picture
This module which handles taking pictures with a camera. It gets as input camera type, ISO, exposure and aperture and returns (non-blocking) a picture (.png, .tiff, ...), wrapped in a struct. Optionally receives feedback for camera settings.

### coordinates
This module receives as input the Picture struct, parses it into star coordinates filtering everything else, and determines how the camera is facing. On the initial setup it solves the lost-in-space problem, and then it recieves feedback on its past position so it does not need to solve the whole problem again, just adjust the facing parameters a bit. If this fails(the camera has been moved), it solves the lost-in-space problem again. It returns the rotation matrix which converts from geocentric to camera coordinate system.
