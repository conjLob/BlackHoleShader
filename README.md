# BlackHoleShader

This is a black hole shader for Unity. Because it is performed by ray tracing, it conforms to the actual appearance of the black hole. 
The shader also condsidered the gravity reshift effect which is in theory caused by the differences of gravitational force between the observer and the light source.

## Usage
1. Put this repository in a directory of your choice within Unity project.
1. Attach a material that set `BlackHole/Schwarzschild` as a shader to a quad object.
1. Set a `Textures/T{temperature}.png` texture to the `Redshift Texture` in the shader's inspector.
   Where `{temperature}` is the temperature inside the accretion disc. 
1. In the inspector, you can adjust the apperarance of the black hole by setting the parameters as the following.

## Parameters
| parameter | description |
| --- | --- |
| Quad Scale | If the black hole goes out of the drawing, please enlarge this. |
| Speed Of Light | Speed Of Light (m / s). The default value is 1. |
| Schwarzschild Radius | If this value is changed significantly, fine tuning of other parameters may be necessary. |
| N Steps | Number of steps in solving differential equation. Increasing this value makes the black hole visible from far away. |
| Step Size | Step size in solving differential equation. Reducing this value improves accuracy but slows down the calculation. |
| Escape Velocity | When the speed of light exceeds this speed, light is considered to have escaped from the gravity of the black hole. This is the ratio to the speed of light. |
| Max Winding Of Light | Maximum winding number of light around the black hole. When this is exceeded, the calculation will be stopped. |
| Ring Radius | Outer radius of accretion disc. |
| Redshift Texture | Lookup texture to represent gravity redshift. |
| N Divisions Of r | Number of divisions of disk noise in the r direction. |
| N Divisions Of phi | Number of divisions of disk noise in the phi direction. |
| sigma | Variance of Gaussian blur. |
| Step Width | Sampling step width of Gaussian blur. |
| Threshold | Threshold of bloom. |
| Suppression | Suppression of bloom. |

## Redshift Texture
Redshift texture is a lookup texture for gravitational redshift.
The horizontal axis is 3a / R, and the vertical axis is a / r0.
Where R, r0 are the position (radius) of the light source and the position (radius) of the observer.
If you want to generate this texture, you can use the python notebook `Notebooks/Redshift.ipynb`.
Please install the following python packages, set the parameters, and run the notebook.

### Requirements
* Python3
* Jupyter Notebook
* Numpy
* Scipy
* Matplotlib
* Pillow

### Color Matching Functions
1. Download data of 2-deg XYZ CMFs [here](http://cvrl.ucl.ac.uk/cmfs.htm). Select 0.1mm as step size and csv as format.
1. Put the downloaded file in the same directory as the notebook.

### Parameters
| parameter | description |
| --- | --- |
| size | Size of texture. If you increase this, the calculation will be slower. |
| temperature | Temperature inside the accretion disc. Black holes are actually very hot, but if set as so, they maybe look bad. |
| base_temperature | Base temperature of color temperature. The light is white at this temperature. |

## TODO
* Draw stars.
* Improve the accuracy of the buffer depth near the Schwarzschild radius.
* Code optimization.
* Correction of rotation speed of accretion disc.
* Doppler effect.

## License
MIT License. See the [LICENSE](./LICENSE) file.
