In this homework you will write a simple ray tracer to render a scene composed of two spheres and
a tetrahedron. The two spheres are of different radii and colors (e.g., blue and green). Pick the
centers and the radii of the spheres so that they do not intersect in the 3-D space.
For simplicity, you may assume a directional light and a parallel camera geometry for rendering
purposes. This means that a single LookAt vector can be used to generate all rays. For example,
you may set the LightDir = (0, 1, 0) and LookAt = (0, 0, −1).

• Write a simple OpenGL program to display an image. Fix the window size (e.g., 256 × 256).
You may use the glDrawPixels function to draw the contents of an array of RGB values.
Verify the correctness of your program by displaying a test image that you generate (e.g., color
stripes). 

• Generate parallel rays and solve the ray-object intersection problems to determine the color
of each pixel. Verify the correctness of your ray tracer by rendering the scene from multiple
angles (e.g., change the LookAt vector). 

• Add ambient, diffuse and specular shading to your ray tracer. Pick the ambient, diffuse and
specular constants and light intensities to reflect a realistic rendering of the scene. 

• Add a glazed surface to the scene. You may choose to set one of the objects in the scene or
have a plane below these objects that appears as a glazed surface. 

• Create a small animation by changing the LookAt vector and assembling the rendered images
into a movie. 10pts

• Along with the source code and makefile, submit a report (a PDF file) that describes and
documents your experiments for each step of the homework. 

