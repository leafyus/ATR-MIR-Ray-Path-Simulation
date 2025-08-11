

# **A Computational Framework for Simulating Attenuated Total Reflectance (ATR) Geometries in Python**

---

## **Section 1: Foundational Principles of Ray Optics for ATR Simulation**

A robust simulation of Attenuated Total Reflectance (ATR) spectroscopy requires a firm grounding in the principles of geometric optics. The behavior of a mid-infrared (MIR) ray as it traverses different media is governed by a precise set of physical laws. This section establishes the essential physical and mathematical context, detailing the phenomena at the crystal-sample interface, the nature of the evanescent wave that facilitates spectroscopic measurement, and the specific optical properties of materials relevant to the simulation.

### **1.1 The Physics of the Crystal-Sample Interface: Refraction, Reflection, and the Critical Angle**

The interaction of a light ray with a planar interface between two different, optically isotropic media is governed by the phenomena of reflection and refraction. When a ray travels from a medium with refractive index n1​ to a medium with refractive index n2​, its path changes. The angle of the incoming ray (θ1​) and the angle of the transmitted, or refracted, ray (θ2​) are related by Snell's Law.1 Both angles are measured with respect to the normal, a line perpendicular to the interface. The scalar form of Snell's Law is:

n1​sin(θ1​)=n2​sin(θ2​)  
This equation is the cornerstone of ray tracing through different materials. However, for ATR spectroscopy, a special case of reflection is paramount: Total Internal Reflection (TIR). TIR is the phenomenon where a ray is completely reflected back into the incident medium, with no light being refracted into the second medium.3 This occurs only when two specific conditions are met 5:

1. The ray must be traveling from a medium of higher refractive index to a medium of lower refractive index. This is often described as moving from an optically "denser" to an optically "rarer" medium (i.e., n1​\>n2​).  
2. The angle of incidence, θ1​, must be greater than a specific value known as the critical angle, θc​.

The critical angle represents the threshold at which the angle of refraction, θ2​, would be 90 degrees, meaning the refracted ray would travel exactly along the interface. Any angle of incidence greater than this makes refraction physically impossible. The critical angle can be calculated directly from Snell's Law by setting θ2​=90∘ (where sin(θ2​)=1) 3:

θc​=arcsin(n1​n2​​)  
This formula highlights a crucial point: since the arcsin function is only defined for arguments between \-1 and 1, a real-valued critical angle only exists if n2​≤n1​. This mathematical constraint directly corresponds to the first physical condition for TIR.8

A particularly elegant and computationally efficient approach, which will be central to the proposed software, involves using a vector-based formulation of Snell's Law. The mathematical structure of this equation naturally handles the transition to total internal reflection. Specifically, the calculation for the refracted vector involves a square root term whose argument becomes negative precisely when the physical conditions for TIR are met (i.e., when sin(θ1​)\>n2​/n1​).9 This allows the simulation to determine whether a ray refracts or reflects through a single computational check, avoiding explicit angle calculations and separate logical checks for the critical angle, thereby making the physics engine more robust and efficient.

### **1.2 The Evanescent Wave: The "Attenuated" in ATR**

While the simulation is based on geometric ray tracing, understanding the physical mechanism of ATR provides essential context. During total internal reflection, although the ray is geometrically depicted as reflecting perfectly at the boundary, an electromagnetic field known as the **evanescent wave** actually penetrates a very short distance into the rarer medium (the sample).4 This wave is non-propagating and its intensity decays exponentially with distance from the interface.

If the sample material absorbs light at the incident wavelength, it will absorb energy from this evanescent field. This absorption "attenuates" the intensity of the reflected ray. By measuring this attenuation across a range of infrared wavelengths, an absorption spectrum of the sample can be generated. This is the fundamental principle of ATR spectroscopy.4

The distance the evanescent wave effectively probes into the sample is characterized by the **penetration depth** (dp​). This is defined as the distance from the interface at which the wave's intensity has decayed to 1/e (about 37%) of its value at the surface. The penetration depth is a critical parameter in ATR experiments and depends on the wavelength of light (λ), the angle of incidence (θ1​), and the refractive indices of the crystal (n1​) and the sample (n2​).4 For a given wavelength, a larger incidence angle results in a shallower penetration depth. This simulation tool, by allowing precise control over geometry and incidence angle, enables the user to explore how these parameters affect the conditions for TIR and, by extension, the effective probing of the sample.

### **1.3 Material and Medium Properties at 5 µm**

The user's specification of a fixed MIR wavelength of 5 µm simplifies the simulation significantly, as it removes the need to account for chromatic dispersion (the change in refractive index with wavelength). At this wavelength, we can use fixed values for the refractive indices of the relevant materials.

**ATR Crystal Materials:** The choice of crystal is fundamental to ATR. Common materials possess high refractive indices, a necessary condition for achieving TIR against most samples. Their properties are summarized in Table 1\.

| Material | Refractive Index (n) @ 5 µm | Key Properties / Notes |
| :---- | :---- | :---- |
| Germanium (Ge) | \~4.017 | Very high index, for highly absorbing or high-index samples. Limited spectral range. 12 |
| Zinc Selenide (ZnSe) | \~2.42 | General purpose, low cost. Sensitive to pH \< 5 and \> 9\. 11 |
| Diamond | \~2.40 | Extremely hard and inert. Expensive. Has absorption bands in mid-IR. 12 |
| Silicon (Si) | \~3.42 | High index, hard. Affected by strong acids/alkalis. 11 |

**Sample and Surrounding Media:** The simulation must also define the media surrounding the ATR crystal. These include the sample itself and the medium at the entrance and exit faces, which is typically air. The properties for these are summarized in Table 2\.

| Medium | Refractive Index (n) @ 5 µm | Notes |
| :---- | :---- | :---- |
| Air | 1.000 | Standard reference, assumed for non-sample boundaries. 15 |
| Water (liquid) | 1.325 | Common liquid sample. This is the real part of the index; the imaginary part related to absorption is not used in this geometric simulation. 16 |

These tabulated values will serve as the default parameters within the simulation's GUI, providing a realistic and scientifically relevant starting point for analysis.

## **Section 2: Core Algorithms for 2D Ray Tracing**

Translating the physical principles into a computational model requires a set of robust and efficient algorithms. This section details the vector-based mathematical framework that forms the core of the ray-tracing engine. This approach avoids trigonometric functions where possible, leading to faster and more stable code.

### **2.1 Vector Mechanics for Rays and Boundaries**

To perform calculations, we must first define our geometric entities in a computationally friendly format. All vector mathematics will leverage the NumPy library for its performance and concise syntax.17

* **Ray Representation:** A ray is defined by its origin point and its direction. It is not an infinite line but has a starting point. It will be represented as an object containing:  
  * An origin point P, which is a 2D vector (e.g., a NumPy array \[x, y\]).  
  * A unit direction vector d, also a 2D NumPy array, indicating the direction of travel.  
    Any point along this ray can be found using the parametric equation P+t⋅d, where the scalar parameter t≥0 represents the distance along the ray from its origin.19  
* **Boundary Representation:** Each side of the ATR crystal geometry is a boundary. A boundary is defined as a line segment between two vertices, Q1​ and Q2​. Its parametric form is Q1​+u⋅(Q2​−Q1​), where the scalar parameter u is valid only in the range $$.19 Each boundary will also be associated with an outward-facing unit normal vector  
  n, which is essential for applying the laws of reflection and refraction.

### **2.2 The Ray-Segment Intersection Algorithm**

The fundamental operation of the ray tracer is to determine where a ray will strike a boundary. This is achieved by finding the intersection point of a ray and a line segment. We solve the system of equations formed by setting the parametric ray and segment equations equal:

P+t⋅d=Q1​+u⋅s  
where s=Q2​−Q1​ is the vector representing the boundary segment. This vector equation can be expanded into a system of two linear equations with two unknowns, t and u. Using the 2D vector cross product (defined as v1​×v2​=v1x​v2y​−v1y​v2x​), we can solve for t and u directly 20:

u=s×d(P−Q1​)×d​  
t=d×s(Q1​−P)×s​  
A valid intersection that is relevant to the simulation must satisfy two conditions 19:

1. t≥0: The intersection point must lie in front of the ray's origin.  
2. 0≤u≤1: The intersection point must lie on the line segment that forms the boundary.

The algorithm must also handle the edge case where the denominator, d×s, is zero. This occurs when the ray and the boundary segment are parallel.20 In this scenario, they do not intersect (or are collinear), and this check prevents a division-by-zero error.

The choice of this parametric algorithm is significant because it naturally provides the distance t to the intersection point. To correctly propagate the ray, the simulation must find the *next* boundary it will hit. This is achieved by calculating the intersection with *all* boundaries of the polygon, filtering for only the valid ones (where t and u are in their respective ranges), and then selecting the intersection corresponding to the **smallest positive value of t**. This ensures the ray advances to the nearest surface in its path, forming the fundamental logic of the simulation loop.

### **2.3 Implementing Optical Laws in Vector Form**

Once an intersection is found, the physics at the interface determines the ray's new direction. Vector formulations are used for both reflection and refraction.

* **Reflection (for TIR):** For specular reflection, the outgoing direction vector, dout​, can be calculated from the incoming direction, din​, and the surface normal, n, without any angle calculations:  
  dout​=din​−2(din​⋅n)n

  Here, (din​⋅n) is the dot product of the two vectors. This formula is used whenever the conditions for TIR are met.  
* **Refraction (Snell's Law):** To calculate the direction of a refracted ray, an explicit vector form of Snell's Law is used. A highly effective formula is given by 1:  
  dout​=μdin​+(1−μ2(1−(din​⋅n)2)​−μ(din​⋅n))n

  where μ=n1​/n2​ is the ratio of the refractive indices, and all vectors (din​, n, dout​) are unit vectors.

As previously discussed, this single formula elegantly unifies the check for refraction and TIR. The simulation first calculates the term under the square root, c=1−μ2(1−(din​⋅n)2).

* If c\<0, the condition for TIR is met, and the reflection formula is applied.  
* If c≥0, refraction occurs, and the full vector refraction formula is used to find the new direction.

This unified logic forms the heart of the physics engine, providing a computationally robust and physically accurate method for tracing the ray's path.

## **Section 3: Software Architecture and Project Structure**

A well-defined software architecture is critical for creating a tool that is maintainable, extensible, and easy to debug. The proposed structure separates the core simulation logic from the graphical user interface, a design principle that ensures maximum flexibility.

### **3.1 Recommended Python Environment and Libraries**

The project will be built using standard, well-supported Python libraries, minimizing setup complexity.

* **Python 3.x:** The project should be developed using a modern version of Python (e.g., 3.8 or newer).  
* **NumPy:** This library is essential for all numerical computations, particularly for representing and manipulating 2D vectors and points. Its high-performance ndarray objects are the foundation of the geometry and physics modules.17  
* **Matplotlib:** This will be the visualization engine, used to render the ATR crystal geometry and the calculated ray path. It is a powerful and flexible library for creating static and interactive 2D plots.17  
* **Tkinter:** As Python's built-in GUI toolkit, Tkinter is the ideal choice for the user interface. It requires no external installation, is lightweight, and is more than capable of handling the required input widgets and embedding the Matplotlib plot.23

To ensure a reproducible environment, a requirements.txt file should be created with the following content:

numpy  
matplotlib

Tkinter is part of the Python standard library and does not need to be listed.

### **3.2 Proposed Directory and File Structure**

A modular project structure separates concerns, making the codebase cleaner and easier to manage. The following structure is recommended:

atr\_simulator/  
├── atr\_sim/  
│   ├── \_\_init\_\_.py  
│   ├── core/  
│   │   ├── \_\_init\_\_.py  
│   │   ├── geometry.py   \# Polygon, Ray, vector operations  
│   │   └── physics.py      \# Snell's law, reflection, TIR logic  
│   ├── gui/  
│   │   ├── \_\_init\_\_.py  
│   │   └── app.py        \# Main Tkinter application class, GUI layout  
│   └── utils/  
│       ├── \_\_init\_\_.py  
│       └── constants.py  \# Refractive indices, material properties  
├── main.py               \# Entry point to start the application  
└── Instruction.md        \# Project documentation

This organization enforces a critical separation: the core modules contain the pure physics and geometry logic and depend only on NumPy. The gui module handles all user interaction and visualization, depending on Tkinter and Matplotlib. This design choice is fundamental to the project's long-term health. It allows the core simulation engine to be tested independently of the interface. Furthermore, if a different GUI framework (like PyQt or a web-based one using Plotly Dash) is desired in the future, only the gui module needs to be rewritten; the complex physics engine in core remains untouched and reusable.23

### **3.3 Data Flow and Class Interaction**

The components of the application will interact in a clear, defined sequence:

1. **main.py:** This script serves as the simple entry point. It imports the main Application class from atr\_sim.gui.app and launches the Tkinter event loop.  
2. **gui.app.Application:** This class is the controller of the entire application. It builds the Tkinter window, lays out all the input widgets and the Matplotlib canvas. When the user clicks the "Run Simulation" button, it orchestrates the simulation by:  
   * Gathering all parameters (geometry, materials, ray definition) from the GUI input fields.  
   * Instantiating the necessary Polygon and Ray objects from the core.geometry module.  
   * Executing the main simulation loop, which calls functions from the core.physics module.  
   * Receiving the final ray path (a list of 2D points) as output.  
   * Calling its internal plotting method to clear the Matplotlib canvas and draw the new geometry and ray path.  
3. **core.geometry:** This module provides simple data classes (Polygon, Ray) to hold the state of the simulation objects.  
4. **core.physics:** This module contains stateless functions that perform the core calculations (e.g., calculate\_interaction, find\_closest\_intersection). They take geometric and physical data as input and return calculated results.

The main simulation loop itself, likely a method within the Application class, will proceed as follows: start with the initial ray, find the closest valid intersection, store the intersection point, apply the physics at that point to get a new direction, create a new ray from that point with the new direction, and repeat until the ray exits the geometry or a maximum number of reflections is reached.

## **Section 4: Component-Level Implementation Guide**

This section provides a more detailed blueprint for the classes and functions within each module, serving as a guide for writing the Python code.

### **4.1 The utils.constants Module**

To avoid hardcoding physical values throughout the application, all constants should be centralized in this module. A dictionary is a suitable structure for this.

**Example constants.py:**

Python

\# Refractive indices for a wavelength of 5 µm  
REFRACTIVE\_INDICES \= {  
    "air": 1.0,  
    "water": 1.325,  
    "znse": 2.42,  
    "ge": 4.017,  
    "diamond": 2.40,  
    "si": 3.42  
}

### **4.2 The core.geometry Module**

This module defines the data structures for the geometric objects.

* **Ray Class:** A simple container for a ray's properties.  
  Python  
  import numpy as np

  class Ray:  
      def \_\_init\_\_(self, origin: np.ndarray, direction: np.ndarray):  
          self.origin \= origin  
          \# Ensure direction is always a unit vector  
          self.direction \= direction / np.linalg.norm(direction)

* **Polygon Class:** Represents the ATR crystal.  
  Python  
  import numpy as np

  class Polygon:  
      def \_\_init\_\_(self, vertices: list, material\_n: float, default\_boundary\_n: float):  
          self.vertices \= \[np.array(v) for v in vertices\]  
          self.material\_n \= material\_n  
          \# Map each boundary segment to an external refractive index  
          self.boundary\_media\_n \= \[default\_boundary\_n\] \* len(self.vertices)

      def set\_boundary\_medium(self, boundary\_index: int, medium\_n: float):  
          if 0 \<= boundary\_index \< len(self.boundary\_media\_n):  
              self.boundary\_media\_n\[boundary\_index\] \= medium\_n

      def get\_segments\_and\_normals(self) \-\> list:  
          \# Returns a list of tuples: (Q1, Q2, outward\_normal)  
          \#... calculation logic...

The design of the boundary\_media\_n list is critical. It directly maps the physical concept of having different media on different faces of the crystal to a clean data structure. The GUI can provide input fields for each boundary, and when the simulation finds an intersection with segment i, it can easily retrieve the corresponding external refractive index n\_2 by accessing polygon.boundary\_media\_n\[i\].

### **4.3 The core.physics Module**

This module contains the stateless functions that comprise the physics engine.

* **find\_closest\_intersection(ray: Ray, polygon: Polygon) function:**  
  * This function implements the logic described in Section 2.2.  
  * It iterates through the segments provided by polygon.get\_segments\_and\_normals().  
  * For each segment, it calculates the intersection parameters t and u.  
  * It filters for valid intersections where t\>ϵ (a small positive number to avoid self-intersection) and 0≤u≤1.  
  * It returns a dictionary containing the intersection point, segment index, and normal vector for the intersection with the minimum valid t. If no valid intersections are found, it returns None.  
* **calculate\_interaction(d\_in: np.ndarray, normal: np.ndarray, n1: float, n2: float) function:**  
  * This function implements the unified physics logic from Section 2.3.  
  * It takes the incident direction, surface normal, and the two refractive indices (n1​ is the medium the ray is in, n2​ is the medium it is hitting).  
  * It calculates the term c under the square root of the vector Snell's Law formula.  
  * If c\<0, it returns the result of a helper reflection function.  
  * Otherwise, it returns the calculated refracted direction vector.

### **4.4 The gui.app Module: The Simulation Engine**

The main Application class in this module will contain the method that drives the simulation.

* **run\_simulation(self) method (pseudocode):**  
  Python  
  def run\_simulation(self):  
      \# 1\. Get all parameters from Tkinter widgets (with error handling)  
      vertices, material\_n, boundary\_ns, ray\_origin, ray\_dir \= self.get\_params\_from\_gui()

      \# 2\. Create core objects  
      polygon \= Polygon(vertices, material\_n, 0) \# Placeholder default\_n  
      for i, n in enumerate(boundary\_ns):  
          polygon.set\_boundary\_medium(i, n)  
      current\_ray \= Ray(ray\_origin, ray\_dir)

      \# 3\. Initialize path and loop  
      ray\_path \= \[current\_ray.origin\]  
      MAX\_BOUNCES \= 50  
      for \_ in range(MAX\_BOUNCES):  
          \# 4\. Find next intersection  
          intersection\_data \= physics.find\_closest\_intersection(current\_ray, polygon)

          if intersection\_data is None:  
              \# Ray has exited the polygon, trace it for a fixed distance  
              final\_point \= current\_ray.origin \+ current\_ray.direction \* EXIT\_DISTANCE  
              ray\_path.append(final\_point)  
              break

          \# 5\. Process intersection  
          new\_start \= intersection\_data\["intersection\_point"\]  
          ray\_path.append(new\_start)  
          normal \= intersection\_data\["normal"\]  
          segment\_idx \= intersection\_data\["segment\_index"\]

          n1 \= polygon.material\_n  
          n2 \= polygon.boundary\_media\_n\[segment\_idx\]

          \# 6\. Apply physics  
          new\_direction \= physics.calculate\_interaction(current\_ray.direction, normal, n1, n2)

          \# 7\. Update ray for next iteration  
          current\_ray \= Ray(new\_start, new\_direction)

      \# 8\. Plot the results  
      self.plot\_results(polygon, ray\_path)

## **Section 5: Designing the Graphical User Interface (GUI)**

The GUI is the user's entry point to the simulation. It must be intuitive, providing clear controls for all simulation parameters and a clean visualization of the results. Tkinter, combined with Matplotlib, provides all the necessary components.

### **5.1 GUI Layout and Widget Selection**

The main application window should be organized logically, for instance, with a control panel for inputs and a larger area for the plot. Tkinter's .grid() geometry manager is well-suited for this, offering more layout flexibility than .pack().25

* **Control Panel:** This frame will contain all user inputs, organized with ttk.Label widgets.  
  * **Geometry:** A ttk.Combobox for selecting a shape (e.g., "Trapezoid," "Rectangle"). ttk.Entry widgets for dimensions like length, height, and angle. A "Presets" ttk.Combobox would be a powerful addition, allowing users to select common, commercially available crystal configurations (e.g., "Ge 50x20x2mm 45° Trapezoid").26 Selecting a preset would automatically populate the geometry, material, and dimension fields, greatly enhancing usability.  
  * **Materials:** A ttk.Combobox for crystal material ("ZnSe," "Ge," etc.) that fills a ttk.Entry with the corresponding refractive index from the constants module. The entry should remain editable for custom values.  
  * **Boundaries:** A dynamically generated set of ttk.Entry widgets, one for each boundary of the selected geometry, allowing the user to define the external medium for each face.  
  * **Ray:** ttk.Entry widgets for the ray's starting X/Y coordinates and its initial angle in degrees.  
  * **Actions:** A main ttk.Button to "Run Simulation" and another to "Clear" the plot.

### **5.2 Embedding the Matplotlib Visualization**

Embedding a Matplotlib plot into a Tkinter application is a standard procedure and is well-documented.28

1. **Imports:** The necessary classes are imported: from matplotlib.figure import Figure and from matplotlib.backends.backend\_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk.  
2. **Canvas Creation:** Inside the GUI's \_\_init\_\_ method, a Matplotlib Figure and Axes object are created. A FigureCanvasTkAgg widget is then instantiated, linking the Matplotlib figure to the Tkinter parent frame.  
3. **Toolbar:** The NavigationToolbar2Tk is added to the layout. This provides the user with built-in, interactive zoom and pan controls, which are invaluable for closely inspecting the ray's path at interfaces.  
4. **Plotting Method:** A dedicated method, plot\_results(self, polygon, ray\_path), will handle all drawing. Its first step must be to call self.ax.clear() to remove previous plots. It will then plot the polygon outline and the ray path. A crucial line is self.ax.set\_aspect('equal', adjustable='box'), which ensures that the geometry is not distorted and that angles are represented correctly. Finally, self.canvas.draw() must be called to update the display with the new plot.

### **5.3 Connecting Controls to the Simulation**

The command option of the "Run Simulation" button will be bound to the run\_simulation method. This method will be responsible for reading the string values from all ttk.Entry widgets using their .get() method and converting them to the appropriate numeric types (e.g., float). This data conversion step should be wrapped in try-except blocks to gracefully handle potential ValueError exceptions if the user enters non-numeric text, preventing the application from crashing and allowing it to display a helpful error message instead.

## **Section 6: Project Documentation: Instruction.md**

Clear documentation is essential for any software project. The following Instruction.md file should be included in the project's root directory to serve as a comprehensive guide for both users and future developers.

\`\# ATR-MIR Ray Path Simulator

This project provides a Python-based application for simulating the path of a mid-infrared (MIR) ray through various 2D geometries, intended for testing and visualizing configurations for Attenuated Total Reflectance (ATR) spectroscopy.

## **Project Overview**

The application allows a user to:

* Define a 2D crystal geometry (e.g., trapezoid, rectangle) with a specific refractive index.  
* Define the surrounding media at each boundary of the crystal.  
* Specify the starting point and initial direction of an MIR ray (fixed at 5 µm wavelength).  
* Visualize the crystal and the ray path, which is calculated based on the physical laws of reflection and refraction, including Total Internal Reflection (TIR).

## **Installation and Setup**

### **Prerequisites**

* Python 3.8 or newer.

### **Instructions**

1. It is highly recommended to use a virtual environment to manage dependencies. Create one using:bash  
   python \-m venv venv  
   source venv/bin/activate \# On Windows, use venv\\Scripts\\activate  
2. Install the required libraries from the requirements.txt file:  
   Bash  
   pip install \-r requirements.txt

## **How to Use the Simulator**

1. Run the application from the project's root directory:  
   Bash  
   python main.py

2. **Define Geometry:**  
   * Use the "Shape" dropdown to select the crystal geometry.  
   * Enter the dimensions (e.g., Length, Height, Angle) in the corresponding fields.  
   * Alternatively, select a common commercial crystal from the "Presets" dropdown to auto-fill these parameters.  
3. **Set Materials:**  
   * Select a crystal material (e.g., "ZnSe", "Ge") to auto-fill its refractive index, or enter a custom value.  
   * For each boundary of the shape, enter the refractive index of the surrounding medium (e.g., 1.0 for air, 1.325 for water).  
4. **Define Ray:**  
   * Enter the starting X and Y coordinates for the ray.  
   * Enter the initial angle of the ray in degrees (0 degrees points right along the positive x-axis).  
5. **Run Simulation:**  
   * Click the "Run Simulation" button to trace the ray path.  
   * The visualization will update to show the geometry and the calculated path.  
   * Use the toolbar at the bottom of the plot to zoom and pan for detailed inspection.  
   * Click "Clear" to reset the plot.

## **Technical Reference**

### **Physics Engine (atr\_sim/core/)**

The simulation is built on a 2D vector-based ray tracing engine.

* **Ray-Boundary Intersection:** Intersections are calculated by solving a parametric system of equations for the ray and each boundary segment. The engine finds the closest valid intersection in the direction of ray travel.  
* **Optical Laws:** Ray interactions at interfaces are handled by vector forms of Snell's Law and the law of reflection. Total Internal Reflection (TIR) is not a special case; it is naturally handled by the vector Snell's Law equation when its solution becomes imaginary, at which point the reflection logic is triggered. This avoids explicit angle calculations.

### **Coordinate System**

The simulation uses a standard 2D Cartesian coordinate system. The bottom-left of the plotting area is typically near the origin (0,0).

### **File Structure**

* main.py: Main entry point to launch the application.  
* atr\_sim/core/: Contains the core physics and geometry logic. It is independent of the GUI.  
  * geometry.py: Defines the Polygon and Ray data classes.  
  * physics.py: Contains functions for intersection finding and optical calculations.  
* atr\_sim/gui/: Contains the Tkinter-based GUI.  
  * app.py: The main Application class that builds the UI and orchestrates the simulation.  
* atr\_sim/utils/: Contains utility modules.  
  * constants.py: Stores all physical constants, like refractive indices.

## **Advanced Topics and Future Expansion**

This simulator provides a solid foundation that can be extended for more complex analyses.

### **Simulating Multiple Wavelengths (Chromatic Dispersion)**

The current simulation operates at a fixed wavelength of 5 µm. A major enhancement would be to account for chromatic dispersion, where the refractive index (n) changes with wavelength (λ). This can be implemented using **Sellmeier equations**, which are empirical formulas that model the refractive index of a material across a spectrum.

Sellmeier Equation for Germanium (Ge):  
For a wavelength range of 2.0-14.0 µm, the refractive index of Ge can be modeled as 31:  
n2(λ)−1=λ2−BAλ2​+λ2−DCλ2​

(Note: The source provides a three-term equation. For simplicity, a two-term example is shown here. The implementation should use the full formula from the source.)  
Sellmeier Equation for Zinc Selenide (ZnSe):  
For a wavelength range starting from 2.5 µm, the refractive index of ZnSe can be modeled as 32:  
n2(λ)=A+λ2−CBλ2​

where A=4.00, B=1.90, and C=0.113 (with λ in microns).  
**Implementation Path:**

1. Add a "Wavelength" input field to the GUI.  
2. In utils/constants.py, store the Sellmeier coefficients for each material.  
3. Create a function, e.g., get\_refractive\_index(material, wavelength), that calculates n using the appropriate Sellmeier equation.  
4. Modify the run\_simulation method to call this function to get the dynamic refractive indices before passing them to the physics engine.

### **Other Potential Enhancements**

* **3D Simulation:** Extend the vector logic to 3D to model more complex beam paths.  
* **Absorption Modeling:** Incorporate the imaginary part of the refractive index (k) to model the attenuation of the ray's intensity at each reflection.  
* **Output Metrics:** Calculate and display useful metrics like the total path length inside the crystal and the number of reflections on the sample-facing surface.

#### **Works cited**

1. Vector Form of Snell's Law \- StarkEffects.com, accessed August 11, 2025, [https://www.starkeffects.com/snells-law-vector.shtml](https://www.starkeffects.com/snells-law-vector.shtml)  
2. Snell's law \- Wikipedia, accessed August 11, 2025, [https://en.wikipedia.org/wiki/Snell%27s\_law](https://en.wikipedia.org/wiki/Snell%27s_law)  
3. Total internal reflection \- Wikipedia, accessed August 11, 2025, [https://en.wikipedia.org/wiki/Total\_internal\_reflection](https://en.wikipedia.org/wiki/Total_internal_reflection)  
4. Attenuated total reflectance \- Wikipedia, accessed August 11, 2025, [https://en.wikipedia.org/wiki/Attenuated\_total\_reflectance](https://en.wikipedia.org/wiki/Attenuated_total_reflectance)  
5. byjus.com, accessed August 11, 2025, [https://byjus.com/physics/total-internal-reflection/\#:\~:text=Following%20are%20the%20two%20conditions,greater%20than%20the%20critical%20angle.](https://byjus.com/physics/total-internal-reflection/#:~:text=Following%20are%20the%20two%20conditions,greater%20than%20the%20critical%20angle.)  
6. What are the 2 conditions for total internal reflection? \- SZPHOTON, accessed August 11, 2025, [https://szphoton.com/blogs/articles/what-are-the-2-conditions-for-total-internal-reflection](https://szphoton.com/blogs/articles/what-are-the-2-conditions-for-total-internal-reflection)  
7. Snell's Law \-- The Law of Refraction \- UBC Math, accessed August 11, 2025, [https://www.math.ubc.ca/\~cass/courses/m309-01a/chu/Fundamentals/snell.htm](https://www.math.ubc.ca/~cass/courses/m309-01a/chu/Fundamentals/snell.htm)  
8. total internal reflection – TIR, critical angle, evanescent wave, frustrated, FTIR \- RP Photonics, accessed August 11, 2025, [https://www.rp-photonics.com/total\_internal\_reflection.html](https://www.rp-photonics.com/total_internal_reflection.html)  
9. Snell's law in vector form \- refraction \- Physics Stack Exchange, accessed August 11, 2025, [https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form](https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form)  
10. Optimized Snell's Law Refraction \- Ryan Brucks, accessed August 11, 2025, [https://shaderbits.com/blog/optimized-snell-s-law-refraction](https://shaderbits.com/blog/optimized-snell-s-law-refraction)  
11. Crystal Selection for ATR \- PIKE Technologies, accessed August 11, 2025, [https://www.piketech.com/atr-crystal-selection/](https://www.piketech.com/atr-crystal-selection/)  
12. Attenuated Total Reflectance (ATR) \- Bruker, accessed August 11, 2025, [https://www.bruker.com/en/products-and-solutions/infrared-and-raman/ft-ir-routine-spectrometer/what-is-ft-ir-spectroscopy/atr-attenuated-total-reflectance.html](https://www.bruker.com/en/products-and-solutions/infrared-and-raman/ft-ir-routine-spectrometer/what-is-ft-ir-spectroscopy/atr-attenuated-total-reflectance.html)  
13. Germanium Material | Infrared Optical Materials, accessed August 11, 2025, [https://wavelength-oe.com/germanium/](https://wavelength-oe.com/germanium/)  
14. F.A.Q. – PIKE Technologies, accessed August 11, 2025, [https://www.piketech.com/frequently-asked-questions/](https://www.piketech.com/frequently-asked-questions/)  
15. Refractive index \- Wikipedia, accessed August 11, 2025, [https://en.wikipedia.org/wiki/Refractive\_index](https://en.wikipedia.org/wiki/Refractive_index)  
16. Optical properties of water and ice \- Wikipedia, accessed August 11, 2025, [https://en.wikipedia.org/wiki/Optical\_properties\_of\_water\_and\_ice](https://en.wikipedia.org/wiki/Optical_properties_of_water_and_ice)  
17. 40 Top Python Libraries Every Data Scientist Should Know in 2025, accessed August 11, 2025, [https://www.stxnext.com/blog/most-popular-python-scientific-libraries](https://www.stxnext.com/blog/most-popular-python-scientific-libraries)  
18. Scientific Computing with Python \- GeeksforGeeks, accessed August 11, 2025, [https://www.geeksforgeeks.org/python/scientific-computing-with-python/](https://www.geeksforgeeks.org/python/scientific-computing-with-python/)  
19. graphics \- How do you check for intersection between a line ..., accessed August 11, 2025, [https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin](https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin)  
20. Finding Intersection Points \- A Trig-less Line of Sight Algorithm in Two Dimensions, accessed August 11, 2025, [https://basstabs.github.io/2d-line-of-sight/Final.html](https://basstabs.github.io/2d-line-of-sight/Final.html)  
21. What is the simplest to implement line segment intersect algorithm? \- Reddit, accessed August 11, 2025, [https://www.reddit.com/r/algorithms/comments/9moad4/what\_is\_the\_simplest\_to\_implement\_line\_segment/](https://www.reddit.com/r/algorithms/comments/9moad4/what_is_the_simplest_to_implement_line_segment/)  
22. Matplotlib — Visualization with Python, accessed August 11, 2025, [https://matplotlib.org/](https://matplotlib.org/)  
23. List of Python GUI Library and Packages \- GeeksforGeeks, accessed August 11, 2025, [https://www.geeksforgeeks.org/python/python3-gui-application-overview/](https://www.geeksforgeeks.org/python/python3-gui-application-overview/)  
24. Plotly Python Graphing Library, accessed August 11, 2025, [https://plotly.com/python/](https://plotly.com/python/)  
25. Embedding a MatPlotLib Graph in Tkinter \[.grid method\], and Customizing MatPlotLib's Navigation Toolbar \- Stack Overflow, accessed August 11, 2025, [https://stackoverflow.com/questions/59550783/embedding-a-matplotlib-graph-in-tkinter-grid-method-and-customizing-matplotl](https://stackoverflow.com/questions/59550783/embedding-a-matplotlib-graph-in-tkinter-grid-method-and-customizing-matplotl)  
26. atr optical elements, accessed August 11, 2025, [https://www.spectral-systems.com/wp-content/uploads/2014/06/ATR-OPTICAL-ELEMENTS-6-19-17.pdf](https://www.spectral-systems.com/wp-content/uploads/2014/06/ATR-OPTICAL-ELEMENTS-6-19-17.pdf)  
27. Germanium ATR Trapezoid Prism 50x20x2mmthk \- Knight Optical, accessed August 11, 2025, [https://www.knightoptical.com/stock/germanium-atr-trapezoid-prism-50x20x2mmthk](https://www.knightoptical.com/stock/germanium-atr-trapezoid-prism-50x20x2mmthk)  
28. How to embed a Matplotlib graph to your Tkinter GUI \- PythonProgramming.net, accessed August 11, 2025, [https://pythonprogramming.net/how-to-embed-matplotlib-graph-tkinter-gui/](https://pythonprogramming.net/how-to-embed-matplotlib-graph-tkinter-gui/)  
29. Embed in Tk — Matplotlib 3.10.5 documentation, accessed August 11, 2025, [https://matplotlib.org/stable/gallery/user\_interfaces/embedding\_in\_tk\_sgskip.html](https://matplotlib.org/stable/gallery/user_interfaces/embedding_in_tk_sgskip.html)  
30. How to embed Matplotlib charts in Tkinter GUI? \- GeeksforGeeks, accessed August 11, 2025, [https://www.geeksforgeeks.org/python/how-to-embed-matplotlib-charts-in-tkinter-gui/](https://www.geeksforgeeks.org/python/how-to-embed-matplotlib-charts-in-tkinter-gui/)  
31. Refractive index of Ge (Germanium) \- Burnett \- RefractiveIndex.INFO, accessed August 11, 2025, [https://refractiveindex.info/?shelf=main\&book=Ge\&page=Burnett](https://refractiveindex.info/?shelf=main&book=Ge&page=Burnett)  
32. Refractive Index of ZnSe, ZnTe, and CdTe\*, accessed August 11, 2025, [https://pubs.aip.org/aip/jap/article-pdf/35/3/539/18331474/539\_1\_online.pdf](https://pubs.aip.org/aip/jap/article-pdf/35/3/539/18331474/539_1_online.pdf)