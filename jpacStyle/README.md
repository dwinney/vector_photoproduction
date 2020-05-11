# jpacStyle
Library for making plots in C++/[ROOT](https://root.cern/) for the JPAC Collaboration.

Compile the library with the following command:
```bash
mkdir build
cd build
cmake ..
make
```
This will make `libJpacStyle.a` in the build directory which can be linked to other code to have access to the header files.

### jpacGraph1D
This object allows you to easily make one-dimensional plots according to the style and with minimal ROOT syntax.

Basic usage is:
```c++
// Initialize the object
jpacGraph1D* my_1Dplotter = new jpacGraph1D();

// Add x and f(x) data in a vector and a string for the legend label
my_1Dplotter->AddEntry(vector<double>, vector<double>, string);

// Plot to file
my_1Dplotter->Plot(string);
```
You can add multiple curves (up to 10, because we only have 10 colors defined).

Additional customization can be set up with the previous functions before plotting:
```c++
// Manually place location of the legend with relative coordinate of the bottom left vertex x and y
my_1Dplotter->SetLegend(double, double);

// Label (units) for x and y axes and upper and lower bounds
my_1Dplotter->SetXaxis(string, double, double);
my_1Dplotter->SetYaxis(string, double, double);
```
Axes labels are TLatex objects and thus follow the same syntax for mathematical symbols (see [doc](https://root.cern.ch/doc/master/classTLatex.html)). For an example script using this object see [bessel.cpp](./examples/bessel.cpp).

### jpacGraph1Dc
This is operates nearly identical to the above but allows for plotting complex valued function defined on the real line. All the functions available in `jpacGraph1D` are present here except all with the possibility of accepting complex vectors when adding entries:
```c++
jpacGraph1Dc* my_1Dcplotter = new jpacGraph1Dc();
my_1Dcplotter->AddEntry(vector<double>, vector<complex<double>>, string);
```
Output is the Real and Imaginary parts plotted seperately in the same file (See [hankel.cpp](./examples/hankel.cpp)).

### jpacGraph2D
This object allows yo uto make two-dimensional plots according to JPAC color scheme with minimal ROOT interfacing.

Basic usage is even easier than above since only one function is plottable at once:
```c++
// Initialize
jpacGraph2D* my_2Dplotter = new jpacGraph2D();

// Add in the x, y, and z data as vectors
my_2Dplotter->AddData(vector<double>, vector<double>, vector<double>);

// Plot to file
my_2Dplotter->Plot(string);
```
As above, additional customization of the axes is available through:
```c++
// Axis label, minimum and maximum values
my_2Dplotter->SetXaxis(string, double, double);
my_2Dplotter->SetYaxis(string, double, double);
```
See the example executable [2dgaussian.cpp](./examples/2dgaussian.cpp).

### importStyle.C
Alternatively if you want to make the plots manually through ROOT, this macro imports the jpacStyle and jpacColors. Also adds a function `AddLogo()` which adds the collaboration logo in the upper right corner.
Simply load when opening ROOT:
```bash
root -l importStyle.C
```
