# RIM

This is the code repository for the Python implementation of the Ray Intersection Model (RIM) - for analyzing spatial relationships (e.g. betweenness, intervisibility) of triplets of spatial objects.

RIM is used to describe the spatial relationship of three spatial objects. It evaluates rays cast between two peripheral spatial objects, and their topological relations with the core object to determine its relative position with respect to the peripheral objects. The interpretation of relationships described with RIM (e.g. core object is between/not between the peripheral objects) is left to the user and application context.

## Getting started

### Prerequisites

The package has been tested with the following requirements:
- Python >= 3.5
- numpy >= 1.17.4
- Shapely >= 1.6.4

Note that Shapely requires GDAL to be installed. More information on how to install GDAL can be found at [https://gdal.org/download.html](https://gdal.org/download.html), and at [https://pypi.org/project/GDAL/](https://pypi.org/project/GDAL/).

```bash
$ pip3 install numpy>=1.17.4
$ pip3 install Shapely>=1.6.4
```

### Installing

The **rim** Python package can be installed via pip

```bash
$ pip3 install rim
```

### Usage

Here is an example of a simple Python script where RIM is used to analyze spatial relationships of three polygon geometries.

```python
#!/usr/bin/python3
from rim import RIM

A = 'POLYGON((1 1, 1 5, 3 5, 3 1, 1 1))'
B = 'POLYGON((7 1, 7 5, 9 5, 9 1, 7 1))'
O = 'POLYGON((4 0, 4 6, 6 6, 6 0, 4 0))'

rimobject = RIM(A,B,O)

print(rimobject.rim)
# prints 'RIM 13'
print(rimobject.rim_matrix)
# prints [[1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]]
print(rimobject.rays)
# prints ['ray5']

if 'ray1' in rimobject.rays or 'ray8' in rimobject.rays:
    B_visible_from_A = 'is'
else:
    B_visible_from_A = 'is not'

print('According to the RIM model, B %s visible from A!' %B_visible_from_A)
# prints 'According to the RIM model, B is not visible from A!'
```

Here is a docstring of the RIM class, explaining how the RIM object is created and which attributes it has:

    A RIM object is created from three spatial objects - A, B, and O,
    provided as WKT strings. A and B are peripheral objects, and O is
    a core object which is being analyzed for being in between A and B.

    Attributes
    ----------
    A : str
        a WKT representation of peripheral object A
    B : str
        a WKT representation of peripheral object B
    O : str
        a WKT representation of core object O
    rays : list of str
        a list of all rays that exist between peripheral
        objects A and B in this specific RIM scenario;
        ray in this case is a straight line that shares
        exactly one point with each A and B;
        rays are encoded as ray1-ray8 since there are 8
        distinct rays that can theoretically occur
    extreme_rays : list of str
        a list of extreme rays that exist between peripheral
        objects A and B in this specific RIM scenario,
        extreme rays are encoded as ray1-ray8, and all extreme
        rays are present in the 'rays' attribute as well
    ray_area : list of dict
        a list of dictionary entries where each entry describes a part of the
        ray area with its WKT geometry (key 'ray area') and a list of extreme
        rays WKT geometries (key 'extreme rays'); ray area is the area between
        peripheral objects A and B that is covered by all rays that exist
        between them; in case of the ray area with a single-part geometry
        there will be only 1 element in the 'ray_area'
    rim_matrix : 2D numpy array
        a matrix representation of this specific RIM scenario;
        the 3 columns represent the interior, boundary, and exterior
        of the core object O, while the 4 rows represent interiors
        of all rays, boundaries of all rays, interiors of extreme rays,
        and boundaries of extreme rays; the value 0 stands for 'none of
        the rays have this intersection', 0-1 stands for 'some of the rays
        have this intersection', and 1 stands for 'all of the rays have
        this intersection'
    rim : str
        a string representation of the rim in this specific scenario (A,B,O);
        if the resulting rim is one of the 29 rims that were previously defined
        in studies, the value will be in range 'rim1' - 'rim29', otherwise the
        value will be string representation of the rim_matrix



## Authors

* **Ivan Majic** - *Initial work*

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## References

The RIM model is featured in the following academic publications:

Ivan Majic, Elham Naghizade, Stephan Winter, and Martin Tomko. RIM: A ray intersection model
for the analysis of the between relationship of spatial objects in a 2D plane. International Journal
of Geographical Information Science, pages 1-26, 2020. [doi:10.1080/13658816.2020.1778002](https://doi.org/10.1080/13658816.2020.1778002)

Ivan Majic, Elham Naghizade, Stephan Winter, and Martin Tomko. There is no way! Ternary
qualitative spatial reasoning for error detection in map data. Transactions in GIS, 2021. [doi:10.1111/tgis.12765](https://doi.org/10.1111/tgis.12765)

Ningran Xu, Ivan Majic, and Martin Tomko. Perceptions of Qualitative Spatial Arrangements of
Three Objects. 15th International Conference on Spatial Information Theory (COSIT 2022), pages
1-14, Kobe, Japan, 2022. [doi:10.4230/LIPIcs.COSIT.2022.9](https://doi.org/10.4230/LIPIcs.COSIT.2022.9)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
