#!/usr/bin/env python3
import shapely.wkt, shapely.geometry
from shapely.ops import polygonize, split, unary_union
import numpy as np
from itertools import combinations, product, chain

class RIM:
    """
    Ray Intersection Model (RIM) is used to describe the spatial relationship
    of three spatial objects. It evaluates rays cast between two peripheral
    spaital objects, and their topological relations with the core object to
    determine its relative position with respect to the peripheral objects.
    The interpretation of relationships described with RIM (e.g. core object is
    between/not between the peripheral objects) is left to the user and
    application context.

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

    """

    def __init__(self, A: str, B: str, O: str):
        """
        Parameters
        ----------
        A : str
            a WKT geometry representation of peripheral object A
        B : str
            a WKT geometry representation of peripheral object B
        O : str
            a WKT geometry representation of core object O

        Example
        -------
        Compute RIM of three polygon geometries

            >>> A = 'POLYGON((1 1, 1 5, 3 5, 3 1, 1 1))'
            >>> B = 'POLYGON((7 1, 7 5, 9 5, 9 1, 7 1))'
            >>> O = 'POLYGON((4 0, 4 6, 6 6, 6 0, 4 0))'
            >>> rimobject = RIM(A,B,O)
            >>> rimobject.rim
            'RIM 13'
            >>> rimobject.rim_matrix
            [[1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]]
            >>> rimobject.rays
            ['ray5']
        """
        self.A = A
        self.B = B
        self.O = O
        self.rays = []
        self.extreme_rays = []
        self.rim_matrix = None
        self.ray_area = self._rayArea()
        self.rim = self._RIM()

    def _checkValid(self, *geoms: shapely.geometry) -> None:
        """
        Checks if geometry is valid. If yes, resumes without returning a
        value. If not, raises a ValueError exception.
        """
        for geom in geoms:
            if geom.is_valid:
                pass
            else:
                raise ValueError('Geometry is invalid:\n%s' % geom.wkt)

    def _getAllPoints(self, O: shapely.geometry) -> list:
        """
        Function for getting all the points of a geometry,
        regardless if it is multi or single geometry
        """

        def pointsFromType(O) -> list:
            """
            Function for getting all the points of a geometry
            if it is a Point, MultiPoint, LineString, MultiLineString,
            Polygon, MultiPolygon
            """
            otype = O.type
            if otype == 'MultiPolygon':
                points = [x for p in O for x in list(
                    shapely.geometry.MultiPoint(p.exterior.coords))]
            elif otype == 'Polygon':
                points = list(shapely.geometry.MultiPoint(O.exterior.coords))
            elif otype == 'MultiLineString':
                points = [x for p in O for x in list(
                    shapely.geometry.MultiPoint(p.coords))]
            elif otype == 'LineString':
                points = list(shapely.geometry.MultiPoint(O.coords))
            elif otype == 'MultiPoint':
                points = list(O)
            elif otype == 'Point':
                points = [O]
            return points

        if O.type == 'GeometryCollection':
            # GeometryCollection objects need to be iterated through
            # in order to get the points of every feature within,
            # especially because geometry types can vary
            # (e.g. points and multipoints)
            pointslist = [pointsFromType(geom) for geom in O]
            points = list(chain.from_iterable(pointslist))
        else:
            # just get the points of a single object O
            points = pointsFromType(O)
        # convert list of points to a set 
        # in order to remove duplicates,
        # also prepare (x,y) touples for each point
        # - this is needed because points might be XYZ
        points = set([(p.x, p.y) for p in points])
        points = [shapely.geometry.Point(p) for p in points]
        return points


    def _extendLine(self, p1x, p1y, p2x, p2y, both_ways=False, scale=10) \
        -> shapely.geometry.linestring:
        """
        Creates a line extrapolated in p1->p2 direction,
        or a line extrapolated in both directions if both_ways=True.
        Default extrapolation scale is 10.
        """
        if both_ways:
            # 
            ax, ay = (
                p2x + scale * (p1x - p2x),
                p2y + scale * (p1y - p2y)
                )
        else:
            ax = p1x
            ay = p1y
        bx = p1x + scale * (p2x - p1x)
        by = p1y + scale * (p2y - p1y)
        extended_linestring = shapely.geometry.LineString([(ax, ay), (bx, by)])
        return extended_linestring


    def _rayArea(self) -> tuple:
        """
        Creates the Ray area between two objects A and B.
        The objects have to be in WKT format.
        """

        try:
            A = shapely.wkt.loads(self.A)
            B = shapely.wkt.loads(self.B)

            # Check if peripheral objects A and B are
            # overlapping - which is not allowed
            # if they are, return the text message
            # Overlapping peripheral objects as self.rim
            if A.intersects(B) and not A.touches(B):
                msg = 'Overlapping peripheral objects'
                raise Exception(msg)

            # 1. Make a convex hull of the two objects
            geomcol = shapely.geometry.GeometryCollection([A, B])
            chullAB = geomcol.convex_hull

            # 2. Difference of the chull and the original objects
            # this needs to be handled differently for Polygon and
            # Line objects because if a portion of the line is going
            # through the chull, nothing will happen with the
            # chull_diff.differnce(geom) function
            chull_diff = chullAB
            for geom in (A, B):
                if geom.geom_type == 'Polygon':
                    chull_diff = chull_diff.difference(geom)
                    # will result in a polygon or a multipolygon
                elif geom.geom_type == 'LineString':
                    chull_diff = split(chull_diff, geom)
                    # will always result in a geometrycollection

            self._checkValid(geomcol, chullAB, chull_diff)

            # 3. Keep only the parts that are touching both original objects.
            # If the chull is a multipolygon, this will remove any polygons
            # which are touching only one of the objects A and B
            if hasattr(chull_diff, '__iter__'):
                # hasattr check returns true for both multi-geometries and
                # geometrycollections
                new_chull_diff = []
                for geom in chull_diff:
                    if geom.intersects(A) and geom.intersects(B):
                        new_chull_diff.append(geom)
                chull_diff = shapely.geometry.GeometryCollection(new_chull_diff)
            else:
                chull_diff = shapely.geometry.GeometryCollection([chull_diff])
            
            # get the portion chull_diff that is touching object A
            # this is usually a line
            rayArea_A = chull_diff.intersection(A)
            rayArea_B = chull_diff.intersection(B)

            points_A = self._getAllPoints(rayArea_A)
            points_B = self._getAllPoints(rayArea_B)

            self._checkValid(chull_diff, rayArea_A, rayArea_B)

            # 4. Find parts of the objects A and B that are touching
            # the chull_diff, but are not visible from the other object
            #
            # For this we create lines that connect points of A with points
            # of B, but only keep the lines that are touching both objects,
            # if the line crosses any of the objects then it does not represent
            # a visibility line
            mline = []
            for p1 in points_A:
                for p2 in points_B:
                    line = shapely.geometry.LineString([p1, p2])
                    if line.touches(A) and line.touches(B):
                        mline.append(line)
            mline = shapely.geometry.MultiLineString(mline)

            self._checkValid(mline)
            
            
            new_chull_diff = []
            for chull_diff_geom in chull_diff:
                boundary_diff_A_and_B = (
                    chull_diff_geom.boundary
                    .difference(A)
                    .difference(B)
                    )
                boundary_not_touching_both_A_and_B = []
                if hasattr(boundary_diff_A_and_B, '__iter__'):
                    for geom in boundary_diff_A_and_B:
                        if geom.disjoint(A) or geom.disjoint(B):
                            boundary_not_touching_both_A_and_B.append(geom)
                elif (
                    boundary_diff_A_and_B.disjoint(A)
                    or
                    boundary_diff_A_and_B.disjoint(B)
                ):
                    (
                        boundary_not_touching_both_A_and_B
                        .append(boundary_diff_A_and_B)
                    )
                boundary_not_touching_both_A_and_B = unary_union(
                    boundary_not_touching_both_A_and_B
                    )
                if not boundary_not_touching_both_A_and_B.is_empty:
                    lines_to_polygonize = mline
                    for peripheral_geom in (A,B):
                        if peripheral_geom.type == 'LineString':
                            lines_to_polygonize = (
                                lines_to_polygonize.union(peripheral_geom)
                            )
                        elif peripheral_geom.type == 'Polygon':
                            lines_to_polygonize = (
                                lines_to_polygonize
                                .union(peripheral_geom.boundary)
                            )
                    polygonized_pol = unary_union(
                        list(
                            polygonize(lines_to_polygonize)
                            )
                        )
                    new_chull_diff_geom = shapely.geometry.Polygon(
                        polygonized_pol.exterior,
                        chull_diff_geom.interiors
                        )
                    new_chull_diff.append(new_chull_diff_geom)
                else:
                    new_chull_diff.append(chull_diff_geom)
            chull_diff = shapely.geometry.GeometryCollection(new_chull_diff)

            # Find the points of A and B which are not visible from the other
            # object - i.e. not reached by one of the lines previously analysed
            unreached_points = []
            polygon1 = []
            for pol in chull_diff:
                polygon_exterior_points = []
                for point in shapely.geometry.MultiPoint(pol.exterior.coords):
                    if mline.intersects(point):
                        polygon_exterior_points.append(point)
                    else:
                        unreached_points.append(point)
                polygon_interiors = []
                for interior in pol.interiors:
                    interior_points = []
                    for point in shapely.geometry.MultiPoint(interior.coords):
                        if mline.intersects(point):
                            interior_points.append(point)
                        else:
                            unreached_points.append(point)
                    polygon_interiors.append(interior_points)
                single_polygon = shapely.geometry.Polygon(
                    [(p.x, p.y) for p in polygon_exterior_points],
                    [[(p.x, p.y) for p in polygon_interior_points]
                     for polygon_interior_points in polygon_interiors]
                    )
                polygon1.append(single_polygon)
            polygon1 = shapely.geometry.GeometryCollection(polygon1)
            unreached_points = shapely.geometry.MultiPoint(unreached_points)

            self._checkValid(unreached_points, polygon1)

            if hasattr(chull_diff, '__iter__'):
                polygons_remaining = chull_diff
            else:
                polygons_remaining = [chull_diff]
            polygons_remaining = []
            for cd in chull_diff:
                for pol in polygon1:
                    cd = cd.difference(pol)
                polygons_remaining.append(cd)
            polygons_remaining = shapely.geometry.GeometryCollection(
                polygons_remaining
                )

            areas_to_add = []
            for polygon_remaining in polygons_remaining:
                area_to_add = []
                for line in mline:
                    if line.intersects(polygon_remaining):
                        p1, p2 = shapely.geometry.MultiPoint(line.coords)
                        if p1.intersects(polygon_remaining):
                            # if the first point of the line (which is always
                            # defined with 2 points) is touching the remaining
                            # area, extend the line in the direction p2->p1
                            extended_line = self._extendLine(
                                p2.x, p2.y, p1.x, p1.y)
                        else:
                            # otherwise, extend the line in the direction p1->p2
                            extended_line = self._extendLine(
                                p1.x, p1.y, p2.x, p2.y)
                        # split the polygon in 2 by the extended line
                        polygon_remaining_split = split(
                            polygon_remaining,
                            extended_line
                            )

                        for pol in polygon_remaining_split:
                            if pol.disjoint(unreached_points):
                                areas_to_add.append(pol)

            final_ray_area = []
            for pol in polygon1:
                for area in areas_to_add:
                    if pol.intersects(area):
                        pol = pol.union(area)
                final_ray_area.append(pol)
            
            final_ray_area = shapely.geometry.GeometryCollection(final_ray_area)
            self._checkValid(final_ray_area)


            # Now I have to calculate the extreme rays
            # Finding the extreme rays between the final ray area and objects A
            # and B. Each polygon in a multipolygon ray_area will be considered
            # separately
            ray_area_with_extreme_rays = []
            for geom in final_ray_area:
                extreme_rays = []
                # Scenario 1. extreme rays that are found by taking the Ra
                # boundary and differencing it from A and B
                extreme_rays_1 = (geom.boundary.difference(A)).difference(B)
                if extreme_rays_1.type == 'MultiLineString':
                    for ggeom in extreme_rays_1:
                        if ggeom.intersects(A) and ggeom.intersects(B):
                            extreme_rays.append(ggeom.wkt)
                elif (
                    extreme_rays_1.intersects(A)
                    and extreme_rays_1.intersects(B)):
                    extreme_rays.append(extreme_rays_1.wkt)
                # Allow for extreme rays to be not touching both objects (cases
                # of U shaped objects) - in my cases of polygons inside u-shaped
                # lines, the extreme rays will be found as a 3-point line that
                # is touching only 1 peripheral object -- split this line into
                # two and append
                elif (
                    extreme_rays_1.type == 'LineString'
                    and len(extreme_rays_1.coords) == 3):
                    extreme_rays.append(extreme_rays_1.wkt) 
                # Scenario 2. extreme rays (points)
                # where A, B and Ra all meet each other
                extreme_rays_2 = geom.intersection(A).intersection(B)
                if extreme_rays_2.type == 'MultiLineString':
                    for ggeom in extreme_rays_2:
                        if not ggeom.is_empty:
                            extreme_rays.append(ggeom.wkt)
                elif not extreme_rays_2.is_empty:
                    extreme_rays.append(extreme_rays_2.wkt)
                # Add the part of the ray area together with the related extreme
                # rays to the result
                ray_area_with_extreme_rays.append(
                    {
                    'ray area': geom.wkt,
                    'extreme rays': extreme_rays
                    }
                ) 

            return (True, ray_area_with_extreme_rays)

        except Exception as e:
            raise e


    def _RIM(self) -> str:
        """
        Calculates RIM between two objects A and B.
        The objects have to be in WKT format.
        """

        # there are 8 distinct rays that can occur between A, B, and O
        # they are explained in the paper and should be shown in the readme
        # TODO this list is not complete for LineString and Point geometries
        _ray_types = {
            'ray1' : {
                'Polygon' : np.array([[0,0,1],[0,0,1]]),
                'LineString' : np.array([[0,0,1],[0,0,1]])
                },
            'ray2' : {
                'Polygon' : np.array([[0,1,1],[0,0,1]]),
                'LineString' : np.array([[0,1,1],[0,0,1]]),
                },
            'ray3' : {
                'Polygon' : np.array([[0,1,1],[0,1,1]]),
                'LineString' : np.array([[0,0,1],[0,1,1]]),
                },
            'ray4' : {
                'Polygon' : np.array([[0,1,0],[0,1,0]]),
                # 'LineString' : np.array([[],[]]),
                },
            'ray5' : {
                'Polygon' : np.array([[1,1,1],[0,0,1]]),
                'LineString' : np.array([[1,0,1],[0,0,1]]),
                },
            'ray6' : {
                'Polygon' : np.array([[1,1,1],[0,1,1]]),
                'LineString' : np.array([[0,0,1],[1,0,1]]),
                },
            'ray7' : {
                'Polygon' : np.array([[1,0,0],[0,1,0]]),
                'LineString' : np.array([[1,0,0],[0,1,0]]),
                },
            'ray8' : {
                'Polygon' : np.array([[0,0,1],[0,1,1]]),
                # 'LineString' : np.array([[],[]]),
                },
        }

        # there are 30 distinct RIM matrices I have come across until now in my
        # experiments they are explained in the paper and should be shown in the
        # readme
        _known_rims = {
            '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 1',
            '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 2',
            '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 3',
            '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 4',
            '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]' :
            'RIM 5',
            '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]' :
            'RIM 6',
            '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 7',
            '[0.5, 0.5, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0]' :
            'RIM 8',
            '[0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0]' :
            'RIM 9',
            '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.5, 1.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 10',
            '[0.5, 1.0, 1.0, 0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.0, 1.0, 1.0]' :
            'RIM 11',
            '[0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0]' :
            'RIM 12',
            '[1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 13',
            '[1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0]' :
            'RIM 14',
            '[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]' :
            'RIM 15',
            '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 16',
            '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0]' :
            'RIM 17',
            '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5]' :
            'RIM 18',
            '[0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 19',
            '[0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]' :
            'RIM 20',
            '[0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5]' :
            'RIM 21',
            '[0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0]' :
            'RIM 22',
            '[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 23',
            '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 24',
            '[1.0, 1.0, 1.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.5, 1.0]' :
            'RIM 25',
            '[1.0, 1.0, 1.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 26',
            '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0]' :
            'RIM 27',
            '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 28',
            '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]' :
            'RIM 29', # a line crossing the ray area and both extreme rays
            '[0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5]' :
            'RIM 30', # disjoint when one extreme ray is a point
        }

        try:
            if self.ray_area[0] is False:
                # Ray area could not be calculated for some reason
                # Print the reason out as self.rim
                msg = self.ray_area[1]
                return msg
            else:
                all_ray_areas = self.ray_area[1]
            
            A = shapely.wkt.loads(self.A)
            B = shapely.wkt.loads(self.B)
            O = shapely.wkt.loads(self.O)

            final_result = []
            for ray_area in all_ray_areas:
                Ra = shapely.wkt.loads(ray_area['ray area'])
                rays = []
                rays_txt = []
                extreme_rays = []

                if O.disjoint(Ra):
                    rays.append(_ray_types['ray1'])

                else:
                    # I am removing 'Multi' from the geometry type string, so
                    # that I don't have to add a Multi-variant to every entry
                    # in the _ray_types dict, and the procedure should be the
                    # same for Multi and single types
                    Otype = O.geom_type.replace('Multi', '')
                    intersection = O.intersection(Ra)

                    # difference of polygon and line doesn't do anything, so use
                    # split in that case
                    if Otype == 'LineString':
                        difference = split(Ra, O)
                        # will always result in a geometrycollection
                    else:
                        difference = Ra.difference(O)
                    diff_touches_AB = False
                    if hasattr(difference, '__iter__'):
                        for geom in difference:
                            if geom.intersects(A) and geom.intersects(B):
                                diff_touches_AB = True

                    elif difference.intersects(A) and difference.intersects(B):
                        diff_touches_AB = True

                    if Otype == 'MultiPolygon':

                        def not_crosses(extended_line, O) -> bool:
                            """
                            Returns True if the extended_line does not cross O
                            """
                            return extended_line.relate_pattern(O, 'F*T******')

                        def check_lines_for_ray1(p1list, p2list = None) -> bool:
                            if p2list:
                                myiter = product(p1list, p2list)
                                both_ways = False
                            else:
                                myiter = combinations(p1list, 2)
                                both_ways = True
                            for p1,p2 in myiter:
                                extended_line = self._extendLine(
                                    p1.x, p1.y, p2.x, p2.y, both_ways,
                                    scale=100)
                                if (not_crosses(extended_line, O)
                                and extended_line.intersects(A)
                                and extended_line.intersects(B)):
                                    return(True)

                        # visibility test for multipolygon core objects
                        O_points = [
                            p for p in self._getAllPoints(O) if p.intersects(Ra)
                            ]

                        if check_lines_for_ray1(O_points):
                            rays.append(_ray_types['ray1'][Otype])
                            rays_txt.append('ray1')
                            rays.append(_ray_types['ray2'][Otype])
                            rays_txt.append('ray2')
                        
                        elif check_lines_for_ray1(
                            self._getAllPoints(Ra), O_points):
                            rays.append(_ray_types['ray1'][Otype])
                            rays_txt.append('ray1')
                            rays.append(_ray_types['ray2'][Otype])
                            rays_txt.append('ray2')

                        elif check_lines_for_ray1(
                            self._getAllPoints(A.intersection(Ra)), O_points):
                            rays.append(_ray_types['ray1'][Otype])
                            rays_txt.append('ray1')
                            rays.append(_ray_types['ray2'][Otype])
                            rays_txt.append('ray2')

                        elif check_lines_for_ray1(
                            self._getAllPoints(B.intersection(Ra)), O_points):
                            rays.append(_ray_types['ray1'][Otype])
                            rays_txt.append('ray1')
                            rays.append(_ray_types['ray2'][Otype])
                            rays_txt.append('ray2')

                    elif diff_touches_AB:
                        rays.append(_ray_types['ray1'][Otype])
                        rays_txt.append('ray1')
                        if not difference.equals(Ra):
                            rays.append(_ray_types['ray2'][Otype])
                            rays_txt.append('ray2')
                    
                    if O.touches(Ra):
                        rays.append(_ray_types['ray1'][Otype])
                        rays_txt.append('ray1')

                    if (diff_touches_AB and intersection.geom_type != 'Point'
                    and (intersection.intersects(A)
                        or intersection.intersects(B))):
                        rays.append(_ray_types['ray3'][Otype])
                        rays_txt.append('ray3')

                    if (diff_touches_AB
                    and intersection.intersects(A)
                    and intersection.intersects(B)):
                        rays.append(_ray_types['ray4'][Otype])
                        rays_txt.append('ray4')
                    
                    if ((O.within(Ra) or O.overlaps(Ra)) 
                        and (difference.intersects(A) 
                            and difference.intersects(B))):
                        rays.append(_ray_types['ray5'][Otype])
                        rays_txt.append('ray5')

                    if (not difference.equals(Ra)
                    and (difference.intersects(A) or difference.intersects(B))
                    and (intersection.intersects(A)
                        or intersection.intersects(B))):
                        rays.append(_ray_types['ray6'][Otype])
                        rays_txt.append('ray6')

                    if (intersection.intersects(A)
                    and intersection.intersects(B)):
                        rays.append(_ray_types['ray7'][Otype])
                        rays_txt.append('ray7')

                for extreme_ray in ray_area['extreme rays']:
                    extreme_ray = shapely.wkt.loads(extreme_ray)
                    relation = extreme_ray.relate(O)
                    relation = [
                        1 if x in ('0','1','2') else 0 for x in list(relation)
                        ]
                    nprelation = np.array([relation[:3], relation[3:6]])
                    extreme_rays.append(nprelation)
                    rays.append(nprelation)

                extreme_rays_txt = []
                inverse_ray_types_for_otype = {
                    str(_ray_types[ray].get(Otype)) : ray
                    for ray in _ray_types.keys()
                }
                for r in extreme_rays:
                    extreme_rays_txt.append(
                        inverse_ray_types_for_otype.get(str(r))
                    )
                self.extreme_rays = list(set(extreme_rays_txt))

                self.rays = list(set(rays_txt + extreme_rays_txt))
                
                RIM_part = sum(rays + extreme_rays)/len(rays + extreme_rays)
                if len(extreme_rays) == 0:
                    RIM_part2 = np.array([[0,0,0],[0,0,0]])
                else:
                    RIM_part2 = sum(extreme_rays)/len(extreme_rays)
                RIM_all = np.concatenate((RIM_part, RIM_part2))
                # normalize all values between 0 and 1 (triangle) to 0.5
                RIM_all[(RIM_all>0) & (RIM_all<1)] = 0.5
                final_result.append(list(RIM_all.flatten()))

            self.rim_matrix = final_result

            result = [_known_rims.get(str(rim),rim) for rim in final_result]

            if len(result) > 1:
                result = [x for x in result if x != 'RIM 23']
            if len(result) == 1:
                result = result[0]
            return str(result)
        except Exception as e:
            raise e