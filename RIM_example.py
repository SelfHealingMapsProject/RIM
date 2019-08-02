from RIM_translate import RIM_translate
from RIM import RIM


A = 'Polygon ((0.38447265664216346 0.32246093782891094, 0.38447265664216346 6.97011719460953483, 4.52685547336740424 6.97011719460953483, 4.52685547336740424 0.32246093782891094, 0.38447265664216346 0.32246093782891094))'

B = 'Polygon ((8.74365235266854768 0.32246093782891094, 8.74365235266854768 6.99492188213483512, 12.96044923196968668 6.99492188213483512, 12.96044923196968668 0.32246093782891094, 11.94345704343235326 0.32246093782891094, 11.94345704343235326 3.42304687849151534, 10.0086914164588876 3.42304687849151534, 10.0086914164588876 0.32246093782891094, 8.74365235266854768 0.32246093782891094))'

o1 = 'MultiPolygon (((4.52685547336740246 5.80741262947195125, 4.52685547336740246 6.97011719460953127, 5.82619646979413286 6.97011719460953127, 5.82619646979413286 5.80741262947195125, 4.52685547336740246 5.80741262947195125)))'

rim = RIM(A,B,o1)

translated_rim = RIM_translate(str(rim[0]))

print(translated_rim, rim)