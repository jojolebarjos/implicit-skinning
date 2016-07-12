# Blender script to export the positions of the spheres in csv
import bpy

file = open("centers.csv", "w")
for i in range(1, 10):
    pos = bpy.data.objects["C0" + str(i)].matrix_world.translation
    out = ",".join([("%.6f" % component) for component in pos])
    file.write(out + "\n")
for i in range(10, 31):
    pos = bpy.data.objects["C" + str(i)].matrix_world.translation
    out = ",".join([("%.6f" % component) for component in pos])
    file.write(out + "\n")

file.close()
