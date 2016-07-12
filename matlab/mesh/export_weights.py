# Blender script to export the skinning weights in csv
import bpy

vertices = bpy.data.objects["Hand"].data.vertices
group_names = [g.name for g in bpy.data.objects["Hand"].vertex_groups]
bones_count = len(group_names)

file = open("weights.csv", "w")
str = ",".join([name for name in group_names])
file.write(str + "\n")

for v in vertices:
    weights = [0 for i in range(bones_count)]
    for g in v.groups:
        weights[g.group] = g.weight
    str = ",".join([("%.6f" % w) for w in weights])
    file.write(str + "\n")
file.close()
