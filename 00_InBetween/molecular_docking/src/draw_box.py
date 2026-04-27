from pymol.cgo import *
from pymol import cmd

def draw_box(center_x, center_y, center_z, size_x, size_y, size_z, name="vina_box"):
    """
    Draws a wireframe box in PyMOL representing the AutoDock Vina search space.
    """
    cx, cy, cz = float(center_x), float(center_y), float(center_z)
    sx, sy, sz = float(size_x), float(size_y), float(size_z)
    
    # Calculate the corners of the box
    minX, maxX = cx - sx/2, cx + sx/2
    minY, maxY = cy - sy/2, cy + sy/2
    minZ, maxZ = cz - sz/2, cz + sz/2
    
    # Define the 12 lines of the box using PyMOL's Compiled Graphics Objects (CGO)
    obj = [
        BEGIN, LINES,
        # Bottom face
        VERTEX, minX, minY, minZ,  VERTEX, maxX, minY, minZ,
        VERTEX, minX, minY, maxZ,  VERTEX, maxX, minY, maxZ,
        VERTEX, minX, minY, minZ,  VERTEX, minX, minY, maxZ,
        VERTEX, maxX, minY, minZ,  VERTEX, maxX, minY, maxZ,
        # Top face
        VERTEX, minX, maxY, minZ,  VERTEX, maxX, maxY, minZ,
        VERTEX, minX, maxY, maxZ,  VERTEX, maxX, maxY, maxZ,
        VERTEX, minX, maxY, minZ,  VERTEX, minX, maxY, maxZ,
        VERTEX, maxX, maxY, minZ,  VERTEX, maxX, maxY, maxZ,
        # Vertical connecting lines
        VERTEX, minX, minY, minZ,  VERTEX, minX, maxY, minZ,
        VERTEX, maxX, minY, minZ,  VERTEX, maxX, maxY, minZ,
        VERTEX, minX, minY, maxZ,  VERTEX, minX, maxY, maxZ,
        VERTEX, maxX, minY, maxZ,  VERTEX, maxX, maxY, maxZ,
        END
    ]
    # Load the drawing into PyMOL
    cmd.load_cgo(obj, name)
    print(f"Drew box '{name}' at ({cx}, {cy}, {cz}) with size [{sx}, {sy}, {sz}]")

# Register the command in PyMOL
cmd.extend("draw_box", draw_box)
