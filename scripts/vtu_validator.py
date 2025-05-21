
import xml.etree.ElementTree as ET

file = "fort_frame_0000.vtu"

xml_string = "";

#open file
with open(file, 'rb') as f:
    content = f.read()

#find the position of <AppendedData> tag
start = content.find(b"<AppendedData")
if start != -1:
    xml_string += content[:start].decode('utf-8')
start = content.find(b"_", start) + 1

xml_string += "</VTKFile>"

print(xml_string)
#parse the xml string 
root = ET.fromstring(xml_string)

piece = root.find("UnstructuredGrid/Piece");
num_points = int(piece.get("NumberOfPoints"))
num_cells = int(piece.get("NumberOfCells"))
print("Number of points: ", num_points)
print("Number of cells: ", num_cells)
#for each DataArray tag, get the name and offset
position_offset = -1
for array in root.findall("UnstructuredGrid/Piece/Points/DataArray"):
    position_offset = int(array.get("offset"))

connectivity_offset = -1
offsets_offset = -1
types_offset = -1
for array in root.findall("UnstructuredGrid/Piece/Cells/DataArray"):
    name = array.get("Name")
    offset = array.get("offset")
    if(name == "connectivity"):
        connectivity_offset = offset
    elif(name == "offsets"):
        offsets_offset = offset
    elif(name == "types"):
        types_offset = offset
    else:
        print("Unknown DataArray: ", name)
        exit(1)
    print("Name: ", name, " Offset: ", offset)

def sizeof(type):
    #end of type string ends with bits
    if type.endswith("8"):
        return 1
    elif type.endswith("16"):
        return 2
    elif type.endswith("32"):
        return 4
    elif type.endswith("64"):
        return 8
    else:
        print("Unknown type: ", type)
        exit(1)

for array in root.findall("UnstructuredGrid/Piece/CellData/DataArray"):
    name = array.get("Name")
    offset = int(array.get("offset"))
    type = array.get("type")
    print("Name: ", name, " Offset: ", offset, " Type: ", type)
    epxected_size = num_cells * sizeof(type)
    size = int.from_bytes(content[(start+offset):(start+offset+8)], "little")
    if(size != epxected_size):
        print("!!!!!! Size mismatch, expected size: ", epxected_size, " actual size: ", size)
    print()

# from start interpret as uint64 little endian
position_length = int.from_bytes(content[(position_offset+start):(position_offset+start+8)], "little")
print ("Length: ", position_length)
expected_length = num_points * 3 * 8
print ("Expected length: ", expected_length)