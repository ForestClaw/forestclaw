function tests(filename)
    patches = FCReadAscii.readFile(filename);
    mesh = FCMesh(patches);
    geom = FCMeshGeometry(mesh);
    trisurf(geom.T, geom.x, geom.y, geom.z, "EdgeColor", 'none')
end