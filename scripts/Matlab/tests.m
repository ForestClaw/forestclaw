% function tests(filename)
    filename = 'fort.q0000';
    patches = FCReadAscii.readFile(filename);
    mesh = FCMesh(patches);
    geom = FCMeshGeometry(mesh);
    trisurf(geom.T, geom.x, geom.y, geom.z, ...
        "EdgeColor", 'black', ...
        'FaceColor','interp')
%     tricontour(geom.T, geom.x, geom.y, geom.z, 10)
% end