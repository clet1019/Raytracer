/* Written by Duncan Levear in Spring 2023 for CS3333 at Boston College */

import { default as YAML } from './js-yaml.js'
import { Vector3 } from './Vector3.js'

function parseVectors(json) {
    /*
    Recursively convert json values to vectors if they starts with 'v3_'. 
    */
    for (const k in json) {
        if (k.startsWith('v3_')) {
            json[k] = new Vector3(json[k]);
        }
        else if (k.startsWith('j_')) {
            json[k] = parseVectors(json[k]);
        }
        else if (k.startsWith('a_')) {
            // assume array of JSON
            for (let i=0; i < json[k].length; i++) {
                json[k][i] = parseVectors(json[k][i]);
            }
        }
    }
    return json;
}

export function parseSceneYaml(s_yaml) {
    const raw_parsed = YAML.load(s_yaml);
    return parseVectors(raw_parsed);
}

export function parseSTLarrayBuffer(b, v3_scale, v3_translate) {
    const view = new DataView(b);
    // nTri is an integer contained in bytes at index 80,81,82,83
    const nTri = view.getInt32(80, true);
    const triangles = new Array(nTri); // the return value
    const xyz=['x','y','z']; // trick for using xyz in for-loop

    // helper for easier formulas 
    const bytesPerVertex = 12 + 36 + 2; 
    for (let i=0; i < nTri; i++) {
        const startingByte = 84 + i*bytesPerVertex;
        const [nX, nY, nZ] = [
            view.getFloat32(startingByte, true), 
            view.getFloat32(startingByte+4, true),
            view.getFloat32(startingByte+8, true)
        ];
        const vertices = [[0,0,0],[0,0,0],[0,0,0]];
        for (let j=0; j < 3; j++) {
            // read v_{j}
            const vertexStartingByte = startingByte + 12 + 12*j;
            for (let c=0; c < 3; c++) {
                // c = 0 means x-coordinate, 1 means y-coordinate, 2 means z-coordinate
                vertices[j][c] = view.getFloat32(vertexStartingByte + 4*c, true); 
                vertices[j][c] *= v3_scale[xyz[c]];
                vertices[j][c] += v3_translate[xyz[c]];
            }
        }
        triangles[i] = {
            v3_pt0: new Vector3(vertices[0]),
            v3_pt1: new Vector3(vertices[1]),
            v3_pt2: new Vector3(vertices[2]),
        };
    }
    return triangles;
}
