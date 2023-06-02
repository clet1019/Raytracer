import {Vector3, vectorSum, vectorDifference, vectorScaled} from './Vector3.js'

const EPSILON = 1e-9
let recursion_amount = 0


export class RayTracer {
    constructor(sceneInfo, image) {
        this.scene = sceneInfo;
        this.image = image;
        // clear image all white
        for (let i = 0; i < image.data.length; i++) {
            image.data[i] = 255;
        }
    }

    putPixel(row, col, r, g, b) {
        /*
        Update one pixel in the image array. (r,g,b) are 0-255 color values.
        */
        if (Math.round(row) != row) {
            console.error("Cannot put pixel in fractional row");
            return;
        }
        if (Math.round(col) != col) {
            console.error("Cannot put pixel in fractional col");
            return;
        }
        if (row < 0 || row >= this.image.height) {
            return;
        }
        if (col < 0 || col >= this.image.width) {
            return;
        }

        const index = 4 * (this.image.width * row + col);
        this.image.data[index + 0] = Math.round(r);
        this.image.data[index + 1] = Math.round(g);
        this.image.data[index + 2] = Math.round(b);
        this.image.data[index + 3] = 255;
    }

    render() {
        this.root = new AABBNode(this.scene.a_geometries);
        const [e, u, v, w] = this.setupCamera()
        for (let row=0; row < this.image.height; row++) {
            for (let col=0; col < this.image.width; col++) {
                const ray = this.pixelToRay(row,col, e,u, v, w)
                const color = this.traceRay(ray)
                color.scaleBy(255)
                this.putPixel(row, col, color.x, color.y, color.z)
            }
        }
    }
    
    setupCamera(e, eyeOut, up) {
        e = this.scene.v3_eye;
        eyeOut = this.scene.v3_eyeOut;
        up = this.scene.v3_up;

        const w = vectorScaled(eyeOut, -1);      
        const u = up.crossProduct(w);
        const v = w.crossProduct(u);

        u.normalize(), v.normalize(), w.normalize();

        return [e, u, v, w];
    }
        
    pixelToRay(rowIndex, columnIndex, eye, uVector, vVector, wVector) {
        const eyePosition = eye;
        
        // distance from the camera and image plane
        const imagePlaneDistance = this.scene.f_imageplaneDistance;

        // vector opposite of the view direction
        const wVectorScaled = vectorScaled(wVector, imagePlaneDistance);
    
        // this calculates the top left corner of the plane
        const topLeftCorner = vectorDifference(
            vectorSum(vectorDifference(eyePosition, wVectorScaled), vectorScaled(vVector, this.scene.f_imageplaneHeight / 2)),
            vectorScaled(uVector, this.scene.f_imageplaneWidth / 2)
        );

        // dividing image plane by num of rows and columns to get the height and width of each pixel
        const pixelHeight = this.scene.f_imageplaneHeight / this.scene.i_height;
        const pixelWidth = this.scene.f_imageplaneWidth / this.scene.i_width;
    
        // get the center of the first pixel in the top left
        const firstPixelCenter = vectorSum(
            vectorDifference(topLeftCorner, vectorScaled(vVector, 0.5 * pixelHeight)),
            vectorScaled(uVector, 0.5 * pixelWidth)
        );
    
        // gets the first pixel center and then translates along the u and v vectors
        const targetPixelCenter = firstPixelCenter
            .increaseByMultiple(vectorScaled(vVector, rowIndex * pixelHeight), -1)
            .increaseBy(vectorScaled(uVector, columnIndex * pixelWidth));
    
        // gets the direction of the ray
        const rayDirection = vectorDifference(targetPixelCenter, eyePosition);

        // generate new ray with the distance from the camera and the direction of the ray then return it
        const ray = new Ray(eyePosition, rayDirection);
        return ray;
    }
    
    traceRay(ray) {
        const nodeIntersections = ray.hitsOfNode(this.root);
    
        let nearestIntersection = {t: Infinity};
        for (const intersectionPoint of nodeIntersections) {
            if (intersectionPoint.t < nearestIntersection.t && intersectionPoint.t > EPSILON) {
                nearestIntersection = intersectionPoint;
            }
        }
    
        if (nearestIntersection.t === Infinity) {
            return new Vector3(0, 0, 0);
        }
    
        const resultingColor = this.getColor(nearestIntersection);
        return resultingColor;
        
    }
    
    getColor(record) {     
        let color_added = new Vector3 (0,0,0);
        let reflectedAmount = 0;
        let reflectedColor = new Vector3 (0,0,0);

        
        const lights = this.scene.a_lights;
        for (let i = 0; i <lights.length; i++) {
            const light = lights[i];
            const diffuse_specular = this.whatLight(record, light);
            color_added.increaseBy(diffuse_specular);
        }
        // add if, black scr
            if (record.struckGeometry.j_material.f_reflectance>0 && recursion_amount < 5 ) {
                recursion_amount +=1 ;
                reflectedAmount = record.struckGeometry.j_material.f_reflectance;
                reflectedColor = this.reflected(record);
                recursion_amount -=1;
                color_added.increaseBy(reflectedColor.scaleBy(reflectedAmount));
            }
        return color_added;
    }
 
    
    // To add shading, break it into steps: whatLight(), diffuse(), highlight(), or similar
    whatLight(intersection, lightSource) {
        const intersectionPoint = intersection.pt;
        const shadowRayDirection = vectorDifference(lightSource.v3_position, intersectionPoint);
        const shadowRay = new Ray(intersectionPoint, shadowRayDirection);
        const shadowIntersections = shadowRay.allHits(this.scene.a_geometries);
    
        for (const shadowIntersection of shadowIntersections) {
            if (shadowIntersection.t > 0.0001 && shadowIntersection.t < 1) {
                return new Vector3([0, 0, 0]);
            }
        }
    
        const lightPosition = lightSource.v3_position;
        const toLightVector = vectorDifference(lightPosition, intersectionPoint).normalize();
        const diffuseColor = this.diffuse(intersection, lightSource, toLightVector);
        const specularLight = this.specular(intersection, lightSource, toLightVector);
        const finalColor = new Vector3(0, 0, 0);
        
        finalColor.increaseBy(diffuseColor);
        finalColor.increaseBy(specularLight);
        
        return finalColor;
    }
        
    reflected(hit) {           
        const mirrorRay = this.bounce(hit);
        const mirrorHit = this.traceRay(mirrorRay);
        return mirrorHit;
    }
    
    bounce(hit) {
        const viewingRayDirection = hit.ray.dir;
        const surfaceNormal = hit.normal;
    
        const inverseViewingRayDirection = new Vector3(viewingRayDirection).scaleBy(-1);
        const dotProductInverseNormal = inverseViewingRayDirection.dotProduct(surfaceNormal);
        const dotProductNormalNormal = surfaceNormal.dotProduct(surfaceNormal);
    
        const scaledSurfaceNormal = new Vector3(surfaceNormal).scaleBy(2 * (dotProductInverseNormal / dotProductNormalNormal));
        const bouncedDirection = vectorDifference(scaledSurfaceNormal, inverseViewingRayDirection);
    
        const bouncedRay = new Ray(hit.pt, bouncedDirection);
        return bouncedRay;
    }

    // for diffuse and specular I got help from classmate Anirudh
    // I was having issues with taking the lecture material and coding it but he showed me the logic
    diffuse(hit, light, toLight) { 
        const normal = hit.normal;
        const alignment = toLight.dotProduct(normal);
        if (alignment < 0){
            return new Vector3(0,0,0);
        }
        const color = vectorScaled(hit.struckGeometry.j_material.v3_diffuse, alignment);
        color.scaleBy(light.f_intensity);
        return color;  
    }
    
     // this is for Phong
    specular (hit, light_source, toLight) {      
        const point = hit.pt;
        const specularity_power = hit.struckGeometry.j_material.f_specularity;
        
        if (specularity_power < 0 || specularity_power === undefined){
            return new Vector3(0,0,0)
        }
        
        const to_eye = vectorDifference(this.scene.v3_eye, point);
        const normal = hit.normal;
        
        const alpha = 2 * toLight.dotProduct(normal);
        const outgoingLight = vectorDifference(vectorScaled(normal,alpha), toLight);
        outgoingLight.normalize();
        to_eye.normalize();
            
        let s = outgoingLight.dotProduct(to_eye);
        if (s < 0) {
            s = 0;
        }
        
        s = Math.pow(s,specularity_power);
        const final_spec = new Vector3(1,1,1);
        final_spec.scaleBy(s*light_source.f_intensity);
        return final_spec;
    }
}

class Ray {
    constructor(start, dir) {
        this.start = start; 
        this.dir = dir; 
    }


    tToPt(t) {
        const ret = new Vector3(this.start).increaseByMultiple(this.dir, t);
        return ret;
    }
    
    allHits(geometries) {
        let ret = [];
        for (const g of geometries) {
            const record = this.hit(g);
            if (record.length === undefined) {
                console.error("Return type of hit() should be an array.");
            }
            ret = ret.concat(record);
        }
        return ret;
    }
    
    hit(g) {
        if (g.s_type === 'sphere') {
            return this.hitSphere(g);
        }
        else if (g.s_type === 'sheet') {
            return this.hitSheet(g);
        }
        else if (g.s_type === 'box') {
            return this.hitBox(g);
        }
        else if (g.s_type === 'cylinder') {
            return this.hitCylinder(g);
        }
        else if (g.s_type === 'triangle') {
            return this.hitTriangle(g);
        }
        else {
            console.error("Shape of type " + g.s_type + " is not supported");
        }
    }
    
    hitSheet(g) {
    
        const pt0 = g.v3_pt0;
        const pt1 = g.v3_pt1;
        const pt2 = g.v3_pt2;
        // compute d, normal, edge1, edge2 once only, to save time
        if (g.edge1 === undefined) {
            g.edge1 = vectorDifference(pt0, pt1);
            g.edge2 = vectorDifference(pt2, pt1);

            const unit1 = vectorDifference(pt0, pt1).normalize();
            const unit2 = vectorDifference(pt2, pt1).normalize();
            if (Math.abs(unit1.dotProduct(unit2)) > 0.01) {
                console.error(`Edges ${edge1} and ${edge2} are not orthogonal`);
            }

            g.normal = unit2.crossProduct(unit1);
            g.normal.normalize();

            
            g.d = g.normal.dotProduct(pt1);
        }
        
        const t = (g.d - g.normal.dotProduct(this.start))/g.normal.dotProduct(this.dir);
        const pt = this.tToPt(t);
        
        let alpha = vectorDifference(pt,pt1).dotProduct(g.edge1);
        alpha /= g.edge1.dotProduct(g.edge1);
        let beta = vectorDifference(pt,pt1).dotProduct(g.edge2);
        beta /= g.edge2.dotProduct(g.edge2);

        if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1) {
            return [];
        }
  
        const hit = new HitRecord(this, t, pt, g, g.normal)
        return [hit];
        
        
    }

    hitSphere(geometry) {
        const sphereRadius = geometry.f_radius;
        const startToCenter = vectorDifference(this.start, geometry.v3_center);
        const direction = this.dir;
        const a = direction.dotProduct(direction);
        const b = 2 * (startToCenter.dotProduct(direction));
        const c = startToCenter.dotProduct(startToCenter) - (sphereRadius * sphereRadius);
    
        const discriminant = (b * b) - (4 * a * c);
        if (discriminant < 0) {
            const t = (-b + Math.sqrt(discriminant)) / (2 * a);
            const intersectionPoint = this.tToPt(t);
            const intersectionNormal = vectorDifference(intersectionPoint, geometry.v3_center).normalize();
            
            return [];
        }
    
        const t1 = (-b + Math.sqrt(discriminant)) / (2 * a);
        const t2 = (-b - Math.sqrt(discriminant)) / (2 * a);
        const intersectionPoint1 = this.tToPt(t1);
        const intersectionPoint2 = this.tToPt(t2);
    
        const intersectionNormal1 = vectorDifference(intersectionPoint1, geometry.v3_center).normalize();
        const intersectionNormal2 = vectorDifference(intersectionPoint2, geometry.v3_center).normalize();
    
        const hitRecord1 = new HitRecord(this, t1, intersectionPoint1, geometry, intersectionNormal1);
        const hitRecord2 = new HitRecord(this, t2, intersectionPoint2, geometry, intersectionNormal2);
    
        if (t1 < t2) {
            return [hitRecord1, hitRecord2];
        } else {
            return [hitRecord2, hitRecord1];
        }
    }

    // from hw4 sols
    hitCylinder(g) {
        const center = g.v3_center;
        const height = g.f_height;
        const radius = g.f_radius;
        
        const testSphere = {
            v3_center: new Vector3(center.x, 0, center.z),
            f_radius: radius,
        }
        const testRay = new Ray(new Vector3(this.start.x, 0, this.start.z), new Vector3(this.dir.x, 0, this.dir.z));
        const hits = testRay.hitSphere(testSphere)
        const ret = []
        for (const h of hits) {
            const pt = this.tToPt(h.t)
            if (center.y - height/2 < pt.y && pt.y < center.y + height/2) {
                h.struckGeometry = g
                h.pt = pt
                h.ray = this
                ret.push(h)
            }
        }
        for (const m of [1,-1]) {
            const capCenter = new Vector3(center.x, center.y + m*height/2, center.z)
            const normal = new Vector3(0,m,0)
            const rhs = normal.dotProduct(capCenter)
            const t = (rhs - this.start.dotProduct(normal))/this.dir.dotProduct(normal)
            const hit = new HitRecord(this, t, this.tToPt(t), g, normal)
            const ptToAxis = vectorDifference(capCenter, hit.pt)
            if (ptToAxis.norm() < radius) {
                ret.push(hit)
            }
        }
        return ret
    }

    hitBox(g) {
        const boxMin = g.v3_minPt;
        const boxMax = new Vector3(boxMin.x + g.v3_dim.x, boxMin.y + g.v3_dim.y, boxMin.z + g.v3_dim.z);
        const rayStart = this.start;
        const rayDir = this.dir;
    
        let xInterval, yInterval, zInterval;
    
        xInterval = this.calculateInterval(rayStart.x, rayDir.x, boxMin.x, boxMax.x);
        if (!xInterval) return [];
    
        yInterval = this.calculateInterval(rayStart.y, rayDir.y, boxMin.y, boxMax.y);
        if (!yInterval) return [];
    
        zInterval = this.calculateInterval(rayStart.z, rayDir.z, boxMin.z, boxMax.z);
        if (!zInterval) return [];
    
        const maxStart = Math.max(xInterval[0], yInterval[0], zInterval[0]);
        const minEnd = Math.min(xInterval[1], yInterval[1], zInterval[1]);
    
        if (maxStart <= minEnd) {
            const [t1, t2] = [maxStart, minEnd];
            const [pt1, pt2] = [this.tToPt(t1), this.tToPt(t2)];
    
            const normal1 = this.calculateNormal(pt1, boxMin, boxMax);
            const normal2 = this.calculateNormal(pt2, boxMin, boxMax);
    
            const hit1 = new HitRecord(this, t1, pt1, g, normal1);
            const hit2 = new HitRecord(this, t2, pt2, g, normal2);
    
            return [hit1, hit2];
        } else {
            return [];
        }
    }
    
    // added helpers to help with render time 
    calculateInterval(rayStart, rayDir, boxMin, boxMax) {
        if (rayDir == 0) {
            if (boxMin <= rayStart && rayStart <= boxMax) {
                return [-Infinity, Infinity];
            } else {
                return null;
            }
        } else {
            let interval = [(boxMin - rayStart) / rayDir, (boxMax - rayStart) / rayDir];
            if (interval[0] > interval[1]) {
                interval = [interval[1], interval[0]];
            }
            return interval;
        }
    }
    
    calculateNormal(point, boxMin, boxMax) {
        const EPSILON = 1e-6;
        let normal;
    
        if (Math.abs(point.x - boxMin.x) < EPSILON) {
            normal = new Vector3(-1, 0, 0);
        } else if (Math.abs(point.x - boxMax.x) < EPSILON) {
            normal = new Vector3(1, 0, 0);
        } else if (Math.abs(point.y - boxMin.y) < EPSILON) {
            normal = new Vector3(0, -1, 0);
        } else if (Math.abs(point.y - boxMax.y) < EPSILON) {
            normal = new Vector3(0, 1, 0);
        } else if (Math.abs(point.z - boxMin.z) < EPSILON) {
            normal = new Vector3(0, 0, -1);
        } else if (Math.abs(point.z - boxMax.z) < EPSILON) {
            normal = new Vector3(0, 0, 1);
        }
    
        return normal;
    }

    // I got help on the triangle code from Anirudh. 

    M(a, b, c, d, e, f, g, h, i) {
        return (a * (e * i-h * f)) + (b *(g * f-d * i)) + (c *(d * h - e * g));
    }
            
    beta(a, b, c, d, e, f, g, h, i, j, k, l) {
        return (j * (e * i - h * f)) + (k *(g * f - d * i)) + (l *(d * h-e * g));
                
    }
            
    gamma(a, b, c, d, e, f, g, h, i, j, k, l) {
        return (i * (a * k-j * b)) + (h* (j * c-a * l)) + (g *(b * l-k * c));
    }
            
    getT(a, b, c, d, e, f, g, h, i, j, k, l) {
        let min_1 = f *(a * k-j * b) + e *(j * c-a * l) + d *(b * l-k * c);
        min_1 = min_1 * -1;
            return min_1;
    }
               
    hitTriangle(g) {       
        const a = g.v3_pt0;
        const b = g.v3_pt1;
        const c = g.v3_pt2;
                
        const e = this.start
        const v = this.dir 
                
                
        const edge1 = vectorDifference(b,a)
        const edge2 = vectorDifference(b,c)
                
        const norm = edge2.crossProduct(edge1).normalize()
        const m = this.M(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z)
                
        let beta = this.beta(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z, e.x-a.x, e.y-a.y, e.z-a.z)
                
        beta = beta / m
                 
        let gamma = this.gamma(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z, e.x-a.x, e.y-a.y, e.z-a.z)
                
        gamma = gamma / m
                 
        let t = this.getT(b.x-a.x, b.y-a.y, b.z-a.z, c.x-a.x, c.y-a.y, c.z-a.z, -v.x, -v.y, -v.z, e.x-a.x, e.y-a.y, e.z-a.z)
                
        t = t / m

        const point = this.tToPt(t)
                
        const alpha = 1-beta-gamma
                
            
        if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1) {
             return [];
        }
                 
        const hit =  new HitRecord(this, t, point, g, norm);
        return [hit];
    }
   
    hitsOfNode(node) {
        const hits = this.hitBox(node.box);
        if (hits.length == 0) {
            return []; 
        }
        if (node.leaf) {
            return this.allHits(node.geometries);    
        }
        
        const left= this.hitsOfNode(node.left);
        const right = this.hitsOfNode(node.right);
        const fullHits = left.concat(right);
        
        return fullHits;
    }
}

class HitRecord {
    constructor(ray, t, pt, struckGeometry, normal) {
        this.ray = ray; // ray that was involved
        this.t = t; // t-value of intersection along ray
        this.pt = pt; // vector3, point where the ray hit
        this.struckGeometry = struckGeometry; // object that was hit
        this.normal = normal; // normal vector of struckGeometry at pt
    }
}


// AABBNode class was mostly formulated using the extra help steps from hw5 page
class AABBNode {
    constructor(g) {
        let boundingBox = null;
        for (const geometry of g) {
            if (boundingBox === null) {
                boundingBox = this.boxAround(geometry);
            } else {
                boundingBox = this.boxAround2(this.boxAround(geometry), boundingBox);
            }
        }

        if (g.length <= 3) {
            this.geometries = g;
            this.left = null;
            this.right = null;
            this.box = boundingBox;
            this.leaf = true;
            return;
        }

        this.box = boundingBox;
        this.leaf = false;

        const boxDimensionX = this.box.v3_dim.x;
        const boxDimensionY = this.box.v3_dim.y;
        const boxDimensionZ = this.box.v3_dim.z;

        const maxDimension = Math.max(boxDimensionX, boxDimensionY, boxDimensionZ);

        let sortedGeometries = [];

        if (maxDimension === boxDimensionX) {
            sortedGeometries = g.sort((a, b) => this.boxAround(a).v3_minPt.x - this.boxAround(b).v3_minPt.x);
        } else if (maxDimension === boxDimensionY) {
            sortedGeometries = g.sort((a, b) => this.boxAround(a).v3_minPt.y - this.boxAround(b).v3_minPt.y);
        } else {
            sortedGeometries = g.sort((a, b) => this.boxAround(a).v3_minPt.z - this.boxAround(b).v3_minPt.z);
        }

        const leftGeometries = sortedGeometries.slice(0, sortedGeometries.length / 2);
        const rightGeometries = sortedGeometries.slice(sortedGeometries.length / 2, sortedGeometries.length);
        this.left = new AABBNode(leftGeometries);
        this.right = new AABBNode(rightGeometries);
    }

    boxAround(g) {
        if (g.s_type === 'sphere') {
            return this.boxSphere(g);
        }
        else if (g.s_type === 'sheet') {
            return this.boxSheet(g);
        }
        else if (g.s_type === 'box') {
            return g;
        }
        else if (g.s_type === 'cylinder') {
            return this.boxCylinder(g);
        }
        else if (g.s_type === 'triangle') {
            return this.boxTriangle(g);
        }
        else {
            console.error("Shape of type " + g.s_type + " is not supported");
        }
        
    }
   
    boxAround2(firstBox, secondBox) {
        const minX = Math.min(firstBox.v3_minPt.x, secondBox.v3_minPt.x);   
        const minY = Math.min(firstBox.v3_minPt.y, secondBox.v3_minPt.y);
        const minZ = Math.min(firstBox.v3_minPt.z, secondBox.v3_minPt.z);

        const maxX = Math.max(firstBox.v3_minPt.x + firstBox.v3_dim.x, secondBox.v3_minPt.x + secondBox.v3_dim.x);
        const maxY = Math.max(firstBox.v3_minPt.y + firstBox.v3_dim.y, secondBox.v3_minPt.y + secondBox.v3_dim.y);
        const maxZ = Math.max(firstBox.v3_minPt.z + firstBox.v3_dim.z, secondBox.v3_minPt.z + secondBox.v3_dim.z);
        
         const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(minX, minY, minZ),
            v3_dim: new Vector3(maxX-minX, maxY-minY, maxZ-minZ),
        };
        
        return ret
    }

   boxPoints(pts) {
        let [minX, minY, minZ] = [Infinity, Infinity, Infinity];
        let [maxX, maxY, maxZ] = [-Infinity, -Infinity, -Infinity];
        for (const p of pts) {
            minX = Math.min(minX, p.x);
            minY = Math.min(minY, p.y);
            minZ = Math.min(minZ, p.z);
            maxX = Math.max(maxX, p.x);
            maxY = Math.max(maxY, p.y);
            maxZ = Math.max(maxZ, p.z);
        }
        const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(minX, minY, minZ),
            v3_dim: new Vector3(maxX-minX, maxY-minY, maxZ-minZ),
        };
        return ret;
    }

    boxSheet(g) {
        return this.boxPoints([g.v3_pt0, g.v3_pt1, g.v3_pt2, vectorSum(g.v3_pt2, vectorDifference(g.v3_pt0, g.v3_pt1))]);
    }
    
    boxCylinder(g) {
        const ret = {
            s_type: 'box',
            v3_minPt: vectorSum(g.v3_center, new Vector3(-g.f_radius, -g.f_height/2, -g.f_radius)),
            v3_dim: new Vector3(2 * g.f_radius, g.f_height, 2 * g.f_radius),
        };
        return ret;
    }
    
    boxSphere(g) {
        const ret = {
            s_type: 'box',
            v3_minPt: new Vector3(g.v3_center.x - g.f_radius, g.v3_center.y - g.f_radius, g.v3_center.z - g.f_radius),
            v3_dim: new Vector3(2 * g.f_radius, 2 * g.f_radius, 2 * g.f_radius)
        };
        return ret;    
    }
    
    boxTriangle(triangle) {
        const minX = Math.min(triangle.v3_pt0.x, triangle.v3_pt1.x, triangle.v3_pt2.x);
        const minY = Math.min(triangle.v3_pt0.y, triangle.v3_pt1.y, triangle.v3_pt2.y);
        const minZ = Math.min(triangle.v3_pt0.z, triangle.v3_pt1.z, triangle.v3_pt2.z);
    
        const maxX = Math.max(triangle.v3_pt0.x, triangle.v3_pt1.x, triangle.v3_pt2.x);
        const maxY = Math.max(triangle.v3_pt0.y, triangle.v3_pt1.y, triangle.v3_pt2.y);
        const maxZ = Math.max(triangle.v3_pt0.z, triangle.v3_pt1.z, triangle.v3_pt2.z);
    
        return {
            s_type: 'box',
            v3_minPt: new Vector3(minX, minY, minZ),
            v3_dim: new Vector3(maxX - minX, maxY - minY, maxZ - minZ),
        };
    }
}