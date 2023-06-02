/* Written by Duncan Levear in Spring 2023 for CS3333 at Boston College */
export class Vector3 {
    /*
    Can construct with any of the following:
    const v1 = new Vector3();
    const v2 = new Vector3(0, 0, 0);
    const v3 = new Vector3([0, 0, 0]);
    const v4 = new Vector3(v3);
    
    Most methods return `this`, so you can also do e.g.
    const v5 = new Vector3(1,1,1).normalize().scaleBy(-1);
    */
    constructor(x, y, z) {
        // if no arguments provided, default to [0,0,0]
        if (x === undefined) {
            x = [0,0,0];
        }
        
        // check if x is an array as in `new Vector3([1,2,3])`
        if (x.length !== undefined) {
            [this.x, this.y, this.z] = x;
        }
        // check if x is a vector as in `new Vector3(v3)`
        else if (x.x !== undefined) {
            [this.x, this.y, this.z] = [x.x, x.y, x.z];
        }
        // assume we are in the `new Vector3(1, 2, 3)` case
        else {
            [this.x, this.y, this.z] = [x, y, z];
        }
        if (this.x === undefined || this.y === undefined || this.z === undefined) {
            console.error("Warning, Vector3 constructor failed");
            console.error(this);
        }
    }
    increaseBy(w) {
        [this.x, this.y, this.z] = [this.x + w.x, this.y + w.y, this.z + w.z];
        return this;
    }
    increaseByMultiple(w, t) {
        if (!(typeof(t) === "number")) {
            console.error(`Warning, invalid argument ${t} to increaseByMultiple()`);
        }
        [this.x, this.y, this.z] = [this.x + t*w.x, this.y + t*w.y, this.z + t*w.z];
        return this;
    }
    sub(w) {
        this.increaseByMultiple(w, -1);
    }
    setAsWeightedSum(v, alpha, w, beta) {
        if (!(typeof(alpha) === "number") || !(typeof(beta) === "number")) {
            console.error(`Warning, invalid argument  ${alpha} or ${beta} to setAsWeightedSum()`);
        }
        this.x = v.x*alpha + w.x*beta;
        this.y = v.y*alpha + w.y*beta;
        this.z = v.z*alpha + w.z*beta;
        return this;
    }
    scaleBy(alpha) {
        if (!(typeof(alpha) === "number")) {
            console.error(`Warning, invalid argument  ${alpha} or ${beta} to scaleBy()`);
        }
        [this.x, this.y, this.z] = [alpha*this.x, alpha*this.y, alpha*this.z];
        return this;
    }
    normalize() {
        const length = Math.sqrt(this.dotProduct(this));
        this.scaleBy(1/length);
        return this;
    }
    dotProduct(w) {
        return this.x*w.x + this.y*w.y + this.z*w.z;
    }
    norm() {
        return Math.sqrt(this.dotProduct(this));
    }
    crossProduct(b) {
        const ret = new Vector3();
        const a = this;
        ret.x = a.y*b.z - a.z*b.y;
        ret.y = a.z*b.x - a.x*b.z;
        ret.z = a.x*b.y - a.y*b.x;
        return ret;
    }
    clampAll(min, max) {
        this.x = Math.max(min, Math.min(this.x, max));
        this.y = Math.max(min, Math.min(this.y, max));
        this.z = Math.max(min, Math.min(this.z, max));
    }
    roundAll() {
        this.x = Math.round(this.x);
        this.y = Math.round(this.y);
        this.z = Math.round(this.z);
    }
    toString() {
        return `[${this.x}, ${this.y}, ${this.z}]`
    }
}

export function vectorSum(a, b) {
    return new Vector3().setAsWeightedSum(a, 1, b, 1);
}

export function vectorDifference(a, b) {
    return new Vector3().setAsWeightedSum(a, 1, b, -1);
}

export function vectorScaled(v, alpha) {
        if (!(typeof(alpha) === "number")) {
            console.error(`Warning, invalid argument ${alpha} to vectorScaled()`);
        }
    return new Vector3(v).scaleBy(alpha);
}
